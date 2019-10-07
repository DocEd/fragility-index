library(tidyverse)

ate <- function(t, y) {
  (sum(y[t == 0]) - sum(y[t == 1]))/(length(y)/2)
}

simulate_fragility <- function(intervention_mort,
                               control_mort,
                               power = 0.8,
                               N = NULL,
                               sims = 1000,
                               alpha = 0.05,
                               permit_negative = FALSE) {
  
  if (intervention_mort > control_mort) {
    rlang::abort("Your intervention mortality is greater than your control
                 mortality. The Fragility Index is not defined under these
                 conditions")
  }
  
  if (is.null(N)) {
    N <- as.integer(
      power.prop.test(
        p1 = intervention_mort,
        p2 = control_mort, power = power)$n)*2
  }
  
  if (N%%2 == 1) {
    # Ensure study uses even number of patients for balance.
    N <- N + 1L
  }
  
  # Set up objects to store values we are interested in
  p_value <- rep(NA, sims)
  fragility <- rep(NA, sims)

  df <- tibble::tibble(
    t = as.list(NULL),
    y = as.list(NULL)
  )
  
  # Loop over simulations
  for (i in seq_len(sims)){
    # Randomly assign half to the treatment, and half to control
    # I want exactly half in each group, so will use sample to
    # guarentee this behaviour
    treatment <- sample(1:N, size = N/2, replace = FALSE)
    
    # Vector for the treatment assignment
    t <- rep(as.integer(NA), N)
    
    # Randomly assign half to treatment, half to control
    t[treatment] <- 1L
    t[-treatment] <- 0L
    
    # vector for the outcome
    y <- rep(NA, N)
    
    # Outcome with treatment, note 1 = death
    # Bernoilli trial with p = intervention_mortality
    y[treatment] <- rbinom(N/2, 1, intervention_mort)
    
    # Outcome WITHOUT treatment, note 1 = death
    # Bernoilli trial with p = control_mortality
    y[-treatment] <- rbinom(N/2, 1, control_mort)
    
    # Test the difference in groups
    test <- chisq.test(y, t)
    p <- p_value[i] <- test$p.value
    
    # Start the fragility index counter at 0
    fi <- 0
    
    # Grab the salient parts of the contingency table
    # (to make code easier to read)
    con_tbl <- as.matrix(table(t, y))
    intervention_death <- con_tbl[2, 2]
    intervention_survive <- con_tbl[2, 1]
    control_death <- con_tbl[1, 2]
    control_survive <- con_tbl[1, 1]
    
    if (p <= alpha) {
      while(p < alpha) {
        # Change the outcome for someone in the intervention arm
        intervention_death <- intervention_death + 1
        intervention_survive <- intervention_survive - 1
        if (intervention_survive < 0) {
          rlang::abort(
          message = "negative people are being generated, you probably have too
          few people in your intervention arm for this to work.")
        }
        fi <- fi + 1
        mat <- matrix(
          c(control_survive,
            control_death,
            intervention_survive,
            intervention_death),
          nrow=2, byrow = TRUE)
        test <- fisher.test(mat)
        p <- test$p.value
      }
    } else {
      if (permit_negative) {
        while(p >= alpha) {
          # Change the outcome for someone in the intervention arm
          # Technically the fragility index is not defined for "non-significant"
          # Studies... but since it is just a transformation of the p value
          # It is no less valid in this direction.
          intervention_death <- intervention_death - 1
          intervention_survive <- intervention_survive + 1
          if (intervention_death < 0) {
            rlang::abort(
              message = "negative people are being generated, you probably have too
            few people in your intervention arm for this to work.")
          }
          fi <- fi - 1
          mat <- matrix(
            c(control_survive,
              control_death,
              intervention_survive,
              intervention_death),
            nrow=2, byrow = TRUE)
          test <- fisher.test(mat)
          p <- test$p.value
        }
      } else {
        fi <- as.integer(NA)
      }
    }
    fragility[i] <- fi
    
    # Capture the study values for later
    df <- bind_rows(
      df, tibble::tibble(t = list(t), y = list(y))
      )
  }
  
  # Bind the simulated study values to p and f values.
  df <- bind_cols(
    df, tibble::tibble(p_value, fragility)
  )
  
  df <- df %>%
    mutate(
      ratio = map2_dbl(.x = t, .y = y, .f = function(a, b) {
        con_tbl <- as.matrix(table(a, b))
        intervention_survive <- con_tbl[2, 1]
        control_survive <- con_tbl[1, 1]
        return(intervention_survive/control_survive)
        }),
      reverse_effect = if_else(
        ratio < 1, TRUE, FALSE
      )
    ) %>%
    select(-ratio)
  
  df <- df %>%
    mutate(power = power,
           arr = control_mort-intervention_mort,
           n = N)
  
  return(df)
}

df_70_2 <- simulate_fragility(intervention_mort = 0.38, control_mort = 0.4, power = 0.7)
df_80_2 <- simulate_fragility(intervention_mort = 0.38, control_mort = 0.4, power = 0.8)
df_90_2 <- simulate_fragility(intervention_mort = 0.38, control_mort = 0.4, power = 0.9)
df_70_10 <- simulate_fragility(intervention_mort = 0.3, control_mort = 0.4, power = 0.7)
df_80_10 <- simulate_fragility(intervention_mort = 0.3, control_mort = 0.4, power = 0.8)
df_90_10 <- simulate_fragility(intervention_mort = 0.3, control_mort = 0.4, power = 0.9)
df_70_20 <- simulate_fragility(intervention_mort = 0.2, control_mort = 0.4, power = 0.7)
df_80_20 <- simulate_fragility(intervention_mort = 0.2, control_mort = 0.4, power = 0.8)
df_90_20 <- simulate_fragility(intervention_mort = 0.2, control_mort = 0.4, power = 0.9)

df <- bind_rows(df_70_2, df_80_2, df_90_2,
          df_70_10, df_80_10, df_90_10,
          df_70_20, df_80_20, df_90_20)

load("./data/sim.RData")
save(df, file = "./data/sim.RData")

df %>%
  filter(arr == 0.2)

# Multiple plot function
#
# ggplot objects can be passed in to plots (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(plots = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  n_plots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(n_plots/cols)),
                     ncol = cols, nrow = ceiling(n_plots/cols))
  }
  
  if (n_plots == 1) {
    rlang::abort("There is only 1 plot, don't use this function")
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in seq_len(n_plots)) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

grid_fra <- df %>%
  mutate(ARR = factor(arr), power = factor(power)) %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(y = fragility)) +
  geom_jitter(shape = 1, width = 0.2, height = 0, alpha = 0.5, aes(x = 0)) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
  facet_grid(rows = vars(power), cols = vars(ARR), labeller = label_both,
             scales = "free_x") +
  #scale_x_continuous(breaks = NULL) +
  #scale_y_continuous(breaks = seq(-25, 50, 5), labels = seq(-25, 50, 5)) +
  inspectEHR::theme_cchic() +
  theme(text = element_text(family = "Helvetica"),
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("Fragility Index") +
  coord_flip()

ggsave(filename = "./frailty_grid.svg", plot = grid_fra, height = 3, width = 9)
ggsave(filename = "./frailty_grid.png", plot = grid_fra, height = 3, width = 9)

df %>%
  mutate(ARR = factor(arr), power = factor(power)) %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(x = p_value, y = fragility)) +
  facet_grid(cols = vars(power), labeller = label_both) +
  geom_point(shape = 1, alpha = 0.25, aes(
    colour = factor(arr))) +
  inspectEHR::theme_cchic() +
  labs(colour = "Abolute Risk Reduction") +
  ylab("Fragility Index") +
  xlab("P Value") +
  theme(text = element_text(family = "Helvetica"))

df %>%
  mutate(
    arr = factor(arr),
    power = factor(power)
    ) %>%
  filter(reverse_effect == FALSE,
         power == "0.8") %>%
  ggplot(aes(x = p_value, y = fragility)) +
  geom_point(shape = 1, alpha = 0.5, aes(
    colour = arr)) +
  inspectEHR::theme_cchic() +
  labs(colour = "Abolute Risk Reduction") +
  ylab("Fragility Index") +
  xlab("P Value") +
  theme(text = element_text(family = "Helvetica"))

ggsave("./power_fragility.png", height = 3, width = 9)

setEPS()
postscript("./frailty_grid.eps", height = 3, width = 9)
grid_fra
dev.off()

fig_1_a <- ardsnet %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(y = fragility)) +
  geom_jitter(shape = 1, aes(x = 0)) +
  geom_text(aes(x = 0.3, y = -25, label = "A"), size = 8) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = seq(-25, 50, 5), labels = seq(-25, 50, 5)) +
  coord_flip() +
  inspectEHR::theme_cchic() +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

fig_1_b <- ardsnet %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(y = fragility)) +
  geom_text(aes(x = 0.3, y = -25, label = "B"), size = 8) +
  geom_boxplot() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = seq(-25, 50, 5), labels = seq(-25, 50, 5)) +
  coord_flip() +
  inspectEHR::theme_cchic() +
  ylab("Fragility Index") +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_blank())

setEPS()
postscript("./sim_frailty.eps", height = 3, width = 9)
multiplot(plots = list(fig_1_a, fig_1_b))
dev.off()

ggsave("./sim_frailty.eps", height = 2, width = 6)

fig_2 <- ardsnet %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(x = p_value, y = fragility)) +
  geom_point(shape = 1, size = 1.5) +
  ggformula::geom_spline(nknots = 10, colour = "red", size = 1.2) +
  inspectEHR::theme_cchic() +
  ylab("Fragility Index") +
  xlab("P Value") +
  theme(text = element_text(family = "Arial"))

fig_2

ggsave("./frailty_p.eps", plot = fig_2, height = 4, width = 6)

library(splines)
mod1 <- df %>%
  filter(reverse_effect == FALSE) %>%
  lm(fragility ~ n * ns(log(p_value), 3), data = .)

summary(mod1)
AIC(mod1)

df %>%
  filter(reverse_effect == FALSE) %>%
  mutate(pred = predict(mod1)) %>%
  ggplot(aes(x = fragility, y = pred)) + geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ylab("Predicted Fragility Index") +
  xlab("Calculated Fragility Index") +
  inspectEHR::theme_cchic()

ggsave("./predict.png", height = 4, width = 6)

## I also want to see what what happen for lots of indeividual studies where I
## don't necessarily know whats going on.

for (i in seq_len(1000)) {
  cm <- rbeta(1, 2, 8)
  rrr <- rbeta(1, 2, 8)
  arr <- cm * rrr
  im <- cm - arr
  pwr <- sample(seq(from = 0.7, to = 0.9, by = 0.05), size = 1)
  if (i == 1) {
    trial <- simulate_fragility(
      intervention_mort = im,
      control_mort = cm,
      power = pwr,
      sims = 1,
      permit_negative = TRUE
      )
  } else {
    new_trial <- simulate_fragility(
      intervention_mort = im,
      control_mort = cm,
      power = pwr,
      sims = 1,
      permit_negative = TRUE
    )
    trial <- bind_rows(
      trial, new_trial
      )
  }
}

trial %>%
  ggplot(aes(x = log(p_value), y = fragility, colour = arr)) +
  geom_point(alpha = 0.5)

mod2 <- trial %>%
  filter(reverse_effect == FALSE) %>%
  mutate(log_p = log(p_value)) %>%
  lm(fragility ~ n + n:arr + n:log_p + arr:log_p + n:arr:log_p, data = .)

summary(mod2)
AIC(mod2)

## From thge calculated characteristics

trial %>%
  filter(reverse_effect == FALSE) %>%
  mutate(pred = predict(mod2)) %>%
  ggplot(aes(x = fragility, y = pred)) + geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ylab("Predicted Fragility Index") +
  xlab("Calculated Fragility Index") +
  inspectEHR::theme_cchic()

equatiomatic::extract_eq(mod2)

## From the obserbed characteristics
