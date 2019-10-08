library(tidyverse)
library(splines)
# Note: need inspectEHR for custom ggplot theme
# can install using:
# devtools::install_github("cchic/inspectEHR")
# Or just comment out `inspectEHR::theme_cchic()` when plotting

set.seed(2019)

# Helper functions ====
ate_calc <- function(t, y) {
  (sum(y[t == 0]) - sum(y[t == 1]))/length(y)
}

# Multiple plot function
#
# ggplot objects can be passed in to plots (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
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

# Simulator Function ====
simulate_fragility <- function(intervention_mort,
                               control_mort,
                               power = 0.8,
                               n = NULL,
                               sims = 1000,
                               alpha = 0.05,
                               permit_negative = FALSE) {
  
  if (intervention_mort > control_mort) {
    rlang::abort(
    "Your intervention mortality is greater than your control
     mortality. The Fragility Index is not defined under these
     conditions"
    )
  }
  
  if (is.null(n)) {
    pwr <- power.prop.test(
      p1 = intervention_mort,
      p2 = control_mort,
      power = power)
    n <- as.integer(pwr$n)*2
  }
  
  # Ensure study uses even number of patients for balance.
  if (n%%2 == 1) {
    n <- n + 1L
  }
  
  # Set up objects to store values we are interested in capturing
  p_value <- rep(NA, sims)
  fragility <- rep(NA, sims)

  # Set up some list columns to store treatment allocation and outcomes
  df <- tibble::tibble(
    t = as.list(NULL),
    y = as.list(NULL)
  )
  
  # Loop over simulations
  for (i in seq_len(sims)){
    
    # Randomly assign half to the treatment, and half to control
    # I want exactly half in each group, so will use sample to
    # guarantee this behaviour
    treatment <- sample(1:n, size = n/2, replace = FALSE)
    
    # Vector for the treatment assignment
    t <- rep(as.integer(NA), n)
    
    # Randomly assign half to treatment, half to control
    t[treatment] <- 1L
    t[-treatment] <- 0L
    
    # vector for the outcome
    y <- rep(NA, n)
    
    # Outcome WITH treatment, note 1 = death, 0 = survive
    # Bernoilli trial with p = intervention_mortality
    y[treatment] <- rbinom(n = n/2, size = 1, prob = intervention_mort)
    
    # Outcome WITHOUT treatment, note 1 = death, 0 = survive
    # Bernoilli trial with p = control_mortality
    y[-treatment] <- rbinom(n = n/2, size = 1, prob = control_mort)
    
    # Test the difference in groups
    test <- chisq.test(y, t)
    p <- p_value[i] <- test$p.value
    
    # Start the fragility index counter at 0
    fi <- 0
    
    # Grab the salient parts of the contingency table
    # (convenience to make code easier to read)
    con_tbl <- as.matrix(table(t, y))
    intervention_death <- con_tbl[2, 2]
    intervention_survive <- con_tbl[2, 1]
    control_death <- con_tbl[1, 2]
    control_survive <- con_tbl[1, 1]
    
    # If "significant" start FI procedure.
    if (p <= alpha) {
      while(p < alpha) {
        # Change the outcome for someone in the intervention arm
        intervention_death <- intervention_death + 1
        intervention_survive <- intervention_survive - 1
        if (intervention_survive < 0) {
          rlang::abort(
          "negative people are being generated, you probably have too
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
          # The fragility index was not originally defined for "non-significant"
          # Studies... but since it is just a transformation of the p value
          # It is no less valid in this direction.
          intervention_death <- intervention_death - 1
          intervention_survive <- intervention_survive + 1
          if (intervention_death < 0) {
            rlang::abort(
            "negative people are being generated, you probably have too
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
    
    # Capture the treatment allocation and outcome for later.
    df <- bind_rows(
      df, tibble::tibble(t = list(t), y = list(y))
      )
  }
  
  # Bind the simulated study values to p and fragility values.
  df <- bind_cols(
    df, tibble::tibble(p_value, fragility)
  )
  
  # Occationally, the reverse effect is seen (randomness). We need
  # to add a flag for this. There is nothing wrong with this, but it
  # creates a symmetrical discontinuity in the FI.
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
           c_y = control_mort,
           t_y = intervention_mort,
           arr = control_mort-intervention_mort,
           rrr = (control_mort-intervention_mort)/control_mort,
           n = n,
           ate = purrr::map2_dbl(t, y, ate_calc))
  
  return(df)
}

# Core Scenarios ====
scenarios <- expand.grid(c(0.02, 0.1, 0.2), c(0.7, 0.8, 0.9))

for (i in seq_len(nrow(scenarios))) {
  if (i == 1L) {
    df <- simulate_fragility(
      intervention_mort = 0.4-scenarios[i, 1],
      control_mort = 0.4,
      power = scenarios[i, 2],
      permit_negative = TRUE)
  } else {
    df <- bind_rows(
      df,
      simulate_fragility(
        intervention_mort = 0.4-scenarios[i, 1],
        control_mort = 0.4,
        power = scenarios[i, 2],
        permit_negative = TRUE))
  }
}


# Simulations of single RCTs ====

## We also want to recreate what would happen should we see lots of RCTs
## each with their own effect size, power and sample size. These represent
## the "single draws" of a frequency statistic that we get to observe
## when performing an RCT
for (i in seq_len(1000)) {
  cm <- rbeta(1, 2, 8)
  rrr <- rbeta(1, 2, 8)
  arr <- cm * rrr
  im <- cm - arr
  pwr <- sample(c(0.7, 0.8, 0.9), size = 1)
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

# Convert fixed parameters to factors
df <- df %>%
  mutate_at(
    vars(power, c_y, t_y, arr, rrr, n),
    factor
  )

trial <- trial %>%
  mutate_at(
    vars(power),
    factor
  )

save(df, trial, file = "./data/sim.RData")

# Figures
figure_1 <- df %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(y = fragility)) +
  geom_jitter(shape = 1, width = 0.2, height = 0, alpha = 0.5, aes(x = 0)) +
  geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
  facet_grid(rows = vars(power),
             cols = vars(arr),
             labeller = label_both,
             scales = "free_x") +
  inspectEHR::theme_cchic() +
  theme(text = element_text(family = "Helvetica"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("Fragility Index") +
  coord_flip()

ggsave(filename = "./figures/figure1.svg", plot = figure_1, height = 3, width = 9)
ggsave(filename = "./figures/figure1.png", plot = figure_1, height = 3, width = 9)

# df %>%
#   filter(reverse_effect == FALSE) %>%
#   ggplot(aes(x = p_value, y = fragility)) +
#   facet_grid(cols = vars(power), labeller = label_both) +
#   geom_point(
#     shape = 1, alpha = 0.25,
#     aes(colour = arr)) +
#   inspectEHR::theme_cchic() +
#   labs(colour = "Abolute Risk Reduction") +
#   ylab("Fragility Index") +
#   xlab("P Value") +
#   theme(text = element_text(family = "Helvetica"))

# Use 80% power for a clean exemplar
figure_2 <- df %>%
  filter(reverse_effect == FALSE,
         power == "0.8") %>%
  ggplot(aes(x = p_value, y = fragility)) +
  geom_point(shape = 1, alpha = 0.5, aes(colour = arr)) +
  inspectEHR::theme_cchic() +
  labs(colour = "Abolute Risk Reduction") +
  ylab("Fragility Index") +
  xlab("P Value") +
  theme(text = element_text(family = "Helvetica"))

ggsave(filename = "./figures/figure2.svg", plot = figure_2, height = 3, width = 9)
ggsave(filename = "./figures/figure2.png", plot = figure_2, height = 3, width = 9)

# Predicting the FI ====

# I want to build the case that the FI is a transformation of the p value
# So we should be able to predict the FI accurately from information in
# the trial directly, without having to do the fishers exact test procedure

# Niave analysis ----
niave_p <- df %>%
  filter(reverse_effect == FALSE) %>%
  lm(fragility ~ p_value, data = .)

niave_n <- df %>%
  filter(reverse_effect == FALSE) %>%
  lm(fragility ~ n, data = .)

summary(niave_p)
summary(niave_n)

# Obviously there is no relationship here. This is for 2 reasons:
# - the p value needs to be transformed
# - this is a multivariable problem (because multiple study designs map onto the
#   same p value) and so we need to include some of the study design elements

good_mod <- trial %>%
  filter(reverse_effect == FALSE) %>%
  mutate(logp = log(p_value)) %>%
  lm(fragility ~ n + n:arr + n:logp + arr:logp + n:arr:logp, data = .)

left_plot <- trial %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(x = fragility, y = n)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_smooth(method = "lm") +
  geom_point(alpha = 0.2) +
  ylab("Study Sample Size") +
  xlab("Calculated Fragility Index") +
  inspectEHR::theme_cchic()

right_plot <- trial %>%
  filter(reverse_effect == FALSE) %>%
  mutate(prediction = predict(good_mod)) %>%
  ggplot(aes(x = fragility, y = prediction)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.2) +
  ylab("Predicted Fragility Index") +
  xlab("Calculated Fragility Index") +
  inspectEHR::theme_cchic()

multiplot(plots = list(left_plot, right_plot), cols = 2, file = "./figures/figure3.svg")
multiplot(plots = list(left_plot, right_plot), cols = 2, file = "./figures/figure3.png")

## Extract Coefficients
equatiomatic::extract_eq(good_mod)