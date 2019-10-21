## This is the main script for the manuscript on the fragility index
## Authors:
## - Code - Edward Palmer, @DocEd (github)
## - Manuscript - Edward Palmer, Giampiero Marra, Mervyn Singer
## If you find any erros or would like to discuss the paper in more
## detail, please do get in touch.

library(tidyverse)
# Note: need inspectEHR for custom ggplot theme
# can install using:
# devtools::install_github("cchic/inspectEHR")
# Or just comment out `inspectEHR::theme_cchic()` when plotting

set.seed(2019)

#' Calculate Average Treatment Effect
#' 
#' returns the average treatment effect given treatment assignment in a
#' randomised scenario
#'
#' @param t treatment assignment (must be 0 or 1)
#' @param y outcome (must be 0 or 1)
ate_calc <- function(t, y) {
  mean(y[t == 1]) - mean(y[t == 0])
}

#' Simulate Fragility Index
#' 
#' This is the workhorse function of the script. It calculates the appropriate
#' number of subjects for a given effect size (specified as a control and
#' intervention mortality), power and alpha. Simulations are repeated a default
#' 1000 times (set by sims). Alternately, one can override this behaviour by
#' specifiying the number of subjects directly (set by N).
#' 
#' The `permit_negative` refers to whether or not you want the FI to return
#' values for when the starting p value is alreadyv > 0.05.
#' 
#' A data frame is returned with list columns containing the original treatment
#' allocations and outcomes. A number of other metrics are returned, some
#' observed, some defined apriori. A column specifying `reverse_effect`
#' identifies cases where the reverse effect was present (i.e. harm observed
#' in the treatment arm). These cases should either be excluded, or flipped for
#' them to be valid.
#'
#' @param intervention_mort treatment arm mortality (probability scalar)
#' @param power desired statistical power (probability scalar)
#' @param control_mort control arm mortality (probability scalar)
#' @param n manual selection of subjects in total (overrides other settings)
#' @param sims number of sims to run (default 1000)
#' @param alpha pre-fixed alpha level to run fragility index
#' @param permit_negative allow negative fragility index
simulate_fragility <- function(intervention_mort,
                               power = 0.8,
                               control_mort,
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
    oyt <- mean(y[t == 1])
    # Outcome WITHOUT treatment, note 1 = death, 0 = survive
    # Bernoilli trial with p = control_mortality
    y[-treatment] <- rbinom(n = n/2, size = 1, prob = control_mort)
    oyc <- mean(y[t == 0])
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
           oc_y = oyc,
           t_y = intervention_mort,
           ot_y = oyt,
           arr = control_mort-intervention_mort,
           rrr = (control_mort-intervention_mort)/control_mort,
           n = n,
           ate = purrr::map2_dbl(t, y, ate_calc))
  
  return(df)
}

# Core Scenarios ====
scenarios <- as_tibble(expand.grid(c(0.25, 0.2, 0.15), c(0.7, 0.8, 0.9)))
names(scenarios) <- c("int_mort", "power")

dt <- map2_dfr(
  .x = scenarios$int_mort,
  .y = scenarios$power,
  .f = ~ simulate_fragility(.x, .y, control_mort = 0.3, permit_negative = TRUE)
  )

# Simulations of single RCTs ====

## We also want to recreate what would happen should we see lots of RCTs
## each with their own effect size, power and sample size. These represent
## the "single draws" of a frequency statistic (i.e. p value) that we get to
## observe when performing an RCT

df <- tibble(
  cm = rbeta(1000, 2, 15),
  rrr = runif(1000, 0.1, 0.5),
  arr = cm * rrr,
  im = cm - arr,
  pwr = sample(c(0.7, 0.8, 0.9), size = 1000, replace = TRUE)
  )

dx <- pmap_dfr(
  .l = list(
    df$im,
    df$pwr,
    df$cm
  ),
  .f = ~ simulate_fragility(..., permit_negative = TRUE, sims = 1)
)

# Convert fixed parameters to factors
dt <- dt %>%
  mutate_at(
    vars(power, c_y, t_y, arr, rrr, n),
    factor
  )

dt <- rename(dt, `absolute risk reduction` = arr)

dx <- dx %>%
  mutate_at(
    vars(power),
    factor
  )

save(dt, dx, file = "./data/sim.RData")

# Figures
# Note, figure 1 is the detail of the FI procedure.
figure_2 <- dt %>%
  filter(reverse_effect == FALSE) %>%
  ggplot(aes(y = fragility)) +
  #geom_jitter(shape = 1, width = 0.2, height = 0, alpha = 0.5, aes(x = 0)) +
  geom_boxplot(outlier.alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(rows = vars(power),
             cols = vars(`absolute risk reduction`),
             labeller = label_both) +
  inspectEHR::theme_cchic() +
  theme(
    panel.background = element_rect(fill = "grey80", colour = "white"),
    text = element_text(family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Fragility Index") +
  coord_flip()

ggsave(filename = "./figures/figure2.svg", plot = figure_2, height = 3, width = 9)
ggsave(filename = "./figures/figure2.png", plot = figure_2, height = 3, width = 9)

# Use 80% power for a clean exemplar
figure_3 <- dt %>%
  filter(reverse_effect == FALSE,
         power == "0.8") %>%
  ggplot(aes(x = p_value, y = fragility)) +
  geom_point(shape = 1, alpha = 0.5, aes(colour = `absolute risk reduction`)) +
  inspectEHR::theme_cchic() +
  labs(colour = "Abolute Risk Reduction") +
  ylab("Fragility Index") +
  xlab("P Value") +
  theme(text = element_text(family = "Helvetica"))

ggsave(filename = "./figures/figure3.svg", plot = figure_3, height = 3, width = 9)
ggsave(filename = "./figures/figure3.png", plot = figure_3, height = 3, width = 9)

# Predicting the FI ====

# I want to build the case that the FI is a transformation of the p value
# So we should be able to predict the FI accurately from information in
# the trial directly, without having to do the fishers exact test procedure

# Niave analysis ----
niave_p <- dx %>%
  filter(reverse_effect == FALSE) %>%
  lm(fragility ~ p_value, data = .)

niave_n <- dx %>%
  filter(reverse_effect == FALSE) %>%
  lm(fragility ~ n, data = .)

summary(niave_p)
summary(niave_n)

# There is "no relationship" here. This is for 2 reasons:
# - the p value needs to be transformed
# - this is a multivariable problem (because multiple study designs map onto the
#   same p value) and so we need to include some of the study design elements

good_mod <- dx %>%
  filter(reverse_effect == FALSE) %>%
  filter(fragility > 0) %>%
  mutate(logp = log(p_value)) %>%
  lm(fragility ~ ate*logp*n*oc_y, data = .)

good_mod %>%
  broom::tidy() %>%
  mutate(l95 = confint(good_mod)[,1],
         u95 = confint(good_mod)[,2])

summary(good_mod)

## Not used in study, but its clear there is a perfect reclaiming of the FI here

dx %>%
  filter(
    reverse_effect == FALSE,
    fragility > 0) %>%
  mutate(prediction = predict(good_mod)) %>%
  ggplot(aes(x = fragility, y = prediction)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.2) +
  ylab("Predicted Fragility Index") +
  xlab("Calculated Fragility Index") +
  inspectEHR::theme_cchic()

## And we're done.
