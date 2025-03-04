# Analyse the Rat hormone data from Fahrmeir et al (2013) Regression: Models,
#  Methods and Applications. Springer
pkgs <- c("mgcv", "lme4", "ggplot2", "readr", "dplyr", "forcats", "tidyr",
          "gratia")

# load packages
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load data
rats_url <- "https://bit.ly/rat-hormone"
rats <- read_table(rats_url, col_types = "dddddddddddd-")
# ignore the warning - it"s due to trailing white space at the ends of each
#   row in the file

rats <- rats %>%
    mutate(treatment = fct_recode(factor(group, levels = c(2,1,3)),
                                  Low = "1",
                                  High = "3",
                                  Control = "2"),
           subject = factor(subject))

rats %>%
    na.omit() %>%
    count(subject) %>%
    count(n, name = "n_rats")

plt_labs <- labs(y = "Head height (distance in pixels)",
                 x = "Age in days",
                 colour = "Treatment")

ggplot(rats, aes(x = time, y = response,
                 group = subject, colour = treatment)) +
    geom_line() +
    facet_wrap(~ treatment, ncol = 3) +
    plt_labs

m1_lmer <- lmer(response ~ treatment:transf_time +
                    (1 | subject) + (0 + transf_time | subject),
                data = rats)

m1_gam <- gam(response ~ treatment:transf_time +
                  s(subject, bs = "re") +
                  s(subject, transf_time, bs = "re"),
              data = rats, method = "REML")

fixef(m1_lmer)

coef(m1_gam)[1:4]

summary(m1_lmer)$varcor

variance_comp(m1_gam)

summary(m1_lmer)

summary(m1_gam)

edf(m1_gam)

new_data <- tidyr::expand(rats, nesting(subject, treatment),
                          transf_time = unique(transf_time))

m1_pred <- fitted_values(m1_gam, data = new_data, scale = "response")

ggplot(m1_pred, aes(x = transf_time, y = fitted, group = subject,
                    colour = treatment)) +
    geom_line() +
    facet_wrap(~ treatment) +
    plt_labs

m2_lmer <- lmer(response ~ treatment:transf_time +
                    (1 | subject),
                data = rats)
library(performance)
library(see)
check_model(m2_lmer)
m3_lmer <- lmer(response ~ treatment*transf_time +
                    (1 | subject),
                data = rats)
emmeans::emtrends(m3_lmer, pairwise~treatment, var = "transf_time")
library(ggeffects)
ggemmeans(m3_lmer, c("transf_time",
                     "treatment")) %>% plot()
m2_gam <- gam(response ~ treatment:transf_time +
                  s(subject, bs = "re"),
              data = rats, method = "REML")
m3_gam <- gam(
    response ~  treatment + s(transf_time,
                              by = treatment,
                              bs = "tp",
                              k = 7) +
        s(subject, bs = "re"),
    select = TRUE,
    data = rats,
    method = "REML"
)

m3_gam0 <- gam( response ~  treatment + s(transf_time,
                              bs = "tp",
                              k = 7) +
        s(subject, bs = "re"),
    select = TRUE,
    data = rats,
    method = "REML"
)

ggemmeans(m3_gam, c("transf_time",
                     "treatment")) %>% plot()

library(emmeans)
em1 = emmeans::emmeans(m3_gam,
                 ~ treatment)
em2 =ref_grid(m3_gam, at = list(transf_time = c(seq(.4,2,.1))),
              cov.keep = TRUE)
confint(emmeans(em2, pairwise~treatment | transf_time))$contrast %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    ggplot(aes(x=transf_time,
               y=estimate,
               ymin = lower_cl,
               ymax = upper_cl,
               color = contrast)) +
    geom_line(size = 1.2) +
    geom_ribbon(alpha = .1) +
    scale_color_viridis_d()

summary(m2_lmer)$varcor

variance_comp(m2_gam)

AIC(m2_gam, m1_gam)

draw(m2_gam, parametric = FALSE)
