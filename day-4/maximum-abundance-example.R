# Analyse the simulated species abundance data

# load packages
pkgs <- c("mgcv", "ggplot2", "readr", "dplyr", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)

# load the data
spp_url <- "https://bit.ly/spp-gradient"
gradient <- read_csv(spp_url, col_types = "dd")
# or if you have the data cloned locally
# gradient <- read_csv(here("data", "simulated-gradient.csv"),
#                      col_types = "dd")
gradient


mod1 = gam(log(abundance+1) ~ s(environment,k=50),
           family = gaussian(link="identity"),
           data = gradient)
summary(mod1)
gratia::appraise(mod1)
gratia::draw(mod1)

library(ggeffects)

ggpredict(mod1, "environment") %>% plot()


mod1 = gam(log(abundance+1) ~ s(environment,k=50),
           family = gaussian(link="identity"),
           data = gradient)
summary(mod1)
gratia::appraise(mod1)
gratia::draw(mod1)

library(ggeffects)

ggpredict(mod1, "environment") %>% plot()

# Answers

max_nb = gam(abundance~ s(environment),
             data = gradient,
             method = "REML",
             family = nb())

max_p = update(max_nb, . ~ ., family = poisson())
ggpredict(max_p, "environment") %>% plot()
AIC(max_nb, max_p)

new_df = data.frame(environment = seq(0,100,.1))

fs = fitted_samples(max_p,
                    n = 10000,
                    newdata = new_df,
                    seed = 42)

fs2 = fitted_samples(max_p,
                    n = 100,
                    newdata = new_df,
                    seed = 42)

max_loc <- function(row,abund,env){
  r = row[which.max(abund)]
  env[r]
}

max_post = fs %>%
  group_by(draw) %>%
  summarize(env = max_loc(row,fitted,new_df$environment))


summ = max_post %>%
  summarize(mean = mean(env),
            median = median(env),
            q2_5 = quantile(env, prob = 0.025),
            q97_5 = quantile(env, prob = 0.975))


summ = summ %>%
  mutate(abundance = 0)

gradient %>%
  ggplot(aes(x=environment, y=abundance)) +
  geom_point() +
  #geom_line(data = fs2,
  #          aes(x=row,y=fitted,group=draw)) +
  geom_pointrange(data=summ,
                  aes(y=abundance, x = median,
                  xmax = q97_5, xmin=q2_5), col = "green")

fs2 %>%
  ggplot(aes(x=row,y=fitted,group=draw)) +
  geom_line()
