# Bone Mineral Densitometry
# Dan Wickins - using data from Barton and Peat
# 2026-06-11


# Using longitudinal data from Barton and Peat.

# Load packages ----
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse,
  janitor, 
  here,
  labelled,
  rstatix,
  GGally,
  car,
  lme4,
  lmerTest,
  emmeans,
  multcomp,
  geepack,
  ggeffects,
  gt,
  Epi,           
  tidy 
)


theme_set(theme_minimal() + theme(legend.position = "bottom")) 

# Read in the data
bmd_long <- read.table("data/BMD_long.txt", header = TRUE)
head(bmd_long)

# Now convert it to a tibble
bmd_long <- as_tibble(bmd_long)
head(bmd_long)
glimpse(bmd_long)

# Change 'Time' and 'Group' to factor variables.
bmd_long <- bmd_long |>
  mutate(Group = as_factor(Group), Time = as_factor(Time))

#? Where are the nas?
bmd_long |> 
  filter(if_any(everything(), is.na))

# Count the rows
bmd_long |> 
  filter(if_any(everything(), is.na)) |>
  nrow()

# Now drop the nas
bmd_long <- drop_na(bmd_long)

# Now check
bmd_long |> nrow()

glimpse(bmd_long)

# Now we will do some exploratory data analysis.

### MEAN BMD over time ####
# Compare BMD by group over time
group_by(bmd_long, Group) |>
  get_summary_stats(BMD)



# Now by Time
group_by(bmd_long, Time) |>
  get_summary_stats(BMD)

# Now by Time and Group
group_by(bmd_long, Time, Group) |>
  get_summary_stats(BMD)

# We can see that, although the intervention group starts with a lower mean and 
# IQR, by the time we get to time 3 and even 2, it has essentially surpassed the 
# Control group.

# Now let's compare visually
jittered_plot <- ggplot(bmd_long, aes(x = interaction(Group, Time), y = BMD, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  guides(fill = "none") +
  labs(x = "", y = "BMD") +
  theme_minimal()

jittered_plot

jittered_plot

ggsave(jittered_plot, filename = "outputs/jittered_plot.png", bg = "white")


# What if we remove the interaction term and just put Time, then leave the fill 
# for differentiating Group?
jittered_two <- ggplot(bmd_long, aes(x = Time, y = BMD, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  #guides(fill = "none") 
  labs(x = "", y = "BMD")

jittered_two

ggsave(jittered_two, filename = "outputs/jittered_two.png", bg = "white")

# Could also do this without the jitter
ggplot(bmd_long, aes(x = Time, y = BMD, fill = Group)) +
  geom_boxplot() +
  labs(x = "", y = "BMD")

# Now let's try a conventional clustered bar plot
group_by(bmd_long, Group, Time) |>
  summarise(mean_BMD = mean(BMD), .groups = "drop")|>
  ggplot(aes(Group, mean_BMD, fill = Time, label = mean_BMD))+
           geom_col(position = "dodge")+
           geom_text(position = position_dodge(width = 0.9), vjust= 0.5) +
           coord_flip() +
           labs(x = "", y = "Mean BMD", fill = "")



##### CORRELATIONS ######

# Will need to read the long file back wide.
bmd_wide <- bmd_long |>
  pivot_wider(
    names_from = Time,
    values_from = BMD
  )

head(bmd_wide)

# But now need to re-order columns
bmd_wide <- bmd_wide |>
  relocate(`2`, .before = `3`)

head(bmd_wide)
bmd_wide |>
  filter(if_any(everything(), is.na))

bmd_wide |>
  nrow()

bmd_wide <-drop_na(bmd_wide)
bmd_wide
glimpse(bmd_wide)

# Good now we have complete cases

# Let's obtain the Covariance matrix
cov_obs <- bmd_wide |> 
  dplyr::select(`1`,`2`,`3`)

# covariance matrix
cov_obs <- cov_obs |> cov()
cov_obs

# Now the correlation matrix
cov2cor(cov_obs)

# Now use ggpairs to illustrate positive (linear) correlation among
# consecutive measurements.

graphic <- bmd_wide |> 
  dplyr::select(`1`,`2`,`3`)

ggpairs(graphic, lower = list(continuous = "smooth"))

ggpairs(bmd_wide, mapping = aes(colour = Group), columns = 3:5,
        lower = list(continuous = "smooth"))

# The review page by Alessio Crippa discusses many further plots which
# can be employed, here we seek to move on to the research question.

##RESEARCH QUESTION --
# Does bone density differ between those receiving treatment and those 
# not receiving treatment? If so, does the density differ between the two 
# groups over time? i.e is there an interaction between 'Group' and 'Time'?

# In essence, we seek to understand:
# 1. Time effect
# 2. Group effect.
# 3. Interaction between time and group.


# Start with the empty model. 
# In this case we model the fixed effect, comprising Beta0, which is the
# intercept of the whole model, u0i which is the random effect of each 
# individual subject ('id'), and the error term. The random effect is 
# independent of the error term.

# Use the lmer function from within the lme4 package, referencing Crippa 2022.
lin_0 <- lmer(BMD ~ 1 + (1 | id), data = bmd_long)
summary(lin_0)


# Now the marginal mean with confidence intervals.
ci.lin(lin_0)

# The estimated marginal mean of bmd is 0.8559264.
# This has an estimated between-subject variability of 0.00016330.
# The estimated variance of the error term is 0.0008778. 
# Thus, the correlation between any two repeated measures (ICC) is equal
# to 0.0016330/(0.0016330 + 0.0008778) = 0.6503903138.

# We check variance components using 'ranova' - which is within the lmer test.
ranova(lin_0)
# See that the Chi-square statistic from the LRT on 1 DF is 70.3, which far 
# exceeds significance. This supports the alternative hypothesis of between-
# subject heterogeneity - i.e. subject ID is an important explanatory 
# variable.

# Graphical representation of the model--
ggplot(bmd_long, aes(id, BMD)) +
  geom_point(aes(col = Time, shape = Time)) +
  geom_point(data = group_by(bmd_long, id) |>
               summarise(BMD = mean(BMD), .groups = "drop"),
             aes(col = "Mean", shape = "Mean"), size = 2.5) +
  geom_hline(yintercept = mean(bmd_long$BMD)) +
  labs(x = "Subject id", y = "BMD", col = "Time", shape = "Time")

# Actually a fairly hard to read plot.


# Time effect: Is the mean BMD varying with time?
lin_1 <- lmer(BMD ~ Time + (1 | id), data = bmd_long)
summary(lin_1)

# The mean estimate of BMD at measure 1 is 0.839. The mean difference between 
# BMD measure at TIme 2 comparted to time 1 is 0.003005, then between Time 3 
# and Time 1 is 0.02341

# The 'emmeans' function can here be used for computing marginal means over
# measurement times with corresponding confidence intervals.

tidy(emmeans(lin_1, "Time"), conf_int = TRUE)