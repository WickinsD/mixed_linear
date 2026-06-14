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
  gt
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
ggplot(bmd_long, aes(x = interaction(Group, Time), y = BMD, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  guides(fill = "none") +
  labs(x = "", y = "BMD")


# What if we remove the interaction term and just put Time, then leave the fill 
# for differentiating Group?
ggplot(bmd_long, aes(x = Time, y = BMD, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  #guides(fill = "none") 
  labs(x = "", y = "BMD")

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
