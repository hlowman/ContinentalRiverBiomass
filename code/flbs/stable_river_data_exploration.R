# Data Exploration for "stable" river data
# October 19, 2021
# Heili Lowman

# Load packages.
library(tidyverse)
library(here)
library(patchwork)

# Load data.
dat <- read_csv("data_raw/flbs/selected_autotrophic_rivers_daily.csv")

# Plot of all data.
f1 <- ggplot(dat, aes(x = Date, y = GPP)) +
  geom_point() +
  facet_wrap(year~sitecode, scales = "free")

f1

# Compare with discharge, light, and temp.
dat_long <- dat %>%
  select(sitecode, Date, year, GPP, discharge_m3s, 
         light_PAR, temp.water) %>%
  mutate(light_1000 = light_PAR/1000) %>%
  pivot_longer(cols = c(GPP, discharge_m3s, 
                        light_1000, temp.water), names_to = "measure",
               values_to = "value")

# Plot of all data colored by external force.
f2 <- ggplot(dat_long, aes(x = Date, y = value, color = measure)) +
  geom_line() +
  facet_wrap(year~sitecode, scales = "free")

f2

# Canyon Creek specific
# 20 km upstream of dam
f3 <- dat_long %>%
  filter(sitecode == "nwis_10133650" & year == 2011) %>%
  ggplot(aes(x = Date, y = value, color = measure)) +
  geom_line() +
  labs(title = "E Canyon Creek") +
  facet_wrap(.~year, scales = "free")

f3

# Canyon Creek specific
# just upstream of dam
f4 <- dat_long %>%
  filter(sitecode == "nwis_10133980") %>%
  ggplot(aes(x = Date, y = value, color = measure)) +
  geom_line() +
  labs(title = "E Canyon Creek - above dam") +
  facet_wrap(.~year, scales = "free")

f4

# create paneled figure
p1.1 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012) %>%
  ggplot(aes(x = Date, y = GPP)) +
  geom_line(color = "#3A5D3D") +
  labs(title = "E Canyon Creek above reservoir (2012)") +
  theme_bw()

p1.2 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012) %>%
  ggplot(aes(x = Date, y = discharge_m3s)) +
  geom_line(color = "#3B4F8E") +
  theme_bw()

p1.3 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012) %>%
  ggplot(aes(x = Date, y = temp.water)) +
  geom_line(color = "#D3105C") +
  theme_bw()

p1.4 <- dat %>%
  filter(sitecode == "nwis_10133650" & year == 2012) %>%
  ggplot(aes(x = Date, y = light_PAR)) +
  geom_line(color = "#E69512") +
  theme_bw()

(site1 <- p1.1 /
  p1.2 /
  p1.4 /
  p1.3)

p2.1 <- dat %>%
  filter(sitecode == "nwis_10133980" & year == 2012) %>%
  ggplot(aes(x = Date, y = GPP)) +
  geom_line(color = "#3A5D3D") +
  labs(title = "E Canyon Creek 20km upstream (2012)") +
  theme_bw()

p2.2 <- dat %>%
  filter(sitecode == "nwis_10133980" & year == 2012) %>%
  ggplot(aes(x = Date, y = discharge_m3s)) +
  geom_line(color = "#3B4F8E") +
  theme_bw()

p2.3 <- dat %>%
  filter(sitecode == "nwis_10133980" & year == 2012) %>%
  ggplot(aes(x = Date, y = temp.water)) +
  geom_line(color = "#D3105C") +
  theme_bw()

p2.4 <- dat %>%
  filter(sitecode == "nwis_10133980" & year == 2012) %>%
  ggplot(aes(x = Date, y = light_PAR)) +
  geom_line(color = "#E69512") +
  theme_bw()

(site2 <- p2.1 /
  p2.2 /
  p2.4 /
  p2.3)

# Combine to create full figure.
(both <- (p1.1 | p2.1) /
    (p1.2 | p2.2) /
    (p1.4 | p2.4) /
    (p1.3 | p2.3))

# ggsave(plot = both,
#        filename = "figures/flbs/canyon_creek_2012.jpg",
#        width = 15,
#        height = 10)

# End of script.
