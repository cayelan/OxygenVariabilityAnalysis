#Analyzing variability in dissolved oxygen dynamics in surface and bottom waters 
# of lakes to accompany SIL Kilman Lecture paper
# Cayelan Carey, created 1 April 2022
# last modified 11 September 2023

#load packages
pacman::p_load(tidyverse, lubridate, zoo, goeveg, magrittr, patchwork, trend)

#download GLEON dissolved oxygen data from EDI repository (Stetler et al. 2021)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/698/3/c0aab240966ee30c6a4248c437b53b0d"
infile1 <- paste0(getwd(),"/Raw_Data.csv")
download.file(inUrl1,infile1,method="curl", extra = '-k')

#read in raw data file
raw_data <- read.csv("Raw_Data.csv", header=T) %>% 
  dplyr::mutate(date = as.POSIXct(date, format = "%Y-%m-%d"),
         year = lubridate::year(date),
         day = lubridate::day(date),
         month = lubridate::month(date))

no_obs <- 3 #minimum number of observations within a year to be included in analysis


####inter-annual surface####
# Also referred to as the multi-annual or among-year surface DO variability analysis
# Let's first analyze variability in surface, inter-annual oxygen (1 m depth)
# following methods of Cusser et al. 2021 Ecol Letters which used rolling windows of 3 yrs
# and found that temporal ecological trends took on average 9.66 years to reach a
# critical temporal threshold and achieve consistent results
# We are calculating the CV in 3 year rolling windows over the time series
surface <- raw_data %>% 
  filter(depth==1) %>% 
  group_by(lake_id, year) %>%  
  summarise(median_do_con = median(do_con),
            n_withinyear = n()) %>% 
  filter(n_withinyear > no_obs) %>% 
  ungroup() %>% 
  group_by(lake_id) %>% 
  mutate(roll = rollapply(median_do_con, width=3, FUN=cv, align="right", fill=NA, na.rm=T)) %>% 
  drop_na() %>% 
  summarise(slope = ifelse(n() >= 10, lm(roll ~ year)$coefficients[2], NA),
            p_value = ifelse(n() >= 10, summary(lm(roll ~ year))$coefficients[2,4], NA)) %>% 
  drop_na() 
#226 lake ids total had at least 10 years of oxygen data at 1m depth 
#with consistently >=3 observations/yr
#178 lake ids total if subset to having >=5 observations/yr

length(which(surface$slope>0)) #94/226 increasing variability for n>=3 obs
length(which(surface$slope<0)) #132/226 decreasing variability for n>=3 obs

hist(surface$slope)
hist(surface$slope[which(surface$p_value<0.05)])

summary_annual_surface_among <- surface %>% 
  filter(p_value<0.05) %>% #76 out of 226 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_annual_surface_among$trend) 
# for n>=3 obs, 27/76 (36%) lakes have positive trend in multi-annual variance in median DO (variance is increasing); 
# 49/76 (64%) lakes have decreasing variance)
# for >=5 obs, 16/55 (29%) have positive trend in variance; 71% have decreasing variance


####inter-annual bottom####
#let's repeat this analysis with bottom, inter-annual oxygen (bottom depth)
#among-year bottom analysis
bottom_among <- raw_data %>% 
  group_by(lake_id, year) %>% 
  filter(depth==max(depth)) %>% 
  summarise(median_do_con = median(do_con),
            n_withinyear = n(),
            depth = median(depth)) %>% 
  filter(n_withinyear> no_obs) %>% 
  ungroup() %>% 
  group_by(lake_id) %>% 
  mutate(roll = rollapply(median_do_con, width=3, FUN=cv, align="right", fill=NA, na.rm=T)) %>% 
  drop_na() %>% 
  summarise(slope = ifelse(n() >= 10, lm(roll ~ year)$coefficients[2], NA),
            p_value = ifelse(n() >= 10, summary(lm(roll ~ year))$coefficients[2,4], NA)) %>% 
  drop_na() 
#76 lakes total that had at least 10 years of oxygen data at the deepest lake sampling depth 
#with consistently >=3 observations/yr
#for >=5 obs/yr, 41 lakes total with data that spanned 10 yrs

length(which(bottom_among$slope>0)) #36/76 increasing variability for n>=3 obs
length(which(bottom_among$slope<0)) #40/76 decreasing variability for n>=3 obs

hist(bottom_among$slope)
hist(bottom_among$slope[which(bottom_among$p_value<0.05)])

summary_annual_bottom_among <- bottom_among %>% 
  filter(p_value<0.05) %>% #19 out of 76 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_annual_bottom_among$trend) 
# for n>=3 obs, 10/19 lakes have positive trend in multi-annual variance in median DO (variance is increasing); 
# 9 lakes have decreasing variance)
# for n>=5 obs, 7/15 have positive trend; 8 have decreasing variance


####intra-annual surface####
#now let's do analyze surface, within-year oxygen (1 m depth)
#calculate intra-annual (within-year) CV, then look at slope over time
surface_within <- raw_data %>% 
  group_by(lake_id, year) %>% 
  filter(depth==1) %>% 
  summarise(median_do_con = median(do_con),
            n_withinyear = n(),
            cv_do = cv(do_con)) %>% 
  filter(n_withinyear> no_obs) %>% 
  ungroup() %>% 
  group_by(lake_id) %>% 
  drop_na() %>% 
  summarise(slope = ifelse(n() >= 10, lm(cv_do ~ year)$coefficients[2], NA),
            p_value = ifelse(n() >= 10, summary(lm(cv_do ~ year))$coefficients[2,4], NA)) %>% 
  drop_na() #247 lakes total for n>=3 obs, 200 lakes for n>=5 obs

hist(surface_within$slope)
hist(surface_within$slope[which(surface_within$p_value<0.05)])

summary_annual_surface_within <- surface_within %>% 
  filter(p_value<0.05) %>% #for N>=3 obs, 37 out of 247 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero;
#for n>=5 obs, 32 lakes have significant changes over time

sum(summary_annual_surface_within$trend) #for N>=3 obs, 18/37 lakes have positive trend in sub-annual 
# variance in median DO (variance is increasing); 19 lakes have decreasing variance)
# for N>=5 obs, 14/32 lakes have increasing variance over time


####intra-annual bottom####
#let's repeat with bottom, intra-annual oxygen (bottom depth)
#referred to as within-year analysis
bottom_within <- raw_data %>% 
  group_by(lake_id, year) %>% 
  filter(depth==max(depth)) %>% 
  summarise(median_do_con = median(do_con),
            n_withinyear = n(),
            cv_do = cv(do_con),
            depth = median(depth)) %>% 
  filter(n_withinyear>no_obs) %>% 
  ungroup() %>% 
  group_by(lake_id) %>% 
  drop_na() %>% 
  summarise(slope = ifelse(n() >= 10, lm(cv_do ~ year)$coefficients[2], NA),
            p_value = ifelse(n() >= 10, summary(lm(cv_do ~ year))$coefficients[2,4], NA)) %>% 
  drop_na() #93 lakes total

hist(bottom_within$slope)
hist(bottom_within$slope[which(bottom_within$p_value<0.05)])

summary_annual_bottom_within <- bottom_within %>% 
  filter(p_value<0.05) %>% #For N>=3 obs, 14 out of 93 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero
#For N>=5 obs, 8 lakes have significant changes over time

sum(summary_annual_bottom_within$trend) #14/14 lakes have positive trend in sub-annual variance 
# in median DO (variance is increasing); 0 lakes have decreasing variance
#For N>=5 obs, 8/8 lakes have increasing variance over time


############################################
####make histogram figures####
#For Figure 1, generate 2 panels: 1 inter-annual and 1 intra-annual, 
# with 2 density plots superimposed for surface/bottom layers
#Figure 1 includes n=3 lakes (line 20 above) which exhibit only significant trends over time

# turn data to long form
surface1 <- surface %>% 
  add_column(Depth = "Surface_Interannual")
surface2 <- surface_within %>% 
  add_column(Depth = "Surface_Intraannual")

bottom1 <- bottom_among %>% 
  add_column(Depth = "Bottom_Interannual")
bottom2 <- bottom_within %>% 
  add_column(Depth = "Bottom_Intraannual")

compare <- rbind(surface1, surface2)
compare1 <- rbind(bottom1, bottom2)

#now, let's filter out the lakes to focus on those that show a significant trend over time
compare_a <- rbind(surface1, surface2) %>% 
  filter(p_value<0.05) %>% 
  mutate(depth = "surface",
         split = ifelse(Depth == "Surface_Interannual", "Interannual", "Intraannual"))
compare1_a <- rbind(bottom1, bottom2) %>% 
  filter(p_value<0.05) %>% 
  mutate(depth = "bottom",
         split = ifelse(Depth == "Bottom_Interannual", "Interannual", "Intraannual"))

combined <- bind_rows(compare_a, compare1_a) %>% 
  select(-Depth)

# Overlaying histograms of surface patterns
p1 <- combined %>% filter(depth == "surface") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(surface1$slope[which(surface1$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(surface2$slope[which(surface2$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.022,0.022)) +
  labs(x = "", y = "Density", title = "a) Surface") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.18, 0.83),
        legend.title =element_blank())
  
xlabel = expression("Trend in dissolved oxygen variablity (yr^-1)")

# Overlaying histograms of surface patterns
p2 <- combined %>% filter(depth == "bottom") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(bottom1$slope[which(bottom1$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = median(bottom2$slope[which(bottom2$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.12,0.12)) +
  labs(x = bquote('Trend in dissolved oxygen variability ('*yr^-1*')'), y = "Density", title = "b) Bottom") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

library(patchwork)

p = p1 / p2

ggsave("Fig1.jpg", plot = p, device = "jpg", dpi = 600, width = 4, height = 6)
#ggsave("Fig1.pdf", plot = p, device = "pdf", width = 4, height = 6)


#####Make SI figures######
#now, let's rerun the code above - but instead of using n=3 observations, we will 
# use n=5 observations to check the robustness of our analysis to the number of lakes
# included. To do this, we will rerun all code from lines 20 to 160 above, but 
# on line 20 we will instead using n=5 observations for the analysis to make Fig S1
compare_a <- rbind(surface1, surface2) %>% 
  filter(p_value<0.05) %>% 
  mutate(depth = "surface",
         split = ifelse(Depth == "Surface_Interannual", "Interannual", "Intraannual"))
compare1_a <- rbind(bottom1, bottom2) %>% 
  filter(p_value<0.05) %>% 
  mutate(depth = "bottom",
         split = ifelse(Depth == "Bottom_Interannual", "Interannual", "Intraannual"))

combined <- bind_rows(compare_a, compare1_a) %>% 
  select(-Depth)

# Overlaying histograms of surface patterns
p1 <- combined %>% filter(depth == "surface") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(surface1$slope[which(surface1$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(surface2$slope[which(surface2$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.022,0.022)) +
  labs(x = "", y = "Density", title = "a) Surface") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.18, 0.83),
        legend.title =element_blank())

xlabel = expression("Trend in dissolved oxygen variability (yr^-1)")

# Overlaying histograms of surface patterns
p2 <- combined %>% filter(depth == "bottom") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(bottom1$slope[which(bottom1$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = median(bottom2$slope[which(bottom2$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.12,0.12)) +
  labs(x = bquote('Trend in dissolved oxygen variability ('*yr^-1*')'), y = "Density", title = "b) Bottom") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

library(patchwork)

p = p1 / p2

ggsave("SI_Fig1.jpg", plot = p, device = "jpg", dpi = 600, width = 4, height = 6)
#ggsave("SI_Fig1.pdf", plot = p, device = "pdf", width = 4, height = 6)
#this is a repeat of Fig 1 except using no_obs=5 on line 20


################################################################
#now, to make Fig S2, we'll use n=3 lakes again, but this time, include all lakes
# not just those who had significant trends in oxygen. To do this, we will rerun 
# all code from lines 20 to 160 above, using n=3 observations for the analysis
compare_a <- rbind(surface1, surface2) %>% 
  #filter(p_value<0.05) %>% 
  mutate(depth = "surface",
         split = ifelse(Depth == "Surface_Interannual", "Interannual", "Intraannual"))
compare1_a <- rbind(bottom1, bottom2) %>% 
  #filter(p_value<0.05) %>% 
  mutate(depth = "bottom",
         split = ifelse(Depth == "Bottom_Interannual", "Interannual", "Intraannual"))
combined <- bind_rows(compare_a, compare1_a) %>% 
  select(-Depth)

# Overlaying histograms of surface patterns
p1 <- combined %>% filter(depth == "surface") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(surface1$slope), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(surface2$slope), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.01,0.01)) +
  labs(x = "", y = "Density", title = "a) Surface") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.18, 0.83),
        legend.title =element_blank())

xlabel = expression("Trend in dissolved oxygen variablity (year^-1)")

# Overlaying histograms of surface patterns
p2 <- combined %>% filter(depth == "bottom") %>% 
  ggplot(aes(x = slope, fill = split)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(bottom1$slope), color = "red", linetype = "dashed") +
  geom_vline(xintercept = median(bottom2$slope), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.08,0.08)) +
  labs(x = bquote('Trend in dissolved oxygen variability ('*year^-1*')'), y = "Density", title = "b) Bottom") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

library(patchwork)

p = p1 / p2

ggsave("SI_Fig2.pdf", plot = p, device = "pdf", width = 4, height = 6)



#####comparing lakes with diverging surface & bottom patterns for text references in the manuscript
summary_annual_surface_among #lakes with increasing among-year surface water variability
summary_annual_surface_within #lakes with increasing within-year surface water variability
summary_annual_bottom_among #lakes with increasing among-year bottom water variability
summary_annual_bottom_within #lakes with increasing within-year bottom water variability


lakesDecliningSurfaceVar <- c(summary_annual_surface_among$lake_id[which(summary_annual_surface_among$trend<1)],
                              summary_annual_surface_within$lake_id[which(summary_annual_surface_within$trend<1)]) 
#create a vector of lakes with declining within OR among year variability in surface waters
lakesDecliningSurfaceVar <- sort(unique(lakesDecliningSurfaceVar)) #removed duplicated elements


lakesIncreasingSurfaceVar <- c(summary_annual_surface_among$lake_id[which(summary_annual_surface_among$trend==1)],
                              summary_annual_surface_within$lake_id[which(summary_annual_surface_within$trend==1)]) 
#create a vector of lakes with increasing within OR among year variability in surface waters
lakesIncreasingSurfaceVar <- sort(unique(lakesIncreasingSurfaceVar)) #removed duplicated elements


lakesIncreasingBottomVar <- c(summary_annual_bottom_among$lake_id[which(summary_annual_bottom_among$trend==1)],
                              summary_annual_bottom_within$lake_id[which(summary_annual_bottom_within$trend==1)]) 
lakesIncreasingBottomVar <- sort(unique(lakesIncreasingBottomVar))
#create a vector of lakes with increasing within OR among year variability in bottom waters

lakesDecliningBottomVar <- c(summary_annual_bottom_among$lake_id[which(summary_annual_bottom_among$trend<1)],
                              summary_annual_bottom_within$lake_id[which(summary_annual_bottom_within$trend<1)]) 
lakesDecliningBottomVar <- sort(unique(lakesDecliningBottomVar))
#create a vector of lakes with increasing within OR among year variability in bottom waters


intersect(lakesIncreasingBottomVar,lakesDecliningSurfaceVar)
#6 lakes exhibited both increasing bottom variability AND declining surface variability
#(144 158 169 216 221 415)

intersect(lakesIncreasingBottomVar,lakesIncreasingSurfaceVar)
intersect(lakesDecliningBottomVar,lakesDecliningSurfaceVar)
intersect(lakesDecliningBottomVar,lakesIncreasingSurfaceVar)

#for analyzing lakes with significant change in one layer but not the other
surface
surface_within
bottom_among
bottom_within

#bottom lakes with no significant change
lakesNoBottomChange <- c(bottom_among$lake_id[which(bottom_among$p_value>=0.05)],
                         bottom_within$lake_id[which(bottom_within$p_value>=0.05)])
lakesNoBottomChange <- sort(unique(lakesNoBottomChange))

#surface lakes with no significant change
lakesNoSurfaceChange <- c(surface$lake_id[which(surface$p_value>=0.05)],
                         surface_within$lake_id[which(surface_within$p_value>=0.05)])
lakesNoSurfaceChange <- sort(unique(lakesNoSurfaceChange))

#one layer has significant change, the other doesn't
OneLayerChange <- c(intersect(lakesNoBottomChange,lakesIncreasingSurfaceVar), 
                    intersect(lakesNoBottomChange,lakesDecliningSurfaceVar),
                    intersect(lakesNoSurfaceChange,lakesIncreasingBottomVar),
                    intersect(lakesNoSurfaceChange,lakesDecliningBottomVar))
OneLayerChange <- sort(unique(OneLayerChange))
length(OneLayerChange) #41 lakes exhibited significant change in one layer but not the other layer

#Bravo! You get a cookie! 
