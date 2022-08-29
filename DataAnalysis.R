#Analyzing variability in dissolved oxygen dynamics in surface and bottom waters 
# of lakes to accompany SIL Kilman Lecture paper
# Cayelan Carey, created 1 April 2022

#load packages
pacman::p_load(tidyverse, lubridate, zoo, goeveg, magrittr)

#download GLEON dissolved oxygen data from EDI repository
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/698/3/ca7001cba78a64c0ad0193580f478353" 
infile1 <- paste0(getwd(),"/Raw_Data.csv")
download.file(inUrl1,infile1,method="curl", extra = '-k')

#read in raw data file
raw_data <- read.csv("Raw_Data.csv", header=T) %>% 
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d"),
         year = lubridate::year(date),
         day = lubridate::day(date),
         month = lubridate::month(date))

no_obs <- 3 #minimum number of observations within a year to be included in analysis

####multi-annual surface####
# Let's first analyze variability in surface, multi-annual oxygen (1 m depth)
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

summary_annual <- surface %>% 
  filter(p_value<0.05) %>% #76 out of 226 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_annual$trend) #27/76 (36%) lakes have positive trend in multi-annual variance in median DO (variance is increasing); 
# 49/76 (64%) lakes have decreasing variance)

####multi-annual bottom####
#let's repeat this analysis with bottom, multi-annual oxygen (bottom depth)
bottom <- raw_data %>% 
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
#with consistently >3 observations/yr

length(which(bottom$slope>0)) #36/76 increasing variability
length(which(bottom$slope<0)) #40/76 decreasing variability

hist(bottom$slope)
hist(bottom$slope[which(surface$p_value<0.05)])

summary_annual <- bottom %>% 
  filter(p_value<0.05) %>% #19 out of 76 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_annual$trend) #10/19 lakes have positive trend in multi-annual variance in median DO (variance is increasing); 
# 9 lakes have decreasing variance)



####within-year surface####
#now let's do analyze surface, within-year oxygen (1 m depth)
#calculate within-year CV, then look at slope over time
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
  drop_na() #247 lakes total

hist(surface_within$slope)
hist(surface_within$slope[which(surface_within$p_value<0.05)])

summary_notannual <- surface_within %>% 
  filter(p_value<0.05) %>% #37 out of 247 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_notannual$trend) #18/37 lakes have positive trend in sub-annual variance in median DO (variance is increasing); 
# 19 lakes have decreasing variance)


####within-year bottom####
#let's repeat with bottom, within-annual oxygen (bottom depth)
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

summary_notannual <- bottom_within %>% 
  filter(p_value<0.05) %>% #14 out of 93 lakes have significant changes over time 
  mutate(trend = ifelse(slope>0, 1, 0)) #if slope is greater than 1, value of 1, otherwise zero

sum(summary_notannual$trend) #14/14 lakes have positive trend in sub-annual variance in median DO (variance is increasing); 
# 0 lakes have decreasing variance



####make histogram plots####
#generate 2 figures: 1 among year and 1 within-year, 2 density plots superimposed for surface/bottom

# turn data to long form
surface1 <- surface %>% 
  add_column(Depth = "Surface_AmongYears")
surface2 <- surface_within %>% 
  add_column(Depth = "Surface_WithinYears")

bottom1 <- bottom %>% 
  add_column(Depth = "Bottom_AmongYears")
bottom2 <- bottom_within %>% 
  add_column(Depth = "Bottom_WithinYears")

compare <- rbind(surface1, surface2)
compare1 <- rbind(bottom1, bottom2)

#Using all data, including lakes with no change in data
# Overlaying histograms of surface patterns
ggplot(compare, aes(x = slope, fill = Depth)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(surface1$slope), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(surface2$slope), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.012,0.012)) +
  xlab("Trend in surface dissolved oxygen variability (/year)") +
  ylab("Density")

# Overlaying histograms of surface patterns
ggplot(compare1, aes(x = slope, fill = Depth)) +
  geom_density(alpha = .5) +
   geom_vline(xintercept = median(bottom1$slope), color = "red", linetype = "dashed") +
    geom_vline(xintercept = median(bottom2$slope), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.1,0.1)) +
  xlab("Trend in bottom dissolved oxygen variability (/year)") +
  ylab("Density")


#now, let's filter out the lakes that show a significant trend over time
compare_a <- rbind(surface1, surface2) %>% 
  filter(p_value<0.05)
compare1_a <- rbind(bottom1, bottom2) %>% 
  filter(p_value<0.05)

#Using all data, including lakes with no change in data
# Overlaying histograms of surface patterns
ggplot(compare_a, aes(x = slope, fill = Depth)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(surface1$slope[which(surface1$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = median(surface2$slope[which(surface2$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.022,0.022)) +
  xlab("Trend in surface dissolved oxygen variability (/year)") +
  ylab("Density")

# Overlaying histograms of surface patterns
ggplot(compare1_a, aes(x = slope, fill = Depth)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = median(bottom1$slope[which(bottom1$p_value<0.05)]), color = "red", linetype = "dashed") +
  geom_vline(xintercept = median(bottom2$slope[which(bottom2$p_value<0.05)]), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_x_continuous(limits = c(-0.12,0.12)) +
  xlab("Trend in bottom dissolved oxygen variability (/year)") +
  ylab("Density")




# #just a check that bottom within-year data look ok
# bottom_raw <- raw_data %>%
#   group_by(name, meas_location, year) %>%
#   filter(depth==max(depth)) %>%
#   summarise(median_do_con = median(do_con),
#             n_withinyear = n(),
#             cv_do = cv(do_con),
#             depth = median(depth)) %>%
#   filter(n_withinyear>no_obs) %>%
#   #unique(surface$name) = 293
#   ungroup() %>%
#   group_by(name,meas_location) %>%
#   filter(n()>10) %>%
#   drop_na() #%>%
#   #summarise(slope = ifelse(n() > 1, lm(cv_do ~ year)$coefficients[2], NA),
#   #          p_value = ifelse(n() > 11, summary(lm(cv_do ~ year))$coefficients[2,4], NA))
# 
# par(mfrow = c(5,5))
# for(i in length(unique(bottom_raw$name))){
#   temp <-subset(bottom_raw, bottom_raw$name == unique(bottom_raw$name)[i])
#   plot(temp$year, temp$cv_do, main=paste0(temp$name[1]," p-value= ",summary(lm(temp$cv_do ~ temp$year))$coefficients[2,4]))
# }

#Take-home messages:
# multi-annual variation appears to be decreasing in surface (but not bottom)
# within-year variation appears to be increasing in bottom (but not surface)
