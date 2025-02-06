rm(list = ls())

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(emmeans)
library(ggpmisc)

# Set working directory
setwd("C:/Users/Madison Kuczek/OneDrive - Louisiana State University/Mussel Project/Data/Excreta") 

# Read in Mussel excretion data
df <- fread("Excretion_in_ppb.csv")%>%
  group_by(Species) %>%
  filter(n() >= 5) %>%
  ungroup()

# Calculate site-specific mean NH₄⁺ and SRP for controls
site_control_means <- df %>%
  filter(Species == "Control") %>%
  filter(NH.ug.L < 190, SRP.ug.L < 90) %>%
  group_by(Site) %>%
  summarise(
    Mean_NH_control = mean(NH.ug.L, na.rm = TRUE),
    Mean_SRP_control = mean(SRP.ug.L, na.rm = TRUE)
  )

# Merge site-specific control values back into the main dataset
df <- left_join(df, site_control_means, by = "Site")

# Correct for volume
df <- df %>%
  mutate(
    NH.ug = (NH.ug.L - Mean_NH_control) * volume..mL. / 1000,
    SRP.ug = (SRP.ug.L - Mean_SRP_control) * volume..mL. / 1000
  )

# Correct for time and convert to umol/hr (per capita rate)
df <- df %>%
  mutate(
    NH.umol.h = (NH.ug / (incubation.time..mins. / 60)) / 14,   # 14=molar mass of nitrogen in ammonia(ug to umol)
    SRP.umol.h = (SRP.ug / (incubation.time..mins. / 60)) / 31  # 31=molar mass of nitrogen in ammonia(ug to umol)
  )

# Calculate molar N:P ratio and remove negatives
df <- df %>%
  filter(NH.umol.h >= 0, SRP.umol.h >= 0) %>%  # Remove rows with negative NH4 or SRP
  mutate(NP = NH.umol.h / SRP.umol.h,          # Calculate molar N:P ratio only for valid data
         NH.umol.g.h = NH.umol.h / STDM.g,
         SRP.umol.g.h = SRP.umol.h / STDM.g)  
   # NH.umol.h = ifelse(NH.umol.h < 0, 0, NH.umol.h),       # Replace negative NH4 values with 0
   # SRP.umol.h = ifelse(SRP.umol.h < 0, 0, SRP.umol.h),    # Replace negative SRP values with 0
   # NP = ifelse(NP < 0, 0, NP)                             # Replace negative N:P ratios with 0
  
df <- df[-77, ]
# Filter out control samples for final dataset
#excretion <- filter(df, Species != "Control")

# Save the per capita excretion rates
write.csv(df, file = "excretion_controls_adjusted.csv", row.names = FALSE)

############################# Visualization ##################################
##Plot ordered by tribe and escalating for volumetric excretion

ggplot(data = df) +
  geom_boxplot(aes(x =fct_reorder(Species, SRP.ug.L, .fun=mean, .desc = F), 
                   y =  NH.ug.L, fill = Tribe), outlier.shape = NA) +
  geom_point(aes(x = fct_reorder(Species, SRP.ug.L, .fun = mean, .desc = F), 
                 y =  NH.ug.L), position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic() +
  labs(x = "", y = "SRP(ug/L)")

ggplot(data = df) +
  geom_boxplot(aes(x =fct_reorder(Species, NH.ug.L, .fun=mean, .desc = F), 
                   y =  NH.ug.L, fill = Tribe), outlier.shape = NA) +
  geom_point(aes(x = fct_reorder(Species, NH.ug.L, .fun = mean, .desc = F), 
                 y =  NH.ug.L), position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic() +
  labs(x = "", y = " NH (ug/L)")


### PLOT OF per capita excreta rates adjusted for mass
ggplot(data = df, aes(x =fct_reorder(Species, SRP.umol.g.h, .fun= median), y = SRP.umol.g.h)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "red", color = "red") +  # Mean points in red
  scale_x_discrete(labels = c("Ecra" = "E.crassidens",
                               "Fcer" = "F.cerina",
                               "Oref" = "O.reflexia",
                               "Ouni" = "O.unicolor",
                               "Pbea" = "P.beadlianum",
                               "Pgra" = "P.grandis",
                               "Ppur" = "P.purpuratus",
                               "Ppus" = "P.pustulosa",
                               "Qapi" = "Q.quadrula",
                               "Qver" = "Q.verrucosa",
                               "Uhar" = "U.hartfieldorum")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +  # Rotate x-axis labels
  labs(x = "", y = "SRP (umol/g/h)")

ggplot(data = df, aes(x =fct_reorder(Species, NH.umol.g.h, .fun= median), y = NH.umol.g.h)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "red", color = "red") +  # Mean points in red
  scale_x_discrete(labels = c("Ecra" = "E.crassidens",
                              "Fcer" = "F.cerina",
                              "Oref" = "O.reflexia",
                              "Ouni" = "O.unicolor",
                              "Pbea" = "P.beadlianum",
                              "Pgra" = "P.grandis",
                              "Ppur" = "P.purpuratus",
                              "Ppus" = "P.pustulosa",
                              "Qapi" = "Q.quadrula",
                              "Qver" = "Q.verrucosa",
                              "Uhar" = "U.hartfieldorum")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +  # Rotate x-axis labels
  labs(x = "", y = "NH (umol/g/h)")

##Plot ordered by tribe and escalating for capita rate
ggplot(data = df, aes(x =fct_reorder(Species, SRP.umol.h, .fun=mean, .desc = F), 
                    y = SRP.umol.h)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic()+
  labs(x = "", y = "SRP (umol/h)")

ggplot(data = df, aes(x =fct_reorder(Species, NH.umol.h, .fun=mean, .desc = F), 
                      y = NH.umol.h)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic()+
  labs(x = "", y = "NH (umol/h)")

ggplot(data = df, aes(x = log10(STDM.g), y = log10(SRP.umol.h))) +
  geom_point(alpha = 0.7, size = 3) +  # Adjust transparency and point size
  geom_smooth(method = "lm", se = F, color = "black", linetype = "solid") +
  stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 5) +  # Add R² and p-value
  theme_classic() +
  labs(x = "STDM (g)", y = "SRP (umol/h)", colour = "Species") +
  scale_color_manual(values = c("Ecra" = "blue", "Fcer" = "red", "Oref" = "green",
                                "Ouni" = "purple", "Pbea" = "orange", "Pgra" = "cyan",
                                "Ppur" = "pink", "Ppus" = "brown", "Qapi" = "yellow",
                                "Qver" = "black", "Uhar" = "darkgreen")) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

ggplot(data = df, aes(x = log10(STDM.g), y = log10(NH.umol.h))) +
  geom_point(alpha = 0.7, size = 3) +  # Adjust transparency and point size
  geom_smooth(method = "lm", se = F, color = "black", linetype = "solid") +
  stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, size = 5) +  # Add R² and p-value
  theme_classic() +
  labs(x = "STDM (g)", y = "NH (umol/h)", colour = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

##Plot NP ratio
ggplot(data = df) +
  geom_boxplot(aes(x =fct_reorder(Species, NP, .fun=mean, .desc = F), 
                   y =  NP, fill = Tribe), outlier.shape = NA) +
  geom_point(aes(x = fct_reorder(Species, NP, .fun = mean, .desc = F), 
                 y =  NP), position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic() +
  labs(x = "", y = "ln(NP)")


# ANOVA - log transformed NP for species
lm1 <- lm(log(NP)~Species, data = df); summary(lm1)
em_lm1 <- emmeans(lm1, pairwise ~ Species); summary(em_lm1)

# ANOVA - log transformed NP for tribe
lm7 <- lm(log(NP)~ Species, data = df); summary(lm7)
em_lm7 <- emmeans(lm7, pairwise ~ Species); summary(em_lm7)

lm2 <- lm(log10(NH.umol.g.h) ~ Species, data = df); summary(lm2)
em_lm2 <- emmeans(lm2, pairwise ~ Species, adjust = "tukey"); summary(em_lm2)

lm3 <- lm(log10(SRP.umol.g.h) ~ Tribe, data = df); summary(lm3)
em_lm3 <- emmeans(lm3, pairwise ~ Tribe, adjust = "tukey"); summary(em_lm3)

###### ANCOVA - log transformed
lm4 <- lm(log10(NH.umol.h)~log10(STDM.g) * Species, data = df); summary(lm4); anova(lm4)
em_lm4 <- emmeans(lm4, pairwise ~ log10(STDM.g) * Species);summary(em_lm4)
plot(em_lm4)
contrast(em_lm4, "pairwise")

lmP <- lm(log10(SRP.umol.h)~log10(STDM.g) * Species, data = df); summary(lmP); anova(lmP)
em_lmP <- emmeans(lmP, pairwise ~ log10(STDM.g) * Species);summary(em_lmP)


contrast(em_lm4, "pairwise")
lm5 <- lm(log10(SRP.umol.h)~log10(STDM.g) * Species, data = df); summary(lm5) ##Interactive model
lm6 <- lm(NH.umol.h ~ STDM.g + Species, data = df); summary(lm6)  ##Additive model

######
emm_Species <- emmeans(lm4, pairwise ~ Species)
emm_Species$contrasts

##ANCOVA
lm5 <- lm(NH.umol.h ~ STDM.g * Species, data = df); summary(lm5)  ##Interactive model
lm6 <- lm(NH.umol.h ~ STDM.g + Species, data = df); summary(lm6)  ##Additive model


AIC(lm5, lm6)
BIC(lm5, lm6)

# Plot of STDM on excreta rates
ggplot(df, aes(x = log10(STDM.g), y = log10(NH.umol.h), color = Species)) +
  geom_point() +  # Scatter plot of data points
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Regression lines per species
  theme_classic() +
  labs(x = "STDM.g", y = "NH.umol.h")

ggplot(df, aes(x = log10(STDM.g), y = log10(SRP.umol.h), color = Species)) +
  geom_point() +  # Scatter plot of data points
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Regression lines per species
  theme_classic() +
  labs(x = "STDM.g", y = "SRP.umol.h")

em_lm4_df <- as.data.frame(em_lm4)

# Plot
ggplot(em_lm4, aes(x = log10.STDM.g.trend, y = Species)) +
  geom_point(size = 4, color = "blue") +  # Points for slopes
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +  # Horizontal error bars
  labs(x = "Slope of log10(STDM.g) on log10(NH.umol.h)", y = "Species",
       title = "Species-Specific Slopes for log10(STDM.g)") +
  theme_minimal()

