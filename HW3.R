library(readxl)
library(vegan)
library(ggplot2)
library(plotly)
library(impute)
library(factoextra)
library(psych)
library(ggforce)
library(rstatix)

# Reading data
df <- read_excel('Sleepy lizard.xlsx')

# Selecting only necessary variables
df <- df %>% 
  select(Treatment, Habitat, Connectivity, Tot_WBC, Het_ABS, Lym_ABS, `H:L Ratio`, Mon_ABS, OthG_ABS, LBSI)

# Converting variables to appropriate types
df$Treatment <- as.factor(df$Treatment)
df$Habitat <- as.factor(df$Habitat)
df$Connectivity <- as.factor(df$Connectivity)

# Exploratory Data Analysis
summary(df)

# Creating a subset for blood-related variables
df_blood <- df %>% 
  select(-Treatment, -Habitat, -Connectivity)

# Visualizing the distribution of blood variables
boxplot(df_blood, las = 2)

# Transformation and normalization of blood variables
df_blood_lognorm <- scale(log(df_blood + 1), scale = FALSE)
boxplot(df_blood_lognorm, las = 2)

# Performing PCA
pca_df <- prcomp(df_blood_lognorm)$x %>% 
  as.data.frame() %>%
  select(PC1, PC2) %>% 
  mutate(Treatment = df$Treatment,
         Habitat = df$Habitat)

# Visualizing PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, shape = Treatment, col = Habitat)) +
  geom_point(size = 3) +
  geom_mark_ellipse(aes(col = Treatment, fill = Treatment))

# perMANOVA
# Checking dispersion equality in groups

# Dispersion is equal in groups by Habitat, Connectivity, and Treatment
dist_blood <- vegdist(df_blood_lognorm, method = "euclidean")
PCO_blood <- betadisper(dist_blood, df$Habitat)
plot(PCO_blood)
anova(PCO_blood)

dist_blood <- vegdist(df_blood_lognorm, method = "euclidean")
PCO_blood <- betadisper(dist_blood, df$Treatment)
plot(PCO_blood)
anova(PCO_blood)

dist_blood <- vegdist(df_blood_lognorm, method = "euclidean")
PCO_blood <- betadisper(dist_blood, df$Connectivity)
plot(PCO_blood)
anova(PCO_blood)

# Statistical Analysis

# Task 1

# a: Blood composition ~ Treatment
adonis2(df_blood_lognorm ~ df$Treatment, method = "euclidean")

# b: Blood composition ~ Habitat
df_modified <- df_blood_lognorm %>% as.data.frame() %>% filter(df$Treatment == 2)
hab_modified <- subset(df, Treatment == 2)$Habitat
adonis2(df_modified ~ hab_modified, method = "euclidean")

# c: Blood composition ~ Connectivity
conn_modified <- subset(df, Treatment == 2)$Connectivity
adonis2(df_modified ~ conn_modified, method = "euclidean")

# Task 2

adonis2(df_modified ~ conn_modified + hab_modified, method = "euclidean")

# Conclusion

# The perMANOVA results indicate a significant difference in blood composition between lizards from severely modified and unmodified landscapes.
# However, no significant difference is observed in blood composition of lizards from severely modified landscapes, both by habitat type and connectivity.
# This outcome is as expected since PCA analysis revealed observable distances between lizards from different "Treatment" groups, 
# while lizards from all habitat types are clustered together (in the modified landscape group), suggesting small differences between them.
