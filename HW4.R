# Import required libraries
library(readxl)
library(corrplot)
library(VIM)
library(vegan)
library(ggvegan)
library(dplyr)

# Load data from Excel sheets
fauna_data <- read_excel('Grazing_Magierowski_et_al_2015.xls', sheet = 'fauna')
environment_data <- read_excel('Grazing_Magierowski_et_al_2015.xls', sheet = 'env')
coordinates_data <- read_excel('Grazing_Magierowski_et_al_2015.xls', sheet = 'coord')
raw_data <- read_excel('Grazing_Magierowski_et_al_2015.xls', sheet = 'raw', skip=1)

# Rename column names for clarity
colnames(environment_data) <- c("SITE","Abstraction","Regulation","Grazing","Fines","Temperature","Conductivity","AvgTurbidity","pH",
                                "Alkalinity","NitrateNitrite","DRP","N_total","P_total","AvgShading","AvgAlgae","Chl","GrazingRank")

# Merge relevant data
merged_data <- merge(environment_data, fauna_data)

# Remove missing values
clean_data <- na.omit(merged_data)

# Visualize correlations
corrplot(cor(clean_data[,2:17]))

# Subset environmental and fauna data
environment <- clean_data[,1:18]
fauna <- clean_data[, 19:length(clean_data)]

# Plot boxplot of environmental data
boxplot(environment[,2:17], las=2)

# Perform log normalization of environmental data
log_normalized_environment <- scale(log(environment[,2:17]+1), scale = FALSE)
boxplot(log_normalized_environment, las=2)
clean_data[,2:17] <- log_normalized_environment

# Perform Canonical Correspondence Analysis (CCA)
grazing_cca <- cca(clean_data[,19:length(clean_data)] ~ Abstraction + Grazing + Fines + Temperature + Conductivity + 
                     AvgTurbidity + pH + N_total + P_total + AvgShading + Chl, data = clean_data)

# Calculate Variance Inflation Factor (VIF) of CCA
vif_cca <- vif.cca(grazing_cca)

# Summarize results of CCA
summary_results <- summary(grazing_cca)

# Visualize Screeplot of CCA
screeplot(grazing_cca,  bstick = TRUE)

# Generate autoplot of CCA
autoplot(grazing_cca, scaling = "sites")

# Plot scaling 1 of CCA
plot(grazing_cca, scaling = "sites", 
     main = "Scaling 1, 'Sites' ")

# Generate biplot of CCA with scaling 2
plot(grazing_cca, scaling = 2, 
     display = c("species", "cn"), 
     main = "Biplot CCA, Scaling 2")

# Perform ANOVA of CCA
anova_results <- anova(grazing_cca)

# Perform ANOVA of CCA by term
anova_by_term_results <- anova(grazing_cca, by="term")

# Fit second CCA model
grazing_cca_2 <- cca(fauna ~ Grazing*Abstraction + Fines + Temperature + Conductivity + 
                       AvgTurbidity + pH + N_total*P_total + AvgShading + Chl + GrazingRank, data = environment)

# Perform ANOVA of second CCA model by term
anova_cca_2_by_term_results <- anova(grazing_cca_2, by="term")

# Perform ANOVA of CCA by marginal effects
anova_cca_by_mar_results <- anova(grazing_cca, by="mar")

