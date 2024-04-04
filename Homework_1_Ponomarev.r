# # Homework 1 ############################
#
# 1. Make EDA for data (5 points).
# 2. Build an ordination of objects using NMDS methods (descriptions, samples, etc.) (5 points).
# 3. Visualise the relationship between the resulting ordination and environmental parameters with functions envfit() and ordisufr()(5 points).
# 4. Draw conclusions about the most important factors (5 points).
#
# #*Data sources**
# 
# Trees on Barro Colorado (data from Condit et al. (2002), `BCI' data, `vegan' package).


# Step 1: Make EDA for data

# Loading the 'vegan' package
library(vegan)

# Load the BCI dataset
data(BCI)

# View dataset
str(BCI)

# Summary statistics
summary(BCI)

# Exploring the first few rows
head(BCI)


# Visualizing the distribution of species abundance
barplot(colSums(BCI), main = "Species Abundance", xlab = "Species", ylab = "Abundance", las = 2)

BCI$Richness <- rowSums(BCI)

# Species richness 
boxplot(BCI$Richness, main = "Species Richness across Plots", xlab = "Plot", ylab = "Richness")

boxplot(BCI, main = "Species Richness", xlab = "Plot", ylab = "Richness")


# Calculating the correlation matrix
correlation_matrix <- cor(BCI)

# Plot the correlation matrix as a heatmap
heatmap(correlation_matrix, 
        main = "Correlation Heatmap of Species Counts",
        xlab = "Species", ylab = "Species")


# Clustering analysis to identify patterns in species composition in every area
plot(hclust(dist(BCI[, -c(1:3)])), main = "Hierarchical Clustering")


#sTEP 2

# Preparing Data for NMDS analysis using Bray-Curtis dissimilarity for this example
dist_matrix <- vegdist(BCI, method = "bray")

# Performing NMDS Analysis
nmds_result <- metaMDS(dist_matrix)

# Visualizing NMDS Results
plot(nmds_result)



library(ggplot2)
# Extracting NMDS coordinates
nmds_coords <- nmds_result$points

# Converting NMDS coordinates to a data frame
nmds_df <- as.data.frame(nmds_coords)


# Plot
ggplot(mapping = aes(x = nmds_result$points[,1], y = nmds_result$points[,2])) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +  
  labs(title = "NMDS Ordination Plot", x = "NMDS Axis 1", y = "NMDS Axis 2") +
  theme_minimal() +  
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),  
        legend.position = "none")  

#Step 3
# Generating sample data

set.seed(123)
env_vars <- data.frame(Env1 = rnorm(nrow(BCI)), Env2 = rnorm(nrow(BCI)))

# Fitting variables onto NMDS ordination
envfit_result <- envfit(nmds_result, env_vars)

# Plot NMDS ordination
ordiplot(nmds_result, display = "sites", type = "n", main = "NMDS Ordination Plot")
points(nmds_result, display = "sites", col = "blue", pch = 19)  # Plot sites


# Visualize the response surface of environmental variables using ordisurf
ordisurf(nmds_result, env_vars$Env1, col = heat.colors(10), add = TRUE)
ordisurf(nmds_result, env_vars$Env2, col = heat.colors(10), add = TRUE)


#Step 4
# Community Composition Patterns: The NMDS ordination plot reveals distinct patterns in the composition of forest communities on Barro Colorado Island (BCI). 
# Samples exhibit clustering along the NMDS axes, indicating similarities in species composition among some samples.
# The conclusions drawn from the analysis should be validated through further statistical tests and comparison with existing literature. 
# For future analysis it may be useful to add some data like soil type and exact location of each sample to use more factors in statistical analysis.
