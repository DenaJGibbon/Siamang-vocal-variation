# Load required libraries
library(randomForest)  # For random forest modeling
library(cluster)       # For clustering algorithms
library(ggpubr)
library(umap)# For data visualization
set.seed(3)

# Load data
# (Assuming that the data required for the following operations is already loaded into 'siamangUS2df')

# Siamang unsupervised clustering

# Check the number of rows in the 'siamangUS2df' data frame
nrow(siamangUS2df)

# Convert certain columns to numeric using dplyr mutate_if
# Columns 1 to 3 and 27 to 30 are excluded from conversion
siamangUS2df[, -c(1:3, 18, 27, 28, 29)] <-
  siamangUS2df[, -c(1:3, 18, 27, 28, 29)] %>% mutate_if(is.character, as.numeric)

# Perform cluster validation with clValid to find the optimal number of clusters
intern <-
  clValid(
    siamangUS2df[, -c(1:3, 18, 27, 28, 29)],
    3:10,
    clMethods = c("hierarchical", "kmeans", "pam", "som", "Mclust"),
    maxitems = nrow(siamangUS2df),
    validation = "internal"
  )

# Get the order of clustering method used by hierarchical clustering
intern@clusterObjs$hierarchical$order

# Summary of cluster validation results
summary(intern)

# Get the measures used for cluster validation
intern@measures

# Get the optimal scores for the clustering methods
optimalScores(intern)

# Retrieve the cluster assignments for the selected clustering method (kmeans)
clusterpam <- clusters(intern, "pam")
nclusters <-
  clusterpam$`4`$cluster  # 'nclusters' contains the cluster assignments
#nclusters <- cutree(clusterpam,2)

# Perform UMAP (Uniform Manifold Approximation and Projection) for dimensionality reduction
SiamangUS2.umap <-
  umap::umap(
    siamangUS2df[, -c(1:3, 18, 27, 28, 29)],
    labels = nclusters,
    controlscale = TRUE,
    scale = 3
  )

# Prepare data for plotting with UMAP results and additional information
plot.for.SiamangUS2s <-
  cbind.data.frame(SiamangUS2.umap$layout[, 1:2], siamangUS2df$site)
colnames(plot.for.SiamangUS2s) <- c("Dim.1", "Dim.2", "site")
plot.for.SiamangUS2s$Clusters <- as.factor(nclusters)

# Create a scatter plot to visualize US-II sites colored by 'site'
Siamang.site <-
  ggplot(plot.for.SiamangUS2s, aes(x = Dim.1, y = Dim.2, col = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors(length(unique(
    plot.for.SiamangUS2s$site
  )))) +
  labs(color = 'Site') +
  ggtitle('A. US-II sites') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Create a scatter plot to visualize US-II pairs colored by 'pair'
pair <- siamangUS2df$individual

Siamang.pair <-
  ggplot(plot.for.SiamangUS2s, aes(x = Dim.1, y = Dim.2, col = pair)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors(length(unique(pair)))) +
  labs(color = 'Pair') +
  ggtitle('B. US-II pairs') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Create a scatter plot to visualize US-II unsupervised clustering colored by 'Clusters'
Siamang.cluster <-
  ggplot(plot.for.SiamangUS2s, aes(x = Dim.1, y = Dim.2, col = Clusters)) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors(length(unique(
    plot.for.SiamangUS2s$Clusters
  )))) +
  ggtitle('C. US-II unsupervised clustering') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Arrange the scatter plots in a grid using cowplot
cowplot::plot_grid(Siamang.site, Siamang.pair, Siamang.cluster, nrow = 2)

