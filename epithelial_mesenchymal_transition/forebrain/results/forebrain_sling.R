setwd("~/Documents/SURF 2020/RNA_seq_files/epithelial_mesenchymal_transition/forebrain/results")
library(slingshot)
library(scales)
library(viridis)
df <- read.csv("forebrain_slingshot_data.csv", header=TRUE, row.names="X")
head(df)
UMAPs <- c("X_umap_1", "X_umap_2")
UMAP_coords <- df[UMAPs]
head(UMAP_coords)
cluster_labels = df$leiden
cluster_labels
sling_output<-slingshot(UMAP_coords, 
                        clusterLabels = cluster_labels, 
                        start.clus = 5, 
                        end.clus = c(2, 0),
                        approx_points = 200)
sling_output
# cluster_labels is an array of length N (cells) with their corresponding cluster
# assuming we have k clusters, we can create a new array of length N from a palette of k colors

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
cell_colors_clust <- cell_pal(cluster_labels, hue_pal())

plot(reducedDim(sling_output), 
     col = cell_colors_clust, # A color array: one color for each cell (color corresponds to cluster)
     pch = 16, # graphical parameters for scatterplot: type of dot
     cex = 0.5, # size of dots
     cex.lab = 1.5, # label font size 
     cex.axis = 1) # axis font size
lines(sling_output, lwd = 2, col = 'black')
pt <- slingPseudotime(sling_output)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(5, 6))

for (i in nms) {
  # This command subsets the trajectory i 
  # And bins the numeric values using 100 breaks
  # Using our palette array pal, each bin gets a color that represents their progression in pseudotime		
  colors <- pal[cut(pt[,i], breaks = 100)]
  
  # Cells that were not mapped to the current trajectory i 
  # will have NA values, which the colors array will keep as NA 
  # and the plotting function ignores by default. So we can do: 		
  plot(UMAP_coords, col = colors, pch = 16, cex = 0.5, main = i)
  
  # This plots the min spanning tree (analogous to PAGA in scanpy)
  #lines(sling_output, lwd = 1, col = 'black', type = 'lineages')
}
write.table(pt, 
            file = './pseudotime_trajectories_slingshot.csv', 
            quote = F, #removes annoying quotes when writing to file,
            sep=","
)

