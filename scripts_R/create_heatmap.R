library( RColorBrewer );
library( BoutrosLab.plotting.general );
library( scales );
library( reshape2 );
library( gridExtra );

#Load RData with LFC data of DEG
load("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/lfc_DEG.RData")

# Create colour pallet list
colour.pallet <- list();

# Define a colour pallete to samples
colour.pallet$samples <- alpha( brewer.pal( n=8, name = "Set2" )[ c(1,2,3,4,5,6) ], 0.8 );
names( colour.pallet$samples ) <- c("Y1","Y2","Y3","F1","F2","F3");

# Define a colour pallete for heatmap
colour.pallet$heatmap.cells = brewer.pal( n=11 , name = "RdBu")[ c(10,6,2) ];

### Create a heatmap of log2FC from RNA-seq results

# Initialize object 'fc.heatmap' that will hold data from heatmap
fc.heatmap <- list();

# Create covariates to the columns (Samples)
fc.heatmap$samples.covariate <- list(
  rect = list(
    col = 'black',
    fill = colour.pallet$samples[ colnames(lfc_DEG) ],
    lwd = 0.5
  )
);

# Join covariate legends
fc.heatmap$covariates.legend <- list(
  legend = list(
    colours = colour.pallet$samples[ colnames(lfc_DEG) ],
    labels = colnames(lfc_DEG),
    title = 'Sample',
    border = 'white'
  ),
  legend = list(
    colours = colour.pallet$heatmap.cells,
    labels = round( max( abs( range(lfc_DEG) ) ), 2) * c(-1,0,1),
    title = 'Fold-change (log2)',
    continuous = TRUE,
    height = 9,
    pos.x = 0.21,
    tck = 1,
    at = c(0,50,100)
  )
);

# Dendrogram provided
fc.heatmap$dendrogram <- create.dendrogram(
  x = lfc_DEG,
  cluster.dimension = 'col'
);

plot(fc.heatmap$dendrogram)

# Create heatmap
fc.heatmap$plot <- create.heatmap(
  x = t(lfc_DEG),
  
  #cluster.dimensions = "rows",
  colour.scheme = colour.pallet$heatmap.cells,
  
  xlab.label = 'Samples',
  xaxis.lab = NULL,
  axis.xlab.padding = 1,
  xlab.cex = 1.5,
  
  ylab.label = 'Genes',
  yaxis.lab = NULL,
  ylab.cex = 1.5,
  
  print.colour.key = FALSE,
  
  covariates.top = fc.heatmap$samples.covariate,
  covariates.top.grid.border = list(col = 'white', lwd = 1.5),
  covariates.top.grid.row = list(col = 'white', lwd = 1.5),
  covariates.top.grid.col = list(col = 'white', lwd = 1.5),
  covariate.legends = fc.heatmap$covariates.legend,
  
  # legend customization
  legend.side = 'right',
  legend.title.just = 'left',
  legend.title.cex = 1.1,
  
  style = "Nature",
  
  description = 'Heatmap created using BoutrosLab.plotting.general'
);

# Export heatmap
setwd("./plots/")
tiff(
  filename = "fold_change_heatmap.tiff", 
  type = "cairo", 
  units = "cm",
  width = 14, 
  height = 25, 
  res = 500, 
  compression = 'lzw'
);

plot(fc.heatmap$plot);

dev.off();
