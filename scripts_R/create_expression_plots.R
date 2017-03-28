### create_expression_plots.R #####################################################################
# This script loads the output from run_DE_analysis.R and creates a series of plots for
# expression analysis and comparison (PCA, heatmap and scatterplots). Figure 2 from paper

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd( "~/angiogenesis_transcriptome/scripts/" );

# Add some locally installed libraries to path
.libPaths( "~/R/lib/" );

# Load required libraries and datasets
library( DESeq2 );
# library( org.Mm.eg.db );
library( RColorBrewer );
library( BoutrosLab.plotting.general );
library( scales );
library( reshape2 );
library( gridExtra );

### FUNCTIONS #####################################################################################

### Function - factor2numeric() #######################################################################
# Input variables: A factor with values to be converted
#
# Output variables: A vector with numeric values
#
# Description: Convert a factor to a vector of numeric values
factor2numeric <- function( x ) {
  
  stopifnot( inherits( x, "factor" ) );
  
  return( as.numeric( levels(x) )[x] );
  }

### DATA ANALYSIS #################################################################################

# Load expression data info
load( file = "aux_files/expression_data.Rdata" );

### Define colour pallet to be used on all plots

# Define heatmap pallet to be used
# display.brewer.all()

# Create colour pallet list
colour.pallet <- list();

# Define a colour pallete to samples
colour.pallet$samples <- alpha( brewer.pal( n=11, name = "RdYlGn" )[ c(7,9,11,5,4,2,1) ], 0.8 );
names( colour.pallet$samples ) <- c("P12", "P15", "P17", "R12", "R12.5", "R15", "R17" );

# Define a colour pallete for heatmap
colour.pallet$heatmap.cells = brewer.pal( n=11 , name = "RdBu")[ c(10,6,2) ];

### Create PCA from samples

# Filter only genes with > 500 reads on all samples for PCA
selected <- rowSums( counts( ddsFull ) >= 500 ) == ncol(ddsFull);
pca <- plotPCA( rlog( ddsFull[selected, ] ), intgroup=c("Treat", "Time", "Batch") );
# plot(pca)

# PCA plot
pca.plot <- list();

# Save PCA plot
pca.plot$final.plot <- create.scatterplot(
  formula = PC2 ~ PC1,
  data = pca$data,
  main = "\n",
  
  # Control points colours and shape
  col = colour.pallet$samples[ paste0( pca$data$Treat, pca$data$Time) ],
  pch = 16,
  cex = 2.5,
  
  # Configure axis labels
  # xlab.label = pca$labels$x,
  xlab.label = paste0( gsub(": ", " (", pca$labels$x ), ")" ),
  xaxis.fontface = 1,

  # ylab.label = pca$labels$y,
  ylab.label = paste0( gsub(": ", " (", pca$labels$y ), ")" ),
  yaxis.fontface = 1,

  key = list(
    corner = c(0.45,1.15),
    cex   = 1.3,
    columns = 4,
    padding.text = 1,
    between = 0.5,
    between.columns = 1,
    rectangles = list(
      col = colour.pallet$samples[ c("P12", "R12", "", "R12.5", "P15", "R15", "P17", "R17" ) ],
      size = 4.375,
      border = "white" ),
    text  = list(
      label = c("P12", "R12", "", "R12.5", "P15", "R15", "P17", "R17" )
      )
    ),

  style = "Nature"
  
  );

# Export PCA figure
tiff(
  filename = "~/Desktop/PCA.tiff", 
  type = "cairo", 
  width = 20, 
  height = 20, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

# grid.draw( pca.plot$final.plot );
plot( pca.plot$final.plot );

dev.off();


### Create a heatmap of log2FC from RNA-seq results

# Initialize object 'fc.heatmap' that will hold data from heatmap
fc.heatmap <- list();

# Create covariates to the columns (Samples)
fc.heatmap$samples.covariate <- list(
  rect = list(
    col = 'black',
    fill = colour.pallet$samples[ colnames(de$P12.fold.change) ],
    lwd = 0.5
    )
  );

# Join covariate legends
fc.heatmap$covariates.legend <- list(
  legend = list(
    colours = colour.pallet$samples[ colnames(de$P12.fold.change) ],
    labels = colnames(de$P12.fold.change),
    title = 'Sample',
    border = 'white'
    ),
  legend = list(
    colours = colour.pallet$heatmap.cells,
    labels = round( max( abs( range(de$P12.fold.change) ) ), 2) * c(-1,0,1),
    title = 'Fold-change (log2)',
    continuous = TRUE,
    height = 9,
    pos.x = 0.21,
    tck = 1,
    at = c(0,50,100)
    )
  );

# Create heatmap
fc.heatmap$plot <- create.heatmap(
  #
  x = t(de$P12.fold.change),
  
  cluster.dimensions = "rows",
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
tiff(
  filename = "~/Desktop/fold_change_heatmap.tiff", 
  type = "cairo", 
  units = "cm",
  width = 14, 
  height = 25, 
  res = 500, 
  compression = 'lzw'
  );

plot( fc.heatmap$plot );

dev.off();

### Create a scatter plot of TLDA ~ RNA-seq results

# Create list to hold info about plot
rnaseq.tlda.scatter <- list();

# Read TLDA results file by biogroups
rnaseq.tlda.scatter$tlda.results <- read.csv(
  file = "aux_files/TLDA/Biogroup Results.csv",
  header = TRUE, 
  comment.char = "#"
  );

# Leave only Gene Symbol as target name 
rnaseq.tlda.scatter$tlda.results$Target.Name <- sub(
  pattern = "-(Mm|Hs).*$",
  replacement = "",
  x = rnaseq.tlda.scatter$tlda.results$Target.Name
  );

# Subset TLDA results to columns of interest and cast to appropriated data type (Warnings to Endogenous Control genes will appear, ok!)
rnaseq.tlda.scatter$tlda.results <- data.frame(
  condition = factor( rnaseq.tlda.scatter$tlda.results$Bio.Group.Name ),
  gene = factor( rnaseq.tlda.scatter$tlda.results$Target.Name ),
  mean.CT = rnaseq.tlda.scatter$tlda.results$Mean.Equivalent.Cq,
  rq = factor2numeric( droplevels( rnaseq.tlda.scatter$tlda.results$Rq ) ),
  rq.max = factor2numeric( droplevels( rnaseq.tlda.scatter$tlda.results$Rq.Max ) ),
  rq.min = factor2numeric( droplevels( rnaseq.tlda.scatter$tlda.results$Rq.Min ) )
  );

# Create factor that will be used in contrast
rnaseq.tlda.scatter$dds <- ddsFull;
rnaseq.tlda.scatter$dds$Group <- factor( paste0( rnaseq.tlda.scatter$dds$Treat, rnaseq.tlda.scatter$dds$Time) );
design(rnaseq.tlda.scatter$dds) <- ~ Group;
rnaseq.tlda.scatter$dds <- DESeq(rnaseq.tlda.scatter$dds);

# Select genes that will be used on plot
# 1) Remove genes not evaluated by RNAseq (rRNA 18S)
# 2) Endogenous control for TLDA (Sdha and Tbp).
# *Muc2 gene was splitted in two on new ENSEMBL annotation so we decided to remove it from analysis as well
rnaseq.tlda.scatter$comparable.genes <- as.character( unique( rnaseq.tlda.scatter$tlda.results$gene ) );
rnaseq.tlda.scatter$comparable.genes <- rnaseq.tlda.scatter$comparable.genes[! rnaseq.tlda.scatter$comparable.genes %in% c("18S","Rn18s", "Sdha", "Tbp", "Muc2") ];

# Create a data.frame with geneSymbols and respective ensemblIDs
rnaseq.tlda.scatter$comparable.genes <- data.frame(
  symbol     = rnaseq.tlda.scatter$comparable.genes,
  ensembl.id = mapIds(
                  x = org.Mm.eg.db, 
                  keys = rnaseq.tlda.scatter$comparable.genes, 
                  keytype = "SYMBOL", 
                  column = "ENSEMBL", 
                  multiVals="asNA"
                  ),
  row.names  = NULL
  );

# Build data.frame for RNAseq results
rnaseq.tlda.scatter$rnaseq.results <- sapply( 
  c("P15", "P17", "R12", "R12.5", "R15", "R17"),
  function(sample.name) results(
    rnaseq.tlda.scatter$dds, 
    contrast = c("Group", sample.name, "P12")
    )[ rnaseq.tlda.scatter$comparable.genes$ensembl.id, "log2FoldChange"]
  );
rownames(rnaseq.tlda.scatter$rnaseq.results) <- rnaseq.tlda.scatter$comparable.genes$symbol;

# P12
rnaseq.tlda.scatter$rnaseq.results <- as.matrix( cbind( data.frame(P12 = 0), rnaseq.tlda.scatter$rnaseq.results ) );

# Melt the data frame into data.frame  
# Work with only genes in gene.list and remove reference condition (P12, Rq=1)
rnaseq.tlda.scatter$scatter.data <- melt( rnaseq.tlda.scatter$rnaseq.results );
colnames(rnaseq.tlda.scatter$scatter.data) <- c("gene", "condition", "RNASeq");

# Identify treat and time columns
rnaseq.tlda.scatter$scatter.data$treat <- sub("(R|P).*", "\\1", rnaseq.tlda.scatter$scatter.data$condition );
rnaseq.tlda.scatter$scatter.data$time <- as.numeric( sub("(R|P)", "", rnaseq.tlda.scatter$scatter.data$condition ) );

# Obtain TLDA results from selected experiments. Transform them into log2 (same as RNA-seq)
rnaseq.tlda.scatter$scatter.tlda.results <- do.call(
  what = "rbind",
  args = apply(
    rnaseq.tlda.scatter$scatter.data,
    1,
    function(x) subset( rnaseq.tlda.scatter$tlda.results, gene == x["gene"] & condition == x["condition"] )[, c("rq", "rq.max", "rq.min") ]
    )
  );
rnaseq.tlda.scatter$scatter.tlda.results <- log2(rnaseq.tlda.scatter$scatter.tlda.results);

# Append results into data.frame
rnaseq.tlda.scatter$scatter.data <- cbind(rnaseq.tlda.scatter$scatter.data, rnaseq.tlda.scatter$scatter.tlda.results);

# # Read TLDA results by well
# rnaseq.tlda.scatter$tlda.wells <- read.csv(
#   file = "aux_files/TLDA/Well Results.csv",
#   header = TRUE, 
#   comment.char = "#"
#   );
# 
# # Leave only Gene Symbol as target name 
# rnaseq.tlda.scatter$tlda.wells$Target.Name <- sub(
#   "-(Mm|Hs).*$",
#   "",
#   rnaseq.tlda.scatter$tlda.wells$Target.Name
#   );
# 
# # Filter wells that were ommited in the preprocessing steps
# rnaseq.tlda.scatter$tlda.wells <- rnaseq.tlda.scatter$tlda.wells[rnaseq.tlda.scatter$tlda.wells$Omit == "false", ];
# 
# # Subset columns of interest
# rnaseq.tlda.scatter$tlda.wells <- data.frame(
#   Gene = rnaseq.tlda.scatter$tlda.wells$Target.Name,
#   CT = factor2numeric( droplevels( rnaseq.tlda.scatter$tlda.wells$Eq..Cq ) ),
#   condition = rnaseq.tlda.scatter$tlda.wells$Biological.Group.Name,
#   Sample = factor( rnaseq.tlda.scatter$tlda.wells$Sample.Name )
#   );
# 
# # Get list of experiments that are above threshold
# rnaseq.tlda.scatter$high.CT.experiments <- as.character(
#   apply(
#     rnaseq.tlda.scatter$tlda.wells[ rnaseq.tlda.scatter$tlda.wells$CT > 35, c("gene", "condition") ],
#     1,
#     paste,
#     collapse = '-'
#     )
#   );

# Identify experiments with CT above safety threshold (35 cycles)
rnaseq.tlda.scatter$high.CT.experiments <- apply(
  rnaseq.tlda.scatter$tlda.results[ rnaseq.tlda.scatter$tlda.results$mean.CT > 35, c("gene", "condition") ],
  1,
  paste,
  collapse = '-'
  );

# Mark experiments with high CT
rnaseq.tlda.scatter$scatter.data$safe.experiment <- ! paste( rnaseq.tlda.scatter$scatter.data$gene, rnaseq.tlda.scatter$scatter.data$condition, sep = "-" ) %in% rnaseq.tlda.scatter$high.CT.experiments;

# Make plot
rnaseq.tlda.scatter$scatter.plot <- create.scatterplot(
  #
  formula = rq ~ RNASeq,
  # data = rnaseq.tlda.scatter$scatter.data,
  data = subset(rnaseq.tlda.scatter$scatter.data, safe.experiment),
  type = c("p", "r"),
  
  # Customize points style  
  pch = 19,
  cex = 1,
  col = "black",
  lwd = 3,

  # Format axis
  xlimits = c(-5,10),
  xlab.label = "RNASeq (log2 FC)",
  xat = seq(-5,10,5),

  ylimits = c(-5,10),
  ylab.label = "qRT-PCR (log2 FC)",
  yat = seq(-5,10,5),

  # Adding legend
  print.new.legend = TRUE,
  legend = list (

  # Adding correlation key to plot 1
  inside = list(
    fun = draw.key,
    args = list(
      key = get.corr.key(
        # x = rnaseq.tlda.scatter$scatter.data$RNASeq,
        # y = rnaseq.tlda.scatter$scatter.data$rq,
        x = subset(rnaseq.tlda.scatter$scatter.data, safe.experiment)$RNASeq,
        y = subset(rnaseq.tlda.scatter$scatter.data, safe.experiment)$rq,
        label.items = c('pearson', "pearson.p"),
        alpha.background = 0,
        key.cex = 2.5
        )
      ),
      x = 0.02,
      y = 0.97,
      corner = c(0,1)
      )
    ),

  style = "Nature",
  
  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Export heatmap
tiff(
  filename = "~/Desktop/rnaseq_tlda_scatter_plot.tiff", 
  type = "cairo", 
  units = "cm",
  width = 20, 
  height = 20, 
  res = 500, 
  compression = 'lzw'
  );

plot( rnaseq.tlda.scatter$scatter.plot );

dev.off();

### Create scatter plots for each gene showing fold-change along with time course

# Add legends to samples
rnaseq.tlda.scatter$timed.fold.plot$legend <- list(
  legend = list(
    colours = colour.pallet$samples[ c("P12", "P15", "P17", "R12", "R12.5", "R15", "R17") ],
    labels = c("P12", "P15", "P17", "R12", "R12.5", "R15", "R17"),
    title = 'Sample',
    border = 'white'
    ),
  # Insert text but leave space for 'key' command (legend.grob can't work with pch)
  legend = list(
    colours = c("transparent", "transparent"),
    labels = c("RNAseq", "TLDA"),
    title = 'Experiment',
    border = 'transparent'
    )
  );

for( c.gene in unique( as.character( rnaseq.tlda.scatter$scatter.data$gene ) ) ){

  # Subset results relative to that gene
  rnaseq.tlda.scatter$timed.fold.plot$gene.subset <- subset(rnaseq.tlda.scatter$scatter.data, gene == c.gene);
  
  # Create rq ~ time plot (TLDA)
  rnaseq.tlda.scatter$timed.fold.plot[[c.gene]] <- create.scatterplot(
    #
    formula = rq ~ time,
    data = rnaseq.tlda.scatter$timed.fold.plot$gene.subset,
    main = c.gene,
  
    pch = 15,
    cex = 2,
    col = colour.pallet$samples[ as.character( rnaseq.tlda.scatter$timed.fold.plot$gene.subset$condition ) ],
    fill = 'transparent',
  
    # Format axes
    xat = c(12, 12.5, 15, 17),
    xaxis.lab = c("12   ", "   12.5", "15", "17"),
    xlab.label = 'Time',
    xaxis.fontface = 1,
  
    ylab.label = 'Fold-change (log2)',
    ylimits = floor( range( rnaseq.tlda.scatter$timed.fold.plot$gene.subset[, c("RNASeq", "rq.max", "rq.min") ] ) ) + c(-0.5, 1.5),
    yaxis.fontface = 1,
  
    # Specify error bars
    error.bar.lwd = 1,
    y.error.up = rnaseq.tlda.scatter$timed.fold.plot$gene.subset$rq.max - rnaseq.tlda.scatter$timed.fold.plot$gene.subset$rq,
    y.error.down = rnaseq.tlda.scatter$timed.fold.plot$gene.subset$rq - rnaseq.tlda.scatter$timed.fold.plot$gene.subset$rq.min,
    y.error.bar.col = "darkgray",
    error.bar.length = 0.05,
    
    legend = list(
      right = list(fun = legend.grob( legends = rnaseq.tlda.scatter$timed.fold.plot$legend ) )
      ),
    
    # Insert pch symbols that couldn't be added by legend command
    key = list(
      columns = 1,
      padding.text = 1,
      between = 0.1,
      between.columns = 0.8,
      x = 1.105,
      y = 0.30,
      points = list(
        pch = c(18,15),
        cex = 2,
        border = "white" )
      ),
    
    style = "Nature",
    
    description = 'Scatter plot created by BoutrosLab.plotting.general'
    ) + as.layer(
      
      # Add RNASeq ~ time plot as a layer
      create.scatterplot(
        #
        formula = RNASeq ~ time,
        data = rnaseq.tlda.scatter$timed.fold.plot$gene.subset,
        pch = 18,
        cex = 2.2,
        col = colour.pallet$samples[ as.character( rnaseq.tlda.scatter$timed.fold.plot$gene.subset$condition ) ],
        
        description = 'Scatter plot created by BoutrosLab.plotting.general'
        )
      
      );
  
    # Export heatmap
    tiff(
      filename = paste0( "~/Desktop/figure_2/" , c.gene , ".tiff" ), 
      type = "cairo", 
      units = "cm",
      width = 15, 
      height = 12.5, 
      res = 500, 
      compression = 'lzw'
      );
  
    plot( rnaseq.tlda.scatter$timed.fold.plot[[c.gene]] );
  
    dev.off();

  }

# ### Aggregate all plot on a multiplot (manually done with GIMP)
# tiff(
#   filename = "~/Desktop/figure2.tiff",
#   type = "cairo",
#   units = "cm",
#   width = 35,
#   height = 45,
#   res = 500,
#   compression = 'lzw'
#   );
# 
# # Could not find a better way to insert spacing between plots
# blank.panel <- grid.rect( gp = gpar( col="white") );
# 
# # Create merged images
# grid.arrange(
#   blank.panel,
#   pca.plot$final.plot,
#   fc.heatmap$plot,
#   rnaseq.tlda.scatter$scatter.plot,
#   rnaseq.tlda.scatter$timed.fold.plot$Vegfa,
#   rnaseq.tlda.scatter$timed.fold.plot$Edn2,
#   layout_matrix = matrix( c(1,1,3,2,1,3,1,1,3,4,1,3,1,1,1,5,1,6), nrow = 6, ncol = 3, byrow = TRUE),
#   widths = c(0.45, 0.1, 0.45),
#   heights = c(0.05, 0.2, 0.05, 0.3, 0.05, 0.25)
#   );
# 
# dev.off();

### Save Rdata object for later use
save.image(file = "~/Desktop/figure2_expression/raw_individual/session.Rdata");
load("~/Desktop/figure2_expression/raw_individual/session.Rdata");
