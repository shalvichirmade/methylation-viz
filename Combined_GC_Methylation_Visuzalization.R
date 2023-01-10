#### Visualize Methylation - SLC9A3
## HG002 CpG methylation annotated file used to produce analysis of CpG/5mC data using pg-cpg-tools.

## Shalvi Chirmade
## December 1, 2022

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly = T)


#### 1. Load required libraries ----

# if (!require(package)) install.packages('package')
# library(package)

#install.packages(c("ggbio", "magrittr", "tidyverse"))
#BiocManager::install(c("ggbio", "GenomicRanges", "GenomicAlignments", "Rsamtools", "TxDb.Hsapiens.UCSC.hg38.knownGene", "karyoploteR", "biomaRt", "regioneR", "Homo.sapiens"))

library(biomaRt)
library(GenomicAlignments)
library(GenomicRanges)
library(ggbio)
library(Homo.sapiens)
library(karyoploteR)
library(magrittr)
library(regioneR)
library(Rsamtools)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


#### 2. Queries to Change ----

# Change directory location.
setwd("~/Documents/pb-cpg-tools & Visualizations/Script/samtools_pb_cpg_tools_files/")

# Change name of BED file in directory. This comes from the output of pb-cpg-tools.
bed <- "cpg_chr7_140686707_140716643.combined.reference.bed"

# Chromosome for visualization.
chromosome <- "chr7"


# Can run the rest of the script without edits.

#### 2. Input BED file ----

dfBED <- read.table(bed,
                    header = F,
                    sep = "\t",
                    stringsAsFactors = F,
                    quote = "")

# Explicitly name the columns
colnames(dfBED) <- c("Chromosome",
                     "Start_Coordinate",
                     "End_Coordinate",
                     "Modification_Probability",
                     "Haplotype",
                     "Coverage",
                     "Estimated_Mod_Site_Count",
                     "Estimated_Unmod_Site_Count",
                     "Discretized_Modification_Probability")

# Find minimum and maximum coordinates.
min <- min(dfBED$Start_Coordinate)
max <- max(dfBED$End_Coordinate)

# Start coordinate for visualization.
start_coord <- format(min, scientific = F)

# End coordinate for visualization.
end_coord <- format(max, scientific = F)


#### 3. Create visualization ----

## Using karyoploteR, a multi-tiered visualization is created to showcase the genomic region in study along with the genes located, methylation probability and GC locations.

# Divide prob by 100 to have values from 0 to 1
dfBED$Discretized_Modification_Probability_100 <- (dfBED$Discretized_Modification_Probability) / 100


## Reduce window size
# dfBED %<>% filter((Start_Coordinate >= start_coord) & 
#                     (End_Coordinate <= end_coord))
#TODO need to figure out the threshold of points that this Loess fucntion can handle

gr <- makeGRangesFromDataFrame(dfBED[,c(1,2,3,10)],
                                keep.extra.columns = T,
                                ignore.strand = T,
                                seqinfo = NULL,
                                seqnames.field = "Chromosome",
                                start.field = "Start_Coordinate",
                                end.field = "End_Coordinate")


# Invisible color https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# Invisible black for points on confident interval on methylation bar graph
color <- t_col("black", percent = 99)


# Using makeGenes

#TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# - - - add
#TxDb <- makeTxDbFromUCSC("hg38", "ncbiRefSeqSelect")
TxDb <- makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.0.refseq_genomic.gff.gz")


plot.params <- getDefaultPlotParams(plot.type = 1)
plot.params$ideogramheight <- 5

{
kp <- plotKaryotype(chromosomes = chromosome,
                      zoom = paste0(chromosome,
                                    ":", 
                                    start_coord,
                                    "-", 
                                    end_coord),
                      plot.type = 1,
                      plot.params = plot.params)
  
  
  kpAddBaseNumbers(kp, tick.dist = 10000, 
                   tick.len = 10, 
                   tick.col = "red", 
                   cex = 0.5,
                   minor.tick.dist = 1000, 
                   minor.tick.len = 5, 
                   minor.tick.col = "gray")
  
  genes.data <- makeGenesDataFromTxDb(TxDb,
                                      karyoplot = kp,
                                      plot.transcripts = T,
                                      plot.transcripts.structure = T)
  #genes.data <- addGeneNames(genes.data)
  #genes.data <- mergeTranscripts(genes.data)
  
  kpPlotGenes(kp, data = genes.data,
              data.panel = 2,
              r0 = 0, 
              r1 = 0.2,
              coding.exons.col = "turquoise3",
              coding.exons.border.col = "turquoise3",
              non.coding.exons.col = "coral",
              non.coding.exons.border.col	= "coral",
              gene.name.cex = 0.5,
              gene.name.position = "top")
  
  kpAddLabels(kp,
              labels = "Genes",
              r0 = 0, 
              r1 = 0.2,
              cex = 0.8,
              label.margin = 0.035)
  
  kpAxis(kp,
         side = 1,
         r0 = 0.3, 
         r1 = 0.7,
         tick.pos = c(0, 25, 50, 75, 100),
         labels = c(0.0, 0.25, 0.5, 0.75, 1.0),
         cex = 0.5,
         ymin = 0,
         ymax = 100)
  
  kpBars(kp, chr = chromosome,
         r0 = 0.3, 
         r1 = 0.7,
         x0 = start(gr),
         x1 = end(gr),
         y1 = gr$Discretized_Modification_Probability_100,
         border = "lightslategrey",
         clipping = T)
  
  kpPoints(kp, chr = chromosome,
           data = gr,
           r0 = 0.3, 
           r1 = 0.7,
           y = gr$Discretized_Modification_Probability_100,
           col = color,
           pch = 20)
  
  kpPlotLoess(kp, chr = chromosome,
              data = gr,
              r0 = 0.3, 
              r1 = 0.7,
              y = gr$Discretized_Modification_Probability_100,
              col = "coral",
              lwd = 2,
              ci.col = color,
              ci.border = color,
              span = 0.3)
  
  kpAddLabels(kp,
              labels = "Methylation\nProbability",
              r0 = 0.3, 
              r1 = 0.7,
              cex = 0.8,
              label.margin = 0.035)
  
  kpPlotRegions(kp,
                data = gr, col = alpha("aquamarine3", 0.8),
                r0 = 0.8,
                r1 = 1)
  
  kpAddLabels(kp,
              labels = "GC\nLocation",
              r0 = 0.8,
              r1 = 1,
              cex = 0.8,
              label.margin = 0.035)

legend("bottomright",
       legend = c("Each GC Location",
                  "Per GC Region",
                  "Smooth Spline",
                  "Gene Coding Region",
                  "Gene Non-coding Region"),
       col = c("aquamarine3",
               "lightslategrey",
               "coral",
               "turquoise3",
               "coral"),
       pch = c(15, 15, NA, 15, 15),
       lty = c(NA, NA, 1, NA, NA),
       bty = "n",
       pt.cex = 1,
       cex = 0.5,
       text.col = "black",
       horiz = F,
       inset = c(-0.1, -0.05))
}
