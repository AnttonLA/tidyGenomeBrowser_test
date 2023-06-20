#!/usr/bin/env Rscript

library(optparse)

# This script uses the tidyGenomeBrowser package to plot genomic regions. It is mainly intended to be used to visualize
# genetic variants discovered by GWAS. The script takes a .gwas file containing genomic positions as input, and plots the
# variants together with additional data (genes, PCHi-C, ATAC-seq). The script is intended to be run from the command line.

option_list <- list(
  make_option(c("-g", "--gwas"), type="character", default=NULL, help=".gwas file containing the genomic positions to be plotted"),
  make_option(c("-p", "--padding"), type="integer", default=100000, help="Padding around the genomic region to be plotted (default: 100000)"),
  make_option(c("-c", "--hic"), type="character", default=NULL, help="Hi-C data file name"),
  make_option(c("-a", "--atac"), type="character", default=NULL, help="ATAC-seq data file name"),
  make_option(c("-t", "--celltypes"), type="character", default=NULL, help="Cell types to be plotted"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file name")
)

help_txt <- "This script uses the tidyGenomeBrowser package to plot genomic regions. It is mainly intended to be used to visualize
genetic variants discovered by GWAS. The script takes a .gwas file containing genomic positions as input, and plots the
variants together with additional data (genes, PCHi-C, ATAC-seq). The script is intended to be run from the command line."
opt_parser <- OptionParser(option_list=option_list, description=help_txt)
args <- parse_args(opt_parser)

if (is.null(args$gwas)) {
  stop(".gwas file not specified! Use -g or --gwas to specify the .gwas file.")
}

if (is.null(args$output)) {
  stop("Output file name not specified! Use -o or --output to specify the output file name.")
}

if (!is.null(args$atac) & is.null(args$celltypes)) {
  print("ATAC-seq file specified, but no cell types given. Plotting all cell types.")
}

library(GenomicRanges)
library(rtracklayer)
library(tidyGenomeBrowser)
library(magrittr)
library(tidyverse)
library(ggplot2)

########################################################################################################################
# Plotting
# Set theme for the plotting. See all themes: https://ggplot2.tidyverse.org/reference/ggtheme.html
theme_set(theme_classic())


compose_and_save_plot <- function(args, plot_stack, output_file) {
  # This function will assing heights to the tracks in the plot, as well as the size of the final plot itself.
  # It will create a 'heights_list' variable based on the number of tracks in the plot_stack.

  if (length(plot_stack) == 2) {  # Bare GWAS + genes plot
    heights_list <- c(3,3)
  } else if (length(plot_stack) == 3) {
    if (is.null(args$atac)) {  # GWAS + genes + Hi-C plot
      heights_list <- c(3,3,3)
    } else if (is.null(args$hic)) {  # GWAS + genes + ATAC-seq plot
      heights_list <- c(3,3,12)
    } else { stop("Something went wrong with the plotting. 3 tracks requested but non-null Hi-C and ATAC-seq files!") }
  } else if (length(plot_stack) == 4) {  # Full plot
    heights_list <- c(3,3,3,12)
  } else { stop("Something went wrong with the plotting. Invalid number of tracks!") }

  all_tracks <- browseStack(plot_stack, squeeze="internal", heights = heights_list)
  # Save output
  ggsave(all_tracks, filename=output_file, height=10, width=10)
  # Finish script execution
  q()
}

########################################################################################################################
# Load GWAS data from .gwas file
gwas_file <- args$gwas

gwas_gp <- gwas_file |>
    read.table(header=TRUE) |>
    makeGRangesFromDataFrame(start.field = "position",
                             end.field = "position",
                             keep.extra.columns = TRUE) |>
    as("GPos") |>
    transform(score=-log10(pval),
              color=phenotype,
              facet=phenotype)

seqlevelsStyle(gwas_gp) <- "UCSC"

# TODO: this processing is outside the scope of this script. It should get a .gwas file with simplified names already.
# use str_replace_all() with regular expressions to replace multiple patterns in the 'phenotype' strings
# mcols(gwas_gp)$facet <- str_replace_all(mcols(gwas_gp)$facet,
#                                             pattern = c("Bpanel_" = "",
#                                                         "Tpanel_" = "",
#                                                         "_div_" = "/",
#                                                         "_Frequency$" = " Freq.",
#                                                         "classical_monocytes" = "Classic Mono.",
#                                                         "doublePosNK" = "PosPosNK",
#                                                         "pos" = "+",
#                                                         "neg" = "-",
#                                                         "pl" = "&",
#                                                         "_" = " "))

# TODO: include a check for the a max number of phenotypes and a max width of the genomic region


pad <- args$padding
w <- GenomicRanges::reduce(as(gwas_gp, "GRanges"), min.gapwidth=1e6) + pad

########################################################################################################################
# GWAS Track
gwas_track <- browsePositions(gwas_gp, region=w) +

    scale_color_brewer(palette="Set1", guide="none") +
    geom_point(aes(x=5,y=5), size=10) +
    theme(panel.grid.major = element_line(color = "red",
                                          size = 0.05,
                                          linetype = 2)) +
    theme(strip.text.y = element_text(size=8, angle=0),
          strip.background = element_rect(colour=NA, fill="white"),
          axis.title.y=element_text(size=8)) +
    ylab("-log10(PValue)")

########################################################################################################################
# Gene Track
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC gene "knownGene" database
library(org.Hs.eg.db)  # Gene names, used to map IDs to symbols

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_models <- exonsBy(txdb, by="gene")  # function from the TxDb.Hsapiens.UCSC.hg38.knownGene package to extract exon info
seqlevelsStyle(gene_models) <- "UCSC"  # Set the style of the sequence levels in meta_models to "UCSC"

# Change the names of the meta_models to gene symbols
names(gene_models) <- mapIds(org.Hs.eg.db, keys = names(gene_models), keytype = "ENTREZID", column = "SYMBOL")
gene_track <- gene_models |>
    browseTranscripts(region=w) +
    ylab("Genes") +


if (is.null(args$hic) & is.null(args$atac)) {
  # If no other arguments were give, plot only GWAS and gene tracks
    plot_stack <- list(gwas_track,
                     gene_track)
  # Add marker lines across all plots
  plot_stack <- map(plot_stack, ~ .x + geom_vline(xintercept=pos(gwas_gp), linetype="dotted", alpha=0.25))
  compose_and_save_plot(args, plot_stack, args$output)
}

########################################################################################################################
# Hi-C Track
if (!is.null(args$hic)) {
hic_df <- args$hic |>
    read.table(sep="\t", header=TRUE) |>
    # Rename columns because 'start' and 'end' are reserved words
    rename("chr_other"="chr","start_other"="start","end_other"="end","score"="log.observed.expected.")

bait_gr <- GRanges(hic_df$chr_bait, IRanges(hic_df$start_bait, hic_df$end_bait))
other_gr <- GRanges(hic_df$chr_other, IRanges(hic_df$start_other, hic_df$end_other))

hic_gi <- GInteractions(bait_gr, other_gr, score=hic_df$score)
hic_track <- browseInteractions(hic_gi, region=w)
}

if (is.null(args$atac)) {
  # If no ATAC-seq file was given, plot GWAS + gene + Hi-C tracks
    plot_stack <- list(gwas_track,
                       gene_track,
                       hic_track)
  # Add marker lines across all plots
  plot_stack <- map(plot_stack, ~ .x + geom_vline(xintercept=pos(gwas_gp), linetype="dotted", alpha=0.25))
  compose_and_save_plot(args, plot_stack, args$output)
}

########################################################################################################################
# ATAC-seq Track

# Read in atac-seq data
atac_file <- read.csv(args$atac, sep="\t")

# Format as GRanges
atac_gr <- atac_file %>%  # convert data from "wide" to "long" format
    gather(key = "facet", value = "score", -chr, -pos, factor_key=TRUE) |>
    mutate(end=pos+99, color=facet) |>
    makeGRangesFromDataFrame(start.field="pos",
                             keep.extra.columns=TRUE) |>
    sort()

# Give each cell type a color
ATAC_cols <- c(HSC="black",
               MPP="grey",
               CMP="lightcoral",
               Ery="firebrick",
               MEP="coral",
               Mega="red",
               GMP.A="hotpink",
               GMP.B="hotpink2",
               GMP.C="hotpink4",
               Mono="blueviolet",
               mDC="darkmagenta",
               LMPP="cyan",
               CLP="darkturquoise",
               NK="dodgerblue",
               CD4="navy",
               CD8="royalblue",
               Bcell="slateblue",
               pDC="steelblue",
               Plasma="skyblue")

# TODO: change this to use the input instead
if (is.null(args$celltypes)) {
  # If no celltypes were given, plot all of them
  pops_to_plot <- unique(atac_gr$color)
} else {
  # Otherwise, plot only the ones given
# TODO: safety check. What happens if file has celltypes that do not exist?
  pops_to_plot <- readLines(args$celltypes)
}
# Subset to only the populations we want to plot
atac_gr <- atac_gr[atac_gr$color %in% pops_to_plot]
ATAC_cols <- ATAC_cols[pops_to_plot]

# ATAC track
atac_track <- browseSignal(atac_gr, region=w) +
    coord_cartesian(ylim=c(0,0.005), xlim = c(start(w), end(w))) + # This can fixate the height of the tracks
    scale_fill_manual("Cell", values=ATAC_cols, guide="none") +
    #theme_void() +
    theme(strip.text.y = element_text(size=8, angle=0),
          strip.background = element_rect(colour=NA, fill="white"),
          axis.title.y = element_text(size=8),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(0.001, "mm"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major = element_blank()) +
    ylab("Accessability") +
    guides(y = "none")

if (is.null(args$hic)) {
  # If no Hi-C file was given, plot GWAS + gene + ATAC-seq tracks
    plot_stack <- list(gwas_track,
                       gene_track,
                       atac_track)
  # Add marker lines across all plots
  plot_stack <- map(plot_stack, ~ .x + geom_vline(xintercept=pos(gwas_gp), linetype="dotted", alpha=0.25))
  compose_and_save_plot(args, plot_stack, args$output)
} else {
  # Plot full plot
  plot_stack <- list(gwas_track,
                       gene_track,
                       hic_track,
                       atac_track)
  # Add marker lines across all plots
  plot_stack <- map(plot_stack, ~ .x + geom_vline(xintercept=pos(gwas_gp), linetype="dotted", alpha=0.25))
  compose_and_save_plot(args, plot_stack, args$output)
}

stop("Ooops, script end reached!!! This should not have happened...")
