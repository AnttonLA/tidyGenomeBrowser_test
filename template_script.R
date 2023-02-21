library(tidyGenomeBrowser)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)  # Gene names
library(AnnotationDbi)
library(magrittr)
library(tidyverse)

folder = "data/GWAS/"
filename = paste(folder,"region_8_chr1:24370664-25439985.bed", sep="")
bed_file = read.csv(filename, sep="\t")

grt = with(bed_file, GPos(chromosome, pos=position, pval=pval, phenotype=phenotype))

# Extract the necessary information to construct GRanges objects for the plotting.
chr_num <- seqnames(grt)[1]
chr_w_prefix <- paste0("chr", seqnames(grt))[1]
start_pos <- pos(grt)[1]
end_pos <- pos(grt)[length(pos(grt))]
buffer <- 100000  # Extra distance on both ends of the range

plotting_region = GRanges(seqnames=chr_num, ranges=IRanges(start=start_pos, end=end_pos) + buffer)
plotting_region_chr = GRanges(seqnames=chr_w_prefix, ranges=IRanges(start=start_pos, end=end_pos) + buffer)

# Reference gene and transcript models
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_models <- genes(txdb)

# Create a data frame called gene_info that contains the gene ID and gene name information for all the genes in the TxDb object.
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(txdb, keytype = "GENEID"), keytype="ENTREZID", columns=c("SYMBOL", "ENTREZID"))

# Create a vector called gene_names using the mapIds() function from the org.Hs.eg.db package to map the gene IDs to gene names
gene_names <- mapIds(org.Hs.eg.db, keys = gene_info$ENTREZID, keytype = "ENTREZID", column = "SYMBOL")

# Add column 'gene_names' to 'gene_models' GRanges object
mcols(gene_models)$gene_name <- gene_names[match(gene_models$gene_id, gene_info$ENTREZID)]

# Read ATAC-seq file
folder = "data/ATAC_seq/"
filename = paste(folder,"region_8_chr1:24370664-25439985_ATAC.tsv", sep="")
atac_file = read.csv(filename, sep="\t")
atac_grt = with(atac_file, GPos(chr, pos=pos, HSC=HSC, MPP=MPP, LMPP=LMPP, CLP=CLP, CMP=CMP, GMP.A=GMP.A, GMP.B=GMP.B, GMP.C=GMP.C, MEP=MEP, Mono=Mono, mDC=mDC, Ery=Ery, Mega=Mega, CD4=CD4, CD8=CD8, Bcell=Bcell, Plasma=Plasma, NK=NK, pDC=pDC))


# Create tracks for plotting
gwas_track <- grt %>%
    transform(score=-log10(pval), color=phenotype)%>%
    browsePositions(plotting_region)

gene_track <- gene_models %>%
    transform(name=gene_name) %>%
    browseIntervals(region=plotting_region_chr)

# Broken !!
atac_track <- atac_grt %>%
    browsePositions(plotting_region)

plot_list <- list(gwas_track + ylab("GWAS (-log10 P value)"),
                 gene_track + ylab("Genes"),
                 atac_track)

print(browseStack(plot_list))

