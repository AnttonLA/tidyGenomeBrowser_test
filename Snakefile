
configfile: "config.yaml"

input_bed_file = "region_12_chr1:42763250-45199832.bed"

rule all:
    input:
        "/home/antton/Desktop/tGB_test_6.png"

rule curate_bed:
    input:
        bed_file = input_bed_file
    output:
        bed_file = "curated.bed"
    shell:
        ""

rule create_HiC_file:
    # TODO: edit so that chr & pos are taken from bed file
    params:
        chromosome = "chr1",
        start = 44612259,
        end = 44719832
    output:
        "output/chr1_44612259_44719832.txt"
    shell:
        ". /home/antton/Resources/Mifsud_Hi-C/extract_region_HiC.sh {params.chromosome} {params.start} {params.end} {output}"

# TODO: generate ATAC-seq file

rule plot_region:
    output:
        "/home/antton/Desktop/tGB_test_6.png"
    # TODO: edit so bed file & arguments are dynamic
    shell:
        "Rscript plot_tracks.R -b data/GWAS/region_12_chr1:42763250-45199832.bed -c output/chr1_44612259_44719832.txt \
        -a output/chr1:44612259-44719832_ATAC.tsv -o {output}"
