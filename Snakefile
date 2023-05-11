import os

configfile: "config.yaml"

input_gwas_file = "region_12_chr1:42763250-45199832.gwas"
gwas_file_path = os.path.join(config["gwas_dir"], input_gwas_file)

# If no directory 'resources' exists, create it
if not os.path.exists("resources"):
    os.makedirs("resources")

########################################################################################################################
# Pipeline starts here

rule all:
    input:
        f"{os.path.splitext(gwas_file_path)[0]}_curated.gwas"

rule curate_gwas_file:
    input:
        gwas_file_path
    params:
        start = None,
        end = None,
        pval_threshold = 5e-8,
        phenotype_file = f"resources/{input_gwas_file}_phenotype_list.txt",
        phenotype_map_file = f"resources/{input_gwas_file}_phenotype_map.tsv",
    output:
        f"{os.path.splitext(gwas_file_path)[0]}_curated.gwas"
    shell:
        "python curate_gwas_file.py {input} -r resources -o {output}"

rule create_HiC_file:
    # TODO: edit so that chr & pos are taken from gwas file
    input:
        f"{os.path.splitext(gwas_file_path)[0]}_curated.gwas"
    params:
        chromosome = f"{os.path.splitext(gwas_file_path)[0]}_curated.gwas".split("_")[1],
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
