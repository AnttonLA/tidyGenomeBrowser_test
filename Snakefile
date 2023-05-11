import os

from curate_gwas_file import curate_gwas_file
from extract_window import extract_window

configfile: "config.yaml"


input_gwas_file = "region_12_chr1:42763250-45199832.gwas"
global_gwasfile_path = os.path.join(config["gwas_dir"], input_gwas_file)

# If no directory 'resources' exists, create it
if not os.path.exists("resources"):
    os.makedirs("resources")

# If no directory 'output' exists, create it
if not os.path.exists("output"):
    os.makedirs("output")

# Curate .gwas file before the rules are applied
argument_list = f"{global_gwasfile_path} -r resources -o data/intermediate/{os.path.splitext(input_gwas_file)[0]}_curated.gwas".split()
curate_gwas_file(argument_list)

window_file = f"resources/{os.path.splitext(input_gwas_file)[0]}_window_file.tsv"
print("window_file:", window_file)
argument_list = f"data/intermediate/{os.path.splitext(input_gwas_file)[0]}_curated.gwas -o {window_file}".split()
print("argument_list", argument_list)
extract_window(argument_list)

with open(window_file, "r") as f:
    line = f.readline()
    chromosome, start, end = line.split("\t")
    start = int(start)
    end = int(end)

os.environ["chromosome"] = chromosome
os.environ["start"] = str(start)
os.environ["end"] = str(end)

########################################################################################################################
# Pipeline starts here

rule all:
    input:
        f"data/intermediate/{os.path.splitext(input_gwas_file)[0]}_curated.gwas"


rule create_HiC_file:
    input:
        f"resources/{os.path.splitext(input_gwas_file)[0]}_window_file.tsv"
    params:
        chromosome = os.environ["chromosome"],
        start = os.environ["start"],
        end = os.environ["end"]
    output:
        f"resources/{os.environ['chromosome']}_{os.environ['start']}_{os.environ['end']}_HiC.tsv"
    shell:
        "bash /home/antton/Resources/Mifsud_Hi-C/extract_region_HiC.sh {params.chromosome} {params.start} {params.end}"

# TODO: generate ATAC-seq file

rule plot_region:
    output:
        "/home/antton/Desktop/tGB_test_6.png"
    # TODO: edit so bed file & arguments are dynamic
    shell:
        "Rscript plot_tracks.R -b data/GWAS/region_12_chr1:42763250-45199832.bed -c output/chr1_44612259_44719832.txt \
        -a output/chr1:44612259-44719832_ATAC.tsv -o {output}"
