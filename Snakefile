import os

from curate_gwas_file import curate_gwas_file
from extract_window import extract_window

configfile: "config.yaml"


input_gwas_file = "region_12_chr1:42763250-45199832.gwas" #"region_2_chr1:4168331-5168331.gwas"#"region_36_chr1:165441996-166819395.gwas"
# NOTE: when changing the file name, remember to change the filters specified for curate_gwas_file() further down !

# If you want to narrow the genomic window, you can use these variables. Otherwise, set them to None
window_start = 44_600_000  # None
window_end = 44_700_000  # None

global_gwasfile_path = os.path.join(config["gwas_dir"], input_gwas_file)

# Check that the input file exists
if not os.path.exists(global_gwasfile_path):
    raise FileNotFoundError(f"File {global_gwasfile_path} does not exist!")

run_name = os.path.splitext(input_gwas_file)[0]  # Same as input file name without the extension

# If no directory 'output' exists, create it
if not os.path.exists("output"):
    os.makedirs("output")

# Create a run-specific directory in 'output'
if not os.path.exists(f"output/{run_name}"):
    os.makedirs(f"output/{run_name}")

run_dir = f"output/{run_name}"

# Curate .gwas file before the Snakemake rules are applied. This uses 'curate_gwas_file()'
curated_gwas_file = f"{run_dir}/{os.path.splitext(input_gwas_file)[0]}_curated.gwas"
argument_str = f"{global_gwasfile_path} "\
               f"-r {run_dir} "\
               f"-s {window_start} "\
               f"-e {window_end} "\
               f"-o {curated_gwas_file}"

argument_list = argument_str.split()
curate_gwas_file(argument_list)

# Extract the genomic window from the curated .gwas file. This uses 'extract_window()'
window_file = f"{run_dir}/{os.path.splitext(input_gwas_file)[0]}_window_file.tsv"
argument_list = f"{curated_gwas_file} -o {window_file}".split()

extract_window(argument_list)

# Read the window file and set the environment variables
with open(window_file, "r") as f:
    line = f.readline()
    chromosome, start, end = line.split("\t")
    start = int(start)
    end = int(end)

os.environ["chromosome"] = chromosome
os.environ["start"] = str(start)
os.environ["end"] = str(end)
os.environ["run_dir"] = run_dir

print(f"Run for file {input_gwas_file} is underway!\nGenomic Window: chr{chromosome}:{start}-{end}")

########################################################################################################################
# Pipeline starts here

rule all:
    input:
        f"{run_dir}/{run_name}_plot.png"


rule create_HiC_file:
    params:
        run_dir = os.environ["run_dir"],
        chromosome = os.environ["chromosome"],
        start = os.environ["start"],
        end = os.environ["end"]
    output:
        f"{os.environ['run_dir']}/chr{os.environ['chromosome']}_{os.environ['start']}_{os.environ['end']}_HiC.tsv"
    shell:
        "bash /home/antton/Resources/Mifsud_Hi-C/extract_region_HiC.sh {params.chromosome} {params.start} {params.end} {params.run_dir}"


rule create_ATACseq_file:
    params:
        run_dir = os.environ["run_dir"],
        chromosome = os.environ["chromosome"],
        start = os.environ["start"],
        end = os.environ["end"]
    output:
        f"{os.environ['run_dir']}/chr{os.environ['chromosome']}_{os.environ['start']}_{os.environ['end']}_ATAC.tsv"
    shell:
        "mamba run -n integrative-scrna-scatac-human-foetal_env python /home/antton/Resources/Ulirsch_PBMC_ATAC-seq/sliding_window_on_bigwig_files.py -c {params.chromosome} -s {params.start} -e {params.end} -o {output}"


rule create_cell_types_file:
    output:
        f"{run_dir}/cell_types.txt"
    run:
        all_cell_types_list = ["HSC", "MPP", "CMP", "Ery", "MEP", "Mega", "GMP.A", "GMP.B", "GMP.C", "Mono", "mDC",
                               "LMPP", "CLP", "NK", "CD4", "CD8", "Bcell", "pDC", "Plasma"]
        selected_cell_types_list = ["HSC", "NK", "CD4", "CD8", "Bcell", "pDC", "Plasma"]
        with open(output[0],"w") as f:
            for cell_type in selected_cell_types_list:
                f.write(f"{cell_type}\n")


rule draw_plot:
    input:
        gwas_file = curated_gwas_file,
        hic_file = f"{run_dir}/chr{chromosome}_{start}_{end}_HiC.tsv",
        atac_file = f"{run_dir}/chr{chromosome}_{start}_{end}_ATAC.tsv",
        celltypes_file = f"{run_dir}/cell_types.txt"
    output:
        f"{run_dir}/{run_name}_plot.png"
    shell:
        "mamba run -n tidygenomebrowser_env Rscript plot_tracks.R -g {input.gwas_file} -c {input.hic_file} -a {input.atac_file} -t {input.celltypes_file} -p 10000 -o {output}"
