# tidyGenomeBrowser_test

This is a test project using the tidyGenomeBrowser package. The goal is to create a tidyGenomeBrowser plot for a
certain genomic region specified in an IGV .gwas file, and optionally integrate a gene track, a Hi-C track and an
ATAC-seq track into the same figure.

The plotting itself is carried out by an R script 'plot_tracks.R', which can be called from the command line or a bash
script (including the Snakefile).

A Snakemake pipeline ensures that the necessary data is extracted and processed, and some utility scripts allow to have
better control over this process.

## Usage

Before runnning the Snakemake pipeline, make sure that the paths in the 'config.yaml' file are correct.
You should also edit the snakefile itself so that the variables `input_gwas_file` and `argument_str` are set to the
correct values.

The Snakemake pipeline can be run with the following command:

```bash
snakemake --cores=all
```

After running the pipeline once, you might want to tweak the resulting plots. This can be accomplished by editing the
`argument_str` variable in the snakefile, as well as the files in the 'resources/' folder for our .gwas file.

## Project structure
The 'data/' directory contains .gwas files. When running the pipeline, a folder will be created inside 'output/' with
the same name as the .gwas file, and all the intermediary and output files will be generated there.
Inside the output folder, there will be a 'resources/' folder for the files that are expected to be edited by the user. 

## References

The tidyGenomeBrowser package can be found here: 
[https://github.com/MalteThodberg/tidyGenomeBrowser](https://github.com/MalteThodberg/tidyGenomeBrowser)

