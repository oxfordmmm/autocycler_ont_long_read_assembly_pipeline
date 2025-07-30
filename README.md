# Autocycler ONT long-read assembly pipeline
A Nextflow pipeline for generating reference-grade consensus assemblies of clinical Enterobacterales isolates using Oxford Nanopore Technologies (ONT) long-read DNA sequences and the consensus assembler Autocycler, with downstream quality control (CheckM2, seqkit stats) and annotation (MLST, Kraken2, Bakta, AMRFinder PLUS, MOB-suite).

# Databases Download
Databases for 
- Hybracter
- CheckM2
- Kraken2
- AMRFinder [PLUS]
- MOB-suite
need to be downloaded separately, and specified using the `--databasesDir` option. Please refer to each tool's own documentation on instruction on how to download the required databases.

# Usage
``` 
nextflow run main.nf --inputFastq_ONT <input_long_reads_dir/> \
--inputmodel <r1041_e82_400bps_bacterial_methylation> \
--manual_assemblies <dir_of_manual_assemblies_with_names_matching_reads> \
--databasesDir <path_to_downloaded_databases> \
--projDir <path_to_project_directory> \
--conda.cacheDir <dir_to_conda_cache_dir> \
--work-dir <path_to_dir_to_store_work_directory> \
-profile conda \
-resume \
-with-dag \ #optional extra run reports
-with-timeline \ #optional extra run reports
-with-report \ #optional extra run reports
-with-trace  #optional extra run reports
```
