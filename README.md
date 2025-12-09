# gaudi_prep_wdl

Workflow to prepare inputs for GAUDI WDL

## gaudi_prep

This workflow is used to create input files to run make_fbm and gaudi. vcf_to_plink and FLARE are run on the same VCF files, and the outputs are used to run convert_lanc. 

Inputs:

input | description
--- | ---
vcf_files | VCF Array to both convert to PLINK2 pgen format and use as target VCF files for FLARE
samples_keep | Optional file with IDs to keep in both VCF to PLINK and FLARE
ref_file_list | VCF Array of reference VCF files for FLARE
out_prefix_list | Required FLARE input, string array of output prefixes
genetic_map_file | Required FLARE input file
reference_map_file | Required FLARE input file