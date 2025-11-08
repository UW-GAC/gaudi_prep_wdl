version 1.0 

import "https://raw.githubusercontent.com/frankp-0/HAUDI_workflow/refs/heads/main/vcf_to_plink2.wdl" as vcf_to_plink2
import "https://raw.githubusercontent.com/frankp-0/HAUDI_workflow/refs/heads/main/convert_lanc.wdl" as convert_lanc
import "https://raw.githubusercontent.com/frankp-0/HAUDI_workflow/refs/heads/main/make_fbm.wdl" as make_fbm
import "run_flare.wdl" as run_flare

workflow gaudi_prep {
    input {
        # vcf to plink input 
        Array[File] vcf_files # also use as target_file_list in FLARE
        File? samples_keep # also use as gt_samples in FLARE

        # convert lanc input

        # FLARE input 
        Array[File] ref_file_list
        Array[String] out_prefix_list
        File genetic_map_file
        File reference_map_file
    }

    call vcf_to_plink2.vcf_to_plink2{
        input:
            vcf_files = vcf_files,
            samples_keep = samples_keep

        # Output files: pgen, pvar, psam
    }

    call run_flare.run_flare {
        input:
            ref_file_list = ref_file_list,
            target_file_list  = vcf_files,
            out_prefix_list = out_prefix_list,
            genetic_map_file = genetic_map_file,
            reference_map_file = reference_map_file, 
            gt_samples = samples_keep
        # Output fils: log_array, model_array, anc_vcf_array, global_anc_array
    }

    call convert_lanc.convert_lanc{
        input: 
            ancestry_files = run_flare.anc_vcf_array,
            ancestry_file_fmt = "FLARE",
            pgen_files = vcf_to_plink2.pgen,
            pvar_files = vcf_to_plink2.pvar,
            psam_files = vcf_to_plink2.psam

        # Output files: lanc_files
    }

    output {
        Array[File] lanc_files = convert_lanc.lanc_files
        Array[File] pgen = vcf_to_plink2.pgen 
        Array[File] pvar = vcf_to_plink2.pvar
        Array[File] psam = vcf_to_plink2.psam
    }
}