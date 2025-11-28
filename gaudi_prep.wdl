version 1.0 

import "https://raw.githubusercontent.com/frankp-0/HAUDI_workflow/refs/heads/main/vcf_to_plink2.wdl" as vcf_to_plink2
import "https://raw.githubusercontent.com/frankp-0/HAUDI_workflow/refs/heads/main/convert_lanc.wdl" as convert_lanc
import "https://raw.githubusercontent.com/bdchen/FLARE_workflow/refs/heads/main/FLARE.wdl" as run_flare

workflow gaudi_prep {
    input {
        # vcf to plink input 
        Array[File] vcf_files # also use as target_file_list in FLARE
        File? samples_keep # also use as gt_samples in FLARE

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
        # Output files: log_array, model_array, anc_vcf_array, global_anc_array
    }

    call combine_flare {
        input:
            flare_files = run_flare.global_anc_array
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
        File merged_global_ancestry = combine_flare.merged_global_ancestry
        File global_ancestry_plot = combine_flare.global_ancestry_plot
        Array[File] lanc_files = convert_lanc.lanc_files
        Array[File] pgen = vcf_to_plink2.pgen 
        Array[File] pvar = vcf_to_plink2.pvar
        Array[File] psam = vcf_to_plink2.psam
    }
}


task combine_flare {
    input {
        Array[File] flare_files
    }

    Int N = length(flare_files)

    command <<<

    for f in ~{sep=' ' flare_files}; do
        echo "$f" >> flare_files.txt
    done

    R << RSCRIPT
    library(tidyverse)
    library(RColorBrewer)

    # get chromosome sizes
    chr_sizes <- read_tsv(
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
        col_names=c('chrom','size')) %>% 
    filter(chrom %in% paste0('chr',1:22)) %>%
    mutate(chr_num = as.integer(sub('chr','',chrom))) %>%
    select(-chrom) %>%
    arrange(chr_num)
    chr_sizes <- chr_sizes[1:N, ]

    # read flare files
    flare_files <- readLines('flare_files.txt')

    df <- tibble(flare_file = flare_files) %>%
      mutate(
        base = basename(flare_file),
        chrom = str_extract(base, '^chr[0-9]+'),
        chr_num = as.integer(str_remove(chrom, 'chr'))
      ) %>%
      inner_join(chr_sizes, by='chr_num') %>%
      arrange(chr_num)
    flare_files <- df[['flare_file']]
    chr_weights <- df[['size']] / sum(df[['size']])

    # combine chromosomes
    combine_chrs <- function(flare_files, chr_weights) {
        tmp <- read_tsv(flare_files[1], show_col_types=FALSE)
        samples <- tmp[['SAMPLE']]
        tmp_numeric <- as.matrix(tmp[, setdiff(names(tmp), "SAMPLE")])
        fracs <- tmp_numeric * chr_weights[1]
        for (i in 2:length(flare_files)) {
            tmp <- read_tsv(flare_files[i], show_col_types=FALSE)
            tmp_numeric <- as.matrix(tmp[, setdiff(names(tmp), "SAMPLE")])
            fracs <- fracs + tmp_numeric * chr_weights[i]
        }
        fracs <- fracs / rowSums(fracs)
        flr <- bind_cols(samples = samples, as_tibble(fracs))
        return(flr)
    }

    flr <- combine_chrs(flare_files, chr_weights)
    write_tsv(flr, 'global_ancestry.tsv')
    # plotting
    flr_long <- flr %>%
      mutate(n=row_number()) %>%
      pivot_longer(-c(samples, n), names_to='Cluster', values_to='Value')
    K <- length(unique(flr_long[['Cluster']]))
    colormap <- setNames(c(brewer.pal(8,'Dark2'), brewer.pal(8,'Set2'))[1:K],
                         unique(flr_long[['Cluster']]))
    p <- ggplot(flr_long, aes(x=n, y=Value, fill=Cluster, color=Cluster)) +
      geom_bar(stat='identity') +
      scale_fill_manual(values=colormap, breaks=rev(names(colormap))) +
      scale_color_manual(values=colormap, breaks=rev(names(colormap))) +
      theme_classic() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    ggsave('global_ancestry.png', p, width=12, height=4)
    RSCRIPT

    >>>

    output {
        File merged_global_ancestry = "global_ancestry.tsv"
        File global_ancestry_plot = "global_ancestry.png"
    }

    runtime {
        docker: "rocker/tidyverse:4.3.1"  # R + tidyverse
        memory: "8G"
        cpu: 1
    }
}