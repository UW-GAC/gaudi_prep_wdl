version 1.0

workflow run_flare {
  input{
      Array[File] ref_file_list
      Array[File] target_file_list
    Array[String] out_prefix_list
             File genetic_map_file
             File reference_map_file
            File? gt_samples
  }
    
  call flare as flare1 {
    input:
      ref = ref_file_list[0],
      gt  = target_file_list[0],
      out = out_prefix_list[0],
      map = genetic_map_file,
      ref_panel = reference_map_file,
      gt_samples = gt_samples
  }

  scatter(i in range(length(ref_file_list))) {
    if(i!=0) {
      call flare as flares {
        input: 
          ref = ref_file_list[i],
          gt  = target_file_list[i],
          out = out_prefix_list[i],
          map = genetic_map_file,
          ref_panel = reference_map_file,
          gt_samples = gt_samples,
          em = false,
          model = flare1.model
      }
    }
  }

  output {
    Array[File] log_array        = select_all(flatten([[flare1.log    ], flares.log    ]))
    Array[File] model_array      = select_all(flatten([[flare1.model  ], flares.model  ]))
    Array[File] anc_vcf_array    = select_all(flatten([[flare1.anc_vcf], flares.anc_vcf]))
    Array[File] global_anc_array = select_all(flatten([[flare1.global_anc], flares.global_anc]))
  }

}

task flare {
  input {
    # Required inputs
     String  out
       File  ref
       File  gt
       File  map
       File  ref_panel

    # Optional inputs (FLARE's defaults restated for clarity)
    Boolean  em      = true
    Boolean  array   = false
    Boolean  probs   = false
      Float  min_maf = 0.005
        Int  min_mac = 50
        Int  gen     = 10
        Int  seed    = -99999
        Int? nthreads
       File? model
       File? gt_samples
       File? gt_ancestries
       File? excludemarkers

    # Runtime specs
    Int gb_disk = 20
    Int gb_mem = 10
    Int n_cpu = 1
    Int preemptible = 0
  }

  command <<< 
    java ~{"-Xmx"+gb_mem+"G"} -jar /flare.jar \
      ~{"out="            + out           } \
      ~{"ref="            + ref           } \
      ~{"gt="             + gt            } \
      ~{"map="            + map           } \
      ~{"ref-panel="      + ref_panel     } \
        "em=~{em}"                          \
        "array=~{array}"                    \
        "probs=~{probs}"                    \
      ~{"min-maf="        + min_maf       } \
      ~{"min-mac="        + min_mac       } \
      ~{"gen="            + gen           } \
      ~{"seed="           + seed          } \
      ~{"model="          + model         } \
      ~{"gt-samples="     + gt_samples    } \
      ~{"gt-ancestries="  + gt_ancestries } \
      ~{"excludemarkers=" + excludemarkers} \
      ~{"nthreads=" + nthreads}
  >>>

  output {
    File log        = "${out}.log"
    File model      = "${out}.model"
    File anc_vcf    = "${out}.anc.vcf.gz"
    File global_anc = "${out}.global.anc.gz"
  }

  runtime {
         docker: "bdchen/run_flare:0.0.2"
          disks: "local-disk ${gb_disk} HDD"
         memory: "${gb_mem} GB"
            cpu: "${n_cpu}"
    preemptible: "${preemptible}"
  }
}