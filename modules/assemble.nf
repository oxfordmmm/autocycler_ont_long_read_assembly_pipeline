process RAVEN_ASSEMBLE_INITIAL {
    tag {sample}
    publishDir "assemblies/raven_initial", pattern: "*.fasta", mode: 'copy', saveAs: { filename -> "${sample}_raven_initial.fasta"}
    cpus 4
    
    input:
    tuple val(sample), val(source), path('reads.fastq')

    output:
    tuple val(sample), val(source), path('reads.fastq'), val('raven_assembly_initial'), path("raven_assembly.fasta"), emit: fasta
   

    script: 

    """
    raven --threads $task.cpus --disable-checkpoints reads.fastq > raven_assembly.fasta
    """


}

process GENOME_SIZE_AND_CONTAMINATION {

    tag {sample}
    publishDir "checkm2_initial/", pattern: "checkm2_out/quality_report.tsv", mode: 'copy', saveAs: {filename -> "${sample}_checkm2_initial.tsv"}
    
    cpus 4

    input:
    tuple val(sample), val(source), path('reads.fastq'), val('raven_assembly_initial'), path("raven_assembly.fasta")

    output:
    tuple val(sample), val(source), path('reads.fastq'), val('raven_assembly_initial'), path("raven_assembly.fasta"),  path('checkm2_out/quality_report.tsv'), path("contamination.txt"), path("genome_size.txt"), emit: data
    tuple val(sample), val(source), path('checkm2_out/quality_report.tsv'), emit: tsv

    script:

    """
    checkm2 predict --threads $task.cpus --input raven_assembly.fasta --output-directory checkm2_out/ --stdout -x "fasta" --database_path  "${params.databasesDir}/bin/checkm2_db/uniref100.KO.1.dmnd"
   

    # Extract Genome_Size (9th column) and Contamination (3rd column)
    contamination=\$(awk 'NR==2 {print \$3}' checkm2_out/quality_report.tsv)
    genome_size=\$(awk 'NR==2 {print \$16}' checkm2_out/quality_report.tsv)

    echo "\$contamination" > contamination.txt
    echo "\$genome_size" > genome_size.txt

    """
}



process AUTOCYCLER_SUBSAMPLE {
    tag {sample}
    cpus 4
    publishDir "subsamples", pattern: "subsampled_reads", mode: "copy", saveAs: { filename -> "${sample}_subsampled_reads" }
    
    input:
    tuple val(sample), val(source), path('reads.fastq'), val('raven_assembly_initial'), path("raven_assembly.fasta"),  path('checkm2_out/quality_report.tsv'), path("contamination.txt"), path("genome_size.txt")
    

    output:
    tuple val(sample), val('subsampled_ont'), path("subsampled_reads"), emit: fol
    //tuple val(sample), val('subsampled_ont'), path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), path("subsampled_reads/sample_02.fastq"), path("subsampled_reads/sample_03.fastq"), path("subsampled_reads/sample_04.fastq"), emit: fqs
    tuple val(sample), val('subsampled_ont'), val("01"), path("genome_size.txt"), path("subsampled_reads/sample_01.fastq"), emit: fq01
    tuple val(sample), val('subsampled_ont'), val("02"), path("genome_size.txt"), path("subsampled_reads/sample_02.fastq"), emit: fq02
    tuple val(sample), val('subsampled_ont'), val("03"), path("genome_size.txt"), path("subsampled_reads/sample_03.fastq"), emit: fq03
    tuple val(sample), val('subsampled_ont'), val("04"), path("genome_size.txt"), path("subsampled_reads/sample_04.fastq"), emit: fq04
   

    script:
    
    """
    genome_size=\$(<genome_size.txt)
        
    autocycler subsample --reads reads.fastq --out_dir subsampled_reads --genome_size \${genome_size}
    """


}


process RAWQC_ONT{

    tag {sample}
    publishDir "raw_QC_CSVs/${source}", mode: 'copy'

    input:
    tuple val(sample), val(source), path('reads.fastq.gz')

    stageInMode 'copy'

    output:
    tuple val(sample), val(source), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all reads.fastq.gz > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
    }


process PUBLISH_CONTAMINATED {
    publishDir "contamination", pattern: "*raven_assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_assembly.fasta"}
    publishDir "contamination", pattern: "*quality_report.tsv", mode: "copy", saveAs: { filename -> "${sample}_quality_report.tsv"}
 
    input:
    tuple val(sample), val(source), path('reads.fastq'), val('raven_assembly_initial'), path("raven_assembly.fasta"),  path('quality_report.tsv'), path("contamination.txt"), path("genome_size.txt")
    
    output:
    path("${sample}_raven_assembly.fasta")
    path("${sample}_quality_report.tsv")

    script:
    """
    cp raven_assembly.fasta   ${sample}_raven_assembly.fasta
    cp quality_report.tsv ${sample}_quality_report.tsv
    """
}




process RAWQC_ONT_SUBSAMPLED {

    tag {sample}
    publishDir "raw_QC_CSVs/${source}_${subsample_number}", mode: 'copy'

    input:
    tuple val(sample), val('source'), val('subsample_number'), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val(source), path("${sample}_${subsample_number}.tsv")

    script:
    """
    seqkit stats -T --all subsampled_reads.fastq > ${sample}_${subsample_number}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
    }




process AUTOCYCLER_ASSEMBLE_CANU {
    tag {sample + ' ' + number}
    cpus 16

    //publishDir "autocycler_assemblies/canu/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/canu_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/canu_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/canu_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}





    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('canu'), val(number), path("canu_${number}.fasta"), emit: fa

    script:
    """
    genome_size=\$(<genome_size.txt)

    ${projectDir}/bin/autocycler_scripts/canu.sh subsampled_reads.fastq canu_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_FLYE {
    tag {sample + ' ' + number}
    cpus 8

    //publishDir "autocycler_assemblies/flye/${sample}", mode: "copy"
    publishDir "assemblies/flye_01", pattern: "*01/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_02", pattern: "*02/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_03", pattern: "*03/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_04", pattern: "*04/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/flye_01", pattern: "*01/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_02", pattern: "*02/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_03", pattern: "*03/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/flye_04", pattern: "*04/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('flye'), val(number), path("flye_out_${number}/assembly.fasta"), emit: fa
    tuple val(sample), val('flye'), val(number), path("flye_out_${number}/assembly_graph.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)
    
    flye --nano-hq subsampled_reads.fastq --threads "${task.cpus}" --out-dir flye_out_${number}

    ##${projectDir}/bin/autocycler_scripts/flye.sh subsampled_reads.fastq flye_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_MINIASM {
    tag {sample + ' ' + number}
    cpus 8

    //publishDir "autocycler_assemblies/miniasm/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/miniasm_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/miniasm_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/miniasm_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.fasta"), emit: fa
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.gfa"), optional: true, emit: gfa
    script:
    """
    genome_size=\$(<genome_size.txt)

    ${projectDir}/bin/autocycler_scripts/miniasm.sh subsampled_reads.fastq miniasm_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_RAVEN {
    tag {sample + ' ' + number}
    cpus 8

    //publishDir "autocycler_assemblies/raven/${sample}", pattern: "*.fasta", mode: "copy"
    publishDir "assemblies/raven_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta"}
    publishDir "assemblies/raven_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}
    publishDir "assemblies/raven_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('raven'), val(number), path("raven_${number}.fasta"), emit: fa
    tuple val(sample), val('raven'), val(number), path("raven_${number}.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)

    ${projectDir}/bin/autocycler_scripts/raven.sh subsampled_reads.fastq raven_${number} "${task.cpus}" "\$genome_size"
    """
}



process HYBRACTER_LONG_INDIVIDUAL {
    scratch true
    errorStrategy 'retry'
    cpus 4
    tag {sample + ' ' + number}

    publishDir "assemblies/hybracter_long_${number}", pattern: "hybracter_long_out*", mode: 'copy', saveAs: { filename -> "${sample}_hybracter_out" }
    publishDir "assemblies/hybracter_long_${number}", pattern: "hybracter_long_${number}_final.fasta", mode: "copy", saveAs: { filename -> "${sample}.fasta" }  

    input:
    tuple val(sample), val(source), val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_${number}_final.fasta"), emit: fasta
    tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out"), emit: fol
    //tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_${number}_plasmid.fasta"), optional: true, emit: plasmids_candidate
    //tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/incomplete/${sample}_${number}_final.fasta"), optional: true, emit: plasmids_fallback
    //tuple val(sample), val('hybracter_long'), val(number), path("hybracter_long_out/FINAL_OUTPUT/{complete,incomplete}/${sample}_${number}_chromosome.fasta"), optional: true, emit: chromosomes
 
    script:
    """
    hybracter long-single -l subsampled_reads.fastq -s hybracter_long_${number} -c 4000000 -o hybracter_long_out -t ${task.cpus} --databases "${params.databasesDir}/bin/hybracter_databases" #--medakaModel ${params.inputmodel}
    
    # Safely copy the final output fasta if it exists
    final_fasta=\$(find hybracter_long_out/FINAL_OUTPUT -name 'hybracter_long_${number}_final.fasta' | head -n 1)
    if [ -f "\$final_fasta" ]; then
        cp "\$final_fasta" "hybracter_long_${number}_final.fasta"
    else
        echo "WARNING: Final fasta not found for ${sample}_${number}" >&2
        touch "hybracter_long_${number}_final.fasta"  # optional placeholder
    fi


    """
}


process CHECKM2_INPUT_ASSEMBLIES {
    tag {sample}

    cpus 32

    publishDir "assembly_QC/input_assemblies", pattern: "*.tsv", mode: 'copy', saveAs: { filename -> "${sample}_input_checkm2.tsv" }
    


    input:
    tuple val(sample), path('assembly_01.fasta'), path('assembly_02.fasta'), path('assembly_03.fasta'), path('assembly_04.fasta'), path('assembly_05.fasta'), path('assembly_06.fasta'), path('assembly_07.fasta'), path('assembly_08.fasta'), path('assembly_09.fasta'), path('assembly_10.fasta'), path('assembly_11.fasta'), path('assembly_12.fasta'), path('assembly_13.fasta'), path('assembly_14.fasta'), path('assembly_15.fasta'), path('assembly_16.fasta'),  path('assembly_17.fasta'), path('assembly_18.fasta'), path('assembly_19.fasta'), path('assembly_20.fasta')
    
    output:
    tuple val(sample), val("input_assemblies"), path("checkm2_input_assemblies_out"), emit: fol
    tuple val(sample), val("autocycler"), path("quality_report.tsv"), emit: tsv 
    


    script:
    """
    # Create the output directory for renamed assemblies
    mkdir -p input_assemblies

    # Copy each symlinked input assembly to the folder using the desired naming logic
    for assembly_file in assembly_*.fasta; do
        real_path=\$(readlink -f "\$assembly_file")
        base_name=\$(basename "\$real_path")

        if [[ "\$base_name" == "assembly.fasta" ]]; then
            # Get parent directory name (e.g., flye_out_01), strip _out if present
            parent_dir=\$(basename \$(dirname "\$real_path"))
            prefix=\$(echo "\$parent_dir" | sed 's/_out//')
            new_name="\${prefix}.fasta"
        else
            new_name="\$base_name"
        fi

        cp "\$assembly_file" "input_assemblies/\$new_name"
    done

    # Run checkm2 on all assemblies
    checkm2 predict -t 32 -i input_assemblies  -o checkm2_input_assemblies_out/ --database_path  "${params.databasesDir}/bin/checkm2_db/uniref100.KO.1.dmnd" -x "fasta" --force
    cp checkm2_input_assemblies_out/quality_report.tsv ./quality_report.tsv
    """


}

process AUTOCYCLER_COMPRESS { 
    tag {sample} 
    cpus 4 
    publishDir "assemblies/autocycler/${sample}", mode: "copy" 
    publishDir "assemblies/autocycler", mode: "copy", pattern: "consensus_assembly.fasta", saveAs: { filename -> "${sample}.fasta" } 
    publishDir "assemblies/autocycler", mode: "copy", pattern: "consensus_assembly.gfa", saveAs: { filename -> "${sample}.gfa" }
    publishDir "assemblies/autocycler", mode: "copy", pattern: "metrics.tsv", saveAs: { filename -> "${sample}.tsv" }
  
  
    input:
    tuple val(sample), path('assembly_01.fasta'), path('assembly_02.fasta'), path('assembly_03.fasta'), path('assembly_04.fasta'), path('assembly_05.fasta'), path('assembly_06.fasta'), path('assembly_07.fasta'), path('assembly_08.fasta'), path('assembly_09.fasta'), path('assembly_10.fasta'), path('assembly_11.fasta'), path('assembly_12.fasta'), path('assembly_13.fasta'), path('assembly_14.fasta'), path('assembly_15.fasta'), path('assembly_16.fasta'), path('assembly_17.fasta'), path('assembly_18.fasta'), path('assembly_19.fasta'), path('assembly_20.fasta')

    output:
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.fasta"), path("unitigs_num.txt"), emit: fa
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.gfa"), path("unitigs_num.txt"), emit: gfa
    tuple val(sample), val("autocycler"), path("autocycler_out"), emit: fol
    tuple val(sample), val("autocycler"), path("metrics.tsv"), emit: metrics
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.fasta"), path("autocycler_out/consensus_assembly.gfa"), path("autocycler_out"), path("metrics.tsv"), path("unitigs_num.txt"), emit: all    

    script:
    """
    # Make a directory for autocycler input
    mkdir -p "autocycler_assemblies"

    # Copy each symlinked input assembly to the autocycler folder using the original filename
    for assembly_file in assembly_*.fasta; do
       real_path=\$(readlink -f "\$assembly_file")
       real_name=\$(basename "\$real_path")
       cp "\$assembly_file" "autocycler_assemblies/\${real_name}"
    done

    
    # Compress the input assemblies into a unitig graph
    autocycler compress -i autocycler_assemblies -a autocycler_out --max_contigs 25

    # Cluster the input contigs into putative genomic sequences
    autocycler cluster -a autocycler_out --max_contigs 25 

    # Trim and resolve each QC-pass clusterfailed_consensus_assemblies_ch
    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c"
        autocycler resolve -c "\$c"
    done

    # Combine resolved clusters into a final assembly
    autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

    # Optional extra assess assemblies - same info as in .yml file, but more easiyl readable
    autocycler table > metrics.tsv
    autocycler table -a "autocycler_out" -n "${sample}" >> metrics.tsv 

     # Move consensus outputs to predictable path
    cp autocycler_out/consensus_assembly.fasta consensus_assembly.fasta
    cp autocycler_out/consensus_assembly.gfa consensus_assembly.gfa

    # Extract number of unitigs from the metrics block to filter failed assemblies more easily
    tail -n 1 metrics.tsv | awk -F '\t' '{ print \$(NF-1) }' > unitigs_num.txt
    
    """
}


process PUBLISH_FAILED_CONSENSUS {
    tag {sample}
    
    cpus 1
    publishDir "failed_assemblies", pattern: "consensus_assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_consensus_assembly.fasta"}
    publishDir "failed_assemblies", pattern: "metrics.tsv", mode: "copy", saveAs: { filename -> "${sample}_metrics.tsv"}

    input:
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.fasta"), path("autocycler_out/consensus_assembly.gfa"), path("autocycler_out"), path("metrics.tsv"), path("unitigs_num.txt")

    output:
    tuple val(sample), val("autocycler"), path("consensus_assembly.fasta"), path("autocycler_out/consensus_assembly.gfa"), path("autocycler_out"), path("metrics.tsv"), path("unitigs_num.txt")
   
    
    script:
    """
    cp "autocycler_out/consensus_assembly.fasta" "./consensus_assembly.fasta"
    echo "metrics.tsv" 
    """

}


process MEDAKA_FULL {

    label 'medaka'
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "assemblies/${source_contigs}_medaka_full", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), val(source_contigs), path('contigs.fasta'), val(source_reads), path('reads.fastq.gz')
  
    output:
    tuple val(sample), val("${source_contigs}_medaka_full"), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus -i reads.fastq.gz -d contigs.fasta -o output -t ${task.cpus}  --bacteria -m "${params.inputmodel}"
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}


process SEQKIT {
    tag {sample + ' ' + assembler}

    publishDir "assembly_QC/${assembler}", pattern: "*.tsv", mode: 'copy', saveAs: { filename -> "${sample}_seqkit_stats.tsv" }

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all assembly.fasta > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
}


process CHECKM2_CONSENSUS {
    tag {sample + ' ' + assembler}
    cpus 4
 
    publishDir "assembly_QC/${assembler}", pattern: "checkm2_out/quality_report.tsv", mode: 'copy', saveAs: { filename -> "${sample}_checkm2.tsv" }

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("checkm2_out/quality_report.tsv")

    script:
    """
    checkm2 predict --threads $task.cpus --input assembly.fasta --output-directory checkm2_out/ --stdout -x "fasta" --database_path  "${params.databasesDir}/bin/checkm2_db/uniref100.KO.1.dmnd"
    """

    stub:
    """
    mkdir checkm2_out
    touch checkm2_out/quality_report.tsv
    """

}


process MLST {
    tag {sample}
    cpus 4

    label 'short'

    publishDir "MLST/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path("assembly.fasta")
  
    output:
    tuple val(sample), val(assembler), path("${sample}.tsv"), optional: true

    script:
    """
    mlst assembly.fasta > ${sample}.tsv
    """
    
    stub:
    """
    touch ${sample}.tsv
    """
}


process KRAKEN2 {
    scratch true
    tag {sample + ' ' + assembler}
    cpus 32

    publishDir "kraken2/${assembler}", pattern: "*_k2_report.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_report.tsv"}
    publishDir "kraken2/${assembler}", pattern: "*_k2_output.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_output.tsv"}

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path('assembly.fasta'), path("${sample}_k2_output.tsv"), path("${sample}_k2_report.tsv"), emit: fasta

    script:
    """
    kraken2 --db "${params.databasesDir}/bin/kraken2_db" assembly.fasta  --use-names --report ${sample}_k2_report.tsv --output ${sample}_k2_output.tsv  --threads ${task.cpus}

    """

    stub:
    """
    touch ${sample}_k2_report.tsv
    touch ${sample}_k2_output.tsv
    """
}


process SPECIES_FROM_KRAKEN2 {
    tag {sample + ' ' + assembler}
    publishDir "kraken2/${assembler}", pattern: "*.txt", mode: 'copy', saveAs: { filename -> "${sample}_chromosome_species.txt"}
    

    input:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('output.tsv'), path('report.tsv')

    output:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('main_chromosome_species.txt')

    script:
    """
    awk -F'\t' '\$4 > 1000000' output.tsv | sort -t\$'\t' -k4,4nr > filtered.tsv

    if [ -s filtered.tsv ]; then
         awk -F'\t' '{print \$3}' filtered.tsv | awk '{print \$1, \$2}' | paste -sd, - > main_chromosome_species.txt

    else
        awk '\$4 == "S"' report.tsv > filtered.tsv
        awk -v max="\$max_value" '\$1 == max {print substr(\$0, index(\$0,\$6))}' filtered.tsv | paste -sd, - > main_chromosome_species.txt
    fi
    """

}

process AMRFINDERPLUS {
    errorStrategy 'ignore'    
    tag {sample + ' ' + assembler}
    cpus 32

    label 'short'

    publishDir "AMRFinderPlus/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('species.txt')

    output:
    tuple val(sample), val(assembler), path("${sample}.tsv")

    script:
    """
    species=\$(<species.txt)

    if [[ \$species == Escherichia* ]]; then
        amrfinder -n assembly.fasta -O Escherichia -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1

    elif [[ \$species == "Klebsiella pneumoniae" ]]; then
        amrfinder -n assembly.fasta -O Klebsiella_pneumoniae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Klebsiella oxytoca" ]]; then
        amrfinder -n assembly.fasta -O Klebsiella_oxytoca -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Enterobacter cloacae" ]]; then
        amrfinder -n assembly.fasta -O Enterobacter_cloacae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Enterobacter asburiae" ]]; then
        amrfinder -n assembly.fasta -O Enterobacter_absuriae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Citrobacter freundii" ]]; then
        amrfinder -n assembly.fasta -O Citrobacter_freundii -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    elif [[ \$species == "Serratia marcescens" ]]; then
        amrfinder -n assembly.fasta -O Serratia_marcescens -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    else
        amrfinder -n assembly.fasta -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1
    fi
    """
}


process BAKTA {
    tag {sample + ' ' + assembler}
    errorStrategy 'ignore' 

    publishDir "bakta/${assembler}", mode: 'copy'
    publishDir "bakta_gff/${assembler}", pattern: "*.gff3", mode: 'copy', saveAs: { filename -> "${sample}.gff3"} 

    cpus 16

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.gff3"), emit: gff3
    tuple val(sample), val(assembler), path("${sample}_bakta"), emit: fol
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.tsv"), emit: tsv    

    script:
    """
    bakta --db ${params.databasesDir}/bin/bakta_db_v6/db -t ${task.cpus} --prefix ${sample} -o ${sample}_bakta/ assembly.fasta
    """
    stub:
    """
    mkdir ${sample}_bakta
    touch ${sample}_bakta/${sample}.gbk
    """
}

process MOB_SUITE {
    tag {sample + ' ' + assembler}
    cpus 8

    label 'short'

    publishDir "MOB_suite/${assembler}", mode: 'copy', pattern: "${sample}_mobsuite_out"
    

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out"), emit: fol
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out/biomarkers.blast.txt"), path("${sample}_mobsuite_out/contig_report.txt"), path("${sample}_mobsuite_out/mge.report.txt"), path("${sample}_mobsuite_out/mobtyper_results.txt"), emit: all
    
    script:
    """
    # NOTE: Mob-typer information is automatically generated for all plasmids reconstructed by mob-recon
    mob_recon --infile assembly.fasta  --outdir "${sample}_mobsuite_out" -g "${params.databasesDir}/bin/mob_suite_db/2019-11-NCBI-Enterobacteriacea-Chromosomes.fasta" --force
    """

    stub:
    """
    mkdir -p "${sample}_mobsuite_out"
    touch ${sample}_mobsuite_out/biomarkers.blast.txt
    touch ${sample}_mobsuite_out/contig_report.txt
    touch ${sample}_mobsuite_out/mge.report.txt
    touch ${sample}_mobsuite_out/mobtyper_results.txt
    """

}


process SUMMARISE_RAWQC {

    publishDir "summaries", mode: 'copy', pattern: "read_qc_merged.tsv"

    input:
    val(qc_items)

    output:
    path("read_qc_merged.tsv")

    script:
    def num_items = qc_items.size()
    if (num_items % 3 != 0)
        throw new IllegalArgumentException("Channel input must be a flat list with triples of (sample, source, path)")

    def tmp_lines = []
    for (int i = 0; i < num_items; i += 3) {
        def sample = qc_items[i]
        def source = qc_items[i + 1]
        def path = qc_items[i + 2]
        def filename = path.getName()
        def quoted_path = path.toString().replace('$','\\$')

        tmp_lines << """tail -n +2 "${quoted_path}" | awk -v s="${sample}" -v src="${source}" -v fn="${filename}" 'BEGIN{OFS="\\t"} {print s, src, fn, \$0}'"""
    }

    def script = """
    head -n 1 "${qc_items[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "source", "filename", \$0}' > read_qc_merged.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> read_qc_merged.tsv
    """

    return script
}


process SUMMARISE_CONTAMINATED {

    tag { "contaminated_summary" }
    publishDir "summaries", mode: 'copy', pattern: "contaminated_reads_checkm2_summary.tsv"

    input:
    val(contam_tuples)

    output:
    path("contaminated_reads_checkm2_summary.tsv")

    script:
    def num_items = contam_tuples.size()
    if (num_items % 8 != 0)
        throw new IllegalArgumentException("Expected flat list of 8-tuples, but got ${num_items} items.")

    def append_lines = []
    for (int i = 0; i < num_items; i += 8) {
        def sample = contam_tuples[i]
        def source = contam_tuples[i + 1]
        def checkm2_report = contam_tuples[i + 5].toString().replace('$','\\$')

        append_lines << """tail -n +2 "${checkm2_report}" | awk -v s="${sample}" -v t="${source}" 'BEGIN{OFS="\\t"} {print s, t, \$0}'"""
    }

    def script = """
    head -n 1 "${contam_tuples[5]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "source", \$0}' > contaminated_reads_checkm2_summary.tsv
    {
    ${append_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> contaminated_reads_checkm2_summary.tsv
    """

    return script
}

process SUMMARISE_INPUT_ASSEMBLY_QC {

    tag { "merge_checkm2_reports" }
    publishDir "summaries", mode: 'copy', pattern: "input_assemblies_checkm2_merged.tsv"

    input:
    val(reports)

    output:
    path("input_assemblies_checkm2_merged.tsv")

    script:
    def num_items = reports.size()
    if (num_items % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of (sample, 'autocycler', path), but got ${num_items} items.")

    def append_lines = []
    for (int i = 0; i < reports.size(); i += 3) {
        def sample = reports[i]
        def path = reports[i + 2]
        def quoted_path = path.toString().replace('$','\\$')

        append_lines << """tail -n +2 "${quoted_path}" | awk -v s="${sample}" 'BEGIN{OFS="\\t"} {print s, \$0}'"""
    }

    def script = """
    head -n 1 "${reports[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", \$0}' > input_assemblies_checkm2_merged.tsv
    {
    ${append_lines.join('\n')}
    } | sort -k1,1 >> input_assemblies_checkm2_merged.tsv
    """

    return script
}




process SUMMARISE_CONSENSUS_ASSEMBLY_SEQKIT {

    tag { "consensus_seqkit_summary" }
    publishDir "summaries", mode: 'copy', pattern: "consensus_assembly_seqkit_summary.tsv"

    input:
    val(input_tuples)

    output:
    path("consensus_assembly_seqkit_summary.tsv")

    script:
    def n = input_tuples.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 3) {
        def sample = input_tuples[i]
        def assembler = input_tuples[i + 1]
        def tsv = input_tuples[i + 2].toString().replace('$', '\\$')
        tmp_lines << """tail -n +2 "${tsv}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """
    head -n 1 "${input_tuples[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > consensus_assembly_seqkit_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> consensus_assembly_seqkit_summary.tsv
    """

    return script
}


process SUMMARISE_CONSENSUS_ASSEMBLY_CHECKM2 {

    tag { "consensus_checkm2_summary" }
    publishDir "summaries", mode: 'copy', pattern: "consensus_assembly_checkm2_summary.tsv"

    input:
    val(input_tuples)

    output:
    path("consensus_assembly_checkm2_summary.tsv")

    script:
    def n = input_tuples.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 3) {
        def sample = input_tuples[i]
        def assembler = input_tuples[i + 1]
        def tsv = input_tuples[i + 2].toString().replace('$', '\\$')
        tmp_lines << """tail -n +2 "${tsv}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """
    head -n 1 "${input_tuples[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > consensus_assembly_checkm2_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> consensus_assembly_checkm2_summary.tsv
    """

    return script
}


process SUMMARISE_MLST {

    tag { "mlst_summary" }
    publishDir "summaries", mode: 'copy', pattern: "mlst_summary.tsv"

    input:
    val(mlst_reports)

    output:
    path("mlst_summary.tsv")

    script:
    def n = mlst_reports.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 3) {
        def sample = mlst_reports[i]
        def assembler = mlst_reports[i + 1]
        def tsv = mlst_reports[i + 2].toString().replace('$', '\\$')
        tmp_lines << """tail "${tsv}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
     }
    
    def script = """
    echo -e "sample\tassembler\tfilename\ttyping_scheme\tmlst\tgene1\tgene2\tgene3\tgene4\tgene5_1\tgene5_2\tgene6" > mlst_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mlst_summary.tsv
    """
    return script
}


process SUMMARISE_KRAKEN2 {
    publishDir "summaries", mode: 'copy', pattern: "kraken2_summary.tsv"

    input:
    val(kraken2_tuples)

    output:
    path("kraken2_summary.tsv")

    script:
    def n = kraken2_tuples.size()
    if (n % 5 != 0) 
        throw new IllegalArgumentException("Expected flat list of 4-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 5) { 
        def sample = kraken2_tuples[i]
        def assembler = kraken2_tuples[i + 1]
        def output_tsv = kraken2_tuples[i + 3].toString().replace('$', '\\$')

        // Properly extract the species name as the full 3rd column using awk
        tmp_lines << '''awk -v s="''' + sample + '''" -v a="''' + assembler + '''" 'BEGIN{OFS="\\t"} {
    completeness=$1;
    contig_id=$2;
    species=$3;
    i=4;
    while(i<=NF && $i !~ "\\)$") {
        species=species " " $i;
        i++
    }
    if(i<=NF) species=species " " $i;
    print s, a, completeness, contig_id, species
}' "''' + output_tsv + '''"'''
    }

    def script = """
    echo -e "sample\tassembler\tcompleteness\tcontig_id\tspecies\tlength" > kraken2_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> kraken2_summary.tsv
    """

    return script
}



process SUMMARISE_AMRFINDER {

    tag { "amrfinder_summary" }
    publishDir "summaries", mode: 'copy', pattern: "amrfinder_summary.tsv"

    input:
    val(amrfinder_input)

    output:
    path("amrfinder_summary.tsv")

    script:
    def n = amrfinder_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 3) {
        def sample = amrfinder_input[i]
        def assembler = amrfinder_input[i + 1]
        def tsv = amrfinder_input[i + 2].toString().replace('$', '\\$')
        tmp_lines << """tail -n +2 "${tsv}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """
    head -n 1 "${amrfinder_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > amrfinder_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> amrfinder_summary.tsv
    """
    return script
}

process SUMMARISE_BAKTA {
    
    tag { "bakta_summary" }
    publishDir "summaries", mode: 'copy', pattern: "bakta_summary.tsv"

    input:
    val(bakta_input)

    output:
    path("bakta_summary.tsv")

    script:
    def n = bakta_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def tmp_lines = []
    for (int i = 0; i < n; i += 3) {
        def sample = bakta_input[i]
        def assembler = bakta_input[i + 1]
        def dir_path = bakta_input[i + 2].toString().replace('$', '\\$')
        def tsv = "${dir_path}/${sample}.tsv"
        tmp_lines << """tail -n +7 "${tsv}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """
    sed -n '6p' "${bakta_input[2]}/${bakta_input[0]}.tsv" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > bakta_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> bakta_summary.tsv
    """
    return script
}

process SUMMARISE_MOBSUITE {
    
    tag { "mobsuite_summary" }
    publishDir "summaries", mode: 'copy', pattern: "*mobsuite*.tsv"

    input:
    val(mobsuite_input)

    output:
    path("mobsuite_biomarkers_summary.tsv")
    path("mobsuite_contigs_summary.tsv")
    path("mobsuite_mge_summary.tsv")
    path("mobsuite_mobtyper_results.tsv")

    script:
    def n = mobsuite_input.size()
    if (n % 6 != 0)
        throw new IllegalArgumentException("Expected flat list of 6-tuples, but got ${n} items.")

    def biomarkers_lines = []
    def contigs_lines = []
    def mge_lines = []
    def mobtyper_lines = []

    for (int i = 0; i < n; i += 6) {
        def sample = mobsuite_input[i]
        def assembler = mobsuite_input[i + 1]
        def biomarker = mobsuite_input[i + 2].toString().replace('$', '\\$')
        def contig = mobsuite_input[i + 3].toString().replace('$', '\\$')
        def mge = mobsuite_input[i + 4].toString().replace('$', '\\$')
        def mobtyper = mobsuite_input[i + 5].toString().replace('$', '\\$')

        biomarkers_lines << """tail -n +2 "${biomarker}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
        contigs_lines    << """tail -n +2 "${contig}"    | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
        mge_lines        << """tail -n +2 "${mge}"       | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
        mobtyper_lines   << """tail -n +2 "${mobtyper}"  | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """
    head -n 1 "${mobsuite_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_biomarkers_summary.tsv
    {
    ${biomarkers_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_biomarkers_summary.tsv

    head -n 1 "${mobsuite_input[3]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_contigs_summary.tsv
    {
    ${contigs_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_contigs_summary.tsv

    head -n 1 "${mobsuite_input[4]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_mge_summary.tsv
    {
    ${mge_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_mge_summary.tsv

    head -n 1 "${mobsuite_input[5]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_mobtyper_results.tsv
    {
    ${mobtyper_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_mobtyper_results.tsv
    """
    return script
}



