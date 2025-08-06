process RAVEN_ASSEMBLE_INITIAL {

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample}
    publishDir "${params.outdir}/assemblies/raven_initial", pattern: "*.fasta", mode: 'copy', saveAs: { filename -> "${sample}_raven_initial.fasta"}
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

    conda = "${projectDir}/envs/checkm2.yml"
    tag {sample}
    publishDir "${params.outdir}/checkm2_initial/", pattern: "checkm2_out/quality_report.tsv", mode: 'copy', saveAs: {filename -> "${sample}_checkm2_initial.tsv"}
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

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample}
    publishDir "${params.outdir}/subsamples", pattern: "subsampled_reads", mode: "copy", saveAs: { filename -> "${sample}_subsampled_reads" }
    //not called in command
    cpus 4
    
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

    conda = "${projectDir}/envs/seqkit.yml"
    tag {sample}
    publishDir "${params.outdir}/raw_QC_CSVs/${source}", mode: 'copy'

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
    publishDir "${params.outdir}/contamination", pattern: "*raven_assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_assembly.fasta"}
    publishDir "${params.outdir}/contamination", pattern: "*quality_report.tsv", mode: "copy", saveAs: { filename -> "${sample}_quality_report.tsv"}
 
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

    conda = "${projectDir}/envs/seqkit.yml"
    tag {sample}
    publishDir "${params.outdir}/raw_QC_CSVs/${source}_${subsample_number}", mode: 'copy'

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

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample + ' ' + number}
    cpus 16

    publishDir "${params.outdir}/assemblies/canu_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}_canu_01.fasta"}
    publishDir "${params.outdir}/assemblies/canu_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}_canu_02.fasta"}
    publishDir "${params.outdir}/assemblies/canu_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}_canu_03.fasta"}
    publishDir "${params.outdir}/assemblies/canu_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}_canu_04.fasta"}
    publishDir "${params.outdir}/assemblies/canu_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}_canu_01.gfa"}
    publishDir "${params.outdir}/assemblies/canu_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}_canu_02.gfa"}
    publishDir "${params.outdir}/assemblies/canu_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}_canu_03.gfa"}
    publishDir "${params.outdir}/assemblies/canu_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}_canu_04.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('canu'), val(number), path("canu_${number}.fasta"), emit: fa
    tuple val(sample), val('canu'), val(number), path("canu_${number}.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)
 
    # Updates autocycler script:
    autocycler helper canu --reads subsampled_reads.fastq --out_prefix canu_${number} --threads ${task.cpus} --genome_size \$genome_size --read_type ont_r10
    # old script
    # ${projectDir}/bin/autocycler_scripts/canu.sh subsampled_reads.fastq canu_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_FLYE {

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample + ' ' + number}
    cpus 8

    publishDir "${params.outdir}/assemblies/flye_01", pattern: "*01/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_flye_01.fasta"}
    publishDir "${params.outdir}/assemblies/flye_02", pattern: "*02/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_flye_01.fasta"}
    publishDir "${params.outdir}/assemblies/flye_03", pattern: "*03/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_flye_01.fasta"}
    publishDir "${params.outdir}/assemblies/flye_04", pattern: "*04/assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_flye_01.fasta"}
    publishDir "${params.outdir}/assemblies/flye_01", pattern: "*01/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}_flye_01.gfa"}
    publishDir "${params.outdir}/assemblies/flye_02", pattern: "*02/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}_flye_01.gfa"}
    publishDir "${params.outdir}/assemblies/flye_03", pattern: "*03/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}_flye_01.gfa"}
    publishDir "${params.outdir}/assemblies/flye_04", pattern: "*04/assembly_graph.gfa", mode: "copy", saveAs: { filename -> "${sample}_flye_01.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('flye'), val(number), path("flye_${number}.fasta"), emit: fa
    tuple val(sample), val('flye'), val(number), path("flye_${number}.gfa"), optional: true, emit: gfa
    
    script:
    """
    genome_size=\$(<genome_size.txt)
    
    # Updates Autocycler script command:
       
    autocycler helper flye --reads subsampled_reads.fastq --out_prefix flye_${number} --threads ${task.cpus} --genome_size \$genome_size --read_type ont_r10

    # old direct Flye command
    #flye --nano-hq subsampled_reads.fastq --threads "${task.cpus}" --out-dir flye_out_${number}
    # old autocycler command
    ##${projectDir}/bin/autocycler_scripts/flye.sh subsampled_reads.fastq flye_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_MINIASM {

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample + ' ' + number}
    cpus 8

    publishDir "${params.outdir}/assemblies/miniasm_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}_miniasm_01.fasta"}
    publishDir "${params.outdir}/assemblies/miniasm_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}_miniasm_02.fasta"}
    publishDir "${params.outdir}/assemblies/miniasm_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}_miniasm_03.fasta"}
    publishDir "${params.outdir}/assemblies/miniasm_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}_miniasm_04.fasta"}
    publishDir "${params.outdir}/assemblies/miniasm_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}_miniasm_01.gfa"}
    publishDir "${params.outdir}/assemblies/miniasm_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}_miniasm_02.gfa"}
    publishDir "${params.outdir}/assemblies/miniasm_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}_miniasm_03.gfa"}
    publishDir "${params.outdir}/assemblies/miniasm_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}_miniasm_04.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.fasta"), emit: fa
    tuple val(sample), val('miniasm'), val(number), path("miniasm_${number}.gfa"), optional: true, emit: gfa
    
    script:
    """
    genome_size=\$(<genome_size.txt)

    # Updates Autocycler script
    autocycler helper miniasm --reads subsampled_reads.fastq --out_prefix miniasm_${number} --threads ${task.cpus} --genome_size \$genome_size --read_type ont_r10

    # old script command
    #${projectDir}/bin/autocycler_scripts/miniasm.sh subsampled_reads.fastq miniasm_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_RAVEN {

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample + ' ' + number}
    cpus 8

    publishDir "${params.outdir}/assemblies/raven_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_01.fasta"}
    publishDir "${params.outdir}/assemblies/raven_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_02.fasta"}
    publishDir "${params.outdir}/assemblies/raven_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_03.fasta"}
    publishDir "${params.outdir}/assemblies/raven_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}_raven_04.fasta"}
    publishDir "${params.outdir}/assemblies/raven_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}_raven_01.gfa"}
    publishDir "${params.outdir}/assemblies/raven_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}_raven_02.gfa"}
    publishDir "${params.outdir}/assemblies/raven_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}_raven_03.gfa"}
    publishDir "${params.outdir}/assemblies/raven_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}_raven_04.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('raven'), val(number), path("raven_${number}.fasta"), emit: fa
    tuple val(sample), val('raven'), val(number), path("raven_${number}.gfa"), optional: true, emit: gfa

    script:
    """
    genome_size=\$(<genome_size.txt)
    
    # Updates Autocycler script:
    autocycler helper raven --reads subsampled_reads.fastq --out_prefix raven_${number} --threads ${task.cpus} --genome_size \$genome_size --read_type ont_r10

    # old script command
    # ${projectDir}/bin/autocycler_scripts/raven.sh subsampled_reads.fastq raven_${number} "${task.cpus}" "\$genome_size"
    """
}

process AUTOCYCLER_ASSEMBLE_PLASSEMBLER {

    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample + ' ' + number}
    cpus 8

    publishDir "${params.outdir}/assemblies/plassembler_01", pattern: "*01.fasta", mode: "copy", saveAs: { filename -> "${sample}_plassembler_01.fasta"}
    publishDir "${params.outdir}/assemblies/plassembler_02", pattern: "*02.fasta", mode: "copy", saveAs: { filename -> "${sample}_plassembler_02.fasta"}
    publishDir "${params.outdir}/assemblies/plassembler_03", pattern: "*03.fasta", mode: "copy", saveAs: { filename -> "${sample}_plassembler_03.fasta"}
    publishDir "${params.outdir}/assemblies/plassembler_04", pattern: "*04.fasta", mode: "copy", saveAs: { filename -> "${sample}_plassembler_04.fasta"}
    publishDir "${params.outdir}/assemblies/plassembler_01", pattern: "*01.gfa", mode: "copy", saveAs: { filename -> "${sample}_plassembler_01.gfa"}
    publishDir "${params.outdir}/assemblies/plassembler_02", pattern: "*02.gfa", mode: "copy", saveAs: { filename -> "${sample}_plassembler_02.gfa"}
    publishDir "${params.outdir}/assemblies/plassembler_03", pattern: "*03.gfa", mode: "copy", saveAs: { filename -> "${sample}_plassembler_03.gfa"}
    publishDir "${params.outdir}/assemblies/plassembler_04", pattern: "*04.gfa", mode: "copy", saveAs: { filename -> "${sample}_plassembler_04.gfa"}

    input:
    tuple val(sample), val(source),  val(number), path("genome_size.txt"), path("subsampled_reads.fastq")

    output:
    tuple val(sample), val('plassembler'), val(number), path("plassembler_${number}.fasta"), emit: fa, optional: true
    tuple val(sample), val('plassembler'), val(number), path("plassembler_${number}.gfa"), emit: gfa, optional: true
    tuple val(sample), val('plassembler'), val(number), path("plassembler_${number}.tsv"), emit: tsv, optional: true

    script:
    """
    # Genome size not used, as plassembler only needs a minimum lower bound chromosome size, and etimated genome size may be to large for this. default: 1000000.
    # genome_size=\$(<genome_size.txt)

    # Plassembler long command directly
    plassembler long -d "${params.databasesDir}/bin/plassembler_db" -l subsampled_reads.fastq -o plassembler_out -t ${task.cpus} --force --skip_qc

    # move outputs to standardised paths
    mv plassembler_out/plassembler_plasmids.fasta plassembler_${number}.fasta
    mv plassembler_out/plassembler_plasmids.gfa plassembler_${number}.gfa
    mv plassembler_out/plassembler_summary.tsv plassembler_${number}.tsv

    # Command directly from automated autocycler script
    # but this uses the standard db path, an no option within Autocycler scripts to specify, so use plassembler long command directly
    #autocycler helper plassembler --reads subsampled_reads/sample_${number}.fastq --out_prefix assemblies/plassembler_${number} --threads ${task.cpus} --genome_size "\$genome_size" --read_type "ont_r10" --min_depth_rel 0.1 
    
    """
}


process CHECKM2_INPUT_ASSEMBLIES {

    conda = "${projectDir}/envs/checkm2.yml"
    tag {sample}
    cpus 32
    publishDir "${params.outdir}/assembly_QC/input_assemblies", pattern: "*.tsv", mode: 'copy', saveAs: { filename -> "${sample}_input_checkm2.tsv" }
    


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
    checkm2 predict  -t ${task.cpus} -i input_assemblies  -o checkm2_input_assemblies_out/ --database_path  "${params.databasesDir}/bin/checkm2_db/uniref100.KO.1.dmnd" -x "fasta" --force
    cp checkm2_input_assemblies_out/quality_report.tsv ./quality_report.tsv
    """


}


process AUTOCYCLER_COMPRESS {
    
    conda = "${projectDir}/envs/autocycler_v5_with_assemblers.yml"
    tag {sample}
    cpus 4 //not called during command

    publishDir "${params.outdir}/assemblies/autocycler/${sample}", mode: "copy"
    publishDir "${params.outdir}/assemblies/autocycler", mode: "copy", pattern: "consensus_assembly.fasta", saveAs: { filename -> "${sample}.fasta" }
    publishDir "${params.outdir}/assemblies/autocycler", mode: "copy", pattern: "consensus_assembly.gfa", saveAs: { filename -> "${sample}.gfa" }
    publishDir "${params.outdir}/assemblies/autocycler", mode: "copy", pattern: "metrics.tsv", saveAs: { filename -> "${sample}.tsv" }    

    input:
    tuple val(sample), path('assembly_01.fasta'), path('assembly_02.fasta'), path('assembly_03.fasta'), path('assembly_04.fasta'), path('assembly_05.fasta'), path('assembly_06.fasta'), path('assembly_07.fasta'), path('assembly_08.fasta'), path('assembly_09.fasta'), path('assembly_10.fasta'), path('assembly_11.fasta'), path('assembly_12.fasta'), path('assembly_13.fasta'),  path('assembly_14.fasta'), path('assembly_15.fasta'), path('assembly_16.fasta'), path('assembly_17.fasta'), path('assembly_18.fasta'), path('assembly_19.fasta'), path('assembly_20.fasta')
    
    output:
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.fasta"), path("unitigs_num.txt"), emit: fa
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.gfa"), path("unitigs_num.txt"), emit: gfa
    tuple val(sample), val("autocycler"), path("autocycler_out"), emit: fol
    tuple val(sample), val("autocycler"), path("metrics.tsv"), emit: metrics
    tuple val(sample), val("autocycler"), path("autocycler_out/consensus_assembly.fasta"), path("autocycler_out/consensus_assembly.gfa"), path("autocycler_out"), path("metrics.tsv"), path("unitigs_num.txt"), emit: all


    script:
    '''
    # Make a directory for autocycler input
    mkdir -p "autocycler_assemblies"

    # Copy each symlinked input assembly to the autocycler folder using the original filename
    for assembly_file in assembly_*.fasta; do
       real_path=\$(readlink -f "\$assembly_file")
       real_name=\$(basename "\$real_path")
       cp "\$assembly_file" "autocycler_assemblies/\${real_name}"
    done

    # Remove any empty files from the autocycler input folder
    find autocycler_assemblies/ -type f -size 0 -print -delete

    # Give circular contigs from Plassembler extra clustering weight
    for f in autocycler_assemblies/plassembler*.fasta; do
        [ -e "$f" ] || continue
        sed -i 's/circular=True/circular=True Autocycler_cluster_weight=2/' "$f"
    done

    # Give contigs from Canu and Flye extra consensus weight
    for f in autocycler_assemblies/canu*.fasta autocycler_assemblies/flye*.fasta; do
        [ -e "$f" ] || continue
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
    done


    # Compress the input assemblies into a unitig graph
    autocycler compress -i "autocycler_assemblies" -a "autocycler_out" --max_contigs 25

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
    autocycler table -a "autocycler_out" >> metrics.tsv

     # Move consensus outputs to predictable path
    cp autocycler_out/consensus_assembly.fasta consensus_assembly.fasta
    cp autocycler_out/consensus_assembly.gfa consensus_assembly.gfa

    # Extract number of unitigs from the metrics block to filter failed assemblies more easily
    tail -n 1 metrics.tsv | awk -F '\t' '{ print \$(NF-1) }' > unitigs_num.txt

    '''
}



process PUBLISH_FAILED_CONSENSUS {
    
    tag {sample}
    cpus 1
    publishDir "${params.outdir}/failed_assemblies", pattern: "consensus_assembly.fasta", mode: "copy", saveAs: { filename -> "${sample}_consensus_assembly.fasta"}
    publishDir "${params.outdir}/failed_assemblies", pattern: "metrics.tsv", mode: "copy", saveAs: { filename -> "${sample}_metrics.tsv"}

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

    conda = "${projectDir}/envs/medaka.yml"
    tag {sample + ' ' + source_contigs}
    cpus 4

    label 'short'

    publishDir "${params.outdir}/assemblies/${source_contigs}_medaka_full", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

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


