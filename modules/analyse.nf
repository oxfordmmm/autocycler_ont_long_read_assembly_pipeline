process SEQKIT {

    conda = "${projectDir}/envs/seqkit.yml"
    tag {sample + ' ' + assembler}
    publishDir "${params.outdir}/assembly_QC/${assembler}", pattern: "*.tsv", mode: 'copy', saveAs: { filename -> "${sample}_seqkit_stats.tsv" }

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

    conda = "${projectDir}/envs/checkm2.yml"
    tag {sample + ' ' + assembler}
    cpus 4
    publishDir "${params.outdir}/assembly_QC/${assembler}", pattern: "checkm2_out/quality_report.tsv", mode: 'copy', saveAs: { filename -> "${sample}_checkm2.tsv" }

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

    conda = "${projectDir}/envs/mlst.yml"
    tag {sample}
    cpus 4 // not called in command
    publishDir "${params.outdir}/MLST/${assembler}", mode: 'copy'

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

    conda = "${projectDir}/envs/kraken2.yml"
    scratch true
    tag {sample + ' ' + assembler}
    cpus 32     //needs a very large amount of RAM to run so need many CPUs
    publishDir "${params.outdir}/kraken2/${assembler}", pattern: "*_k2_report.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_report.tsv"}
    publishDir "${params.outdir}/kraken2/${assembler}", pattern: "*_k2_output.tsv", mode: 'copy', saveAs: { filename -> "${sample}_k2_output.tsv"}

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
    publishDir "${params.outdir}/kraken2/${assembler}", pattern: "*.txt", mode: 'copy', saveAs: { filename -> "${sample}_chromosome_species.txt"}

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


process BAKTA {

    conda = "${projectDir}/envs/bakta.yml"
    tag {sample + ' ' + assembler}
    errorStrategy 'retry'
    publishDir "${params.outdir}/bakta/${assembler}", mode: 'copy'

    cpus 16

    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.gff3"), emit: gff3
    tuple val(sample), val(assembler), path("${sample}_bakta"), emit: fol
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.tsv"), emit: tsv
    tuple val(sample), val(assembler), path("${sample}_bakta/${sample}.faa"), path("${sample}_bakta/${sample}.gff3"), emit: protein

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


process AMRFINDERPLUS {

    conda = "${projectDir}/envs/amrfinder.yml"
    errorStrategy 'retry'
    tag {sample + ' ' + assembler}
    cpus 2     //cpus set to low number as not enough RAM on VMs for blast search if >4 multi-threading and get segmentation fault

    publishDir "${params.outdir}/AMRFinderPlus/${assembler}", mode: 'copy'

    input:
    tuple val(sample), val(assembler), path('assembly.fasta'), path('species.txt'), val(assembler), path('assembly.faa'), path('assembly.gff3')

    output:
    tuple val(sample), val(assembler), path("${sample}.tsv")

    script:
    """
    # Convert assembly header line to be in format '>contig_1' instead of '>1'
    awk '/^>/{print ">contig_" ++i; next} {print}' assembly.fasta > amrfinder_input_assembly.fasta

    species=\$(<species.txt)

    if [[ \$species == Escherichia* ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Escherichia -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1 --annot>

    elif [[ \$species == "Klebsiella pneumoniae" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Klebsiella_pneumoniae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18>
    elif [[ \$species == "Klebsiella oxytoca" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Klebsiella_oxytoca -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1 >
    elif [[ \$species == "Enterobacter cloacae" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Enterobacter_cloacae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.>
    elif [[ \$species == "Enterobacter asburiae" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Enterobacter_absuriae -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18>
    elif [[ \$species == "Citrobacter freundii" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Citrobacter_freundii -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.>
    elif [[ \$species == "Serratia marcescens" ]]; then
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -O Serratia_marcescens -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1>
    else
        amrfinder --nucleotide amrfinder_input_assembly.fasta --protein assembly.faa --gff assembly.gff3 -o ${sample}.tsv --plus -d ${params.databasesDir}/bin/amrfinder_db/2024-12-18.1 --annotation_format ba>
    fi
    """
}



process MOB_SUITE {

    conda = "${projectDir}/envs/mob_suite.yml"
    tag {sample + ' ' + assembler}
    cpus 8  //not called in command
    publishDir "${params.outdir}/MOB_suite/${assembler}", mode: 'copy', pattern: "${sample}_mobsuite_out"


    input:
    tuple val(sample), val(assembler), path('assembly.fasta')

    output:
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out"), emit: fol
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out/biomarkers.blast.txt"), emit: biomarkers, optional: true
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out/contig_report.txt"), emit: contigs, optional: true
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out/mge.report.txt"), emit: mges, optional: true
    tuple val(sample), val(assembler), path("${sample}_mobsuite_out/mobtyper_results.txt"), emit: mobtyper, optional: true
    //tuple val(sample), val(assembler), path("${sample}_mobsuite_out/biomarkers.blast.txt", optional: true), path("${sample}_mobsuite_out/contig_report.txt", optional: true), path("${sample}_mobsuite_out/mg>

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
