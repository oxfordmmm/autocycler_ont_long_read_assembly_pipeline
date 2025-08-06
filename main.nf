#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2


// include modules
include {RAVEN_ASSEMBLE_INITIAL} from './modules/assemble.nf'
include {GENOME_SIZE_AND_CONTAMINATION} from './modules/assemble.nf'
include {PUBLISH_CONTAMINATED} from './modules/assemble.nf'
include {RAWQC_ONT} from './modules/assemble.nf'
include {RAWQC_ONT_SUBSAMPLED} from './modules/assemble.nf'
include {AUTOCYCLER_SUBSAMPLE} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_CANU} from './modules/assemble.nf' 
include {AUTOCYCLER_ASSEMBLE_FLYE} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_MINIASM} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_RAVEN} from './modules/assemble.nf'
include {AUTOCYCLER_ASSEMBLE_PLASSEMBLER} from './modules/assemble.nf'
include {CHECKM2_INPUT_ASSEMBLIES} from './modules/assemble.nf'
include {AUTOCYCLER_COMPRESS} from './modules/assemble.nf'
include {PUBLISH_FAILED_CONSENSUS} from './modules/assemble.nf'
include {MEDAKA_FULL} from './modules/assemble.nf'
include {SEQKIT} from './modules/analyse.nf'
include {CHECKM2_CONSENSUS} from './modules/analyse.nf'
include {BAKTA} from './modules/analyse.nf'
include {MLST} from './modules/analyse.nf'
include {KRAKEN2} from './modules/analyse.nf'
include {SPECIES_FROM_KRAKEN2} from './modules/analyse.nf'
include {AMRFINDERPLUS} from './modules/analyse.nf'
include {MOB_SUITE} from './modules/analyse.nf'
include {SUMMARISE_RAWQC} from './modules/summarise.nf'
include {SUMMARISE_CONTAMINATED} from './modules/summarise.nf'
include {SUMMARISE_INPUT_ASSEMBLY_QC} from './modules/summarise.nf'
include {SUMMARISE_CONSENSUS_ASSEMBLY_CHECKM2} from './modules/summarise.nf'
include {SUMMARISE_CONSENSUS_ASSEMBLY_SEQKIT} from './modules/summarise.nf'
include {SUMMARISE_MLST} from './modules/summarise.nf'
include {SUMMARISE_KRAKEN2} from './modules/summarise.nf'
include {SUMMARISE_AMRFINDER} from './modules/summarise.nf'
include {SUMMARISE_BAKTA} from './modules/summarise.nf'
include {SUMMARISE_MOBSUITE_BIOMARKERS} from './modules/summarise.nf'
include {SUMMARISE_MOBSUITE_CONTIGS} from './modules/summarise.nf'
include {SUMMARISE_MOBSUITE_MGES} from './modules/summarise.nf'
include {SUMMARISE_MOBSUITE_MOBTYPER} from './modules/summarise.nf'

workflow {
    
    // INPUT CHANNELS: 
    Channel.fromPath( "${params.inputFastq_ONT}/*.fastq" )
           .map{ file -> tuple(file.simpleName, 'raw_ont', file) }
           .set{ ont_labeled_ch }
    //ont_labeled_ch.view()
     
         
    
    main:
    
    RAVEN_ASSEMBLE_INITIAL(ont_labeled_ch)
    GENOME_SIZE_AND_CONTAMINATION(RAVEN_ASSEMBLE_INITIAL.out.fasta)
    //GENOME_SIZE_AND_CONTAMINATION.out.data.view()

    //Split genomes based on contamination
    contam_free_reads_ch = GENOME_SIZE_AND_CONTAMINATION.out.data
        .filter { sample, source, reads, raven, raven_assembly, checkm2_report, contamination, genome_size -> contamination.text.trim().toFloat() < 50.0}
    
    contam_reads_ch = GENOME_SIZE_AND_CONTAMINATION.out.data
        .filter { sample, source, reads, raven, raven_assembly, checkm2_report, contamination, genome_size -> contamination.text.trim().toFloat() >= 50.0}
    PUBLISH_CONTAMINATED(contam_reads_ch)    


        
    AUTOCYCLER_SUBSAMPLE(contam_free_reads_ch)
    subsampled_fqs_mixed = AUTOCYCLER_SUBSAMPLE.out.fq01.mix(AUTOCYCLER_SUBSAMPLE.out.fq02, AUTOCYCLER_SUBSAMPLE.out.fq03, AUTOCYCLER_SUBSAMPLE.out.fq04)
  
    RAWQC_ONT(ont_labeled_ch)
    RAWQC_ONT_SUBSAMPLED(subsampled_fqs_mixed)
    

    // Generate input assemblies
    AUTOCYCLER_ASSEMBLE_CANU(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_FLYE(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_MINIASM(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_RAVEN(subsampled_fqs_mixed)
    AUTOCYCLER_ASSEMBLE_PLASSEMBLER(subsampled_fqs_mixed)
    //replaced hybracter with plassembler, and added 2x weight to Flye (and Canu) assemblies instead
    
    assemblies_for_autocycler =  AUTOCYCLER_ASSEMBLE_CANU.out.fa.mix(AUTOCYCLER_ASSEMBLE_FLYE.out.fa, AUTOCYCLER_ASSEMBLE_MINIASM.out.fa, AUTOCYCLER_ASSEMBLE_RAVEN.out.fa, AUTOCYCLER_ASSEMBLE_PLASSEMBLER.out.fa)
    assemblies_for_autocycler_grouped = assemblies_for_autocycler.groupTuple(sort: true)
    assemblies_for_autocycler_grouped = assemblies_for_autocycler_grouped.map {sample, assemblers, numbers, paths ->
         return tuple(sample, paths)
         }
    //assemblies_for_autocycler_grouped.view()

    flattened_assemblies = assemblies_for_autocycler_grouped.map {sample, paths -> tuple(sample, *paths) }
    

    CHECKM2_INPUT_ASSEMBLIES(flattened_assemblies)

    AUTOCYCLER_COMPRESS(flattened_assemblies)

    //FILTER out failed Autocycler consensus assemblies where >100 unitigs
    // use channel with all outputs to make publishing failed assemblies easier
    successful_consensus_assemblies_ch = AUTOCYCLER_COMPRESS.out.all
        .filter { sample, assembler, consensus_assembly_fasta, consensus_assembly_gfa, consensus_assembly_fol, consensus_assembly_metrics, unitigs_num -> unitigs_num.text.trim().toFloat() < 100.0}
        .map { sample, assembler, consensus_assembly_fasta, consensus_assembly_gfa, consensus_assembly_fol, consensus_assembly_metrics, unitigs_num -> tuple(sample, assembler, consensus_assembly_fasta) }

    failed_consensus_assemblies_ch = AUTOCYCLER_COMPRESS.out.all
        .filter { sample, assembler, consensus_assembly_fasta, consensus_assembly_gfa, consensus_assembly_fol, consensus_assembly_metrics, unitigs_num -> unitigs_num.text.trim().toFloat() >= 100.0}
    PUBLISH_FAILED_CONSENSUS(failed_consensus_assemblies_ch)

    // Add in manual assemblies for those where autocycler consensus failed and manual intervention needed
    // MANUAL INPUT CHANNELS:
    Channel.fromPath( "${params.manual_assembies}/*.fasta" )
         .map{ file -> tuple(file.simpleName, 'manual', file) }
         .set{ manual_assembly_ch }

    all_assemblies_ch = successful_consensus_assemblies_ch.mix(manual_assembly_ch)

    medaka_full_polishing_input_ch =  all_assemblies_ch.combine(ont_labeled_ch, by: 0)
    MEDAKA_FULL(medaka_full_polishing_input_ch)    

    SEQKIT(MEDAKA_FULL.out)

    CHECKM2_CONSENSUS(MEDAKA_FULL.out)
 
    MLST(MEDAKA_FULL.out)
     
    KRAKEN2(MEDAKA_FULL.out)
    SPECIES_FROM_KRAKEN2(KRAKEN2.out.fasta)
 
    BAKTA(MEDAKA_FULL.out)

    amrfinder_input_ch = SPECIES_FROM_KRAKEN2.out.combine(BAKTA.out.protein, by: 0)
    AMRFINDERPLUS(amrfinder_input_ch)
        
    //Mob-suite give more comprehensive annotation than plasmidfinder- as does relaxase typer as well as replicon type, and other outputs.
    MOB_SUITE(MEDAKA_FULL.out)
    

    // COLLATE SUMMARIES

    // 1. Raw and subsampled reads QC
    all_reads_ch = RAWQC_ONT.out.mix(RAWQC_ONT_SUBSAMPLED.out).collect()
    SUMMARISE_RAWQC(all_reads_ch)

    // 1b. Contaminated reads
    all_contam_reads_ch = contam_reads_ch.collect()
    SUMMARISE_CONTAMINATED(all_contam_reads_ch)

    // 2. Input assemblies QC
    all_input_assemblies_qc_ch = CHECKM2_INPUT_ASSEMBLIES.out.tsv.collect()
    SUMMARISE_INPUT_ASSEMBLY_QC(all_input_assemblies_qc_ch)
    
     // 3a. Consensus assembly (Medaka polished) QC - Seqkit stats
    all_consensus_seqkit_ch = SEQKIT.out.collect()
    SUMMARISE_CONSENSUS_ASSEMBLY_SEQKIT(all_consensus_seqkit_ch)
 
    // 3b. Consensus assembly (Medaka polished) QC- CheckM2
    all_consensus_checkm2_ch = CHECKM2_CONSENSUS.out.collect()
    SUMMARISE_CONSENSUS_ASSEMBLY_CHECKM2(all_consensus_checkm2_ch)
 
    // 3c. Failed assemblies


    // 4. MLST
    all_mlst_ch = MLST.out.collect()
    SUMMARISE_MLST(all_mlst_ch)

    // 5. Kraken2
    all_kraken2_ch = KRAKEN2.out.collect()
    SUMMARISE_KRAKEN2(all_kraken2_ch)

    // 6. AMRFinder Plus
    all_amrfinder_ch = AMRFINDERPLUS.out.collect()
    SUMMARISE_AMRFINDER(all_amrfinder_ch)

    // 7. Bakta 
    all_bakta_ch = BAKTA.out.fol.collect()
    SUMMARISE_BAKTA(all_bakta_ch)

    // 8. MOB-suite 
    all_mobsuite_biomarkers_ch = MOB_SUITE.out.biomarkers.collect()
    all_mobsuite_contigs_ch = MOB_SUITE.out.contigs.collect()
    all_mobsuite_mges_ch = MOB_SUITE.out.mges.collect()
    all_mobsuite_mobtyper_ch = MOB_SUITE.out.mobtyper.collect()
    SUMMARISE_MOBSUITE_BIOMARKERS(all_mobsuite_biomarkers_ch)
    SUMMARISE_MOBSUITE_CONTIGS(all_mobsuite_contigs_ch)
    SUMMARISE_MOBSUITE_MGES(all_mobsuite_mges_ch)
    SUMMARISE_MOBSUITE_MOBTYPER(all_mobsuite_mobtyper_ch)
   

}
