process SUMMARISE_RAWQC {

    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "read_qc_merged.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "contaminated_reads_checkm2_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "input_assemblies_checkm2_merged.tsv"

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


process SUMMARISE_AUTOCYCLER_METRICS {

    tag { "autocycler_metrics_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "autocycler_metrics_summary.tsv"

    input:
    val(input_tuples)

    output:
    path("autocycler_metrics_summary.tsv")

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
    head -n 1 "${input_tuples[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > autocycler_metrics_summary.tsv
    {
    ${tmp_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> autocycler_metrics_summary.tsv
    """

    return script
}



process SUMMARISE_CONSENSUS_ASSEMBLY_SEQKIT {

    tag { "consensus_seqkit_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "consensus_assembly_seqkit_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "consensus_assembly_checkm2_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "mlst_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "kraken2_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "amrfinder_summary.tsv"

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
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "bakta_summary.tsv"

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



process SUMMARISE_MOBSUITE_BIOMARKERS {

    tag { "mobsuite_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "*mobsuite*.tsv"

    input:
    val(mobsuite_input)

    output:
    path("mobsuite_biomarkers_summary.tsv")

    script:
    def n = mobsuite_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def biomarkers_lines = []

    for (int i = 0; i < n; i += 3) {
        def sample = mobsuite_input[i]
        def assembler = mobsuite_input[i + 1]
        def biomarker = mobsuite_input[i + 2].toString().replace('$', '\\$')


        biomarkers_lines << """tail -n +2 "${biomarker}" | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
     }

    def script = """
    head -n 1 "${mobsuite_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_biomarkers_summary.tsv
    {
    ${biomarkers_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_biomarkers_summary.tsv

    """
    return script
}


process SUMMARISE_MOBSUITE_CONTIGS {

    tag { "mobsuite_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "*mobsuite*.tsv"

    input:
    val(mobsuite_input)

    output:
    path("mobsuite_contigs_summary.tsv")

    script:
    def n = mobsuite_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def contigs_lines = []

    for (int i = 0; i < n; i += 3) {
        def sample = mobsuite_input[i]
        def assembler = mobsuite_input[i + 1]
        def contig = mobsuite_input[i + 2].toString().replace('$', '\\$')

        contigs_lines    << """tail -n +2 "${contig}"    | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """

    head -n 1 "${mobsuite_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_contigs_summary.tsv
    {
    ${contigs_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_contigs_summary.tsv

    """
    return script
}



process SUMMARISE_MOBSUITE_MGES {

    tag { "mobsuite_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "*mobsuite*.tsv"

    input:
    val(mobsuite_input)

    output:
    path("mobsuite_mge_summary.tsv")

    script:
    def n = mobsuite_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def mge_lines = []

    for (int i = 0; i < n; i += 3) {
        def sample = mobsuite_input[i]
        def assembler = mobsuite_input[i + 1]
        def mge = mobsuite_input[i + 2].toString().replace('$', '\\$')

        mge_lines        << """tail -n +2 "${mge}"       | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""

    }

    def script = """

    head -n 1 "${mobsuite_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_mge_summary.tsv
    {
    ${mge_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_mge_summary.tsv

    """
    return script
}


process SUMMARISE_MOBSUITE_MOBTYPER {

    tag { "mobsuite_summary" }
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "*mobsuite*.tsv"

    input:
    val(mobsuite_input)

    output:
    path("mobsuite_mobtyper_results.tsv")

    script:
    def n = mobsuite_input.size()
    if (n % 3 != 0)
        throw new IllegalArgumentException("Expected flat list of 3-tuples, but got ${n} items.")

    def mobtyper_lines = []

    for (int i = 0; i < n; i += 3) {
        def sample = mobsuite_input[i]
        def assembler = mobsuite_input[i + 1]
        def mobtyper = mobsuite_input[i + 2].toString().replace('$', '\\$')

        mobtyper_lines   << """tail -n +2 "${mobtyper}"  | awk -v s="${sample}" -v a="${assembler}" 'BEGIN{OFS="\\t"} {print s, a, \$0}'"""
    }

    def script = """

    head -n 1 "${mobsuite_input[2]}" | awk 'BEGIN{OFS="\\t"} {print "sample", "assembler", \$0}' > mobsuite_mobtyper_results.tsv
    {
    ${mobtyper_lines.join('\n')}
    } | sort -k2,2 -k1,1 >> mobsuite_mobtyper_results.tsv
    """
    return script
}




























