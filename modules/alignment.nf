
process BWA_INDEX {
    tag "Indexing: ${concatenated_ref}"
    publishDir \
        "${params.outdir}/genome/indices/${concatenated_ref_name}/", \
        mode: 'copy'
    label "bwa"
    label "process_medium"
    label "process_high_memory"

    input:
    tuple \
        val(concatenated_ref_name), \
        path(concatenated_ref)
    
    output:
    tuple \
        val(concatenated_ref_name), \
        path("${concatenated_ref_name}/${concatenated_ref_name}.*"), \
        emit: index

    script:
    """
    mkdir -p ${concatenated_ref_name}

    bwa-mem2 index \
        -p ${concatenated_ref_name}/${concatenated_ref_name} \
        ${concatenated_ref}
    """
}

process BWA_ALIGN {
    tag "Aligning: ${sample_id}"
    publishDir \
        "${params.outdir}/alignments/unfiltered/", \
        mode: 'copy'
    label "bwa"
    label "process_large"

    input:
    tuple \
        val(sample_id), \
        val(concatenated_ref_name), \
        path(ref_index), \
        path(trimmed_r1_fastq), \
        path(trimmed_r2_fastq)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.sam"), \
        emit: unfiltered_sam

    script:
    """
    mkdir bwa_index
    for index_file in ${ref_index}; do
        mv \$index_file bwa_index/
    done

    bwa-mem2 mem \\
        -o ${sample_id}.sam \\
        -t ${task.cpus} \\
        "bwa_index/${concatenated_ref_name}" \\
        ${trimmed_r1_fastq} \\
        ${trimmed_r2_fastq}
    """
}