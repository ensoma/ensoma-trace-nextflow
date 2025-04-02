
process FILTER_ALIGNMENTS {
    tag "Filtering alignments: ${sample_id}"
    publishDir \ 
        "${params.outdir}/alignments/flag_filtered/", \ 
        mode: 'copy', pattern: "*.flag_filtered.*"
    publishDir \
        "${params.outdir}/alignments/deduped/", \ 
        mode: 'copy', pattern: "*.deduped.*"
    publishDir \
        "${params.outdir}/alignments/filtered_final/", \ 
        mode: 'copy', pattern: "*.filtered_final.*"
    publishDir \
        "${params.outdir}/alignments/isatools_discarded/", \ 
        mode: 'copy', pattern: "*_discarded.*"
    label "process_medium"

    input:
    tuple \
        val(sample_id), \
        path(unfiltered_sam)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.flag_filtered.bam*"), \
        path("${sample_id}.deduped.bam*"), \
        path("${sample_id}.filtered_final.bam*"), \
        emit: filtered_bams
    tuple \
        val(sample_id), \
        path("${sample_id}.alt_sup_discarded.bam"), \
        path("${sample_id}.softclip_discarded.bam"), \
        path("${sample_id}.orphans_discarded.bam"), \
        emit: isatools_discarded

    script:
    def max_memory = '2G'
    """
    # 1) Filter based on flags and quality scores.
    # 2) Then position sort.
    samtools view \\
            -@ ${task.cpus} \\
            -u \\
            -F ${params.exclude_flags} \\
            -f ${params.required_flags} \\
            -q ${params.min_mapq} \\
            ${unfiltered_sam} |
        samtools sort \\
            -@ ${task.cpus} \\
            -m ${max_memory} \\
            -O 'BAM' \\
            --write-index \\
            -o "${sample_id}.flag_filtered.bam##idx##"${sample_id}.flag_filtered.bam.bai"

    # 1) Deduplicate with UMI-tools.
    # 2) Index the deduplicated BAM file.
    umi_tools dedup \\
        --paired \\
        -I "${sample_id}.flag_filtered.bam" \\
        -S "${sample_id}.deduped.bam"

    samtools index \\
        -b \\
        -o "${sample_id}.deduped.bam.bai" \\
        "${sample_id}.deduped.bam"

    # 1) Remove supplementary and alternative alignments.
    # 2) Remove too large of 5' R1 softclip.
    # 3) Remove orphaned reads (missing mate pair).
    # 4) Position sort the further filtered bam.
    isatools sam mapping-filter \\
            -i "${sample_id}.deduped.bam" \\
            -d "${sample_id}.alt_sup_discarded.bam" \\
            -u |
        isatools sam fiveprime-filter \\
            -m ${params.softclip_threshold} \\
            -d "${sample_id}.softclip_discarded.bam" \\
            -u |
        isatools sam remove-orphaned \\
            -d "${sample_id}.orphans_discarded.bam" \\
            -u |
        samtools sort \\
            -@ ${task.cpus} \\
            -O 'BAM' \\
            --write-index \\
            -m ${max_memory} \\
            -o "${sample_id}.filtered_final.bam##idx##${sample_id}.filtered_final.bam.bai"
    """
}