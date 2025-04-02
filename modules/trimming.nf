process EXTRACT_UMIS {
    tag "Extracting UMIs: ${sample_id}"
    publishDir \
        "${params.outdir}/trimmed_fastqs/umi_extracted", \
        mode: 'copy'
    label "process_large"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_R1.umi_extracted.fastq.gz"), \
        path("${sample_id}_R2.umi_extracted.fastq.gz"), \
        emit: umi_extracted_fastqs
    path "${sample_id}.umi_extracted.log"

    script:
    // Use a regex to capture three parts:
    //  1. Everything before the first group of N’s (discard1)
    //  2. The continuous N’s (umi)
    //  3. Everything after the N’s (discard2)
    def matcher = (params.R2_adapter =~ /^([ATGC]+)(N+)([ATGC]+)$/)
    if (!matcher.matches()) {
        throw new IllegalArgumentException("Input string format is not as expected.")
    }
    def discard1 = matcher[0][1]
    def umi = matcher[0][2]
    def discard2 = matcher[0][3]
    // Build the final regex pattern for UMI-tools, using the length of the UMI (number of N’s)
    def umitools_regexp = \
        "(?P<discard_1>" + \
        discard1 + \
        "){e<=3}(?P<umi_1>.{" + \
        umi.length() + \
        "})(?P<discard_2>" + \
        discard2 + \
        "){e<=3}"
    """
    umi_tools extract \\
        --log="${sample_id}.umi_extracted.log" \\
        --extract-method='regex' \\
        --bc-pattern2="${umitools_regexp}" \\
        -I "${fastq_files[0]}" \\
        -S "${sample_id}_R1.umi_extracted.fastq.gz" \\
        --read2-in="${fastq_files[1]}" \\
        --read2-out="${sample_id}_R2.umi_extracted.fastq.gz"
    """
}

process TRIM_FASTQS {
    tag "Trimming: ${sample_id}"
    publishDir \
        "${params.outdir}/trimmed_fastqs/r1_adapter_trimmed", \
        mode: 'copy', pattern: "*.r1_adapter_*trimmed.*"
    publishDir \
        "${params.outdir}/trimmed_fastqs/ir_trimmed", \
        mode: 'copy', pattern: "*.ir_*trimmed.*"
    publishDir \
        "${params.outdir}/trimmed_fastqs/adapter_readthrough_trimmed", \
        mode: 'copy', pattern: "*.adapter_readthrough_trimmed.*"
    publishDir \
        "${params.outdir}/trimmed_fastqs/ir_readthrough_trimmed", \
        mode: 'copy', pattern: "*.ir_readthrough_trimmed.*"
    label "process_medium"

    input:
    tuple val(sample_id), path(r1_fastq), path(r2_fastq)

    output:
    // R1 adapter trimmed
    tuple \
        val(sample_id), \
        path("${sample_id}_R1.r1_adapter_trimmed.fastq.gz"), \
        path("${sample_id}_R2.r1_adapter_trimmed.fastq.gz"), \
        emit: r1_trimmed_fastqs
    path "${sample_id}.r1_adapter_trimmed.log"
    path "${sample_id}_R*.r1_adapter_untrimmed.fastq.gz"

    // IR trimmed
    tuple \
        val(sample_id), \
        path("${sample_id}_R1.ir_trimmed.fastq.gz"), \
        path("${sample_id}_R2.ir_trimmed.fastq.gz"), \
        emit: ir_trimmed_fastqs
    path "${sample_id}.ir_trimmed.log"
    path "${sample_id}_R*.ir_untrimmed.fastq.gz"

    // Adapter readthrough trimmed
    tuple \
        val(sample_id), \
        path("${sample_id}_R1.adapter_readthrough_trimmed.fastq.gz"), \
        path("${sample_id}_R2.adapter_readthrough_trimmed.fastq.gz"), \
        emit: adapter_readthrough_trimmed_fastqs
    path "${sample_id}.adapter_readthrough_trimmed.log"

    // IR readthrough trimmed
    tuple \
        val(sample_id), \
        path("${sample_id}_R1.ir_readthrough_trimmed.fastq.gz"), \
        path("${sample_id}_R2.ir_readthrough_trimmed.fastq.gz"), \
        emit: ir_readthrough_trimmed_fastqs

    script:
    """
    # R1 adapter trimming
    cutadapt \\
        -j "${task.cpus}" \\
        -g "^${params.R1_adapter}" \\
        -e "${params.allowed_errors_r1_adapter}" \\
        -m "${params.min_length}" \\
        -o "${sample_id}_R1.r1_adapter_trimmed.fastq.gz" \\
        -p "${sample_id}_R2.r1_adapter_trimmed.fastq.gz" \\
        --untrimmed-output "${sample_id}_R1.r1_adapter_untrimmed.fastq.gz" \\
        --untrimmed-paired-output "${sample_id}_R2.r1_adapter_untrimmed.fastq.gz" \\
        "${r1_fastq}" \\
        "${r2_fastq}" \\
        > "${sample_id}.r1_adapter_trimmed.log"
    
    # IR trimming
    cutadapt \\
        -j "${task.cpus}" \\
        -g "^${params.IR_seq}" \\
        -e "${params.allowed_errors_ir}" \\
        -m "${params.min_length}" \\
        -o "${sample_id}_R1.ir_trimmed.fastq.gz" \\
        -p "${sample_id}_R2.ir_trimmed.fastq.gz" \\
        --untrimmed-output "${sample_id}_R1.ir_untrimmed.fastq.gz" \\
        --untrimmed-paired-output "${sample_id}_R2.ir_untrimmed.fastq.gz" \\
        "${sample_id}_R1.r1_adapter_trimmed.fastq.gz" \\
        "${sample_id}_R2.r1_adapter_trimmed.fastq.gz" \\
        > "${sample_id}.ir_trimmed.log"
    
    # Adapter readthrough trimming
    cutadapt \\
        -j "${task.cpus}" \\
        -a "\$(echo ${params.R2_adapter} | tr 'ATGCN' 'TACGN' | rev)" \\
        -A "\$(echo ${params.R1_adapter} | tr 'ATGCN' 'TACGN' | rev)" \\
        -m "${params.min_length}" \\
        -e "${params.allowed_errors_adapter_readthrough}" \\
        -o "${sample_id}_R1.adapter_readthrough_trimmed.fastq.gz" \\
        -p "${sample_id}_R2.adapter_readthrough_trimmed.fastq.gz" \\
        "${sample_id}_R1.ir_trimmed.fastq.gz" \\
        "${sample_id}_R2.ir_trimmed.fastq.gz" \\
        > "${sample_id}.adapter_readthrough_trimmed.log"

    # IR readthrough trimming
    cutadapt \
        -j "${task.cpus}" \
        -A "\$(echo ${params.IR_seq} | tr 'ATGCN' 'TACGN' | rev)" \
        -m "${params.min_length}" \
        -e "${params.allowed_errors_ir_readthrough}" \
        -o "${sample_id}_R1.ir_readthrough_trimmed.fastq.gz" \
        -p "${sample_id}_R2.ir_readthrough_trimmed.fastq.gz" \
        "${sample_id}_R1.adapter_readthrough_trimmed.fastq.gz" \
        "${sample_id}_R2.adapter_readthrough_trimmed.fastq.gz" \
        > "${sample_id}.ir_readthrough_trimmed.log"
    """
}