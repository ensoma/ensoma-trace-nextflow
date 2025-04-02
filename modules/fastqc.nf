
process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'
    label "process_medium"

    input:
    path fastq_files

    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    script:
    fastq_files_arg = fastq_files.join(' ')
    """
    fastqc \
        --outdir . \
        --threads ${task.cpus - 1} \
        --quiet \
        ${fastq_files_arg}
    """
}