
process CONCAT_REFERENCE_FASTAS {
    tag "Concatenating References: ${payload_fa.name} ${integration_fa.name} ${genome_fa.name}"
    publishDir "${params.outdir}/genome/concatenated", mode: 'copy'
    label "process_low"
    label "ubuntu"

    input:
    tuple path(payload_fa), path(integration_fa), path(genome_fa)

    output:
    tuple \
        val("${payload_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}_${integration_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}_${genome_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}"), \
        path("${payload_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}_${integration_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}_${genome_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")}.fa"),
        emit: concatenated_refs

    script:
    def payload_name = payload_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
    def integration_name = integration_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
    def genome_name = genome_fa.name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
    """
    # Unzip FASTAs if required, and standardize the names

    ## If the payload vector FASTA is gzipped, gunzip it
    if [[ "${payload_fa.name}" = *.gz ]]; then
        gunzip -c ${payload_fa} > ${payload_name}.fa
        payload_file=${payload_name}.fa
    else
        payload_file=${payload_fa}
    fi

    ## If the integration FASTA is gzipped, gunzip it
    if [[ "${integration_fa.name}" = *.gz ]]; then
        gunzip -c ${integration_fa} > ${integration_name}.fa
        integration_file=${integration_name}.fa
    else
        integration_file=${integration_fa}
    fi

    ## If the genome FASTA is gzipped, gunzip it
    if [[ "${genome_fa.name}" = *.gz ]]; then
        gunzip -c ${genome_fa} > ${genome_name}.fa
        genome_file=${genome_name}.fa
    else
        genome_file=${genome_fa}
    fi

    # Ensure the FASTA files end in a newline

    ## Ensure the payload vector FASTA ends in a newline
    sed -i '\$a\\' "\${payload_file}"

    ## Ensure the integration FASTA ends in a newline
    sed -i '\$a\\' "\${integration_file}"

    ## Ensure the genome FASTA ends in a newline
    sed -i '\$a\\' "\${genome_file}"

    # Concatenate the FASTA files
    cat \
        "\${payload_file}" \
        "\${integration_file}" \
        "\${genome_file}" \
        > ${payload_name}_${integration_name}_${genome_name}.fa
    """
}