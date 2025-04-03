
process CALL_INTEGRATION_SITES {
    tag "Calling integration sites: ${sample_id}"
    publishDir \
        "${params.outdir}/integration_sites", \
        mode: 'copy'
    label "process_low"
    label "isatoolkit"

    input:
    tuple \
        val(sample_id), \
        path(filtered_bam), \
        path(filtered_bai)
    
    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.integration_sites.bed"), \
        emit: integration_sites_bed
    
    script:
    """
    isatools sam count \\
            -i ${filtered_bam} |
        isatools bed merge \\
            -d ${params.grouping_distance} \\
            --mode ${params.grouping_mode} |
        isatools bed sort \\
            -o ${sample_id}.integration_sites.bed
    """
}

process AGGREGATE_INTEGRATION_SITES {
    tag "Aggregating integration sites"
    publishDir \
        "${params.outdir}/integration_sites", \
        mode: 'copy'
    label "process_low"

    input:
    path(integration_site_beds)

    output:
    path("aggregated_integration_sites.bed")

    script:
    def file_suffix = '.integration_sites.bed'
    """
    mkdir 'integration_site_csvs'

    for f in ${integration_site_beds}; do
        SAMPLE_NAME=\$(basename \${f} ${file_suffix})
        awk \\
            -F'\\t' -v OFS=',' -v SAMPLE_NAME=\$SAMPLE_NAME \\
            '{\$1=\$1; print SAMPLE_NAME,\$0}' \\
            "\$f" \\
            > "integration_site_csvs/\${SAMPLE_NAME}.csv"
    done

    find \\
        integration_site_csvs -name '*.csv' \\
        -exec csvtk concat \\
            -H \\
            -o 'aggregated_integration_sites.bed' \\
            {} +
    """
}