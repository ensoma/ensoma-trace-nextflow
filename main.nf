nextflow.enable.dsl=2

// Parameter validation
include { 
    validateParameters;
    paramsSummaryLog;
    samplesheetToList
} from 'plugin/nf-schema'

validateParameters()
log.info paramsSummaryLog(workflow)

// Create a FASTQ file channel
Channel
    .fromFilePairs( "${params.fastq_dir}/*_R{1,2}_*.fastq.gz" )
    // Remove the remaining Illumina part of the filename
    .map { it -> 
        def sample_id = it[0].split('_')[0]
        return [sample_id, it[1]]
    }
    .set { fastq_ch }

// Create a sample sheet channel
Channel
    .fromList(
        samplesheetToList(
            file(params.sample_sheet),
            './assets/schemas/sample_sheet_schema.json'
        )
    )
    .map { it -> 
        it[3] = file("${params.genome_dir}/${it[3]}")
        it[4] = file("${params.genome_dir}/${it[4]}")
        it[5] = file("${params.genome_dir}/${it[5]}")
        return it
    }
    .set { sample_sheet_ch }

include { FASTQC } from './modules/fastqc.nf'
include { CONCAT_REFERENCE_FASTAS } from './modules/reference_preprocessing.nf'
include { EXTRACT_UMIS; TRIM_FASTQS } from './modules/trimming.nf'
include { BWA_INDEX; BWA_ALIGN } from './modules/alignment.nf'
include { FILTER_ALIGNMENTS } from './modules/alignment_filtering.nf'
include { CALL_INTEGRATION_SITES; AGGREGATE_INTEGRATION_SITES } from './modules/integration_sites.nf'

// Run the workflow
workflow {

    // // FastQC of the raw FASTQ files
    // fastq_ch
    //     .map { it[1] }
    //     .collect() |
    //     FASTQC

    // Extract the UMIs from the FASTQ files
    EXTRACT_UMIS( fastq_ch )

    // Trim the FASTQ files
    TRIM_FASTQS( EXTRACT_UMIS.out.umi_extracted_fastqs )

    // Concatenate the reference FASTA files
    sample_sheet_ch
        .map { it[ 3..5 ] }
        .unique() |
        CONCAT_REFERENCE_FASTAS
    
    // Generate the genome indices
    BWA_INDEX( CONCAT_REFERENCE_FASTAS.out.concatenated_refs )

    // Align the FASTQ files to the reference genome
    sample_sheet_ch
        .map { it -> 
            def sample_id = it[0]
            
            def payload_name = it[3].name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
            def integration_name = it[4].name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
            def genome_name = it[5].name.replaceAll(/\.f(ast)?a(\.gz)?$/, "")
            def concatenated_ref_name = "${payload_name}_${integration_name}_${genome_name}"

            return [ concatenated_ref_name, sample_id ]
        }
        .combine( BWA_INDEX.out.index, by: 0 )
        .map { it.swap(0, 1) }
        .join( TRIM_FASTQS.out.ir_readthrough_trimmed_fastqs ) |
        BWA_ALIGN
    
    // Filter the alignments
    FILTER_ALIGNMENTS( BWA_ALIGN.out.unfiltered_sam )

    // Call the integration sites
    FILTER_ALIGNMENTS.out.filtered_bams
        .map { it ->
            def sample_id = it[0]
            def filtered_bam = it[3].find { it.name.endsWith('.bam') }
            def filtered_bai = it[3].find { it.name.endsWith('.bai') }

            return [ sample_id, filtered_bam, filtered_bai ]
        } |
        CALL_INTEGRATION_SITES

    // Aggregate the integration sites
    CALL_INTEGRATION_SITES.out.integration_sites_bed
        .map { it[1] }
        .collect() |
        AGGREGATE_INTEGRATION_SITES
}