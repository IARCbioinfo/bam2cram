#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run iarcbioinfo/bam2cram-nf --bams input_folder --fasta ref.fa --fai ref.fa.fai -profile singularity
    Mandatory arguments:
      --bams [directory]                input folder with BAM files and indexes
    References
      --fasta [file]                  Path to fasta reference to encode the CRAM file
      --fai    [file]                 Path to fasta reference index.
    Input alternatives:
      --bam_csv                    file with tabular data for each sample to process [label bam index ]
      -profile [str]              Configuration profile to use.
                                  Available: singularity
    Output:
    --output_folder                ouput folder [def: ./results]
    """.stripIndent()
}

// we star coding the pipeline
params.help = null
params.bams = null
params.bam_csv = null
params.fasta= null
params.fai = null

// Show help message
if (params.help) exit 0, show_help()

//print header and tool name
log.info IARC_Header()
log.info tool_header()

//check fasta and fasta index
if (!params.fasta || !params.fai) exit 1, "The reference fasta file with its index need to be provided!"

//we check the reference
ch_fasta = Channel.value(file(params.fasta)).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}
ch_fai =   Channel.value(file(params.fai)).ifEmpty{exit 1, "fai index file not found: ${params.fai}"}
//check that BAM input is provided!
if(!params.bam_csv && !params.bams) exit 1, "No --bam_csv nor --bams options provided!"

//we print the parameters
log.info "\n"
log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
log.info "-\033[2m--------------------------------------------------------\033[0m-"
log.info "\n"

if(params.bam_csv) {
      Channel.fromPath(file(params.bam_csv)).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row.label, file(row.bam), file(row.index)]}
                      .ifEmpty{exit 1, "params.bams_csv was empty - no input files supplied" }
                      .into { inputbams; bamstats; sizebams }

}else{
  if(file(params.bams).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        bams = Channel.fromPath( params.bams+'/*.bam' )
                .map {path -> [ path.name.replace(".bam",""),path]}
        bams_index = Channel.fromPath( params.bams+'/*.bam.bai')
                .map {  path -> [ path.name.replace(".bam.bai",""), path ] }
        //we create the chanel
        files = bams.join(bams_index)
        files.into {inputbams; bamstats; sizebams}
  }else{
    println "ERROR: input folder ${params.bams} contains no BAM files"; System.exit(1)
  }
}



//map the rna-seq reads to the genome
process bam2cram{
  tag "${sample}-2cram"
  //we can remove this to don't keep the bam files
  publishDir "${params.output_folder}/CRAM", mode: 'copy'

  input:
      set val(sample), file(bam), file(bam_index) from inputbams
      file(ref) from ch_fasta
      file(ref_index) from ch_fai
  output:
      //output CRAM files
      set val(sample), file("${sample}.cram"), file("${sample}.cram.crai") into cramfiles, cramstats, sizecrams

  script:
  """
  samtools view -C  -T ${ref} ${bam} -o ${sample}.cram
  samtools index ${sample}.cram
  """
}

//we compute some stat for original bam files
process stats_bams{
  tag "${sample}-bam_stats"
  publishDir "${params.output_folder}/qc/bam", mode: 'copy'
  input:
    set val(sample), file(bam), file(bam_index) from bamstats
  output:
    set val(sample), file("${sample}.bam.flagstat"), file("${sample}.bam.stats") into bam_qc
  script:
  """
  samtools flagstat ${bam} > ${sample}.bam.flagstat
  samtools stats ${bam} > ${sample}.bam.stats
  """
}

//we compute some stat for original bam files
process stats_crams{
  tag "${sample}-cram_stats"

  publishDir "${params.output_folder}/qc/cram", mode: 'copy'
  input:
    set val(sample), file(cram), file(cram_index) from cramstats
  output:
    set val(sample), file("${sample}.cram.flagstat"), file("${sample}.cram.stats") into cram_qc
  script:
  """
  samtools flagstat ${cram} > ${sample}.cram.flagstat
  samtools stats ${cram} > ${sample}.cram.stats
  """
}

//we join bam and cram qc output to compare them
cram_bam_size=sizecrams.join(sizebams)
cram_bam_qc=cram_qc.join(bam_qc)

//we merge both cram_qc and bam_qc to check if they are identical
process check_conversion{
  tag "${sample}-check"
  publishDir "${params.output_folder}/qc/check", mode: 'copy'

  input:
    set val(sample), file(c_fs), file(c_stat), file(b_fs), file (b_stat) from cram_bam_qc
    set val(sample), file(cram), file(cram_index), file(bam), file (bam_index) from cram_bam_size
    //file qc_check from cram_bam_qc.collect()
  output:
    file("${sample}_check.report.txt") into report_qc
   script:
   """
   if diff ${c_fs} ${b_fs} > /dev/null
   then
         fs_test="OK"
   else
         fs_test="fail"
    fi
    #we remove the line of stat that is diferent because of file extension
    cat ${c_stat} | grep -v "# The command line was:" > ${c_stat}.no_cmd_line
    cat ${b_stat} | grep -v "# The command line was:" > ${b_stat}.no_cmd_line
    if diff ${c_stat}.no_cmd_line ${b_stat}.no_cmd_line > /dev/null
    then
          s_test="OK"
    else
          s_test="fail"
     fi

     #get file size
     c_size=`du -Hh ${cram} |cut -f1`
     b_size=`du -Hh ${bam}| cut -f1`
   echo "ID\tflagstat\tstats\tBAM_size\tCRAM_size" > ${sample}_check.report.txt
   echo "${sample}\t\$fs_test\t\$s_test\t\$b_size\t\$c_size" >> ${sample}_check.report.txt
   """
}

report_qc.collectFile(name: 'bam2cram_summary.txt', storeDir: params.output_folder, seed: 'ID\tflagstat\tstats\tBAM_size\tCRAM_size\n', newLine: false, skip: 1)


//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        BAM\u001b[32;1m 2\u001b[33;1m CRAM\u001b[31;1m (${workflow.manifest.version})
        """
}
