
//map the rna-seq reads to the genome

process bam2cram{
  
  tag "${sample}-2cram"
  //we can remove this to don't keep the bam files
  publishDir "${params.outdir}/CRAM", mode: 'copy'

  input:
      tuple val(sample), path(bam), path(bam_index)
      path ref
      path ref_index
  
  output:
      tuple val(sample), file("${sample}.cram"), file("${sample}.cram.crai"), emit: crams

  script:

  if(params.debug == true){
    """
     echo samtools view -C  -T ${ref} ${bam} -o ${sample}.cram
     echo samtools index ${sample}.cram
     touch ${sample}.cram ${sample}.cram.crai
    """
  }else{
  """
      samtools view -C  -T ${ref} ${bam} -o ${sample}.cram
      samtools index ${sample}.cram
  """
  }
}


//we compute some stat for original bam files
process stats_bams{
  tag "${sample}-bam_stats"
  publishDir "${params.outdir}/qc/bam", mode: 'copy'
  input:
    tuple val(sample), file(bam), file(bam_index)
  output:
    tuple val(sample), file("${sample}.bam.flagstat"), file("${sample}.bam.stats"), emit: bam_qc
  script:
  if(params.debug == true){
  """
  echo samtools flagstat ${bam} > ${sample}.bam.flagstat
  echo samtools stats ${bam} > ${sample}.bam.stats
  touch ${sample}.bam.flagstat ${sample}.bam.stats
  """
  }else{
  """
  samtools flagstat ${bam} > ${sample}.bam.flagstat
  samtools stats ${bam} > ${sample}.bam.stats
  """
  }
}


//we compute some stat for original bam files
process stats_crams{
  tag "${sample}-cram_stats"

  publishDir "${params.outdir}/qc/cram", mode: 'copy'
  input:
    tuple val(sample), file(cram), file(cram_index) 
  output:
    tuple val(sample), file("${sample}.cram.flagstat"), file("${sample}.cram.stats"), emit: cram_qc
  script:
  if(params.debug ==  true){
  """
  echo samtools flagstat ${cram} > ${sample}.cram.flagstat
  echo samtools stats ${cram} > ${sample}.cram.stats
  touch ${sample}.cram.flagstat
  touch ${sample}.cram.stats
  """
    
    }else{
  """
  samtools flagstat ${cram} > ${sample}.cram.flagstat
  samtools stats ${cram} > ${sample}.cram.stats
  """
 }
}


//we join bam and cram qc output to compare them
//cram_bam_size=sizecrams.join(sizebams)
//cram_bam_qc=cram_qc.join(bam_qc)

//we merge both cram_qc and bam_qc to check if they are identical
process check_conversion{
  tag "${sample}-check"
  publishDir "${params.outdir}/qc/check", mode: 'copy'

  input:
    tuple val(sample), file(c_fs), file(c_stat), file(b_fs), file (b_stat)
    tuple val(sample), file(cram), file(cram_index), file(bam), file (bam_index)
    
  output:
    
    path("${sample}_check.report.txt") , emit: report_qc
   
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

//report_qc.collectFile(name: 'bam2cram_summary.txt', storeDir: params.output_folder, seed: 'ID\tflagstat\tstats\tBAM_size\tCRAM_size\n', newLine: false, skip: 1)

workflow{


// Show help message
if (params.help) exit 0, show_help()


//print header and tool name
 DGL_Header()
 tool_header()

//check fasta and fasta index
if (!params.fasta || !params.fai) exit 1, "The reference fasta file with its fai index have to be provided!"

//we check the reference
//ch_fasta = Channel.fromPath(params.fasta).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}
//ch_fai =   Channel.fromPath(params.fai).ifEmpty{exit 1, "fai index file not found: ${params.fai}"}

ch_fasta = file(params.fasta)
ch_fai =   file(params.fai)
//check that BAM input is provided!
if(!params.bam_csv && !params.bams) exit 1, "No --bam_csv nor --bams options provided!"

//we print the parameters
log.info "\n"
log.info "-------------------Calling PARAMETERS--------------------"
log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
log.info "--------------------------------------------------------"
log.info "\n"
if(params.bam_csv){
read_bams_ch=Channel.fromPath(params.bam_csv) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, file(row.bam), file(row.bai))}
}
  B2C=bam2cram(read_bams_ch,ch_fasta,ch_fai)
  BMS=stats_bams(read_bams_ch)
  CMS=stats_crams(B2C.crams)
  //we check conversion
  CRC=check_conversion(CMS.join(BMS),B2C.crams.join(read_bams_ch))
  //we collect all stats
  CRC.report_qc.collectFile(name: 'bam2cram_summary.txt', storeDir: params.outdir, seed: 'ID\tflagstat\tstats\tBAM_size\tCRAM_size\n', newLine: false, skip: 1)
}


//help function for the tool
def show_help (){
  DGL_Header()
  tool_header()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run digenomalab/bam2cram --bam_csv bams.csv --fasta ref.fa --fai ref.fa.fai -profile kutral
    Mandatory arguments:
      --bam_csv                    file with tabular data for each sample to process [label bam index ]

    References
      --fasta [file]                  Path to fasta reference to encode the CRAM file
      --fai    [file]                 Path to fasta reference index.
    Input alternatives:
      -profile [str]              Configuration profile to use.
                                  Available: kutral
    Output:
      --outdir                ouput folder [def: ./results]
    """.stripIndent()
}


// header for the Di Genoma Lab tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def DGL_Header (){

log.info"""
#########################################################################
#  _____  _  _____                                  _           _       #
# |  __ \\(_)/ ____|                                | |         | |      #
# | |  | |_| |  __  ___ _ __   ___  _ __ ___   __ _| |     __ _| |__    #
# | |  | | | | |_ |/ _ \\ '_ \\ / _ \\| '_ ` _ \\ / _` | |    / _` | '_ \\   #
# | |__| | | |__| |  __/ | | | (_) | | | | | | (_| | |___| (_| | |_) |  #
# |_____/|_|\\_____|\\___|_| |_|\\___/|_| |_| |_|\\__,_|______\\__,_|_.__/   #
# Workflows for genomics.################################################
""".stripIndent()
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header(){
        log.info"""
        BAM\u001b[32;1m 2\u001b[33;1m CRAM\u001b[31;1m (${workflow.manifest.version})
        """.stripIndent()
}
