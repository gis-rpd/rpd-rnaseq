#!/usr/bin/env nextflow
/*
========================================================================================
                         rpd-rnaseq
========================================================================================
rpd-rnaseq Analysis Pipeline.
----------------------------------------------------------------------------------------
*/

workflow_name = "rpd-rnaseq"
log.info "======================================"
log.info " ${workflow_name}"
log.info "======================================"


def helpMessage() {
    log.info"""
    Usage:
    nextflow main.nf -params-file sample.yaml --publishdir outdir -profile gis 
    -params-file    Sample config file
    --publishDir    copies the process output files to a specified folder
    --keep_workdir  Don't delete workdir 

    Options:
    -profile        config for jobs (use 'gis' to submit jobs to Aquila or the jobs will run
                    locally if not specified)
    --singleEnd     Specifies that the input is single end reads
    --trimming      Trimming is required (Trimgalore step)
    
    Trimming options
    --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
    --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
    --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
    --three_prime_clip_r2 [int]   Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
    --saveTrimmed                 Save trimmed FastQ file intermediates

    star mapping options
    --outFilterMultimapNmax       Int Number of multi mapped reads, default is unique (i.e 1)

    Rsem options
    --single_cell_prior
    --fragment_length_mean        Only for SingleEnd reads
    --fragment_length_sd          Only for SingleEnd reads

    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behaviour

    QC options:
    --skip_qc                     Skip all QC steps aside from MultiQC

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Check that Nextflow version is up to date enough
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException("Nextflow version too old: ${workflow.nextflow.version} < $params.nf_required_version")
    }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/* File setup and parameter checks
 * ----------------------------------------------------------------------
 */
if (params.publishdir == null)
    exit 1, "Missing publishdir param"
if (params.samples == null)
    exit 1, "No samples given"
log.info "List of samples: " +  params.samples.keySet()

if (params.keep_workdir) {
   log.warn "Not cleaning up work automatically"
   cleanup = false
}

/* Input validation
 * ----------------------------------------------------------------------
 */
ref_staridx = file  ( params.references.staridx )
if ( ! ref_staridx.exists())
    exit 1, "Missing star index reference: ${ref_staridx}"

star_gtf = file ( params.references.star_gtf )
if ( ! star_gtf.exists())
    exit 1, "Missing star index reference: ${ref_staridx}"

refflat = file( params.references.refflat )
if ( ! refflat.exists())
    exit 1, "Missing star index reference: ${refflat}"

gtf_bed = file( params.references.gtf_bed )
if ( ! gtf_bed.exists())
    exit 1, "Missing star index reference: ${gtf_bed}"
rsemidx = params.references.rsemidx
if ( rsemidx == null )
    exit 1, "Missing rsem index reference: ${gtf_bed}"
// Define regular variables
//Trim-galore optional values
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
//STAR options
outFilterMultimapNmax = params.outFilterMultimapNmax
//Rsem optional values
fragment_length_mean = params.fragment_length_mean
fragment_length_sd = params.fragment_length_sd
single_cell_prior = params.single_cell_prior
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

/* Channel setup
 * ----------------------------------------------------------------------
 */
sample_keys = params.samples.keySet()
def GetReadUnitKeys = { sk ->
    params.samples[sk].readunits.keySet()
}
def GetSingleRead = { sk, rk ->
    tuple(file(params.samples[sk].readunits[rk]['fq1']))
}
def GetReadPair = { sk, rk ->
    // FIXME if files don't exist, their path might be relative to the input yaml
    // see https://gist.github.com/ysb33r/5804364
    tuple(file(params.samples[sk].readunits[rk]['fq1']),
          file(params.samples[sk].readunits[rk]['fq2']))
}
if(params.singleEnd){
  Channel
      .from(sample_keys)
      .map { sk -> tuple(sk, GetReadUnitKeys(sk).collect{GetSingleRead(sk, it)}.flatten()) }
      .set { raw_fastq_ch }
  //raw_fastq_ch.subscribe { println "$it" }
} else {
  Channel
      .from(sample_keys)
      .map { sk -> tuple(sk, GetReadUnitKeys(sk).collect{GetReadPair(sk, it)}.flatten()) }
      .set { raw_fastq_ch }
  //raw_fastq_ch.subscribe { println "$it" }
}
/*
 * STEP 1 - Merge fastq per sample!
 */
process readunit_merge {
  tag { "Merge readunit fastq for each $sample_id" }
  input:
    set sample_id, file(reads) from raw_fastq_ch
  output:
    set sample_id, file("*fastq.gz") into merged_fastq_ch, fastqc_in_ch
  script:    
    if (params.singleEnd) {
      """
      ls ${reads} | grep "_R1_" | sort | xargs cat > ${sample_id}_R1.fastq.gz;
      """
    }
    else {
      """
      ls ${reads} | grep "_R1_" | sort | xargs cat > ${sample_id}_R1.fastq.gz;
      ls ${reads} | grep "_R2_" | sort | xargs cat > ${sample_id}_R2.fastq.gz;
      """
  }
}

/*
 * STEP 2 - Trim Galore!
 */
if(params.trimming) {
  process trim_galore {
    tag {"Trimgalore from $sample_id"}
    publishDir "${params.publishdir}/${sample_id}/trimgalore", mode: 'copy',
    saveAs: {filename ->
              if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
              else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
              else if (params.saveTrimmed) filename
              else null
          }
    input:
      set sample_id, file(reads) from merged_fastq_ch  
    output:
      set sample_id, file("fq.gz_count"), file("*fq.gz") into fastq_ch
      //set sample_id, file("*trimming_report.txt") into trimgalore_results
      set sample_id, file("*") into fastqc_reports
    script:
      c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
      c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
      tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
      tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
      if (params.singleEnd) {
        """
        trim_galore --fastqc  --gzip $c_r1 $tpc_r1 $reads
        count=\$(zcat *R1_trimmed.fq.gz | wc -l);
        touch fq.gz_count
        echo \${count} >> fq.gz_count
        """
      } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        count=\$(zcat *val_1.fq.gz | wc -l);
        touch fq.gz_count
        echo \${count} >> fq.gz_count
        """
      }
  }
} else {
  process check_fastq_input {
    tag {"checking input fastq from $sample_id"}   
    input:
      set sample_id, file(reads) from merged_fastq_ch  
    output:
      set sample_id, file("fq.gz_count"), file(reads) into fastq_ch
      set sample_id, file("*") into fastqc_reports
    script:
        """
        fastqc -q $reads
        count=\$(zcat *R1.fastq.gz | wc -l);
        touch fq.gz_count
        echo \${count} >> fq.gz_count
        """
  }
}

def check_reads(count){
  def read_count = 0;
  no_of_reads = count.text.trim().toInteger()
  if (no_of_reads > 400) {
      return true
    } else {
      return false
    }
}
//Filter samples with low read counts
fastq_ch
  .filter { sample_id, count, reads -> check_reads(count) }
  .into { fastq_input_ch }
//fastq_input_ch.subscribe { println "$it" }

//Filter STAR aligned bams
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}

/*
 * STEP 3 - Star mapping
 */
process star {
    tag {"star mapping for from $sample_id"}
    publishDir "${params.publishdir}/${sample_id}/STAR", mode: 'copy'
    input:
      set sample_id, file(count), file(reads) from fastq_input_ch
      file (ref_staridx)
      file (star_gtf)
    output:
      set sample_id, file("*Log.final.out"), file ("${sample_id}_Aligned.sortedByCoord.out.bam"), file ("${sample_id}_Aligned.toTranscriptome.out.bam") into star_aligned
      file "*.out" into alignment_logs
      file "*SJ.out.tab"
      file "*Log.out" into star_log
      file "*ReadsPerGene.out.tab"
      file "*Unique.str*.out.bg"
    script:
    outFilterMultimapNmax_val = outFilterMultimapNmax > 0 ? "--outFilterMultimapNmax ${outFilterMultimapNmax}" : ' --outFilterMultimapNmax 1'
      """
      STAR --genomeDir $ref_staridx \\
          --readFilesIn $reads \\
          --runThreadN ${task.cpus} \\
          --readFilesCommand zcat \\
          --runDirPerm All_RWX \\
          --outFilterType BySJout \\
          --outSAMtype BAM SortedByCoordinate \\
          --quantMode TranscriptomeSAM GeneCounts \\
          --outSAMattributes NH HI AS nM NM MD \\
          --outBAMsortingThreadN ${task.cpus} \\
          ${outFilterMultimapNmax_val} \\
          --outFilterMismatchNmax 999 \\
          --outFilterMismatchNoverLmax 0.04 \\
          --alignEndsType EndToEnd \\
          --alignSJoverhangMin 8 \\
          --alignSJDBoverhangMin 1 \\
          --alignIntronMin 20 \\
          --alignIntronMax 1000000 \\
          --alignMatesGapMax 1000000 \\
          --limitBAMsortRAM 2001634664 \\
          --sjdbScore 1 \\
          --outWigType bedGraph \\
          --outFileNamePrefix ${sample_id}_ \\
      """
}
// Filter removes all 'aligned' channels that fail the check
star_aligned
  .filter { sample_id, logs, sorted_bam, ranscriptome_bam -> check_log(logs) }
  .into { collectRNAMetrice; bam_rseqc; bam_rsem; bam_createBigWig}
//collectRNAMetrice.subscribe { println "$it" }

/*
 * STEP 5 - picard collectRnaSeqMetrics analysis
 */ 
process collectRnaSeqMetrics {
  tag {"picard collectRNAMetrice for from $sample_id"}
  publishDir "${params.publishdir}/${sample_id}/RnaSeqMetrics", mode: 'copy'
  input:
    set sample_id, file(log), file (sorted_bam), file(ranscriptome_bam) from collectRNAMetrice
    file (refflat)
  output:
    set sample_id, file("${sample_id}_RNA_Metrics.txt") into picard_metrics
  script:
    """
    picard -Dsamjdk.compression_level=2  -XX:-UseGCOverheadLimit -Xms4000m -Xmx${task.memory.toGiga()}G \
      -XX:ConcGCThreads=${task.cpus} -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=${task.cpus} \
      CollectRnaSeqMetrics \
      I=${sorted_bam} \
      O=${sample_id}_RNA_Metrics.txt \
      REF_FLAT=${refflat} \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    """
}

/*
 * STEP 6 - RSeQC analysis
 */ 
process rseqc {
  tag {"rseqc for from $sample_id"}
  publishDir "${params.publishdir}/${sample_id}/rseqc", mode: 'copy',
  saveAs: {filename ->
            if (filename.indexOf("bam_stat.txt") > 0)    "bam_stat/$filename"
            else if (filename.indexOf("_distribution.txt") > 0)     "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else filename
        }
  input:
    set sample_id, file(log), file (bam_rseqc), file(ranscriptome_bam) from bam_rseqc
    file (gtf_bed)
  output:
    file "*.{txt,pdf,r,xls}" into rseqc_results
  when:
    ! params.skip_qc
  script:
    """
    samtools index $bam_rseqc
    #infer_experiment.py -i $bam_rseqc -r $gtf_bed > ${bam_rseqc.baseName}.infer_experiment.txt
    #junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $gtf_bed
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    #junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $gtf_bed 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    #inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $gtf_bed
    read_distribution.py -i $bam_rseqc -r $gtf_bed > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    """
}

/*
 * Step 6.1 Rseqc create BigWig coverage
 */

process createBigWig {
    tag {"createBigWig for $sample_id"}
    publishDir  "${params.publishdir}/${sample_id}/bigwig", mode: 'copy'

    when:
    !params.skip_qc

    input:
    set sample_id, file(log), file (bam_createBigWig), file(ranscriptome_bam) from bam_createBigWig

    output:
    file "*.bigwig" into bigwig_for_genebody

    script:
    """
    samtools index $bam_createBigWig
    bamCoverage -b $bam_createBigWig -p ${task.cpus} -o ${sample_id}.bigwig
    """
}

/*
 * STEP 7 - RSEM analysis
 */ 
process rsem {
  tag {"rsem for from $sample_id"}
  publishDir "${params.publishdir}/${sample_id}/rsem", mode: 'copy'
input:
    set sample_id, file(log), file (bam_rseqc), file(ranscriptome_bam) from bam_rsem
    val (rsemidx)
    output:
    set sample_id, file("${sample_id}.isoforms.results"), file("${sample_id}.genes.results")into rsem_results
    script:
    if (params.singleEnd) {
      fragment_length_mean_val = fragment_length_mean > 0 ? "--fragment-length-mean ${fragment_length_mean}" : ''
      fragment_length_sd_val = fragment_length_sd > 0 ? "--fragment-length-sd ${fragment_length_sd}" : ''
      } else {
        fragment_length_mean_val = ''
        fragment_length_sd_val = ''
    }
    def st_single_cell_prior = ''
    if (single_cell_prior) {
      st_single_cell_prior = "--single-cell-prior"
    } 
    def rnastrandness = ''
    if (forward_stranded && !unstranded){
        rnastrandness = '--strandedness forward'
    } else if (reverse_stranded && !unstranded){
        rnastrandness = '--strandedness reverse'
    }
    def pairedend = ''
    if ( !params.singleEnd ) {
        pairedend = '--paired-end' 
    }
    """
    rsem-calculate-expression \
      --alignments \
      ${pairedend} \
      ${rnastrandness} \
      -p ${task.cpus} \
      --append-names \
      --seed 12345 \
      ${st_single_cell_prior} \
      --calc-ci \
      --sort-bam-by-coordinate \
      --estimate-rspd \
      --ci-memory 30000 \
      ${fragment_length_mean_val} ${fragment_length_sd_val} \
      ${ranscriptome_bam} \
      ${rsemidx} ${sample_id};
    rsem-plot-model ${sample_id} ${sample_id}.pdf    
    """
}

fastqc_ch_test = fastqc_reports
     .flatten()
     .filter{  it.toString().endsWith("fastqc.zip") || it.toString().endsWith("report.txt")}       

/*
 * STEP 8 - MultiQC analysis
 */
process multiqc {
  tag "Running MultiQC for sample $sample_id"
  publishDir "${params.publishdir}/MultiQC", mode: 'copy'
  input:
    file ('*') from fastqc_ch_test.toList()
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
  output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
  script:
  """
  multiqc . -m rseqc  -m star -m cutadapt -m fastqc
  """
}

/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
    def msg = """\
    Pipeline execution summary
    ---------------------------
    Status:      : ${ workflow.success ? 'COMPLETED' : 'FAILED' }

    Started at   : ${workflow.start}
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}

    Work Dir     : ${workflow.workDir}
    Launch Dir   : ${workflow.launchDir}
    Project Dir  : ${workflow.projectDir}
    """.stripIndent()

    if (! workflow.success) {    
       def errmsg = """\

       Report for task that caused the workflow execution to fail:
       Exit status   : ${workflow.exitStatus}
       Error message : ${workflow.errorMessage}
       Error report  : ${ workflow.errorReport ? workflow.errorReport : '-' }
       """.stripIndent()

       msg = msg + errmsg
    }
    
    status = workflow.success ? 'completed' : 'failed'
    sendMail(from: 'rpd@gis.a-star.edu.sg', to: "${params.mail_to}", 
             subject: "Nextflow execution ${status}: ${workflow_name}", body: msg)
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
