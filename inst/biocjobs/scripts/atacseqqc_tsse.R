## BiocJobs job script: atacseqqc_tsse
## The declared interface lives in ../atacseqqc_tsse.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/atacseqqc_tsse.yaml"), params = list(inputbam = "test-data/test_chr1_chr2.bam", gtf = "test-data/gencode.vM39.annotation_short.gtf.gz", region_cond = "first", output_tsse_values = TRUE, output_tsse_score = TRUE))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . atacseqqc_tsse --out wrappers/atacseqqc_tsse.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . atacseqqc_tsse --out wrappers/atacseqqc_tsse.nf


params <- BiocJobs::jobParams("ATACseqQC", "atacseqqc_tsse")

## -- your analysis ---------------------------------------------------------
## Read inputs from params$<input name>, honour params$<option name>, and
## write every declared output to params$<output name>. Fail loudly with
## stop() on bad input; never prompt or open devices interactively.

suppressPackageStartupMessages({
  library(ATACseqQC)
  library(rtracklayer)
  library(Rsamtools)
})
bamfile <- params$inputbam
gtf <- rtracklayer::import.gff(params$gtf)
seqlengths_from_bam <- scanBamHeader(bamfile)[[1]]$targets

if (params$region_cond == "first") {
  which <- GRanges(
    seqnames= names(seqlengths_from_bam)[1],
    ranges=IRanges(start = 1, end = seqlengths_from_bam[1]))
} else if (params$region_cond == "GW") {
  which <- GRanges(
    seqnames= names(seqlengths_from_bam),
    ranges=IRanges(start = 1, end = seqlengths_from_bam))
} else {
  stopifnot("choosen chromosome is not part of the BAM"=region_cond["chromosome"] %in% names(seqlengths_from_bam))
  which <- GRanges(
    seqnames = region_cond["chromosome"],
    ranges=IRanges(start = 1, end = seqlengths_from_bam[region_cond["chromosome"]]))
}

# Here I consider the BAM file as single-read as the paired info is not used.
# As the bin sizes for the TSSE are relatively large the shift of 5bp is unnecessary.
gal <- readBamFile(bamfile, which=which, bigFile=TRUE)

# Subset the gtf
gtf_sub <- subsetByOverlaps(gtf, which)

if ('transcript' %in% gtf_sub$type) {
    txs <- subset(gtf_sub, type == "transcript")
} else {
    # Aggregate exons to transcripts
    if (! 'transcript_id' %in% colnames(mcols(gtf_sub))) {
        stop("The gtf contains no transcript and the transcript_id is not part of attributes. This gtf cannot be used for this tool.")
    }
    gtf_sub_exons <- subset(gtf_sub, type == "exon")
    txs <- unlist(range(GenomicRanges::split(gtf_sub_exons, gtf_sub_exons$transcript_id)))
    names(txs) <- NULL
}

# Will be exposed one day maybe:
updown <- 1000
endSize <- 100
width <- 100
step <- 100

tsse <- TSSEscore(gal, txs,
                  upstream = updown, downstream = updown,
                  endSize = endSize, width = width, step = step)
x_values <- seq(- updown + width / 2, updown - width / 2, step)

png(params$tsse_plot)
plot(x_values, tsse$values, type = "b",
     xlab = "distance to TSS",
     ylab = "aggregate TSS score",
     main = paste("TSSEscore:", tsse$TSSEscore))
dev.off()

if (params$output_tsse_values) {
  write.table(
    data.frame(x = x_values, v = tsse$values), "tsse_values.txt",
    sep = "\t", row.names = FALSE, col.names = FALSE)
}

if (params$output_tsse_score) {
  cat(tsse$TSSEscore, file = "tsse_score.txt")
}
