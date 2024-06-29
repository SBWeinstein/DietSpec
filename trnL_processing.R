##Sequence processing template
#repeat for each illumina run

# R version (4.2.2)
#17 Jul 2023
#code based on:
#Dada2: https://benjjneb.github.io/dada2/tutorial.html
#For multiple run modifications: https://benjjneb.github.io/dada2/bigdata.html

###Load packages
.cran_packages <- c("ggplot2", "gridExtra", "knitr")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn", "BiocStyle")
# Load packages into session
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

##################################################
###############################################

#set seed for replicatibility between analyses
set.seed(100)

#set working directory
setwd("[WORKING DIRECTORY HERE]") 

###############from FastQ files to (merged) OTU table#####################
###use Cutadapt to first remove primers and filter sequences
#sequence folders
path <- "FILE PATH HERE"  ##  directory containing the fastq files.
head(list.files(path))

#generate matched lists of the forward and reverse read files, and parse out the sample names
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

FWD <- "GGGCAATCCTGAGCCAA"  ##  forward primer sequence, trnl-g
##REV <- "CCWTTGAGTCTCTGCACCTWTC"  ## Reverse primer sequence, trnl-h53  new modified primer
REV <- "CCATTGAGTCTCTGCACCTATC"  ## Reverse primer sequence, trnl-h old trnl primer

###  verify the presence and orientation of these primers in the data
#the function "complement()" is not returning the same object type as others.  removed for now
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult.
#“pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

###count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. 
#Identifying and counting the primers on one set of paired end FASTQ files is sufficient, 
#assuming all the files were created using the same library preparation
library(ShortRead)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))

####tell R where to find cutadapt


cutadapt <-"/LOCATION ON COMPUTER"
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                              "--discard-untrimmed", "--minimum-length", 8,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check,  count the presence of primers in a cutadapt-ed sample (#2)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[3]]))

#primer-free sequence files are now ready to be analyzed through the DADA2 pipeline
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#inspect read quality
plotQualityProfile(cutFs[1:2])

plotQualityProfile(cutRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#filter and trim sequences
#set standard filtering parameters, the most important being the enforcement
# of a maximum of 2 expected errors per-read (Edgar and Flyvbjerg 2015). 
#Trimming and filtering is performed on paired reads jointly, i.e. 
#both reads must pass the filter for the pair to pass.

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  
head(out)

##estimate error rates to parameterize model of substitution errors to distinguish
#sequencing error from biological variation
#set multithread FALSE, avoid crashing current dada2 version
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE) #TRUE crashed R

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#infer sequence variants using DADA2 to infer amplicon sequence variants
#Dereplication: combining identical sequence reads in "unique sequences" with 
#corresponding abundance. 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#dereplicate data
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
#head(track)
track

#save sequence table from run
saveRDS(seqtab, "seqtab.rds") 
saveRDS(track, "read_trackertrnL.rds") #save read tracker
#####################################################