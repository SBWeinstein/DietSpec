## Rules of Herbivory: Large-scale surveys in woodrats (*Neotoma* spp.) reveal constraints on herbivore diet breadth and niche expansion

### R scripts for manuscript are as follows:
trnl_processing.R: Initial processing of plant trnl amplicon sequences. Amplicon sequences are available from the NCBI sequence archive under BioProjects PRJNA722312, PRJNA875083, PRJNA1120941 and Klure and Dearing 2023 (https://doi.org/10.5061/dryad.ghx3ffbss).  T Inputs are fastq.gz files from each sequencing run. Output is a sequence table table.

merge_filter_seqtables.R: Merge sequence tables, remove chimeras, construct a filtered, rarefied phyloseq object with taxonomy and metadata. Taxonomy is assigned using separate custom python script: see https://github.com/robertgreenhalgh/stand. Output is a phyloseq object will all diet samples: ps_rarF_Bigtrnl_29Nov23.rds ( family level ids only) and ps_rar_Bigtrnl_29Nov23.rds (including genera and species level ids).

all_pop_specialization.R: Code to examine diet variation, richness, individual specialization, niche expansion and effects of resource abundance on 57 woodrat populations. Includes contruction of figures 2 and 3, and related supplemental figures. Includes monte carlo resampling null model referenced in figure 3 and maintext.
