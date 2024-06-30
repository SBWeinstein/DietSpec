## Rules of Herbivory: Large-scale surveys in woodrats (*Neotoma* spp.) reveal constraints on herbivore diet breadth and niche expansion

### R scripts for manuscript are as follows:
__trnl_processing.R:__ Initial processing of plant trnl amplicon sequences. Amplicon sequences are available from the NCBI sequence archive under BioProjects PRJNA722312, PRJNA875083, PRJNA1120941 and Klure and Dearing 2023 (https://doi.org/10.5061/dryad.ghx3ffbss).  T Inputs are fastq.gz files from each sequencing run. Output is a sequence table table.

__merge_filter_seqtables.R:__ Merge sequence tables, remove chimeras, construct a filtered, rarefied phyloseq object with taxonomy and metadata. Taxonomy is assigned using separate custom python script: see https://github.com/robertgreenhalgh/stand. Output is a phyloseq object will all diet samples: ps_rarF_Bigtrnl_29Nov23.rds ( family level ids only) and ps_rar_Bigtrnl_29Nov23.rds (including genera and species level ids).

__all_pop_specialization.R:__ Code to examine diet variation, richness, individual specialization, niche expansion and effects of resource abundance on 57 woodrat populations. Includes contruction of figures 2 and 3, and related supplemental figures. Includes monte carlo resampling null model referenced in figure 3 and maintext. Requires ps_rarF_Bigtrnl_29Nov23.rds, diet_pops_30yrnorm_ppt_1980to2010_prism.csv (for precipitation effects on richness and niche width), and Diet_pops_factors_23mar24.csv (for effects of county-level plant family counts on richness and niche width).

__Range_maps.R:__ Code to make maps included in Figure 4. Requires csv with population coordinates (Diet_pops_20Mar24.csv) and downloaded geojson files of range maps (From IUCN, see main text for citation).

__five_species.R:__ Code to make other components of Figure 4 and run analyses on 5 populations with repeated sampling across their range. Requires ps_rarF_Bigtrnl_29Nov23.rds.
