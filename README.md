## Large-scale surveys of woodrats (Neotoma spp.) reveal constraints on diet breadth in herbivorous mammals

### R scripts for manuscript are as follows:
__trnl_processing.R:__ Initial processing of plant trnl amplicon sequences. Amplicon sequences are available from the NCBI sequence archive under BioProjects PRJNA722312, PRJNA875083, PRJNA1120941 and Klure and Dearing 2023 (https://doi.org/10.5061/dryad.ghx3ffbss). Inputs are fastq.gz files from each sequencing run. Output is a sequence table table.

__merge_filter_seqtables.R:__ Merge sequence tables, remove chimeras, construct a filtered, rarefied phyloseq object with taxonomy and metadata. Taxonomy is assigned using separate custom python script: see https://github.com/robertgreenhalgh/stand using a list of plant genera in the regions (plant_list_21Jul23.txt) generated from the USDA plant database. Output is a phyloseq object will all diet samples: ps_rarF_Bigtrnl_29Nov23.rds ( family level ids only) and ps_rar_Bigtrnl_29Nov23.rds (including genera and species level ids).

__all_pop_specialization.R:__ Code to examine diet variation, richness, individual specialization, niche expansion and effects of resource abundance on 57 woodrat populations. Includes contruction of figures 2 and 3, and related supplemental figures. Includes monte carlo resampling null model referenced in figure 3 and maintext. Requires ps_rarF_Bigtrnl_29Nov23.rds, diet_pops_30yrnorm_ppt_1980to2010_prism.csv (for precipitation effects on richness and niche width), and Diet_pops_factors_23mar24.csv (for effects of county-level plant family counts on richness and niche width).

__range_maps.R:__ Code to make maps included in Figure 4. Requires csv with population coordinates (Diet_pops_20Mar24.csv) and downloaded geojson files of range maps (From IUCN, see main text for citation).

__five_species.R:__ Code to make other components of Figure 4 and run analyses on 5 populations with repeated sampling across their range. Requires ps_rarF_Bigtrnl_29Nov23.rds.

__White_water_NDVI.R:__ Code to create z-score scaled NDVI data for White water, requires files hitewater_ndvi_evi_colnames_19891001_20231218.csv and whitewater_precip_colnames_19801001_202311XX.csv, produces WW_scaled_NDVI.csv and supplemental figure with average monthly NDVI and rainfall patterns.

__White_water_multiyear.R:__ Code to examine diet patterns across 13 surveys conducted at one site (White Water),focusing on seasonal differences, producing Figure 5, Figure 6A, 6D, and related figures in supplemental material. Requires ps_rarF_Bigtrnl_29Nov23.rds and vegetation greenness data from WW_scaled_NDVI.csv. 

__White_water_IS.R:__ Code to test for individual specialization in repeated White water surveys using a Monte Carlo simulation based null model, creates Figure 6C. This script requires multiple components from White_water_multiyear.R, easiest if run immediatly after running previous script.

__niche_density.R:__ Code to make Figure 6B. This script requires multiple components from White_water_multiyear.R, easiest if run immediatly after running previous script.

__mark_recap.R:__ Code to test for within individual diet similarity in mark-recapture data from White water. Requires ps_rarF_Bigtrnl_29Nov23.rds

__leave_out.R:__ Code for "leave one out" analysis to test which plants contribute to within individual diet similarity in mark-recapture data from White water. This script requires multiple components from White_water_multiyear.R, easiest if run immediatly after running previous script.

__diet_subsamp.R:__ Code for Monte Carlo subsampling analysis, using 8 large diet datasets to test how smaller sample sizes in the range used in this study influence estimated diet parameters. This script requires ps_rarF_Bigtrnl_29Nov23.rds and ps_rarF_MatocqBL_26Oct24.rds, which comes from sequences downloaded from SRA accession PRJNA887535 (see citations in supporting information).  

