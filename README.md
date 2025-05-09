PROJECT DESCRIPTION
This repository contains data and analysis code for the paper
"INTROGRESSION, GENETIC DIVERSITY AND STRUCTURE IN QUERCUS BICOLOR:   COMPARING THE RANGE CORE AND SOUTHERN MARGIN"
Jesse B Parker, Sean Hoban, Laura Thompson, Scott E Schlarbaum
DOI: ____________


Data processing and analysis for RADseq data of Quercus bicolor
Sequence reads available at: SRR____________


DESCRIPTION OF FILES IN main directory
File: "hybridization_MS_script.qmd"
Description: Quarto file containing code for all analyses. Code chunks that should be implemented in the command line are preceded with "``` {.bash eval= FALSE}", while those to be run in R are preceded with "``` {r eval= FALSE}". This script is not directly executable, as it requires prior generation of intermediate files across multiple platforms and the download of sequence reads from the SRA. 

File: "hybridization_MS_script.html"
Description: HTML rendering of the .qmd file of the same name


DESCRIPTION OF FILES IN ./data/

File: "C1591A_barcodes.tsv"
Description: Text file containing sequences for the 10bp inline index barcodes for each sample for Library A.

File: "C1591B_barcodes.tsv"
Description: Text file containing sequences for the 10bp inline index barcodes for each sample for Library B.

File: "Dsuite_df.csv"
Description: .csv file containing Results of running the Dsuite Dinvestigate module with the specified trios. 

File: "popmap_dataset1.txt"
Description: Contains the population map for STACKS populations module including reference samples (available at NCBI SRA). All populations are pooled into a single population "1". Some samples of low quality were removed.

File: "popmap_gstacks.txt"
Description: Contains the population map for STACKS gstacks module including reference samples (available at NCBI SRA). All populations are pooled into a single population "1". All 142 samples are included.

File: "pure.txt"
Description: Contains the population map for the "pure" Q. bicolor samples with population information.

File: "Qbicolor_161_pops.csv"
Description: Contains the sample information for dataset1 (all samples and reference samples excluding those quality filtered). Includes the following fields: ID, pop (population name), pop2 (EDGE or CORE designation), DBH (diameter breast height of tree), LON (longitude in decimal degrees. coordinates are reported to one decimal at the request of private landowners), LAT (latitude in decimal degrees. coordinates are reported to one decimal at the request of private landowners), sNMF_V1 through sNMF_V5 (ancestry proportion results of sNMF population structure analysis. this is reproducible if following provided code)  

File: "SETS_2.txt"
Description: Required input for Dsuite analysis. Identifies populations for each target individual, outgroup, and parental populations. 

File: "test_trios_2.txt"
Description: Identifies test crosses for introgression analysis for Dsuite Dinvestigate analyses. 


DESCRIPTION OF FILES IN ./geospatial/

File: "Qbicolor_occurences.csv"
Description: Table of Quercus bicolor occurences and background points. Columns are: Species (all Quercus bicolor), pa (meaning "presence/ absence" - presences are "1" and background points are "0"), x (longitude in decimal degrees), y (latitude in decimal degrees). These points were generated from a combination of SERNEC, iNaturalist, and GBIF records and were rigorsouly filtered (details in supplementary materials).

File: "Qlyrata_occurences.csv"
Description: Table of Quercus bicolor occurences obtained from GBIF post-filtering.