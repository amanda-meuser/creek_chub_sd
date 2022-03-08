# Limited evidence of a genetic basis for sex determination in the common creek chub, *Semotilus atromaculatus*

R scripts posted here were used for creating plots in the aforementioned journal article. These data and scripts which support the findings of this study are also publicly available on DataDryad (doi:10.5061/dryad.pnvx0k6nf) and raw fastq files can be found in the NCBI Sequence Read Archive (BioProject ID: PRJNA804573). 

Manhattan_function.R contains the function that was adapted from the Abecasis lab to create the Manhattan plots, while the manhattanplot_...R files contain the actual script for generating the plot itself. The 1132.R file requires the use of a .csv file with the chromosome lengths from the assembly. Suggested code for creating this .csv file from the assembly.fasta file can be found at the top of the 1132.R file. 

venndiagram.R requires FST_DAPC_outliers.txt, which is a tab-delimited file that contains columns of significant markers. Columns must be separated by data set and analysis, ie. fst_16020, dapc_dapc, etc. 

The marker_depth.R file uses .idepth files, created by using the --site-mean-depth command in VCFtools, for each VCF file.

dapc.R takes VCF files to calculate Discriminant analysis of principal components (DAPC).

The FST_GEMMA_plots.R and FST_plots_outliers.R files take the .assoc.txt output files from GEMMA and the .weir.fst files from FST calculation as input. 
