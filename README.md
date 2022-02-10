# Limited evidence of a genetic basis for sex determination in the common creek chub, *Semotilus atromaculatus*

R scripts posted here were used for creating plots in the aforementioned journal article

Manhattan_function.R contains the function that was adapted from the Abecasis lab to create the Manhattan plots, while the manhattanplot_...R files contain the actual script for generating the plot itself. The 1132.R file requires the use of a .csv file with the chromosome lengths from the assembly. Suggested code for creating this .csv file from the assembly.fasta file can be found at the top of the 1132.R file. 

venndiagram.R requires a text file that contains columns of significant markers. Columns must be separated by data set and analysis, ie. 16020_fst, 16020_dapc, 3052_fst, 3052_dapc, 1132_fst, 1132_dapc. 

The marker_depth.R file uses .idepth files, created by using the --site-mean-depth command in VCFtools, for each VCF file.
