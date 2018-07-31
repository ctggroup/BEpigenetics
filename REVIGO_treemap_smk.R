## Taken from Revigo, with slight modifications by Daniel Trejo Banos
## + changed the title
## + changed the output file

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000122","negative regulation of transcription from RNA polymerase II promoter",0.199,10.4131,0.442,0.000,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0050794","regulation of cellular process",18.840,2.2059,0.497,0.698,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0090304","nucleic acid metabolic process",21.449,3.1968,0.644,0.562,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:1901576","organic substance biosynthetic process",30.365,2.9655,0.706,0.658,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0048519","negative regulation of biological process",1.984,3.0199,0.644,0.261,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0051482","positive regulation of cytosolic calcium ion concentration involved in phospholipase C-activating G-protein coupled signaling pathway",0.004,2.0425,0.562,0.659,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0010467","gene expression",19.671,3.0125,0.806,0.143,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",26.547,2.9646,0.678,0.597,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0072507","divalent inorganic cation homeostasis",0.111,2.9799,0.556,0.194,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0019438","aromatic compound biosynthetic process",16.954,3.6514,0.671,0.477,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0018130","heterocycle biosynthetic process",17.388,3.6593,0.669,0.215,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",17.871,3.5665,0.694,0.454,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:1901360","organic cyclic compound metabolic process",30.324,2.8425,0.877,0.159,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0006366","transcription from RNA polymerase II promoter",1.430,5.1227,0.660,0.382,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0016070","RNA metabolic process",15.951,3.3454,0.639,0.577,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0044260","cellular macromolecule metabolic process",34.276,2.1746,0.725,0.417,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0009059","macromolecule biosynthetic process",19.548,3.3174,0.696,0.495,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0044249","cellular biosynthetic process",30.048,3.0131,0.678,0.585,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0043170","macromolecule metabolic process",39.491,2.0871,0.870,0.211,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0006355","regulation of transcription, DNA-templated",9.917,4.0301,0.320,0.695,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",22.502,3.4458,0.651,0.519,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0034654","nucleobase-containing compound biosynthetic process",14.533,3.6768,0.604,0.477,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0034645","cellular macromolecule biosynthetic process",19.291,3.3697,0.604,0.536,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0019222","regulation of metabolic process",11.942,2.4003,0.542,0.438,"negative regulation of transcription from RNA polymerase II promoter"),
                     c("GO:0008150","biological_process",100.000,3.7271,1.000,0.000,"biological_process"),
                     c("GO:0008152","metabolic process",75.387,3.0591,0.995,0.000,"metabolism"),
                     c("GO:0009410","response to xenobiotic stimulus",0.059,60.9739,0.914,0.000,"response to xenobiotic stimulus"),
                     c("GO:0070887","cellular response to chemical stimulus",1.007,5.3943,0.854,0.562,"response to xenobiotic stimulus"),
                     c("GO:0051716","cellular response to stimulus",9.561,2.7996,0.842,0.633,"response to xenobiotic stimulus"),
                     c("GO:0006805","xenobiotic metabolic process",0.051,30.8376,0.846,0.450,"response to xenobiotic stimulus"),
                     c("GO:0007200","phospholipase C-activating G-protein coupled receptor signaling pathway",0.019,3.8621,0.703,0.231,"response to xenobiotic stimulus"),
                     c("GO:0035025","positive regulation of Rho protein signal transduction",0.007,2.0489,0.715,0.279,"response to xenobiotic stimulus"),
                     c("GO:0042221","response to chemical",3.071,7.0381,0.909,0.338,"response to xenobiotic stimulus"),
                     c("GO:0007186","G-protein coupled receptor signaling pathway",0.882,2.2446,0.603,0.449,"response to xenobiotic stimulus"),
                     c("GO:0009987","cellular process",63.780,2.8498,0.992,0.000,"cellular process"),
                     c("GO:0050896","response to stimulus",12.210,3.9816,0.981,0.000,"response to stimulus"),
                     c("GO:0065007","biological regulation",20.498,2.1179,0.983,0.000,"biological regulation"),
                     c("GO:0033363","secretory granule organization",0.006,2.0412,0.937,0.020,"secretory granule organization"),
                     c("GO:0060155","platelet dense granule organization",0.001,2.0634,0.941,0.611,"secretory granule organization"),
                     c("GO:0044237","cellular metabolic process",53.061,3.2214,0.863,0.048,"cellular metabolism"),
                     c("GO:0046483","heterocycle metabolic process",29.664,2.9309,0.828,0.176,"cellular metabolism"),
                     c("GO:0034641","cellular nitrogen compound metabolic process",34.137,2.7768,0.754,0.260,"cellular metabolism"),
                     c("GO:0006725","cellular aromatic compound metabolic process",29.628,2.9198,0.828,0.245,"cellular metabolism"),
                     c("GO:0009058","biosynthetic process",31.611,2.9428,0.937,0.078,"biosynthesis"),
                     c("GO:0006807","nitrogen compound metabolic process",38.744,2.0200,0.935,0.088,"nitrogen compound metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_smk.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "value",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap for smoking",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
