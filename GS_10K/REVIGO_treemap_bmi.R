

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
revigo.data <- rbind(c("GO:0007176","regulation of epidermal growth factor-activated receptor activity",0.005,13.0463,0.610,0.000,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:1902653","secondary alcohol biosynthetic process",0.064,10.0090,0.814,0.223,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0055088","lipid homeostasis",0.041,3.2875,0.841,0.668,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0032787","monocarboxylic acid metabolic process",2.485,2.0467,0.810,0.290,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0055081","anion homeostasis",0.045,4.0594,0.840,0.143,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0048522","positive regulation of cellular process",1.585,1.8265,0.745,0.549,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0070741","response to interleukin-6",0.005,3.9467,0.832,0.644,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0048519","negative regulation of biological process",1.984,1.9973,0.784,0.229,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0048518","positive regulation of biological process",1.744,1.9255,0.787,0.332,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0023051","regulation of signaling",0.934,1.8815,0.758,0.536,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0051716","cellular response to stimulus",9.561,2.1381,0.811,0.633,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006366","transcription from RNA polymerase II promoter",1.430,1.8409,0.837,0.194,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0019216","regulation of lipid metabolic process",0.095,3.2994,0.702,0.482,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0048583","regulation of response to stimulus",1.120,1.8507,0.734,0.463,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0016070","RNA metabolic process",15.951,1.9758,0.819,0.577,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0010646","regulation of cell communication",0.929,1.8448,0.744,0.521,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0044260","cellular macromolecule metabolic process",34.276,2.0086,0.854,0.211,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0009059","macromolecule biosynthetic process",19.548,1.8785,0.855,0.558,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0044249","cellular biosynthetic process",30.048,1.8958,0.844,0.658,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0043170","macromolecule metabolic process",39.491,2.0831,0.928,0.224,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006355","regulation of transcription, DNA-templated",9.917,1.8566,0.593,0.695,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:1902275","regulation of chromatin organization",0.069,2.8366,0.789,0.585,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",22.502,1.9075,0.835,0.587,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0044267","cellular protein metabolic process",14.293,1.7359,0.851,0.373,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0019222","regulation of metabolic process",11.942,2.0752,0.727,0.698,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0050794","regulation of cellular process",18.840,2.3088,0.700,0.474,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0050796","regulation of insulin secretion",0.025,2.4146,0.679,0.510,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:1901576","organic substance biosynthetic process",30.365,1.9419,0.867,0.290,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0046483","heterocycle metabolic process",29.664,1.9811,0.905,0.260,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:1901564","organonitrogen compound metabolic process",17.886,1.8045,0.912,0.415,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0010565","regulation of cellular ketone metabolic process",0.049,2.6151,0.729,0.217,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006139","nucleobase-containing compound metabolic process",26.547,1.9747,0.837,0.597,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0030258","lipid modification",0.377,2.4016,0.789,0.681,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0097006","regulation of plasma lipoprotein particle levels",0.011,4.0583,0.853,0.150,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0019438","aromatic compound biosynthetic process",16.954,1.9097,0.847,0.531,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0070848","response to growth factor",0.150,2.6321,0.811,0.656,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0071398","cellular response to fatty acid",0.008,5.8383,0.815,0.499,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0071402","cellular response to lipoprotein particle stimulus",0.003,5.0248,0.819,0.574,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0032924","activin receptor signaling pathway",0.009,4.6869,0.701,0.676,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0032926","negative regulation of activin receptor signaling pathway",0.002,9.1902,0.659,0.572,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006725","cellular aromatic compound metabolic process",29.628,1.9642,0.906,0.260,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0071772","response to BMP",0.044,5.2432,0.800,0.543,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",0.279,2.9343,0.656,0.590,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0007166","cell surface receptor signaling pathway",0.920,2.0940,0.669,0.372,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0042221","response to chemical",3.071,1.9320,0.857,0.452,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0010033","response to organic substance",0.900,2.0308,0.822,0.543,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0090304","nucleic acid metabolic process",21.449,1.9907,0.820,0.444,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0019538","protein metabolic process",18.489,1.7413,0.899,0.408,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0010467","gene expression",19.671,1.9344,0.903,0.417,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0042987","amyloid precursor protein catabolic process",0.004,5.8125,0.945,0.138,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0042982","amyloid precursor protein metabolic process",0.006,4.3141,0.944,0.141,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0018130","heterocycle biosynthetic process",17.388,1.9414,0.846,0.536,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0031063","regulation of histone deacetylation",0.006,7.4657,0.742,0.525,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0009719","response to endogenous stimulus",0.526,2.1521,0.875,0.244,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0008202","steroid metabolic process",0.161,4.0065,0.812,0.451,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006638","neutral lipid metabolic process",0.042,3.7331,0.819,0.454,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006641","triglyceride metabolic process",0.038,4.4615,0.816,0.354,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006631","fatty acid metabolic process",0.878,2.6187,0.757,0.607,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006635","fatty acid beta-oxidation",0.080,3.7117,0.776,0.476,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0043388","positive regulation of DNA binding",0.013,3.5615,0.843,0.424,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006082","organic acid metabolic process",9.086,2.0342,0.792,0.608,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0010893","positive regulation of steroid biosynthetic process",0.003,10.8955,0.682,0.188,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0098732","macromolecule deacylation",0.085,3.6827,0.930,0.248,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0060334","regulation of interferon-gamma-mediated signaling pathway",0.003,9.6610,0.672,0.462,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0044283","small molecule biosynthetic process",5.677,2.5927,0.798,0.488,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:1901362","organic cyclic compound biosynthetic process",17.871,2.0024,0.862,0.541,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0006968","cellular defense response",0.002,4.5168,0.904,0.479,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0034641","cellular nitrogen compound metabolic process",34.137,1.9730,0.878,0.277,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0001676","long-chain fatty acid metabolic process",0.056,3.1923,0.795,0.642,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0034654","nucleobase-containing compound biosynthetic process",14.533,1.9183,0.809,0.562,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0070102","interleukin-6-mediated signaling pathway",0.003,5.5446,0.689,0.659,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0034645","cellular macromolecule biosynthetic process",19.291,1.8852,0.807,0.555,"regulation of epidermal growth factor-activated receptor activity"),
c("GO:0008150","biological_process",100.000,2.9710,1.000,0.000,"biological_process"),
c("GO:0008152","metabolic process",75.387,2.2920,0.998,0.000,"metabolism"),
c("GO:0009987","cellular process",63.780,2.6680,0.997,0.000,"cellular process"),
c("GO:0023052","signaling",6.765,2.0787,0.993,0.000,"signaling"),
c("GO:0042755","eating behavior",0.006,5.9516,0.988,0.000,"eating behavior"),
c("GO:0050896","response to stimulus",12.210,2.2670,0.994,0.000,"response to stimulus"),
c("GO:0051179","localization",18.495,1.7829,0.994,0.000,"localization"),
c("GO:0065007","biological regulation",20.498,2.4701,0.994,0.000,"biological regulation"),
c("GO:0071704","organic substance metabolic process",58.357,2.2148,0.966,0.013,"organic substance metabolism"),
c("GO:0044237","cellular metabolic process",53.061,2.2115,0.922,0.119,"organic substance metabolism"),
c("GO:0044238","primary metabolic process",53.743,2.2083,0.966,0.120,"organic substance metabolism"),
c("GO:0061077","chaperone-mediated protein folding",0.043,7.2502,0.968,0.018,"chaperone-mediated protein folding"),
c("GO:0007032","endosome organization",0.020,3.9991,0.949,0.020,"endosome organization"),
c("GO:0006457","protein folding",0.903,3.7299,0.961,0.026,"protein folding"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,2.8269,0.958,0.031,"organic hydroxy compound metabolism"),
c("GO:0043585","nose morphogenesis",0.000,10.8345,0.857,0.036,"nose morphogenesis"),
c("GO:0010171","body morphogenesis",0.015,3.7124,0.851,0.377,"nose morphogenesis"),
c("GO:0071827","plasma lipoprotein particle organization",0.007,8.3308,0.821,0.368,"nose morphogenesis"),
c("GO:0071825","protein-lipid complex subunit organization",0.010,6.3322,0.941,0.388,"nose morphogenesis"),
c("GO:0090077","foam cell differentiation",0.005,5.3136,0.833,0.444,"nose morphogenesis"),
c("GO:0070207","protein homotrimerization",0.005,6.7357,0.939,0.374,"nose morphogenesis"),
c("GO:0070206","protein trimerization",0.015,4.4800,0.936,0.667,"nose morphogenesis"),
c("GO:0035019","somatic stem cell population maintenance",0.010,2.5524,0.846,0.440,"nose morphogenesis"),
c("GO:0043584","nose development",0.003,7.6386,0.843,0.641,"nose morphogenesis"),
c("GO:0001843","neural tube closure",0.019,2.7776,0.816,0.597,"nose morphogenesis"),
c("GO:0014010","Schwann cell proliferation",0.001,10.2289,0.797,0.349,"nose morphogenesis"),
c("GO:0045668","negative regulation of osteoblast differentiation",0.008,3.5409,0.664,0.457,"nose morphogenesis"),
c("GO:0060325","face morphogenesis",0.006,4.6400,0.849,0.673,"nose morphogenesis"),
c("GO:0060324","face development",0.010,3.7617,0.849,0.692,"nose morphogenesis"),
c("GO:0002089","lens morphogenesis in camera-type eye",0.005,6.6427,0.838,0.620,"nose morphogenesis"),
c("GO:0048144","fibroblast proliferation",0.014,2.9752,0.911,0.587,"nose morphogenesis"),
c("GO:0048147","negative regulation of fibroblast proliferation",0.005,5.9908,0.748,0.509,"nose morphogenesis"),
c("GO:0060349","bone morphogenesis",0.018,2.7297,0.830,0.637,"nose morphogenesis"),
c("GO:0009948","anterior/posterior axis specification",0.018,4.6151,0.836,0.431,"nose morphogenesis"),
c("GO:0007422","peripheral nervous system development",0.019,2.8304,0.826,0.563,"nose morphogenesis"),
c("GO:0021772","olfactory bulb development",0.007,5.1530,0.824,0.479,"nose morphogenesis"),
c("GO:0010742","macrophage derived foam cell differentiation",0.005,5.3136,0.833,0.408,"nose morphogenesis"),
c("GO:0048741","skeletal muscle fiber development",0.010,4.6388,0.797,0.513,"nose morphogenesis"),
c("GO:0030301","cholesterol transport",0.020,7.7025,0.860,0.044,"cholesterol transport"),
c("GO:0006810","transport",17.616,1.8041,0.942,0.320,"cholesterol transport"),
c("GO:1901264","carbohydrate derivative transport",0.208,3.4724,0.884,0.429,"cholesterol transport"),
c("GO:1905039","carboxylic acid transmembrane transport",0.548,2.6351,0.866,0.636,"cholesterol transport"),
c("GO:0032370","positive regulation of lipid transport",0.011,3.2428,0.740,0.667,"cholesterol transport"),
c("GO:0032365","intracellular lipid transport",0.011,7.0660,0.872,0.670,"cholesterol transport"),
c("GO:1901998","toxin transport",0.008,4.6160,0.965,0.160,"cholesterol transport"),
c("GO:1902001","fatty acid transmembrane transport",0.002,6.9977,0.877,0.602,"cholesterol transport"),
c("GO:0010883","regulation of lipid storage",0.010,3.9427,0.787,0.655,"cholesterol transport"),
c("GO:0015850","organic hydroxy compound transport",0.082,3.9276,0.891,0.365,"cholesterol transport"),
c("GO:0010876","lipid localization",0.296,3.6211,0.931,0.473,"cholesterol transport"),
c("GO:1990542","mitochondrial transmembrane transport",0.079,4.0016,0.958,0.183,"cholesterol transport"),
c("GO:0007154","cell communication",7.219,2.0627,0.952,0.044,"cell communication"),
c("GO:0030212","hyaluronan metabolic process",0.008,4.0520,0.951,0.058,"hyaluronan metabolism"),
c("GO:0006629","lipid metabolic process",3.522,2.4020,0.866,0.069,"lipid metabolism"),
c("GO:0044281","small molecule metabolic process",15.138,2.1477,0.866,0.300,"lipid metabolism"),
c("GO:0009058","biosynthetic process",31.611,1.9349,0.969,0.082,"biosynthesis"),
c("GO:0006807","nitrogen compound metabolic process",38.744,2.1118,0.968,0.094,"nitrogen compound metabolism"),
c("GO:1901360","organic cyclic compound metabolic process",30.324,2.0497,0.932,0.099,"organic cyclic compound metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_bmi.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "value",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap BMI",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
