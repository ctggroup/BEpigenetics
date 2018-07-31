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
revigo.data <- rbind(c("GO:0007623","circadian rhythm",0.057,2.3961,0.994,0.000,"circadian rhythm"),
                     c("GO:0007631","feeding behavior",0.026,3.2659,0.990,0.000,"feeding behavior"),
                     c("GO:0008150","biological_process",100.000,3.2003,1.000,0.000,"biological_process"),
                     c("GO:0008152","metabolic process",75.387,2.4020,0.999,0.000,"metabolism"),
                     c("GO:0009987","cellular process",63.780,2.7267,0.998,0.000,"cellular process"),
                     c("GO:0023052","signaling",6.765,2.0322,0.994,0.000,"signaling"),
                     c("GO:0030301","cholesterol transport",0.020,18.3023,0.789,0.000,"cholesterol transport"),
                     c("GO:1901998","toxin transport",0.008,5.1802,0.941,0.160,"cholesterol transport"),
                     c("GO:1990542","mitochondrial transmembrane transport",0.079,3.5867,0.928,0.183,"cholesterol transport"),
                     c("GO:0033036","macromolecule localization",3.030,2.7532,0.915,0.259,"cholesterol transport"),
                     c("GO:0006839","mitochondrial transport",0.182,2.4598,0.929,0.281,"cholesterol transport"),
                     c("GO:0034436","glycoprotein transport",0.001,5.1080,0.872,0.328,"cholesterol transport"),
                     c("GO:0015850","organic hydroxy compound transport",0.082,15.2261,0.850,0.365,"cholesterol transport"),
                     c("GO:0006820","anion transport",1.956,2.6295,0.900,0.370,"cholesterol transport"),
                     c("GO:0071702","organic substance transport",4.980,2.6727,0.908,0.423,"cholesterol transport"),
                     c("GO:1901264","carbohydrate derivative transport",0.208,4.6802,0.841,0.429,"cholesterol transport"),
                     c("GO:0015748","organophosphate ester transport",0.144,4.2460,0.845,0.448,"cholesterol transport"),
                     c("GO:0010876","lipid localization",0.296,9.3025,0.881,0.473,"cholesterol transport"),
                     c("GO:0050796","regulation of insulin secretion",0.025,2.5626,0.644,0.480,"cholesterol transport"),
                     c("GO:0015849","organic acid transport",1.024,2.7628,0.821,0.530,"cholesterol transport"),
                     c("GO:0006810","transport",17.616,2.3893,0.897,0.600,"cholesterol transport"),
                     c("GO:1902001","fatty acid transmembrane transport",0.002,3.6273,0.817,0.602,"cholesterol transport"),
                     c("GO:0010878","cholesterol storage",0.003,5.0104,0.761,0.614,"cholesterol transport"),
                     c("GO:0033700","phospholipid efflux",0.003,5.3579,0.817,0.619,"cholesterol transport"),
                     c("GO:0006853","carnitine shuttle",0.000,3.6565,0.834,0.621,"cholesterol transport"),
                     c("GO:0046942","carboxylic acid transport",1.022,2.8885,0.802,0.651,"cholesterol transport"),
                     c("GO:0032370","positive regulation of lipid transport",0.011,4.7234,0.693,0.667,"cholesterol transport"),
                     c("GO:0032365","intracellular lipid transport",0.011,4.2907,0.810,0.670,"cholesterol transport"),
                     c("GO:0034367","macromolecular complex remodeling",0.007,10.4433,0.961,0.000,"macromolecular complex remodeling"),
                     c("GO:0071827","plasma lipoprotein particle organization",0.007,10.2953,0.860,0.382,"macromolecular complex remodeling"),
                     c("GO:0071825","protein-lipid complex subunit organization",0.010,9.9843,0.960,0.390,"macromolecular complex remodeling"),
                     c("GO:0044259","multicellular organismal macromolecule metabolic process",0.012,2.3617,0.868,0.423,"macromolecular complex remodeling"),
                     c("GO:0044236","multicellular organism metabolic process",0.016,2.2816,0.879,0.443,"macromolecular complex remodeling"),
                     c("GO:0061138","morphogenesis of a branching epithelium",0.042,2.0569,0.875,0.476,"macromolecular complex remodeling"),
                     c("GO:0042982","amyloid precursor protein metabolic process",0.006,4.6766,0.954,0.000,"amyloid precursor protein metabolism"),
                     c("GO:0042987","amyloid precursor protein catabolic process",0.004,5.0659,0.955,0.139,"amyloid precursor protein metabolism"),
                     c("GO:0048511","rhythmic process",0.077,2.1207,0.994,0.000,"rhythmic process"),
                     c("GO:0050896","response to stimulus",12.210,2.3237,0.995,0.000,"response to stimulus"),
                     c("GO:0051179","localization",18.495,2.1401,0.995,0.000,"localization"),
                     c("GO:0065007","biological regulation",20.498,2.5258,0.995,0.000,"biological regulation"),
                     c("GO:0097006","regulation of plasma lipoprotein particle levels",0.011,9.5697,0.868,0.000,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0010743","regulation of macrophage derived foam cell differentiation",0.004,4.9917,0.788,0.137,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0007259","JAK-STAT cascade",0.033,8.3600,0.705,0.155,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0055088","lipid homeostasis",0.041,7.6891,0.832,0.159,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0050730","regulation of peptidyl-tyrosine phosphorylation",0.045,3.4531,0.754,0.172,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0009726","detection of endogenous stimulus",0.001,5.3857,0.854,0.185,"regulation of plasma lipoprotein particle levels"),
                     c("GO:1990830","cellular response to leukemia inhibitory factor",0.000,2.7906,0.826,0.185,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0008285","negative regulation of cell proliferation",0.128,2.4147,0.731,0.188,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0006968","cellular defense response",0.002,2.2289,0.879,0.195,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0031331","positive regulation of cellular catabolic process",0.058,2.0436,0.753,0.250,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0051606","detection of stimulus",0.351,3.6408,0.859,0.270,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0016311","dephosphorylation",1.250,2.8592,0.922,0.292,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0001101","response to acid chemical",0.124,2.3789,0.814,0.312,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0006357","regulation of transcription from RNA polymerase II promoter",1.273,2.0334,0.717,0.321,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0009719","response to endogenous stimulus",0.526,3.3406,0.855,0.339,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.115,6.2502,0.683,0.343,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0097696","STAT cascade",0.033,8.3600,0.735,0.386,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0042221","response to chemical",3.071,2.6380,0.833,0.423,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0060707","trophoblast giant cell differentiation",0.002,3.3877,0.843,0.426,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0055099","response to high density lipoprotein particle",0.001,5.3850,0.822,0.432,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0050794","regulation of cellular process",18.840,2.2539,0.748,0.439,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0090077","foam cell differentiation",0.005,4.9645,0.876,0.440,"regulation of plasma lipoprotein particle levels"),
                     c("GO:1990823","response to leukemia inhibitory factor",0.000,2.7906,0.843,0.440,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0007166","cell surface receptor signaling pathway",0.920,2.3510,0.674,0.452,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0051235","maintenance of location",0.129,2.7170,0.788,0.452,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0009720","detection of hormone stimulus",0.000,5.3941,0.822,0.459,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0048583","regulation of response to stimulus",1.120,2.0712,0.741,0.463,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0006366","transcription from RNA polymerase II promoter",1.430,2.1487,0.874,0.465,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0055094","response to lipoprotein particle",0.004,5.1816,0.804,0.473,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0034654","nucleobase-containing compound biosynthetic process",14.533,2.1560,0.851,0.477,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0014902","myotube differentiation",0.028,3.0629,0.860,0.488,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0010646","regulation of cell communication",0.929,2.0319,0.770,0.521,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0060334","regulation of interferon-gamma-mediated signaling pathway",0.003,3.0200,0.638,0.535,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0023051","regulation of signaling",0.934,2.0115,0.782,0.536,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0042063","gliogenesis",0.050,2.4482,0.842,0.561,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0090304","nucleic acid metabolic process",21.449,2.0575,0.862,0.562,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0070848","response to growth factor",0.150,4.2274,0.773,0.563,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0055091","phospholipid homeostasis",0.003,5.3497,0.852,0.566,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0016070","RNA metabolic process",15.951,2.1328,0.861,0.577,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0071398","cellular response to fatty acid",0.008,3.2277,0.781,0.584,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0050869","negative regulation of B cell activation",0.006,2.4812,0.740,0.589,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0018212","peptidyl-tyrosine modification",0.215,2.4581,0.896,0.598,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0048878","chemical homeostasis",0.543,2.6528,0.807,0.619,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0009593","detection of chemical stimulus",0.285,4.5273,0.800,0.631,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0051716","cellular response to stimulus",9.561,2.0587,0.790,0.633,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0010033","response to organic substance",0.900,2.9157,0.787,0.639,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0031347","regulation of defense response",0.132,2.5679,0.702,0.659,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0055081","anion homeostasis",0.045,5.0737,0.831,0.668,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0071495","cellular response to endogenous stimulus",0.402,2.8782,0.794,0.669,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0050853","B cell receptor signaling pathway",0.015,2.2866,0.646,0.670,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0045596","negative regulation of cell differentiation",0.128,2.1900,0.685,0.671,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0001890","placenta development",0.024,2.4631,0.867,0.682,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0080090","regulation of primary metabolic process",11.675,2.1606,0.731,0.693,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0019222","regulation of metabolic process",11.942,2.1345,0.769,0.698,"regulation of plasma lipoprotein particle levels"),
                     c("GO:0044238","primary metabolic process",53.743,2.3766,0.973,0.013,"primary metabolism"),
                     c("GO:0044237","cellular metabolic process",53.061,2.3021,0.934,0.111,"primary metabolism"),
                     c("GO:0071704","organic substance metabolic process",58.357,2.3537,0.973,0.120,"primary metabolism"),
                     c("GO:0006793","phosphorus metabolic process",13.507,2.1773,0.928,0.124,"primary metabolism"),
                     c("GO:0006457","protein folding",0.903,2.5609,0.964,0.023,"protein folding"),
                     c("GO:0061077","chaperone-mediated protein folding",0.043,2.0837,0.971,0.026,"chaperone-mediated protein folding"),
                     c("GO:1901615","organic hydroxy compound metabolic process",0.831,4.1427,0.964,0.031,"organic hydroxy compound metabolism"),
                     c("GO:0006577","amino-acid betaine metabolic process",0.041,3.6826,0.946,0.036,"amino-acid betaine metabolism"),
                     c("GO:0009437","carnitine metabolic process",0.010,3.6962,0.950,0.166,"amino-acid betaine metabolism"),
                     c("GO:0007154","cell communication",7.219,2.0791,0.957,0.044,"cell communication"),
                     c("GO:0016125","sterol metabolic process",0.106,6.8135,0.752,0.053,"sterol metabolism"),
                     c("GO:0006629","lipid metabolic process",3.522,3.9330,0.868,0.154,"sterol metabolism"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",17.871,2.1683,0.894,0.194,"sterol metabolism"),
                     c("GO:0032787","monocarboxylic acid metabolic process",2.485,3.3082,0.804,0.223,"sterol metabolism"),
                     c("GO:0010565","regulation of cellular ketone metabolic process",0.049,2.5887,0.748,0.283,"sterol metabolism"),
                     c("GO:0044281","small molecule metabolic process",15.138,2.8203,0.859,0.300,"sterol metabolism"),
                     c("GO:0042180","cellular ketone metabolic process",0.423,2.4276,0.841,0.350,"sterol metabolism"),
                     c("GO:0006006","glucose metabolic process",0.430,2.8323,0.856,0.350,"sterol metabolism"),
                     c("GO:0006641","triglyceride metabolic process",0.038,3.0302,0.800,0.438,"sterol metabolism"),
                     c("GO:0006638","neutral lipid metabolic process",0.042,2.9430,0.803,0.449,"sterol metabolism"),
                     c("GO:0019438","aromatic compound biosynthetic process",16.954,2.1614,0.880,0.450,"sterol metabolism"),
                     c("GO:0019216","regulation of lipid metabolic process",0.095,4.4453,0.713,0.468,"sterol metabolism"),
                     c("GO:0018130","heterocycle biosynthetic process",17.388,2.1555,0.880,0.477,"sterol metabolism"),
                     c("GO:0008202","steroid metabolic process",0.161,4.7655,0.802,0.486,"sterol metabolism"),
                     c("GO:0044283","small molecule biosynthetic process",5.677,3.5783,0.799,0.488,"sterol metabolism"),
                     c("GO:0044249","cellular biosynthetic process",30.048,2.1657,0.874,0.539,"sterol metabolism"),
                     c("GO:0034645","cellular macromolecule biosynthetic process",19.291,2.0283,0.853,0.553,"sterol metabolism"),
                     c("GO:0009059","macromolecule biosynthetic process",19.548,2.0145,0.892,0.556,"sterol metabolism"),
                     c("GO:0006631","fatty acid metabolic process",0.878,3.8331,0.735,0.578,"sterol metabolism"),
                     c("GO:1902930","regulation of alcohol biosynthetic process",0.008,4.7471,0.748,0.584,"sterol metabolism"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",22.502,2.0869,0.874,0.585,"sterol metabolism"),
                     c("GO:0034440","lipid oxidation",0.087,3.2683,0.781,0.596,"sterol metabolism"),
                     c("GO:0006082","organic acid metabolic process",9.086,2.6068,0.788,0.608,"sterol metabolism"),
                     c("GO:0009062","fatty acid catabolic process",0.111,3.3219,0.748,0.625,"sterol metabolism"),
                     c("GO:0016042","lipid catabolic process",0.401,2.4806,0.781,0.630,"sterol metabolism"),
                     c("GO:0001676","long-chain fatty acid metabolic process",0.056,3.3077,0.775,0.657,"sterol metabolism"),
                     c("GO:1901576","organic substance biosynthetic process",30.365,2.1598,0.898,0.658,"sterol metabolism"),
                     c("GO:0030258","lipid modification",0.377,3.0469,0.772,0.681,"sterol metabolism"),
                     c("GO:1902652","secondary alcohol metabolic process",0.095,6.2168,0.838,0.688,"sterol metabolism"),
                     c("GO:1901360","organic cyclic compound metabolic process",30.324,2.1438,0.945,0.077,"organic cyclic compound metabolism"),
                     c("GO:0044260","cellular macromolecule metabolic process",34.276,2.0420,0.887,0.198,"organic cyclic compound metabolism"),
                     c("GO:0043170","macromolecule metabolic process",39.491,2.1093,0.942,0.224,"organic cyclic compound metabolism"),
                     c("GO:0006725","cellular aromatic compound metabolic process",29.628,2.0026,0.919,0.260,"organic cyclic compound metabolism"),
                     c("GO:0046483","heterocycle metabolic process",29.664,2.0033,0.919,0.260,"organic cyclic compound metabolism"),
                     c("GO:0010467","gene expression",19.671,2.0159,0.928,0.417,"organic cyclic compound metabolism"),
                     c("GO:0009058","biosynthetic process",31.611,2.1444,0.975,0.078,"biosynthesis"),
                     c("GO:0097164","ammonium ion metabolic process",0.178,2.6069,0.975,0.082,"ammonium ion metabolism"),
                     c("GO:0006807","nitrogen compound metabolic process",38.744,2.1845,0.974,0.089,"nitrogen compound metabolism"));

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
	title = "REVIGO Gene Ontology treemap for BMI",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
