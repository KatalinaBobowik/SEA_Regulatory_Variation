

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );
library(ggpubr)

outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/GeneSetAnalysis/"

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0007155","cell adhesion", 0.544,-2.123, 5.184, 4.844,-5.0526,0.581,0.000),
c("GO:0022610","biological adhesion", 0.550, 3.837, 3.750, 4.849,-4.9355,0.886,0.000),
c("GO:0050866","negative regulation of cell activation", 0.030,-0.266,-7.443, 3.586,-5.6478,0.393,0.000),
c("GO:0051607","defense response to virus", 0.098, 5.690,-2.958, 4.098,-7.1129,0.241,0.000),
c("GO:0043269","regulation of ion transport", 0.244,-2.549,-4.322, 4.496,-4.3429,0.646,0.187),
c("GO:0015837","amine transport", 0.018,-6.123,-1.025, 3.375,-4.2248,0.824,0.196),
c("GO:0007166","cell surface receptor signaling pathway", 0.920, 2.858,-4.821, 5.072,-4.1681,0.576,0.318),
c("GO:0006613","protein targeting to membrane", 0.139,-5.743,-2.796, 4.253,-4.1694,0.735,0.376));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.4, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 5 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank(), legend.text=element_text(size=10))  +  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title = element_text(size = 14)) + ggtitle("SMB vs KOR");
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(-15,15);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# rename object
SMBvMPI=p1;

# --------------------------------------------------------------------------

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0002376","immune system process", 0.600, 2.157,-3.960, 4.886,-8.7447,0.988,0.000),
c("GO:0007155","cell adhesion", 0.544, 2.251,-6.248, 4.844,-9.0921,0.868,0.000),
c("GO:0007166","cell surface receptor signaling pathway", 0.920,-4.873, 2.534, 5.072,-10.4067,0.630,0.000),
c("GO:0022610","biological adhesion", 0.550, 3.963,-4.044, 4.849,-8.9747,0.988,0.000),
c("GO:0023052","signaling", 6.765,-1.229,-6.090, 5.939,-5.6308,0.989,0.000),
c("GO:0032501","multicellular organismal process", 2.373,-1.648,-6.927, 5.483,-7.5850,0.988,0.000),
c("GO:0032502","developmental process", 2.812,-3.334,-6.332, 5.557,-4.4935,0.988,0.000),
c("GO:0040011","locomotion", 0.997, 0.018,-5.137, 5.107,-4.4647,0.988,0.000),
c("GO:0048856","anatomical structure development", 2.540, 5.866, 2.568, 5.513,-6.0511,0.741,0.000),
c("GO:0050896","response to stimulus",12.210, 1.122,-6.044, 6.195,-4.2933,0.989,0.000),
c("GO:0007154","cell communication", 7.219, 4.379,-4.980, 5.967,-6.1057,0.964,0.045),
c("GO:0043062","extracellular structure organization", 0.061, 4.468,-1.823, 3.894,-3.4762,0.893,0.068),
c("GO:0008283","cell proliferation", 0.394, 5.748,-3.248, 4.704,-5.2487,0.920,0.080),
c("GO:0070234","positive regulation of T cell apoptotic process", 0.002,-2.130, 7.498, 2.484,-4.4949,0.716,0.177),
c("GO:0042391","regulation of membrane potential", 0.135,-3.836, 7.008, 4.238,-3.9100,0.805,0.223),
c("GO:0002683","negative regulation of immune system process", 0.080,-4.107, 4.480, 4.012,-6.6990,0.559,0.224),
c("GO:0043269","regulation of ion transport", 0.244,-2.476, 5.728, 4.496,-5.0022,0.701,0.247),
c("GO:0042127","regulation of cell proliferation", 0.313,-4.435, 6.294, 4.603,-5.6861,0.632,0.259),
c("GO:0051239","regulation of multicellular organismal process", 0.628, 0.560, 5.260, 4.906,-7.3391,0.592,0.271),
c("GO:0032879","regulation of localization", 0.726,-3.001, 5.730, 4.969,-3.7457,0.748,0.275),
c("GO:0051674","localization of cell", 0.633, 1.533, 1.748, 4.910,-3.6454,0.931,0.275),
c("GO:0072376","protein activation cascade", 0.018,-5.575,-2.535, 3.367,-3.6103,0.860,0.278),
c("GO:0006954","inflammatory response", 0.110,-6.067,-1.967, 4.151,-8.5452,0.798,0.321),
c("GO:0032101","regulation of response to external stimulus", 0.160,-5.579, 1.706, 4.313,-7.7670,0.628,0.332),
c("GO:0048583","regulation of response to stimulus", 1.120,-5.603, 2.287, 5.158,-6.4498,0.677,0.402),
c("GO:0009605","response to external stimulus", 1.370,-6.632,-1.091, 5.245,-5.3215,0.816,0.420),
c("GO:0070887","cellular response to chemical stimulus", 1.007,-6.775,-0.686, 5.111,-3.3370,0.779,0.457),
c("GO:0042221","response to chemical", 3.071,-6.241,-1.147, 5.595,-3.5033,0.804,0.475),
c("GO:0032689","negative regulation of interferon-gamma production", 0.006, 1.249, 5.585, 2.852,-4.1475,0.620,0.491),
c("GO:0006952","defense response", 0.568,-6.280,-1.597, 4.863,-7.3372,0.809,0.491),
c("GO:0007267","cell-cell signaling", 0.407,-1.620, 1.707, 4.718,-3.5202,0.830,0.493),
c("GO:0050866","negative regulation of cell activation", 0.030,-1.566, 6.323, 3.586,-6.5560,0.656,0.504),
c("GO:0010646","regulation of cell communication", 0.929,-4.226, 5.141, 5.076,-3.6674,0.731,0.521),
c("GO:0007186","G-protein coupled receptor signaling pathway", 0.882,-4.810, 2.316, 5.054,-7.2457,0.631,0.538),
c("GO:0010469","regulation of receptor activity", 0.025,-5.094, 1.846, 3.502,-4.8210,0.635,0.603),
c("GO:0021953","central nervous system neuron differentiation", 0.040, 4.360, 3.988, 3.706,-3.7955,0.668,0.603),
c("GO:0042088","T-helper 1 type immune response", 0.006,-5.849,-0.619, 2.921,-3.6217,0.676,0.606),
c("GO:0070374","positive regulation of ERK1 and ERK2 cascade", 0.034,-4.037, 2.177, 3.636,-6.4413,0.549,0.615),
c("GO:0007517","muscle organ development", 0.080, 4.820, 4.220, 4.009,-3.4520,0.682,0.642),
c("GO:0045165","cell fate commitment", 0.091, 5.444, 2.132, 4.067,-5.6946,0.725,0.650),
c("GO:0001822","kidney development", 0.058, 4.830, 3.887, 3.870,-3.4152,0.663,0.660),
c("GO:0035239","tube morphogenesis", 0.112, 5.022, 3.839, 4.158,-5.7747,0.688,0.663),
c("GO:0051240","positive regulation of multicellular organismal process", 0.303, 0.292, 5.339, 4.589,-5.9136,0.493,0.672),
c("GO:0006959","humoral immune response", 0.035,-5.368,-0.928, 3.650,-5.2314,0.653,0.674),
c("GO:0030855","epithelial cell differentiation", 0.139, 5.358, 2.373, 4.252,-3.4658,0.715,0.677),
c("GO:0007165","signal transduction", 6.621,-4.609, 2.518, 5.929,-4.6778,0.569,0.679),
c("GO:0072358","cardiovascular system development", 0.148, 4.860, 3.530, 4.278,-5.2336,0.660,0.681),
c("GO:0048663","neuron fate commitment", 0.019, 4.438, 3.820, 3.376,-3.3245,0.685,0.692),
c("GO:0035295","tube development", 0.178, 5.002, 3.488, 4.359,-4.7570,0.679,0.694),
c("GO:0010647","positive regulation of cell communication", 0.359,-3.662, 5.143, 4.663,-3.8134,0.636,0.696),
c("GO:0050900","leukocyte migration", 0.060, 0.019,-0.854, 3.890,-3.9161,0.571,0.700));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.5, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4.5 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank(), legend.text=element_text(size=10))  +  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title = element_text(size = 14)) + ggtitle("MTW vs KOR") ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(-20,15);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

# rename object
MTWvMPI=p1;

# now plot both
pdf(paste0(outputdir,"RevigoPlot_SMBvsMPI_MTWvsMPI_LFC05_Pval01.pdf"), height=7, width=14)
ggarrange(SMBvMPI, MTWvMPI)
dev.off()

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood/GeneSetAnalysis/RevigoPlot_SMBvsMPI_LFC05_Pval01.pdf");
