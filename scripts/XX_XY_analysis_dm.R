#  If you need help, please contact isabelle.stevant@gmail.com 


source("analysis_functions.R")
require(gridExtra)


###########################################
#                                         #
#          Load and prepare files         #
#                                         #
###########################################

# Load female RPKM data
load(file="../data/female_rpkm.Robj")
female_rpkm <- female_rpkm[,!colnames(female_rpkm) %in% grep("rep",colnames(female_rpkm), value=TRUE)]
prot_coding_genes <- read.csv(file="../data/prot_coding.csv", row.names=1)
females <- female_rpkm[rownames(female_rpkm) %in% as.vector(prot_coding_genes$x),]

# Load female read count data
load(file="../data/female_count.Robj")
female_count <- female_count[,!colnames(female_count) %in% grep("rep",colnames(female_count), value=TRUE)]
female_count <- female_count[rownames(female_count) %in% rownames(females),]

# Load female clustering results
female_clustering <- read.csv(file="../data/female_clustering.csv", row.names=1)$x
female_clustering <- paste("XX", female_clustering, sep="_")
names(female_clustering) <- colnames(females)

# Get cell embryonic stages
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)


# Load male RPKM data
load(file="../data/male_rpkm.Robj")
prot_coding_genes <- read.csv(file="../data/prot_coding.csv", row.names=1)
males <- male_rpkm[rownames(male_rpkm) %in% as.vector(prot_coding_genes$x),]
# this cell passed the QC but is actually an artefact and was removed from the analysis
males <- males[,!colnames(males) %in% "E16.5_XY_20150202_C94_150331_8"]

# Load male read count data
load(file="../data/male_count.Robj")
male_count <- male_count[rownames(male_count) %in% rownames(males),]
colnames(male_count) <- colnames(male_rpkm)
# this cell passed the QC but is actually an artefact and was removed from the analysis
male_count <- male_count[,!colnames(male_count) %in% "E16.5_XY_20150202_C94_150331_8"]

# Load male clustering results
male_clustering <- read.csv(file="../data/male_clustering.csv", row.names=1)$x
male_clustering <- paste("XY", male_clustering, sep="_")
names(male_clustering) <- colnames(males)

# Get cell embryonic stages
male_stages <- sapply(strsplit(colnames(males), "_"), `[`, 1)
names(male_stages) <- colnames(males)

# Merge male and female RPKM and count
all <- cbind(males, females)
all_count <- cbind(male_count, female_count)

# Remove genes never expressed
all <- all[rowSums(all)>0,]
all_count <- all_count[rowSums(all_count)>0,]

# Get cell embryonic stages
all_stages <- sapply(strsplit(colnames(all), "_"), `[`, 1)
names(all_stages) <- colnames(all)

# Get cell sex
all_sex <- sapply(strsplit(colnames(all), "_"), `[`, 2)
names(all_sex) <- colnames(all)

# Get cell types
all_clustering <- c(male_clustering, female_clustering)
names(all_clustering) <- c(names(male_clustering), names(female_clustering))

# Define color palettes
all_sexPalette <- c(
	"#fb7072",
	"#77b0f3"
)

all_stagePalette <- c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21",
	"#43f14b"
)

male_clusterPalette <- c(
	"#457cff", 
	"#aeff01", 
	"#00c5ec", 
	"#009900",  
	"#a38cff",
	"#8dfaff"
)

female_clusterPalette <- c(
	"#560047", 
	"#a53bad", 
	"#eb6bac", 
	"#ffa8a0"
)

all_clusterPalette <- c(
	female_clusterPalette, 
	male_clusterPalette
)


###########################################
#                                         #
#            Var Gene Selection           #
#                                         #
###########################################

# Extract genes with a squared coefficient of variation >1.75 times the fit regression (Brennecke et al 2013 method)
all_data <- getMostVarGenes(all, fitThr=1.75)
all_data <- log(all_data+1)

###########################################
#                                         #
#              RtSNE Analysis             #
#                                         #
###########################################

all_sub_pca <- PCA(
	t(all_data), 
	ncp = ncol(all_data), 
	graph=FALSE
)

# Compute how many PCs to include in the t-SNE using Jackstraw
# the result is 19 significant PCs
significant_pcs <- permutationPA(
	all_sub_pca$ind$coord, 
	B = 100, 
	threshold = 0.05, 
	verbose = TRUE, 
	seed = NULL
)$r

# Compute and plot the t-SNE colored by embryonic stages
all_t_sne <- run_plot_tSNE(
	pca=all_pca,
	pc=significant_pcs,
	iter=5000,
	conditions=all_stages,
	colours=all_stagePalette
)

# Plot the t-SNE colored by sex
all_t_sne_sex <- plot_tSNE(
	tsne=all_t_sne, 
	conditions=all_sex, 
	colours= all_sexPalette
)

# Plot the t-SNE colored by sex-respective cell clusters
all_t_sne_all_clusters <- plot_tSNE(
	tsne=all_t_sne, 
	conditions=all_clustering, 
	colours= all_clusterPalette
)

# Recompute PCA with FactomineR package for the clustering
res.pca <- PCA(
	t(all_data),
	ncp = significant_pcs, 
	graph=FALSE
)

# Clustering cells based on the PCA, with a minimum of 5 expected clusters
res.hcpc <- HCPC(
	res.pca, 
	graph = FALSE,
	min=5
)

# Plot the hierarchical clustering result
plot(res.hcpc, choice ="tree", cex = 0.6)

# Extract clustering results
hc_clustering <- res.hcpc$data.clust$clust
hc_clustering <- paste("C", hc_clustering, sep="")
names(hc_clustering) <- colnames(all_data)
all_cluster <- hc_clustering

# Define new cluster colors
all_newClusteringPalette <- c(
	"#6d8c99",
	"#eec584",
	"#17a398",
	"#313b3c",
	"#c29979"
)

# Plot the t-SNE colored by the new cell clusters
all_t_sne_new_clusters <- plot_tSNE(
	tsne=all_t_sne, 
	conditions=hc_clustering, 
	colours= all_newClusteringPalette
)

###########################################
#                                         #
#          Pseudotime Analysis            #
#                                         #
###########################################

# Foetal Leydig cells were removed (XY_C5) as they are only 7 and fail to create a separate branch in the lineage reconstruction
# We are convinced that if we had more cells from this cell type, the reconstruction would have been easier
FLC <- names(all_clustering[all_clustering=="XY_C5"])
all_data_without_FLC <- all_data[,!colnames(all_data) %in% FLC]
all_stages_without_FLC <- all_stages[!names(all_stages) %in% FLC]
all_clustering_without_FLC <- all_clustering[!names(all_clustering) %in% FLC]
hc_clustering_without_FLC <- hc_clustering[!names(hc_clustering) %in% FLC]
all_sex_without_FLC <- all_sex[!names(all_sex) %in% FLC]
all_without_FLC <- all[,!colnames(all) %in% FLC]
all_count_without_FLC <- all_count[,!colnames(all_count) %in% FLC]

# Compute the Diffusion map
all_dm <- run_diffMap(
	all_data_without_FLC, 
	all_clustering_without_FLC,
	sigma=14
)

# Plot the 3D diffusion map using DCs 2 to 4 colored by the new cell clustering
plot_dm_3D(
	dm=all_dm, 
	dc=c(2:4),
	condition=hc_clustering_without_FLC, 
	colour=all_newClusteringPalette
)


# Plot the Eigen values per diffusion component (similar to screeplots for PCAs)
plot_eigenVal(
	dm=all_dm
)

# Compute the pseudotime using Slingshot based on the informative diffusion map components (DCs)
# We see on the eigenVal plot a first elbow between DC 5 and 6, so we consider the first 6 DCs for further analysis
# DC1 is driven by one outlier cell, so we don't consider it in our analysis
# We added DC 8 to help the lineage reconstruction 
# We also help the pseudotime calculation by giving the starting cliuster and the end clusters
lineage <- get_lineage(
	dm=all_dm, 
	dim=c(2:6,8), 
	condition=hc_clustering_without_FLC,
	start="C2", 
	end=c("C1", "C4", "C5"), 
	shrink.method="tricube"
)

# Plot the 3D diffusion map using DCs 2 to 4 colored by the new cell clustering
# The black line is the cell lineage reconstruction by Slinshot
plot_dm_3D(
	dm=all_dm, 
	dc=c(2:4),
	condition=hc_clustering_without_FLC, 
	colour=all_newClusteringPalette
)
plot3d(lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_all_dm_new_clustering_cols.svg"), fmt = "svg", drawText = TRUE )

# Plot the 3D diffusion map using DCs 2 to 4 colored by embryonic stages
plot_dm_3D(
	dm=all_dm, 
	dc=c(2:4),
	condition=all_stages_without_FLC, 
	colour=all_stagePalette
)
plot3d(lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_all_dm_stages_cols.svg"), fmt = "svg", drawText = TRUE )


# Plot the 3D diffusion map using DCs 2 to 4 colored by sex
plot_dm_3D(
	dm=all_dm, 
	dc=c(2:4),
	condition=all_sex_without_FLC, 
	colour=all_sexPalette
)
plot3d(lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_all_dm_sex_cols.svg"), fmt = "svg", drawText = TRUE )

# Plot the 3D diffusion map using DCs 2 to 4 colored by sex-respective cell clusters
plot_dm_3D(
	dm=all_dm, 
	dc=c(2:4),
	condition=all_clustering_without_FLC, 
	colour=all_clusterPalette[-9]
)
plot3d(lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_all_dm_old_clustering_cols.svg"), fmt = "svg", drawText = TRUE )


# Attribute cells to lineages with a distance of 0.8 from the smoothed line and compute the pseudotime
pseudotime <- get_pseudotime(lineage, wthres=0.8, ranked=TRUE)
rownames(pseudotime) <- colnames(all_data_without_FLC)

## Lineage result: 
# Lineage 1 -> progenitors to Sertoli
# Lineage 2 -> progenitors to Granulosa
# Lineage 3 -> progenitors to progenitors

####################################################
#
# Plot gene expression along pseudotime in cell lineages 1 and 2 (Supporting cells)
#
####################################################

# Twisted palette to plot gene expression by lineages along pseudotime
all_stagePalette2 <- c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21",
	"#3272bd",
	"#ca3032",
	"#43f14b",
	"#ca3032",
	"#3272bd"
)

all_clusterPalette2 <- all_clusterPalette[-9]

# Example of plot: Sox9 in the supporting cell lineage (1 and 2) along the pseidotime
# Lines are smoothed expression by sex
p <- plot_smoothed_gene_per_lineage_sex(
	rpkm_matrix=all_without_FLC, 
	pseudotime=pseudotime, 
	lin=c(1:2),
	gene="Sox9", 
	stages=all_stages_without_FLC, 
	clusters=all_clustering_without_FLC, 
	stage_colors=all_stagePalette2,
	cluster_colors=all_clusterPalette2,
	lineage_colors=all_clusterPalette2
)
print(p)

# List o gene of interest
interestingGenes <- c(
	"Gadd45g",
	"Trpc7",
	"Adam8",
	"Pnlip",
	"Runx1",
	"Ifitm1",
	"Podxl",
	"Msx1",
	"Emx2",
	"Wt1",
	"Nr5a1",
	"Lhx9"
)

# Plot all the genes of interest in one file
p <- list()
for (genes in interestingGenes) {
	p[[genes]] <- plot_smoothed_gene_per_lineage_sex(
	rpkm_matrix=all_without_FLC, 
	pseudotime=pseudotime, 
	lin=c(1:2),
	gene=genes, 
	stages=all_stages_without_FLC, 
	clusters=all_clustering_without_FLC, 
	stage_colors=all_stagePalette2,
	cluster_colors=all_clusterPalette2,
	lineage_colors=all_clusterPalette2
)
}
pdf(file=paste("../graph/", format(Sys.time(), "%Y%m%d"), "_supporting_genes.pdf"), width=10, height=12)
do.call(grid.arrange,c(p, ncol=2))
dev.off()

####################################################
#
# Plot gene expression along pseudotime in cell lineages 3 (Progenitor cells)
#
####################################################

# Twisted palette to plot gene expression by lineages along pseudotime
all_clusterPalette2 <- all_clusterPalette[-9]

all_stagePalette3 <- c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21",
	"#43f14b",
	"#ca3032",
	"#3272bd",
	"#3272bd"
)

# List o gene of interest
interestingGenes <- c(
	"Acta2",
	"Pdgfra",
	"Arx",
	"Ptch1",
	"Inhba",
	"Etl4",
	"Foxl2",
	"Bmp3",
	"Cfh",
	"Znrf3",
	"Enc1",
	"Ablim1"
)

# Example of plot: Cfh in the progenitor cell lineage (3) along the pseidotime
# Lines are smoothed expression by sex
p <- plot_smoothed_gene_per_lineage_sex(
	rpkm_matrix=all_without_FLC, 
	pseudotime=pseudotime, 
	lin=c(3),
	gene="Cfh", 
	stages=all_stages_without_FLC, 
	clusters=all_clustering_without_FLC, 
	stage_colors=all_stagePalette3,
	cluster_colors=all_clusterPalette2,
	lineage_colors=all_clusterPalette2
)
print(p)

# Plot all the genes of interest in one file
p <- list()
for (genes in interestingGenes) {
	p[[genes]] <- plot_smoothed_gene_per_lineage_sex(
	rpkm_matrix=all_without_FLC, 
	pseudotime=pseudotime, 
	lin=c(3),
	gene=genes, 
	stages=all_stages_without_FLC, 
	clusters=all_clustering_without_FLC, 
	stage_colors=all_stagePalette3,
	cluster_colors=all_clusterPalette2,
	lineage_colors=all_clusterPalette2
)
}
pdf(file=paste("../graph/", format(Sys.time(), "%Y%m%d"), "_progenitors_genes.pdf"), width=10, height=12)
do.call(grid.arrange,c(p, ncol=2))
dev.off()

#######################################################
#                                                     #
#            Heatmap Supporting lineage               #
#                                                     #
#######################################################

####################################################
#
# Data preparation
#
####################################################

# Extract female supportig cells
grac_pseudotime <- pseudotime[,2]
names(grac_pseudotime) <- rownames(pseudotime)
grac_pseudotime <- grac_pseudotime[!is.na(grac_pseudotime)]
grac_pseudotime <- grac_pseudotime[grep("XX", names(grac_pseudotime), value=TRUE)]

# Extract male supportig cells
sert_pseudotime <- pseudotime[,1]
names(sert_pseudotime) <- rownames(pseudotime)
sert_pseudotime <- sert_pseudotime[!is.na(sert_pseudotime)]
sert_pseudotime <- sert_pseudotime[grep("XY", names(sert_pseudotime), value=TRUE)]

# Get genes that display dynamic expression through pseudotime in their sex-respective study
de_sert <- read.csv(file="../data/male_sertoli_lineage_pseudotime_DE_genes.csv", row.names=1)
de_grac <- read.csv(file="../data/female_granulosa_lineage_pseudotime_DE_genes.csv", row.names=1)

# Get expression data from supporting cells
supporting_lineage_data <- all[,colnames(all) %in% c(names(grac_pseudotime), names(sert_pseudotime))]

# Get dynamic gene names
supporting_male_sig_gene_pseudoT <- as.vector(de_sert[de_sert$qval<0.05,"genes"])
supporting_female_sig_gene_pseudoT <- as.vector(de_grac[de_grac$qval<0.05,"genes"])

# Merge gene lists
gene_list <- unique(c(supporting_male_sig_gene_pseudoT, supporting_female_sig_gene_pseudoT))

# Get expression matrix of the dynamic genes in the supporting cell lineage for both sexes
de_matrix <- log(supporting_lineage_data[rownames(supporting_lineage_data) %in% gene_list,]+1)


####################################################
#
# Double heatmap lineage
#
####################################################

L1_lineage <- grac_pseudotime
L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), L1_rpkm_exp_lineage))]

L2_lineage <- sert_pseudotime
L2_ordered_lineage <- L2_lineage[order(L2_lineage, decreasing = TRUE)]
L2_rpkm_exp_lineage <- de_matrix[,names(L2_ordered_lineage)]
L2_cells <- L2_rpkm_exp_lineage[,order(match(names(L2_ordered_lineage), L2_rpkm_exp_lineage))]

L1_lineage_cells <- names(L1_ordered_lineage)
L2_lineage_cells <- names(L2_ordered_lineage)

comp_list <- comparelists(L1_lineage_cells, L2_lineage_cells)
common_cells <- comp_list$intersect
L1_spe_cells <- L1_lineage_cells[!L1_lineage_cells %in% comp_list$intersect]
L2_spe_cells <- L2_lineage_cells[!L2_lineage_cells %in% comp_list$intersect]

L1_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L1 cells", L1_spe_cells)
)
names(L1_cellLin) <- c(common_cells, L1_spe_cells)
L1_cellLin <- L1_cellLin[order(match(names(L1_cellLin), colnames(L1_cells)))]

L2_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L2 cells", L2_spe_cells)
)
names(L2_cellLin) <- c(common_cells, L2_spe_cells)
L2_cellLin <- L2_cellLin[order(match(names(L2_cellLin), colnames(L2_cells)))]

cellType_L1 <- all_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)

cellType_L2 <- all_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)

cellLin <- c(
	L2_cellLin,
	L1_cellLin
)

data_heatmap <- data.frame(
	L2_cells,
	L1_cells
)

cellType <- c(
	cellType_L2,
	cellType_L1
)

annotation_col <- data.frame(
	cellType=cellType,
	Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)

names(cellLinCol) <- unique(cellLin)

cellTypeCol <- c(
	XY_C1="#457cff",
	XY_C2="#aeff01",
	XY_C3="#00c5ec",
	XY_C4="#009900",  
	XY_C6="#8dfaff",
	XX_C1="#560047",
	XX_C2="#a53bad",
	XX_C3="#eb6bac", 
	XX_C4="#ffa8a0"
)

annotation_colors <- list(
	cellType=cellTypeCol,
	clustering=gene_cluster_colors,
		Stages=c(
			E10.5="#2754b5",
			E11.5="#8a00b0", 
			E12.5="#d20e0f", 
			E13.5="#f77f05", 
			E16.5="#f9db21",
			P6="#43f14b"
		)
)

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026','#800026'))
mypalette <- c(rev(cold(13)), warm(12))
breaksList = seq(-2.5, 2.5, by = 0.2)

scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
    min = apply(x, 1, min, na.rm = T)
    max = apply(x, 1, max, na.rm = T)
    sdev = apply(x, 1, sd, na.rm = T)
    return((x - m) / sdev)
}

scaled_data_heatmap <- scale_rows(data_heatmap)

pdf(file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_supporting_heatmap_DE_genes_granulosa_sertoli_k_25_p-val_0.05_span_0.4.pdf"),  height = 6, width = 7)

set.seed(123)

gene_clustering <- pheatmap(
	scaled_data_heatmap, 
	kmeans_k=25,
	gaps_col=length(cellType_L2),
	show_colnames=FALSE, 
	cluster_cols=FALSE,
	clustering_method="ward.D",
	annotation_col=annotation_col,
	annotation_colors=annotation_colors,
	annotation_names_row=FALSE,
	color=mypalette,
	breaks=breaksList
)

dev.off()

write.csv(gene_clustering$kmeans$cluster, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_supporting_heatmap_DE_genes_granulosa_sertoli_k-mean_25_p-val_0.05.csv"))


##################################################################################################################################################

#######################################################
#                                                     #
#            Progenitor lineage analysis              #
#                                                     #
#######################################################



prog_pseudotime <- pseudotime[,3]
names(prog_pseudotime) <- rownames(pseudotime)
prog_pseudotime <- prog_pseudotime[!is.na(prog_pseudotime)]
fem_prog_pseudotime <- prog_pseudotime[grep("XX", names(prog_pseudotime), value=TRUE)]
mal_prog_pseudotime <- prog_pseudotime[grep("XY", names(prog_pseudotime), value=TRUE)]



de_mal_prog <- read.csv(file="Supplementary_Data_6_Progenitors_DE_genes_pseudotime.csv", row.names=1)
de_fem_prog <- read.csv(file="../graph/female/20180712_female_lineage2_pseudotime_DE_genes.csv", row.names=1)


###########################################
#										  #
#		Double heatmap lineage			  #
#										  #
###########################################

prog_lineage_data <- all[,colnames(all) %in% c(names(fem_prog_pseudotime), names(mal_prog_pseudotime))]
prog_male_sig_gene_pseudoT <- rownames(de_mal_prog)
prog_female_sig_gene_pseudoT <- as.vector(de_fem_prog[de_fem_prog$qval<0.05,"genes"])

gene_list <- unique(c(prog_male_sig_gene_pseudoT, prog_female_sig_gene_pseudoT))

de_matrix <- log(prog_lineage_data[rownames(prog_lineage_data) %in% gene_list,]+1)


L1_lineage <- fem_prog_pseudotime

L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), L1_rpkm_exp_lineage))]


L2_lineage <- mal_prog_pseudotime

L2_ordered_lineage <- L2_lineage[order(L2_lineage, decreasing = TRUE)]
L2_rpkm_exp_lineage <- de_matrix[,names(L2_ordered_lineage)]
L2_cells <- L2_rpkm_exp_lineage[,order(match(names(L2_ordered_lineage), L2_rpkm_exp_lineage))]


L1_lineage_cells <- names(L1_ordered_lineage)
L2_lineage_cells <- names(L2_ordered_lineage)

comp_list <- comparelists(L1_lineage_cells, L2_lineage_cells)


common_cells <- comp_list$intersect
L1_spe_cells <- L1_lineage_cells[!L1_lineage_cells %in% comp_list$intersect]
L2_spe_cells <- L2_lineage_cells[!L2_lineage_cells %in% comp_list$intersect]


L1_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L1 cells", L1_spe_cells)
)
names(L1_cellLin) <- c(common_cells, L1_spe_cells)

L1_cellLin <- L1_cellLin[order(match(names(L1_cellLin), colnames(L1_cells)))]

L2_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L2 cells", L2_spe_cells)
)
names(L2_cellLin) <- c(common_cells, L2_spe_cells)
L2_cellLin <- L2_cellLin[order(match(names(L2_cellLin), colnames(L2_cells)))]


cellType_L1 <- all_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)

cellType_L2 <- all_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)


cellLin <- c(
	L2_cellLin,
	L1_cellLin
)

data_heatmap <- data.frame(
	L2_cells,
	L1_cells
)

cellType <- c(
	cellType_L2,
	cellType_L1
)



gene_clustering <- pheatmap(
	data_heatmap, 
	scale="row", 
	cutree_rows=15,
	clustering_method="ward.D2",
	silent=TRUE
)


clusters <- cutree(gene_clustering$tree_row, k = 15)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"

gene_cluster_palette <- c(
	"#a6cee3",
	"#1f78b4",
	"#b2df8a",
	"#33a02c",
	"#fb9a99",
	"#e31a1c",
	"#fdbf6f",
	"#ff7f00",
	"#cab2d6",
	"#6a3d9a",
	"#ffff99",
	"#b15928", 
	"#49beaa", 
	"#611c35", 
	"#2708a0"
)

gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
names(gene_cluster_colors) <- 1:max(clusters)

annotation_row <- clustering


annotation_col <- data.frame(
	# Cell_Clusters=clusters_in_lineage,
	# cellLineages=cellLin,
	cellType=cellType,
	Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)


names(cellLinCol) <- unique(cellLin)

cellTypeCol <- c(
	XY_C1="#457cff",
	XY_C2="#aeff01",
	XY_C3="#00c5ec",
	XY_C4="#009900",
	XY_C5="#a38cff",	
	XY_C6="#8dfaff",
	XX_C1="#560047",
	XX_C2="#a53bad",
	XX_C3="#eb6bac", 
	XX_C4="#ffa8a0"
)


annotation_colors <- list(
	cellType=cellTypeCol,
	# cellLineages=cellLinCol,
	clustering=gene_cluster_colors,
		Stages=c(
			E10.5="#2754b5",
			E11.5="#8a00b0", 
			E12.5="#d20e0f", 
			E13.5="#f77f05", 
			E16.5="#f9db21",
			P6="#43f14b"
		)
)

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026','#800026'))
mypalette <- c(rev(cold(13)), warm(12))
breaksList = seq(-2.5, 2.5, by = 0.2)






scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
    min = apply(x, 1, min, na.rm = T)
    max = apply(x, 1, max, na.rm = T)
    sdev = apply(x, 1, sd, na.rm = T)
    return((x - m) / sdev)
}

scaled_data_heatmap <- scale_rows(data_heatmap)


# tiff(file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_progenitor_heatmap_DE_genes_k_15_p-val_0.05.tiff"), res = 300, height = 21, width = 18, units = 'cm')

pdf(file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_progenitor_heatmap_DE_genes_male_female_k_20_p-val_0.05_span_0.4.pdf"),  height = 6, width = 7)

gene_clustering <- pheatmap(
	scaled_data_heatmap, 
	# data_heatmap,
	# scale="row",
	kmeans_k=25,
	# breaks=c(-0.5, 0, 0.5),
	gaps_col=length(cellType_L2),
	show_colnames=FALSE, 
	# show_rownames=FALSE, 
	cluster_cols=FALSE,
	# cluster_rows=FALSE,
	# cutree_rows=6,
	clustering_method="ward.D",
	# annotation_row=annotation_row,
	annotation_col=annotation_col,
	annotation_colors=annotation_colors,
	# cutree_rows=5,
	annotation_names_row=FALSE,
	# color=c(viridis(100)),
	color=mypalette,
	# color=seurat_palette,
	breaks=breaksList
)

dev.off()


write.csv(gene_clustering$kmeans$cluster, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_progenitor_heatmap_DE_genes_k_20_p-val_0.05.csv"))

# tf_dynamics <- clustering[clustering$Genes %in% tf_list$V1,]
# write.csv(tf_dynamics, file="male_female_supporting_cells_pseudotime_TF_new.csv")


