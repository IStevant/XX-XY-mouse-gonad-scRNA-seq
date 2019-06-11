#  If you need help, contact me at isabelle.stevant@gmail.com 

source("analysis_functions.R")
library("dendextend")

###########################################
#                                         #
#          Load and prepare files         #
#                                         #
###########################################

# Load the RPKM matrix
load(file="../data/female_rpkm.Robj")
# Remove the duplicated cell (experimental duplicate, i.e. cDNA from one cell were used twice for library prep as a library prep control)
female_rpkm <- female_rpkm[,!colnames(female_rpkm) %in% grep("rep",colnames(female_rpkm), value=TRUE)]

# Load the protein coding gene list
prot_coding_genes <- read.csv(file="../data/prot_coding.csv", row.names=1)
# Subset the protein coding genes from the RPKM matrix
females <- female_rpkm[rownames(female_rpkm) %in% as.vector(prot_coding_genes$x),]

# Load the read count matrix
load(file="../data/female_count.Robj")
# Remove the duplicated cell
female_count <- female_count[rownames(female_count) %in% rownames(females),!colnames(female_count) %in% grep("rep",colnames(female_count), value=TRUE)]

# Extract embryonic stage data from cell names
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
# Extract capture date from cell names
female_captures <- sapply(strsplit(colnames(females), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(females)


# Color palette
female_stagePalette <- 	c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21",
	"#43f14b"
)

# Remove gene never expressed
females <- females[rowSums(females)>0,]

###########################################
#										  #
#			Var Gene Selection		      #
#										  #
###########################################

# Extract genes with a squared coefficient of variation >2 times the fit regression (Brennecke et al 2013 method)
females_data <- getMostVarGenes(females, fitThr=2)
females_data <- log(females_data+1)

###########################################
#										  #
#			  RtSNE Analysis			  #
#										  #
###########################################

# Run a PCA on the highly variable genes
female_sub_pca <- PCA(
	t(females_data), 
	ncp = ncol(females_data), 
	graph=FALSE
)

# Estimate which PCs contain significant information with Jackstraw
# the result is 9 significant PCs
significant_pcs <- permutationPA(
	female_sub_pca$ind$coord, 
	B = 100, 
	threshold = 0.05, 
	verbose = TRUE, 
	seed = NULL
)$r

# Compute and plot the t-SNE using the significant PCs
female_t_sne <- run_plot_tSNE(
	pca=female_sub_pca,
	pc=significant_pcs,
	iter=5000,
	conditions=female_stages,
	colours=female_stagePalette
)

# Recompute PCA with FactomineR package for the clustering
res.pca <- PCA(
	t(females_data), 
	ncp = significant_pcs, 
	graph=FALSE
)

# Clustering cells based on the PCA with a minimum of 4 expected clusters
res.hcpc <- HCPC(
	res.pca, 
	graph = FALSE,
	min=4
	)

# Plot the hierarchical clustering result
plot(res.hcpc, choice ="tree", cex = 0.6)

# Extract clustering results
female_clustering <- res.hcpc$data.clust$clust
female_clustering <- paste("C", female_clustering, sep="")
names(female_clustering) <- rownames(res.hcpc$data.clust)

# Switch cluster names for convinience purpose (not very clean, quick and dirty method...)
female_clustering[female_clustering=="C1"] <- "C11"
female_clustering[female_clustering=="C2"] <- "C22"
female_clustering[female_clustering=="C22"] <- "C1"
female_clustering[female_clustering=="C11"] <- "C2"

# Cluster color palette
female_clusterPalette <- c(
		"#560047", 
		"#a53bad", 
		"#eb6bac", 
		"#ffa8a0"
	)

write.csv(female_clustering, file="../data/female_clustering.csv")

# Plot the t-SNE colored by cell clusters
female_t_sne_new_clusters <- plot_tSNE(
	tsne=female_t_sne, 
	conditions=female_clustering, 
	colours= female_clusterPalette
)

###########################################
#										  #
#	     DE genes between clusters	  	  #
#										  #
###########################################

# Prepare data to be loaded in Monocle
DE_female <- prepare_for_DE (
	female_count, 
	female_clustering, 
	female_stages
)

# Get genes diffet=rentially expressed between clusters (Monocle)
female_DE_genes <- findDEgenes(
	DE_female, 
	qvalue=0.05
)

# Get in which cluster the DE genes ar ethe most expressed (see analysis_functions.R for details)
de_clusters <- get_up_reg_clusters(
	females, 
	female_clustering, 
	female_DE_genes
)

# Write the result in a file
write.csv(
	de_clusters, 
	quote = FALSE, 
	file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_DE_genes_per_clusters_4_groups.csv")
)

###########################################
#										  #
#	  		   GO term analysis	  		  #
#										  #
###########################################

# For details see "GOSemSim" package manual

# Extract DE gene names
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes

# Convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,c("genes", "cluster")]
de_gene_clusters <- data.frame(
	ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
	cluster=de_gene_clusters$cluster
)

list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, de_gene_clusters$cluster)


# Run full GO enrichment test
formula_res <- compareCluster(
	ENTREZID~cluster, 
	data=de_gene_clusters, 
	fun="enrichGO", 
	OrgDb="org.Mm.eg.db",
	ont		   = "BP",
	pAdjustMethod = "BH",
	pvalueCutoff  = 0.01,
	qvalueCutoff  = 0.05
)

# Run GO enrichment test and merge terms that are close to each other to remove result redundancy
lineage1_ego <- simplify(
	formula_res, 
	cutoff=0.5, 
	by="p.adjust", 
	select_fun=min
)

# Plot both analysis results
dotplot(formula_res, showCategory=5)
dotplot(lineage1_ego, showCategory=5)

# Save results
write.csv(formula_res@compareClusterResult, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster.csv"))
write.csv(lineage1_ego@compareClusterResult, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster_simplified.csv"))


###########################################
#										  #
#		  Plot stats about cells		  #
#										  #
###########################################

# PCA using all the genes
female_pca <- prcomp(
	t(log(females+1)), 
	center=TRUE, 
	scale=TRUE
)

# Plot PCA colored by embryonic stages
plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=female_stages, 
	colours=female_stagePalette
)
# Plot PCA coored by cell clusters
plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=female_clustering, 
	colours=female_clusterPalette
)

# Count how many genes are expressed in total
gene_per_cell_females <- females
gene_per_cell_females <- rowSums(gene_per_cell_females)
gene_per_cell_females[gene_per_cell_females>0] <- 1
sum(gene_per_cell_females)



# Count how many genes are expressed per cell
gene_per_cell_females <- females
gene_per_cell_females[gene_per_cell_females>0] <- 1
gene_per_cell_females <- colSums(gene_per_cell_females)
gene_per_cell_females <- data.frame(
	cells=names(gene_per_cell_females),
	geneNb = gene_per_cell_females
)

# Plot the distribution of the number of expressed genes per cell
ggplot(gene_per_cell_females, aes(geneNb)) +
	geom_histogram(color="black", fill="grey") +
	theme_bw() +
	labs(x = "Detected genes", y="Cell count")+
	ggtitle(paste("Median=", median(gene_per_cell_females$geneNb)," genes per cell", sep="")) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)



# Endogenous Sf1 (Nr5a1) expression vs Sf1-GFP transgene expression
sf1_vs_gfp <- as.data.frame(
	t(
		log(females[c("eGFP", "Nr5a1"),]+1)
	)
)
# Correlation between both expression
cor(sf1_vs_gfp, method="spearman") # 0.838
cor(sf1_vs_gfp, method="pearson") # 0.781

# Plot endogenous Sf1 (Nr5a1) expression vs Sf1-GFP transgene expression and the linear regression line
ggplot(sf1_vs_gfp, aes(x=Nr5a1, y=eGFP)) +
	geom_point(color="black") +
	geom_smooth(method=lm, color="red", se = FALSE) +
	theme_bw() +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)


# Plot PCAs of the replicated stages (E11.5 to P6) to verify if any batch effect

# E11.5
e11.5 <- females[,colnames(females) %in% names(female_stages[female_stages=="E11.5"])]
e11.5 <- e11.5[rowSums(e11.5)>0,]
e11.5_captures <- as.factor(sapply(strsplit(colnames(e11.5), "_"), `[`, 3))
levels(e11.5_captures) <- c("Capture 1", "Capture 2")
e11.5_palette <- c("#0000ff", "#8a00b0")

female_pca <- prcomp(
	t(log(e11.5+1)), 
	center=TRUE, 
	scale=TRUE
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=e11.5_captures, 
	colours=e11.5_palette,
	title="E11.5"
)

# E12.5
e12.5 <- females[,colnames(females) %in% names(female_stages[female_stages=="E12.5"])]
e12.5 <- e12.5[rowSums(e12.5)>0,]
e12.5_captures <- as.factor(sapply(strsplit(colnames(e12.5), "_"), `[`, 3))
levels(e12.5_captures) <- c("Capture 1", "Capture 2")
e12.5_palette <- c("#ff0000", "#ff7674")

female_pca <- prcomp(
	t(log(e12.5+1)), 
	center=TRUE, 
	scale=TRUE
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=e12.5_captures, 
	colours=e12.5_palette,
	title="E12.5"
)

#E13.5
e13.5 <- females[,colnames(females) %in% names(female_stages[female_stages=="E13.5"])]
e13.5 <- e13.5[rowSums(e13.5)>0,]
e13.5_captures <- as.factor(sapply(strsplit(colnames(e13.5), "_"), `[`, 3))
levels(e13.5_captures) <- c("Capture 1", "Capture 2")
e13.5_palette <- c("#803c00", "#f77f05")

female_pca <- prcomp(
	t(log(e13.5+1)), 
	center=TRUE, 
	scale=TRUE
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=e13.5_captures, 
	colours=e13.5_palette,
	title="E13.5"
)

# E16.5
e16.5 <- females[,colnames(females) %in% names(female_stages[female_stages=="E16.5"])]
e16.5 <- e16.5[rowSums(e16.5)>0,]
e16.5_captures <- as.factor(sapply(strsplit(colnames(e16.5), "_"), `[`, 3))
levels(e16.5_captures) <- c("Capture 1", "Capture 2")
e16.5_palette <- c("#7ca328", "#f9db21")

female_pca <- prcomp(
	t(log(e16.5+1)), 
	center=TRUE, 
	scale=TRUE
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=e16.5_captures, 
	colours=e16.5_palette,
	title="E16.5"
)

# P6
P6 <- females[,colnames(females) %in% names(female_stages[female_stages=="P6"])]
P6 <- P6[rowSums(P6)>0,]
P6_captures <- as.factor(sapply(strsplit(colnames(P6), "_"), `[`, 3))
levels(P6_captures) <- c("Capture 1", "Capture 2")
P6_palette <- c("#21c934", "#2f707c")

female_pca <- prcomp(
	t(log(P6+1)), 
	center=TRUE, 
	scale=TRUE
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=P6_captures, 
	colours=P6_palette,
	title="P6"
)

###########################################
#										  #
#		  heatmap_marker_genes		      #
#										  #
###########################################

markerGenes <- c(
	"Nr2f1",
	"Nr2f2",
	"Maf",
	"Foxl2",
	"Rspo1",
	"Lgr5",
	"Bmp2",
	"Runx1",
	"Amhr2",
	"Kitl",
	"Fst",
	"Esr2",
	"Amh",
	"Ptges"
)

gene_subset <- as.matrix(log(females[rownames(females) %in% markerGenes,]+1))

cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C3"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C4"])]

heatmap_gene_subset <- cbind(
	cl1_gene_subset, 
	cl2_gene_subset,
	cl3_gene_subset,
	cl4_gene_subset
)

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
rowbreaks <- c(6, 15)

colbreaks <- c(
	ncol(cl1_gene_subset),
	ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
	ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)
)

cluster_color <- c(
		C1="#560047",
		C2="#a53bad", 
		C3="#eb6bac", 
		C4="#ffa8a0"
)

stage_color=c(
	E10.5="#2754b5", 
	E11.5="#8a00b0", 
	E12.5="#d20e0f", 
	E13.5="#f77f05", 
	E16.5="#f9db21",
	P6="#43f14b"
)

tiff(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_de_genes_per_cell_clusters_small.tiff"), res = 300, height = 10, width = 16, units = 'cm')
	plot_heatmap_2(
		heatmap_gene_subset, 
		female_clustering, 
		female_stages, 
		rowbreaks, 
		colbreaks,
		cluster_color,
		stage_color
	)
dev.off()



# Marker gene expression per cluster as violin plot
pdf(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_violin_markers.pdf"), width=10, height=22)

require(gridExtra)

p <- list()
for (genes in markerGenes) {
	p[[genes]] <- violin_gene_exp(
				genes, 
				females, 
				female_clustering, 
				female_clusterPalette
			)
}

do.call(grid.arrange,c(p, ncol=3))

dev.off()

# Marker gene expression on the t-SNE plot
pdf(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_tSNE_markers.pdf"), width=16, height=28)

require(gridExtra)

p <- list()
for (genes in markerGenes) {
	p[[genes]] <- tsne_gene_exp(
				female_t_sne,
				genes, 
				females
			)
}

do.call(grid.arrange,c(p, ncol=4))

dev.off()

###########################################
#										  #
#		  Heatmap top DE genes		      #
#										  #
###########################################

# DE genes
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.0001)
gene_names <- gene_names$genes

gene_names <- get_top_up_reg_clusters(de_clusters, 20)
gene_subset <- as.matrix(log(females[rownames(females) %in% gene_names,]+1))

cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C3"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C4"])]

heatmap_gene_subset <- cbind(
	cl1_gene_subset, 
	cl2_gene_subset,
	cl3_gene_subset,
	cl4_gene_subset
)

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
rowbreaks <- c(6, 15)

colbreaks <- c(
	ncol(cl1_gene_subset),
	ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
	ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)
)

cluster_color <- c(
		C1="#560047",
		C2="#a53bad", 
		C3="#eb6bac", 
		C4="#ffa8a0"
)

stage_color=c(
	E10.5="#2754b5", 
	E11.5="#8a00b0", 
	E12.5="#d20e0f", 
	E13.5="#f77f05", 
	E16.5="#f9db21",
	P6="#43f14b"
)

tiff(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_de_genes_per_cell_clusters_small.tiff"), res = 300, height = 10, width = 16, units = 'cm')
	plot_heatmap_2(
		heatmap_gene_subset, 
		female_clustering, 
		female_stages, 
		rowbreaks, 
		colbreaks,
		cluster_color,
		stage_color
	)
dev.off()

###########################################
#										  #
#		   Pseudotime Analysis		  	  #
#										  #
###########################################

# Compute the Diffusion map
female_dm <- run_diffMap(
	females_data, 
	female_clustering,
	sigma=15
)

# Plot the Eigen values per diffusion component (similar to screeplots for PCAs)
plot_eigenVal(
	dm=female_dm
)

# Plot the 3D diffusion map using DCs 1 to 3 colored by the cell clustering
plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_clustering, 
	colour=female_clusterPalette
)

# Plot the 3D diffusion map using DCs 1 to 3 colored by embryonic stages
plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_stages, 
	colour=female_stagePalette
)

# Compute the pseudotime using Slingshot based on the informative diffusion map components (DCs)
# We see on the eigenVal plot a first elbow between DC 3 and 4, so we consider the first 4 DCs for further analysis
# We also help the pseudotime calculation by giving the starting cliuster and the end clusters
female_lineage <- get_lineage(
	dm=female_dm, 
	dim=c(1:4), 
	condition=factor(female_clustering),
	start="C1",
	end=c("C2", "C4"),
	shrink.method="cosine"
)

# Plot the 3D diffusion map using DCs 1 to 3 colored by embryonic stages
# The black line is the cell lineage reconstruction by Slinshot
plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_stages, 
	colour=female_stagePalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file (first orient the graph as you want, the svg save does a snapshot of the graph the way you orient it)
rgl.postscript( file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_dm_stages_cols.svg"), fmt = "svg", drawText = TRUE )


# Plot the 3D diffusion map using DCs 1 to 3 colored by cell clusters
plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_clustering, 
	colour=female_clusterPalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_dm_clusters_cols.svg"), fmt = "svg", drawText = TRUE )

# Attribute cells to lineages with a distance of 0.9 from the smoothed line and compute the pseudotime
female_pseudotime <- get_pseudotime(female_lineage, wthres=0.9)
rownames(female_pseudotime) <- colnames(females)

write.csv(female_pseudotime, file="../data/female_pseudotime.csv")


# Make the pseudotime of each lineage between 0 and 100 (percentage of cell progression over the pseudotime)
# thos makes the two cell lineage progression comparable
pseudotime_lin <- female_pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime

pseudotime_lin <- female_pseudotime[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime

female_pseudotime[,"curve1"] <- pseudotime_lin1_percent
female_pseudotime[,"curve2"] <- pseudotime_lin2_percent


## Lineage result: 
# Lineage 1 -> progenitors to Granulisa
# Lineage 2 -> progenitors to progenitors

####################################################
#
# Plot gene expression along pseudotime in cell lineages 1 and 2
#
####################################################

# Twisted palette to plot gene expression by lineages along pseudotime
female_clusterPalette2 <- c(
	"#ff6663", 
	"#3b3561"
	)

# Plot one gene expression over pseudotime for the lineage 1
plot_smoothed_gene_per_lineage(
	rpkm_matrix=females, 
	pseudotime=female_pseudotime, 
	lin=c(1),
	gene="Amhr2", 
	stages=female_stages, 
	clusters=female_clustering, 
	stage_colors=female_stagePalette,
	cluster_colors=female_clusterPalette,
	lineage_colors=female_clusterPalette2
)

# Plot one gene expression over pseudotime for the lineage 2
plot_smoothed_gene_per_lineage(
	rpkm_matrix=females, 
	pseudotime=female_pseudotime, 
	lin=c(2),
	gene="Amhr2", 
	stages=female_stages, 
	clusters=female_clustering, 
	stage_colors=female_stagePalette,
	cluster_colors=female_clusterPalette,
	lineage_colors=female_clusterPalette2
)

# Plot one gene expression over pseudotime for the lineage 1 and 2
plot_smoothed_gene_per_lineage(
	rpkm_matrix=females, 
	pseudotime=female_pseudotime, 
	lin=c(1,2),
	gene="Amhr2", 
	stages=female_stages, 
	clusters=female_clustering, 
	stage_colors=female_stagePalette,
	cluster_colors=female_clusterPalette,
	lineage_colors=female_clusterPalette2
)


# List of genes of interest to plot lineage by lineage through pseudotime
gene_list <- c(
	"Sall4",
	"Sox11",
	"Gata4",
	"Lgr5",
	"Runx1",
	"Foxl2",
	"Hey2",
	"Wnt5a",
	"Pdgfra",
	"Nr2f2",
	"Sfrp1",
	"Ifitm3",
	"Ptch1",
	"Wnt4",
	"Rspo1",
	"Cdkn1b",
	"Gli1",
	"Tcf21",
	"Nr0b1",
	"Nr0b2",
	"Nr5a1",
	"Nr6a1"
)

plot_smoothed_genes <- function(genes, lin){
	female_clusterPalette2 <- c("#ff6663", "#3b3561")
	for (gene in genes){
		plot_smoothed_gene_per_lineage(
			rpkm_matrix=females, 
			pseudotime=female_pseudotime, 
			lin=lin,
			gene=gene, 
			stages=female_stages, 
			clusters=female_clustering, 
			stage_colors=female_stagePalette,
			cluster_colors=female_clusterPalette,
			lineage_colors=female_clusterPalette2
		)
	}
}

pdf(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_interesting_profiles.pdf"), width=4, height=4)
	plot_smoothed_genes(gene_list, 1) # plot only lineage 1
	plot_smoothed_genes(gene_list, 2) # plot only lineage 2
	plot_smoothed_genes(gene_list, c(1,2)) # plot the two moleages in the same graph to see the divergence
dev.off()

######################################################
#										             #
#   Genes dynamically expressed through pseudotime   #
#										             #
######################################################

# Using Monocle, calculate which gene display a dynamic expression pattern through the pseudotime for each lineage
# Results are stored n files to avoid to recompute anytime you want to re-run the heatmap

female_lineage1_sig_gene_pseudoT <- get_var_genes_pseudotime(
	females, 
	female_count, 
	female_pseudotime, 
	lineageNb=1, 
	female_clustering
)
female_lineage1_sig_gene_pseudoT <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage1_sig_gene_pseudoT, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_lineage1_pseudotime_DE_genes.csv"))

#############################################################
female_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
	females, 
	female_count, 
	female_pseudotime, 
	lineageNb=2, 
	female_clustering
)
female_lineage2_sig_gene_pseudoT <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage2_sig_gene_pseudoT, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_lineage2_pseudotime_DE_genes.csv"))


###########################################
#										  #
#		Double heatmap lineage			  #
#										  #
###########################################

# Load dynamic gene list per lineage if not in your R session (results from the paper, load your own files intead if you have generated them)
female_lineage1_sig_gene_pseudoT <- read.csv(file="../data/female_lineage1_pseudotime_DE_genes.csv", row.names=1)
female_lineage2_sig_gene_pseudoT <- read.csv(file="../data/female_lineage2_pseudotime_DE_genes.csv", row.names=1)

# Extract the significant genes
female_lineage1_clustering <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
female_lineage2_clustering <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]

# Get gene names
gene_list <- unique(rownames(female_lineage1_clustering), rownames(female_lineage2_clustering))

# Get RPKM matrix of the dynamic genes
de_matrix <- log(females[rownames(females) %in% gene_list,]+1)

# Get cells attributed to the lineage 1
L1_lineage <- female_pseudotime[!is.na(female_pseudotime[,1]),1]
# Order the cells through the pseudotime
L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
# Get expression values of the dynamic genes in the cells from lineage 1
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
# Order the cells (colums of the expression matrix) through the pseudotime
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), L1_rpkm_exp_lineage))]

# Same as previously but for the linage 2 and but cells are ordered decreasingly (to be ploted on the left side)
L2_lineage <- female_pseudotime[!is.na(female_pseudotime[,2]),2]
L2_ordered_lineage <- L2_lineage[order(L2_lineage, decreasing = TRUE)]
L2_rpkm_exp_lineage <- de_matrix[,names(L2_ordered_lineage)]
L2_cells <- L2_rpkm_exp_lineage[,order(match(names(L2_ordered_lineage), L2_rpkm_exp_lineage))]

# Extract cell names
L1_lineage_cells <- names(L1_ordered_lineage)
L2_lineage_cells <- names(L2_ordered_lineage)

# Compare nales in the two lisis
comp_list <- comparelists(L1_lineage_cells, L2_lineage_cells)
# As the two cell lineages share cells before the lineages separate, we want to know which cells are in common and which are lineage specific to label them
common_cells <- comp_list$intersect
L1_spe_cells <- L1_lineage_cells[!L1_lineage_cells %in% comp_list$intersect]
L2_spe_cells <- L2_lineage_cells[!L2_lineage_cells %in% comp_list$intersect]

# label the cells if they are part of the common branch or only in the lineage 1 branch
L1_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L1 cells", L1_spe_cells)
)
names(L1_cellLin) <- c(common_cells, L1_spe_cells)
L1_cellLin <- L1_cellLin[order(match(names(L1_cellLin), colnames(L1_cells)))]

# label the cells if they are part of the common branch or only in the lineage 2 branch
L2_cellLin <- c(
	rep_along("common cells", common_cells), 
	rep_along("L2 cells", L2_spe_cells)
)
names(L2_cellLin) <- c(common_cells, L2_spe_cells)
L2_cellLin <- L2_cellLin[order(match(names(L2_cellLin), colnames(L2_cells)))]

# Get the cell cluster of the cells for each cell lineage
cellType_L1 <- female_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)
cellType_L2 <- female_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)

# Merge cell lineages annotation
cellLin <- c(
	L2_cellLin,
	L1_cellLin
)

# Merge the cell type annotation
cellType <- c(
	cellType_L2,
	cellType_L1
)

# Smooth the gene expression (Loess regression)
L2_cells_smooth <- smooth_gene_exp(
	L2_cells, 
	L2_ordered_lineage, 
	span=0.4
)
L1_cells_smooth <- smooth_gene_exp(
	L1_cells, 
	L1_ordered_lineage, 
	span=0.4
)

# Merge the smoothed expression matrices
data_heatmap <- data.frame(
	L2_cells_smooth,
	L1_cells_smooth
)

# Run pheatmap to cluster the genes per expression pattern
set.seed(123)
gene_clustering <- pheatmap(
	data_heatmap, 
	scale="row", 
	clustering_method="ward.D",
	silent=TRUE
)

# Extract the gene pattern clusters with k-means (17 expected clusters)
clusters <- cutree(gene_clustering$tree_row, k = 17)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"

# Save clustering results
write.csv(clustering, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_mean_norm_female_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k17_scaled.csv"))

# Palette for the heatmap (pheatmap does not take the color palette into account, I don't know why but I kept the code)
gene_cluster_palette <- c(
	'#a6cee3',
	'#1f78b4',
	'#b2df8a',
	'#33a02c',
	'#fb9a99',
	'#e31a1c',
	'#fdbf6f',
	'#ff7f00',
	'#cab2d6',
	'#6a3d9a',
	'#ffff99',
	'#b15928', 
	'#49beaa', 
	'#611c35', 
	'#2708a0',
	'#fccde5',
	'#bc80bd'
)
gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
names(gene_cluster_colors) <- 1:max(clusters)


# Annotate the row of the hearmap with the gene clusters
annotation_row <- data.frame(clustering=clustering)

# Annotation of the colums of the heatmap with the cell lineage, the cell type and the embryonic stage of the cells
annotation_col <- data.frame(
	cellLineages=cellLin,
	cellType=cellType,
	Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)

# Colors for the cell lineages
cellLinCol <- c(
	"#3b3561", 
	"#c8c8c8", 
	"#ff6663"
)
names(cellLinCol) <- unique(cellLin)

# Colors for the cell clusters
cellTypeCol <- c(
	C2="#a53bad", 
	C1="#560047", 
	C3="#eb6bac", 
	C4="#ffa8a0"
)
names(cellTypeCol) <- unique(cellType)

# Legend for all the aannotation to plot on pheatmap
annotation_colors <- list(
	cellType=cellTypeCol,
	cellLineages=cellLinCol,
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

# Color palette for gene expression
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
breaksList = seq(-2.2, 2.5, by = 0.2)

# Plot the heatmap
tiff(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_heatmap_DE_genes_granulosa_progenitors_k_17_pval_005.tiff"), res = 300, height = 21, width = 18, units = 'cm')
gene_clustering <- pheatmap(
	data_heatmap, 
	scale="row",
	gaps_col=length(cellType_L2),
	show_colnames=FALSE, 
	show_rownames=FALSE, 
	cluster_cols=FALSE,
	clustering_method="ward.D",
	annotation_row=annotation_row,
	annotation_col=annotation_col,
	annotation_colors=annotation_colors,
	cutree_rows=17, 
	annotation_names_row=FALSE,
	color=mypalette
)
dev.off()


###########################################
#										  #
#	  		   GO term analysis	  		  #
#										  #
###########################################

# Load manually annotated gene clustering results (gene custers classified as a, b, c... according to their expression pattern redundandy, see paper fig. 2)
dyn_genes <- read.csv(file="../data/female_lineages_DE_gene_pseudotime_clustered_annotated.csv")
gene_names <- dyn_genes$Genes

#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# Subset genes with a correnc entrez ID
gene_clusters <- dyn_genes[dyn_genes$Genes %in% entrez_genes$SYMBOL,,drop=FALSE]

# Generate data frame for GoSemSim
de_gene_clusters <- data.frame(
	ENTREZID=entrez_genes[!duplicated(entrez_genes$SYMBOL),"ENTREZID"],
	Gene_Clusters=gene_clusters$Gene.categories
)

# Run enrichment analysis
formula_res <- compareCluster(
	ENTREZID~Gene_Clusters, 
	data=de_gene_clusters, 
	fun="enrichGO", 
	OrgDb="org.Mm.eg.db",
	ont		   = "BP",
	pAdjustMethod = "BH",
	pvalueCutoff  = 0.05,
	qvalueCutoff  = 0.05
)

# Run simplified GO enrochment analysis
lineage1_ego <- simplify(
	formula_res, 
	cutoff=0.5, 
	by="p.adjust", 
	select_fun=min
)

write.csv(formula_res@compareClusterResult, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster.csv"))
write.csv(lineage1_ego@compareClusterResult, file=paste0("../results/", format(Sys.time(), "%Y%m%d"), "_female_compared_symplified_GO_term_DE_cluster.csv"))

pdf(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_GO_term_DE_genes_clusters.pdf"), width=11, height=8)
	dotplot(formula_res, showCategory=3)+ theme(aspect.ratio=0.8)
	dotplot(lineage1_ego, showCategory=3)+ theme(aspect.ratio=2)
dev.off()