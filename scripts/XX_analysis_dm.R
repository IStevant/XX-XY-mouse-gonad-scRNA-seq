#### WARNING ####
#
#  This script is my working script, thus it is not yet cleaned and contains lots of rubish code and plot tests.
#  If you need help, contact me at isabelle.stevant@gmail.com 

source("analysis_functions.R")
library("dendextend")

###########################################
#                                         #
#          Load and prepare files         #
#                                         #
###########################################

load(file="../data/female_rpkm.Robj")
female_rpkm <- female_rpkm[,!colnames(female_rpkm) %in% grep("rep",colnames(female_rpkm), value=TRUE)]

prot_coding_genes <- read.csv(file="../data/prot_coding.csv", row.names=1)
females <- female_rpkm[rownames(female_rpkm) %in% as.vector(prot_coding_genes$x),]

load(file="../data/female_count.Robj")
female_count <- female_count[rownames(female_count) %in% rownames(females),!colnames(female_count) %in% grep("rep",colnames(female_count), value=TRUE)]

female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
female_captures <- sapply(strsplit(colnames(females), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(females)



female_stagePalette <- 	c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21",
	"#43f14b"
)
females_data <- log(females+1)

females <- females[rowSums(females)>0,]

###########################################
#										  #
#			Var Gene Selection		      #
#										  #
###########################################
pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_clustering.pdf"), width=5, height=5)

females_data <- getMostVarGenes(females, fitThr=2)
females_data <- log(females_data+1)

###########################################
#										  #
#			  RtSNE Analysis			  #
#										  #
###########################################

female_sub_pca <- PCA(
	t(females_data), 
	ncp = ncol(females_data), 
	graph=FALSE
)


######## Long to compute, but the result is 9 significant PCs

# significant_pcs <- permutationPA(
# 	female_sub_pca$ind$coord, 
# 	B = 100, 
# 	threshold = 0.05, 
# 	verbose = TRUE, 
# 	seed = NULL
# )$r
significant_pcs <- 9

female_t_sne <- plot_tSNE(
	pca=female_sub_pca,
	pc=significant_pcs,
	iter=5000,
	conditions=female_stages,
	colours=female_stagePalette
)

write.csv(female_t_sne, "female_t-sne.csv")

res.pca <- PCA(
	t(females_data), 
	ncp = significant_pcs, 
	graph=FALSE
)

res.hcpc <- HCPC(
	res.pca, 
	graph = FALSE,
	min=4
	# consol=FALSE
	)

# plot(res.hcpc, choice ="tree", cex = 0.6)

female_clustering <- res.hcpc$data.clust$clust
female_clustering <- paste("C", female_clustering, sep="")
names(female_clustering) <- rownames(res.hcpc$data.clust)

# Switch cluster names for convinience purpose (not very clean, quick and dirty method...)
female_clustering[female_clustering=="C1"] <- "C11"
female_clustering[female_clustering=="C2"] <- "C22"
female_clustering[female_clustering=="C22"] <- "C1"
female_clustering[female_clustering=="C11"] <- "C2"


female_clusterPalette <- c(
		"#560047", 
		"#a53bad", 
		"#eb6bac", 
		"#ffa8a0"
	)

#write.csv(female_clustering, file="../data/180712_female_clustering_final.csv")

female_t_sne_new_clusters <- plot_tSNE_2(
	tsne=female_t_sne, 
	conditions=female_stages, 
	colours= female_stagePalette
)

female_t_sne_new_clusters <- plot_tSNE_2(
	tsne=female_t_sne, 
	conditions=female_clustering, 
	colours= female_clusterPalette
)

dev.off()


###########################################
#										  #
#			     DE Analysis			  #
#										  #
###########################################


DE_female <- prepare_for_DE (
	female_count, 
	female_clustering, 
	female_stages
)


female_DE_genes <- findDEgenes(
	DE_female, 
	qvalue=0.05
)

de_clusters <- get_up_reg_clusters(
	females, 
	female_clustering, 
	female_DE_genes
)

write.csv(
	de_clusters, 
	quote = FALSE, 
	file=paste0("../data/", format(Sys.time(), "%Y%m%d"), "_female_DE_genes_per_clusters_4_groups.csv")
)

# GO terms

de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes

#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,c("genes", "cluster")]

de_gene_clusters <- data.frame(
	ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
	cluster=de_gene_clusters$cluster
)

list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, de_gene_clusters$cluster)

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

lineage1_ego <- simplify(
	formula_res, 
	cutoff=0.5, 
	by="p.adjust", 
	select_fun=min
)

dotplot(lineage1_ego, showCategory=5)

write.csv(formula_res@compareClusterResult, file=paste0("../data/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster.csv"))

pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_GO_term_DE_genes_clusters.pdf"), width=7, height=6)
dotplot(formula_res, showCategory=5)
dotplot(lineage1_ego, showCategory=5)
dev.off()

write.csv(lineage1_ego@compareClusterResult, file=paste0("../data/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster_simplified.csv"))


###########################################
#										  #
#		  Plot stats about cells		  #
#										  #
###########################################

female_pca <- prcomp(
	t(log(females+1)), 
	center=TRUE, 
	scale=TRUE
)

pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_graphs_qc_pcq.pdf"), width=5)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=female_stages, 
	colours=female_stagePalette
)

plot_pca(
	pca=female_pca, 
	pc=2, 
	conditions=female_clustering, 
	colours=female_clusterPalette
)

dev.off()


# rle(female_captures)
# rle(female_stages)
# rle(paste(female_stages, female_captures))
# rle(female_clustering[order(female_clustering)])


gene_per_cell_females <- females
gene_per_cell_females <- rowSums(gene_per_cell_females)
gene_per_cell_females[gene_per_cell_females>0] <- 1

sum(gene_per_cell_females)



gene_per_cell_females <- data.frame(
	cells=names(gene_per_cell_females),
	geneNb = gene_per_cell_females
)


gene_per_cell_females <- females
gene_per_cell_females[gene_per_cell_females>0] <- 1
gene_per_cell_females <- colSums(gene_per_cell_females)

gene_per_cell_females <- data.frame(
	cells=names(gene_per_cell_females),
	geneNb = gene_per_cell_females
)

median(gene_per_cell_females$geneNb)




sf1_vs_gfp <- as.data.frame(
	t(
		log(females[c("eGFP", "Nr5a1"),]+1)
	)
)

cor(sf1_vs_gfp, method="spearman") # 0.838
cor(sf1_vs_gfp, method="pearson") # 0.781

pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_gene_per_cell_gfp_sf1.pdf"), width=5)

ggplot(gene_per_cell_females, aes(geneNb)) +
	geom_histogram(color="black", fill="grey") +
	theme_bw() +
	labs(x = "Detected genes", y="Cell count")+
	# ggtitle(paste("Median=", median(gene_per_cell_females$geneNb)," genes per cell", sep="")) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)

ggplot(sf1_vs_gfp, aes(x=Nr5a1, y=eGFP)) +
	geom_point(color="black") +
	# geom_smooth(method=lm, color="red", se = FALSE) +
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

dev.off()



pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_graphs_qc_pca_batch.pdf"), width=4, height=5)

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

dev.off()


###########################################
#										  #
#		  heatmap_marker_genes		      #
#										  #
###########################################

markerGenes <- c(
	# "eGFP",
	# "Nr5a1",
	# "Gata4",
	"Nr2f1",
	# "Cbx2",
	# "Actb",
	"Nr2f2",
	"Maf",
	# "Lhx9",
	# "Pdgfra",
	# "Mafb",
	# "Gli1",
	# "Ptch1",
	# "Wnt4",
	"Foxl2",
	"Rspo1",
	"Lgr5",
	"Bmp2",
	# "Gata2",
	# "Gli1",
	# "Wt1",
	# "Cdkn1a",
	# "Cdkn1b",
	# "Cdkn1c",
	# "Axin2",
	"Runx1",
	"Amhr2",
	"Kitl",
	"Fst",
	"Esr2",
	"Amh",
	"Ptges"
	# "Hsd17b1",
	# "Hsd3b1"
	# "Cyp11a1"
	# "Cyp17a1",
	# "Cyp19a1",
	# "Fshr",
	# "Lhcgr"
	# "Insl3"
)

# marker_gene_list <- read.csv(file="../data/marker_gene_list.txt", header=TRUE)
# markerGenes <- marker_gene_list$gene


gene_subset <- as.matrix(log(females[rownames(females) %in% markerGenes,]+1))

# DE genes
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.0001)
gene_names <- gene_names$genes

gene_names <- get_top_up_reg_clusters(de_clusters, 20)



gene_subset <- as.matrix(log(females[rownames(females) %in% gene_names,]+1))

# gene_subset <- gene_subset[order(match(rownames(gene_subset), markerGenes)),]

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

tiff(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_de_genes_per_cell_clusters_small.tiff"), res = 300, height = 10, width = 16, units = 'cm')
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


pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_de_genes_per_cell_clusters.pdf"), height = 900, width = 400)
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


pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_violin_markers.pdf"), width=10, height=22)

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



markerGenes <- c(
	"Ptges"
)

pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_tSNE_markers.pdf"), width=16, height=28)

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




pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_GFP_violin.pdf"), width=6, height=6)
g <- violin_gene_exp(
			"eGFP", 
			females, 
			female_stages, 
			female_stagePalette,
			test=FALSE
		)
print(g)
dev.off()

g <- tsne_gene_exp(
			female_t_sne,
			"eGFP", 
			females
		)
print(g)


g <- tsne_gene_per_cell(
	female_t_sne,
	females
)
print(g)





###########################################
#										  #
#		   Pseudotime Analysis		  	  #
#										  #
###########################################


female_dm <- run_diffMap(
	females_data, 
	female_clustering,
	sigma=15
)


plot_eigenVal(
	dm=female_dm
)


plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_clustering, 
	colour=female_clusterPalette
)


plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_stages, 
	colour=female_stagePalette
)

plot_dm_2D(
	dm=female_dm, 
	dc=2,
	condition=female_stages, 
	colour=female_stagePalette
)

plot_dm_2D(
	dm=female_dm, 
	dc=3,
	condition=female_stages, 
	colour=female_stagePalette
)



female_lineage <- get_lineage(
	dm=female_dm, 
	dim=c(1:4), 
	condition=factor(female_clustering),
	start="C1",
	end=c("C2", "C4")
)


plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_stages, 
	colour=female_stagePalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
rgl.postscript( file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_dm_stages_cols.svg"), fmt = "svg", drawText = TRUE )



plot_dm_3D(
	dm=female_dm, 
	dc=c(1:3),
	condition=female_clustering, 
	colour=female_clusterPalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
rgl.postscript( file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_dm_clusters_cols.svg"), fmt = "svg", drawText = TRUE )


female_pseudotime <- get_pseudotime(female_lineage, wthres=0.9)
rownames(female_pseudotime) <- colnames(females)


pseudotime_lin <- female_pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime

pseudotime_lin <- female_pseudotime[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime

female_pseudotime[,"curve1"] <- pseudotime_lin1_percent
female_pseudotime[,"curve2"] <- pseudotime_lin2_percent



female_clusterPalette2 <- c("#ff6663", "#3b3561")

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

plot_genes_tsne <- function(genes){
	for (gene in genes){
		print(gene)
		tsne_gene_exp(
			female_t_sne,
			gene, 
			females
		)
	}
}

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


pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_interesting_profiles.pdf"), width=4, height=4)
	plot_smoothed_genes(gene_list, c(1,2)) # plot the two moleages in the same graph to see the divergence
	# plot_smoothed_genes(gene_list, 1) # plot only lineage 1
	# plot_smoothed_genes(gene_list, 2) # plot only lineage 2
dev.off()


###########################################
#										  #
#			DE genes lineages			  #
#										  #
###########################################

female_lineage1_sig_gene_pseudoT <- get_var_genes_pseudotime(
	females, 
	female_count, 
	female_pseudotime, 
	lineageNb=1, 
	female_clustering
)
female_lineage1_sig_gene_pseudoT <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage1_sig_gene_pseudoT, file=paste0("../data/", format(Sys.time(), "%Y%m%d"), "_female_lineage1_pseudotime_DE_genes.csv"))

#############################################################
female_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
	females, 
	female_count, 
	female_pseudotime, 
	lineageNb=2, 
	female_clustering
)
female_lineage2_sig_gene_pseudoT <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage2_sig_gene_pseudoT, file=paste0("../data/", format(Sys.time(), "%Y%m%d"), "_female_lineage2_pseudotime_DE_genes.csv"))


###########################################
#										  #
#		Double heatmap lineage			  #
#										  #
###########################################

# female_pseudotime <- read.csv(file="../graph/female/female_pseudotime.csv", row.names=1)

female_lineage1_sig_gene_pseudoT <- read.csv(file="../data/20180720_female_lineage1_pseudotime_DE_genes.csv", row.names=1)
female_lineage2_sig_gene_pseudoT <- read.csv(file="../data/20180720_female_lineage2_pseudotime_DE_genes.csv", row.names=1)


female_lineage1_clustering <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
female_lineage2_clustering <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]


gene_list <- unique(rownames(female_lineage1_clustering), rownames(female_lineage2_clustering))

de_matrix <- log(females[rownames(females) %in% gene_list,]+1)

L1_lineage <- female_pseudotime[!is.na(female_pseudotime[,1]),1]
L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), L1_rpkm_exp_lineage))]


L2_lineage <- female_pseudotime[!is.na(female_pseudotime[,2]),2]
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


cellType_L1 <- female_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)

cellType_L2 <- female_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)


cellLin <- c(
	L2_cellLin,
	L1_cellLin
)


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



data_heatmap <- data.frame(
	L2_cells_smooth,
	L1_cells_smooth
)


cellType <- c(
	cellType_L2,
	cellType_L1
)



scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
    min = apply(x, 1, min, na.rm = T)
    max = apply(x, 1, max, na.rm = T)
    return((x - m) / (max-min))
}


scaled_data_heatmap <-scale_rows(data_heatmap)


gene_clustering <- pheatmap(
	data_heatmap, 
	scale="row", 
	clustering_method="ward.D2",
	silent=TRUE
)


clusters <- cutree(gene_clustering$tree_row, k = 17)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"

gene_cluster_palette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#49beaa', '#611c35', '#2708a0')
gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
names(gene_cluster_colors) <- 1:max(clusters)

annotation_row <- clustering


annotation_col <- data.frame(
	# Cell_Clusters=clusters_in_lineage,
	cellLineages=cellLin,
	cellType=cellType,
	Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)


cellLinCol <- c(
	"#3b3561", 
	"#c8c8c8", 
	"#ff6663"
)
names(cellLinCol) <- unique(cellLin)

cellTypeCol <- c(
	C1="#a53bad", 
	C2="#560047", 
	C3="#eb6bac", 
	C4="#ffa8a0"
)

names(cellTypeCol) <- unique(cellType)


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

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
breaksList = seq(-2.2, 2.5, by = 0.2)


# cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
# warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
# mypalette <- c(rev(cold(20)), warm(21))
# breaksList = seq(-2.2, 2.5, by = 0.2)




tiff(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_heatmap_DE_genes_granulosa_progenitors_k_17_pval_005.tiff"), res = 300, height = 21, width = 18, units = 'cm')

gene_clustering <- pheatmap(
	data_heatmap, 
	scale="row",
	# kmeans_k=15,
	# breaks=c(-0.5, 0, 0.5),
	gaps_col=length(cellType_L2),
	show_colnames=FALSE, 
	show_rownames=FALSE, 
	cluster_cols=FALSE,
	# cluster_rows=FALSE,
	# cutree_rows=6,
	clustering_method="ward.D",
	annotation_row=annotation_row,
	annotation_col=annotation_col,
	annotation_colors=annotation_colors,
	cutree_rows=17, 
	annotation_names_row=FALSE,
	# color=viridis(20),
	color=mypalette
	# breaks=breaksList
)

dev.off()

write.csv(clustering, file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_mean_norm_female_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k17_scaled.csv"))

tf_dynamics <- clustering[rownames(clustering) %in% tf_list$V1,,drop=FALSE]
write.csv(tf_dynamics, file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_mean_norm_female_lineage1_gene_clustering_kmeans_k17_scaled_qval_0.05_TFt.csv"))


# de_genes <- read.csv(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_mean_norm_female_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k10_scaled.csv"))
# gene_names <- de_genes$X

# de_genes <- clustering
# gene_names <- rownames(clustering)

de_genes <- read.csv(file="../graph/female/20180720_female_mean_norm_female_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k17_scaled.csv")
gene_names <- de_genes$Genes


#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

de_gene_clusters <- de_genes[de_genes$Genes %in% entrez_genes$SYMBOL,,drop=FALSE]

de_gene_clusters <- data.frame(
	ENTREZID=entrez_genes[!duplicated(entrez_genes$SYMBOL),"ENTREZID"],
	Gene_Clusters=de_gene_clusters$Gene.categories
)

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

dotplot(formula_res, showCategory=5)

lineage1_ego <- simplify(
	formula_res, 
	cutoff=0.5, 
	by="p.adjust", 
	select_fun=min
)

dotplot(lineage1_ego, showCategory=5) + theme(aspect.ratio=1.1)

write.csv(formula_res@compareClusterResult, file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_compared_GO_term_DE_cluster.csv"))
write.csv(lineage1_ego@compareClusterResult, file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_compared_symplified_GO_term_DE_cluster.csv"))

pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_GO_term_DE_genes_clusters.pdf"), width=11, height=8)
dotplot(formula_res, showCategory=3)+ theme(aspect.ratio=0.8)
dotplot(lineage1_ego, showCategory=3)+ theme(aspect.ratio=2)
dev.off()



###########################################
#										  #
#		Plot Foxl2 vs Lgr5 in GraC  	  #
#										  #
###########################################

pre_graC <- t(log(females[rownames(females) %in% c("Foxl2", "Lgr5"), colnames(females) %in% names(female_clustering[female_clustering=="C3"])]+1))
# graC <- t(log(females[rownames(females) %in% c("Foxl2", "Lgr5"), colnames(females) %in% names(female_clustering[female_clustering=="C4"])]+1))

# pre_graC <- rbind(pre_graC, graC)
pre_graC <- pre_graC[rowSums(pre_graC)>0,]

cor.test(pre_graC[,1],pre_graC[,2])


mean_exp <- rowMeans(pre_graC)
ratio_exp <- log2(pre_graC[,1]/pre_graC[,2])
ratio_exp[ratio_exp=="10"] <- 0
ratio_exp[ratio_exp=="-10"] <- 0
ratio_exp[ratio_exp=="NaN"] <- 0

stage_ordered_cells <- sapply(strsplit(names(ratio_exp[order(ratio_exp)]), "_"), `[`, 1)

pre_graC <- pre_graC[order(ratio_exp),]
# rownames(pre_graC) <- paste(seq_along(rownames(pre_graC)), rownames(pre_graC), sep="_")
rownames(pre_graC) <- seq_along(rownames(pre_graC))


library("reshape")
pre_graC__foxl2_lgr5 <- as.data.frame(melt(pre_graC))
pre_graC__foxl2_lgr5 <- data.frame(
	pre_graC__foxl2_lgr5,
	Stages=stage_ordered_cells
)

colnames(pre_graC__foxl2_lgr5) <- c("Cell", "gene", "Expression", "Stages")


pdf(file=paste0("../graph/female/", format(Sys.time(), "%Y%m%d"), "_female_pre-granulosa_Foxl2_vs_Lgr5.pdf"))
ggplot(pre_graC__foxl2_lgr5, aes(x=Cell, y=Expression, fill=Stages))  + 
	geom_bar(data=subset(pre_graC__foxl2_lgr5 ,gene=="Foxl2"), stat="identity") + 
	geom_bar(data=subset(pre_graC__foxl2_lgr5 ,gene=="Lgr5"),aes(y=Expression*(-1)), stat="identity") +
	geom_hline(yintercept = 0) +
	ggtitle("Foxl2 vs Lgr5 expression") +
	scale_fill_manual(
		values=female_stagePalette[-1],
		name=""
	) +
	theme_bw() +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)
dev.off()
