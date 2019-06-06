
#### WARNING ####
#
#  This script is my working script, thus it is not yet cleaned and contains lots of rubish code and plot tests.
#  If you need help, contact me at isabelle.stevant@gmail.com 



gene_list <- rownames(all)

DE_supporting_lineage_a_b <- function(
	rpkm_matrix=rpkm_matrix, 
	count_matrix=count_matrix, 
	pseudotime=pseudotime, 
	clusters=clusters
	){

	a <- 240
	b1 <- 290
	b2 <- 380
	c <- 410

	early_prog_cells <- c(
		rownames(subset(pseudotime, pseudotime[,1]<=a)),
		rownames(subset(pseudotime, pseudotime[,2]<=a))
	)


	bip_cells1 <- c(
		rownames(subset(pseudotime, pseudotime[,1]>=b1)),
		rownames(subset(pseudotime, pseudotime[,2]>=b1))
	)
	bip_cells2 <- c(
		rownames(subset(pseudotime, pseudotime[,1]<=b2)),
		rownames(subset(pseudotime, pseudotime[,2]<=b2))
	)
	bip_cells <- bip_cells1[bip_cells1 %in% bip_cells2]

	sup_cells <- c(
		rownames(subset(pseudotime, pseudotime[,1]>=c)),
		rownames(subset(pseudotime, pseudotime[,2]>=c))
	)

	phase <- data.frame(
		cells=colnames(all),
		phase=rep_along("0", colnames(all))
	)

	levels(phase$phase) <- c(levels(phase$phase), c("a", "b", "c"))

	phase[phase$cells %in% early_prog_cells, "phase"] <- "a"
	phase[phase$cells %in% bip_cells, "phase"] <- "b"
	phase[phase$cells %in% sup_cells, "phase"] <- "c"


	print("Normalise data prior to DE analysis...")
	prep_de <- prepare_for_DE(
		count_matrix=count_matrix, 
		clustering=phase$phase, 
		stages=clusters
	)

	print("Select early progenitor cells...")
	#Select early prog cells
	early_prog <- prep_de[,colnames(prep_de) %in% c(early_prog_cells, bip_cells, sup_cells)]


	print(dim(early_prog))

	# # Re-detect how many cells express each genes in the subset of cells
	early_prog <- detectGenes(early_prog, min_expr = 5)

	# # Remove genes expressed in less than 10 cells
	early_prog <- early_prog[fData(early_prog)$num_cells_expressed >= 10, ]


	print("Perform DE analysis between a and b...")
	early_prog_DE_genes <- differentialGeneTest(
		early_prog, 
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	early_prog_sig_genes_0.05 <- subset(early_prog_DE_genes, qval < 0.05)
	early_prog_sig_genes_0.01 <- subset(early_prog_DE_genes, qval < 0.01)

	print(paste(nrow(early_prog_sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(early_prog_sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	print("Attribute sex to DE genes...")

	early_prog_cells_sex_0.05 <- phase$phase
	names(early_prog_cells_sex_0.05) <- early_prog_cells
	early_prog_count <- count_matrix[,colnames(count_matrix) %in% early_prog_cells]

	early_prog_DE <- get_up_reg_clusters(
		count=early_prog_count, 
		clustering=early_prog_cells_sex_0.05, 
		DE_genes=early_prog_sig_genes_0.05
	)

	print("Done!")



	return(early_prog_DE)

}


DE_a_vs_b <- DE_supporting_lineage_a_b(
	rpkm_matrix=all,
	count_matrix=all_count,
	pseudotime=pseudotime,
	clusters=all_clustering
	)



write.csv(DE_a_vs_b, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_DE_supporting_lineage_a_vs_b.csv"))




DE_supporting_lineage_bip_sup_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	count_matrix=count_matrix, 
	pseudotime=pseudotime, 
	dyn_genes=dyn_genes,
	clusters=clusters
	){

	early_prog_lim <- 290
	bip_sup_lim <- 380


	sup_cells <- c(
		rownames(subset(pseudotime, pseudotime[,1]>=early_prog_lim)),
		rownames(subset(pseudotime, pseudotime[,2]>=early_prog_lim))
	)

	sup_pseudotime <- pseudotime[rownames(pseudotime) %in% sup_cells,]

	bip_sup_cells <- c(
		rownames(subset(sup_pseudotime, sup_pseudotime[,1]<=bip_sup_lim)),
		rownames(subset(sup_pseudotime, sup_pseudotime[,2]<=bip_sup_lim))
	)


	print("Normalise data prior to DE analysis...")
	prep_de <- prepare_for_DE(
		count_matrix=count_matrix, 
		clustering=sapply(strsplit(colnames(count_matrix), "_"), `[`, 2), 
		stages=clusters
	)

	print("Select early progenitor cells...")
	#Select early prog cells
	bip_sup <- prep_de[,colnames(prep_de) %in% bip_sup_cells]

	# Select dynamic genes
	bip_sup <- bip_sup[rownames(bip_sup) %in% dyn_genes,]

	print(dim(bip_sup))

	# Re-detect how many cells express each genes in the subset of cells
	bip_sup<- detectGenes(bip_sup, min_expr = 5)

	# Remove genes expressed in less than 10 cells
	bip_sup <- bip_sup[fData(bip_sup)$num_cells_expressed >= 10, ]


	print("Perform DE analysis between XY and XX...")
	bip_sup_DE_genes <- differentialGeneTest(
		bip_sup, 
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	bip_sup_sig_genes_0.05 <- subset(bip_sup_DE_genes, qval < 0.05)
	bip_sup_sig_genes_0.01 <- subset(bip_sup_DE_genes, qval < 0.01)

	print(paste(nrow(bip_sup_sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(bip_sup_sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))


	print("Attribute sex to DE genes...")

	bip_sup_cells_sex_0.05 <- sapply(strsplit(bip_sup_cells, "_"), `[`, 2)
	names(bip_sup_cells_sex_0.05) <- bip_sup_cells
	bip_sup_count <- count_matrix[,colnames(count_matrix) %in% bip_sup_cells]

	bip_sup_DE <- get_up_reg_clusters(
		count=bip_sup_count, 
		clustering=bip_sup_cells_sex_0.05, 
		DE_genes=bip_sup_sig_genes_0.05
	)

	print("Done!")



	return(bip_sup_DE)
}




bip_sup_sex_DE <- DE_supporting_lineage_bip_sup_sex(
	rpkm_matrix=all,
	count_matrix=all_count,
	pseudotime=pseudotime,
	dyn_genes=gene_list,
	clusters=all_clustering
	)

write.csv(bip_sup_sex_DE, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_DE_supporting_lineage_bip_sup_sex.csv"))



DE_supporting_lineage_sup_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	count_matrix=count_matrix, 
	pseudotime=pseudotime, 
	dyn_genes=dyn_genes,
	clusters=clusters
	){

	bip_sup_lim <- 410


	sup_cells <- c(
		rownames(subset(pseudotime, pseudotime[,1]>=bip_sup_lim)),
		rownames(subset(pseudotime, pseudotime[,2]>=bip_sup_lim))
	)


	print("Normalise data prior to DE analysis...")
	prep_de <- prepare_for_DE(
		count_matrix=count_matrix, 
		clustering=sapply(strsplit(colnames(count_matrix), "_"), `[`, 2), 
		stages=clusters
	)

	print("Select early progenitor cells...")
	#Select early prog cells
	sup <- prep_de[,colnames(prep_de) %in% sup_cells]

	# Select dynamic genes
	sup <- sup[rownames(sup) %in% dyn_genes,]

	print(dim(sup))

	# Re-detect how many cells express each genes in the subset of cells
	sup<- detectGenes(sup, min_expr = 5)

	# Remove genes expressed in less than 10 cells
	sup <- sup[fData(sup)$num_cells_expressed >= 10, ]


	print("Perform DE analysis between XY and XX...")
	sup_DE_genes <- differentialGeneTest(
		sup, 
		fullModelFormulaStr="~cellType",
		cores = 3
	)


	sup_sig_genes_0.05 <- subset(sup_DE_genes, qval < 0.05)
	sup_sig_genes_0.01 <- subset(sup_DE_genes, qval < 0.01)

	print(paste(nrow(sup_sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(sup_sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	print("Attribute sex to DE genes...")

	sup_cells_sex_0.05 <- sapply(strsplit(sup_cells, "_"), `[`, 2)
	names(sup_cells_sex_0.05) <- sup_cells
	sup_count <- count_matrix[,colnames(count_matrix) %in% sup_cells]

	sup_DE <- get_up_reg_clusters(
		count=sup_count, 
		clustering=sup_cells_sex_0.05, 
		DE_genes=sup_sig_genes_0.05
	)

	print("Done!")



	return(sup_DE)
}




sup_sex_DE <- DE_supporting_lineage_sup_sex(
	rpkm_matrix=all,
	count_matrix=all_count,
	pseudotime=pseudotime,
	dyn_genes=gene_list,
	clusters=all_clustering
	)

write.csv(sup_sex_DE, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_DE_supporting_lineage_sup_sex.csv"))


#######
# Join DE result to gene clustering
#######

library(dplyr)

genes <- as.data.frame(rownames(all))
colnames(genes) <- "genes"


sup_genes <- read.csv(file="../graph/all/supplementary_data_4_supporting_lineage_k-mean_25_p-val_0.05.csv")
colnames(sup_genes) <- c("genes", "clusters")

early_prog <- early_prog_sex_DE[,c("genes", "cluster","qval")]
bip_sup <- bip_sup_sex_DE[,c("genes", "cluster","qval")]
sup <- sup_sex_DE[,c("genes", "cluster","qval")]

join0 <- left_join(genes,sup_genes)
colnames(join0) <- c("genes", "clusters")
join1 <- left_join(join0,early_prog)
colnames(join1) <- c("genes", "clusters", "DE in Early prog", "qvalue1")
join2 <- left_join(join1,bip_sup)
colnames(join2) <- c("genes", "clusters", "DE in Early prog", "qvalue1", "DE in bip. sup.", "qvalue2")
join3 <- left_join(join2,sup)
colnames(join3) <- c("genes", "clusters", "DE in Early prog", "qvalue1", "DE in bip. sup.", "qvalue2", "DE in sup.", "qvalue3")


write.csv(join3, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_test_all_genes.csv"))

########################################################################################

gene_list <- unique(c(prog_male_sig_gene_pseudoT, prog_female_sig_gene_pseudoT))


DE_prog_lineage_early_prog_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	count_matrix=count_matrix, 
	pseudotime=pseudotime, 
	dyn_genes=dyn_genes,
	clusters=clusters
	){

	early_prog_lim <- 320


	early_prog_cells <- c(
		rownames(subset(pseudotime, pseudotime[,3]<=early_prog_lim))
	)

	print("Normalise data prior to DE analysis...")
	prep_de <- prepare_for_DE(
		count_matrix=count_matrix, 
		clustering=sapply(strsplit(colnames(count_matrix), "_"), `[`, 2), 
		stages=clusters
	)

	print("Select early progenitor cells...")
	#Select early prog cells
	early_prog <- prep_de[,colnames(prep_de) %in% early_prog_cells]

	# Select dynamic genes
	early_prog <- early_prog[rownames(early_prog) %in% dyn_genes,]

	print(dim(early_prog))

	# Re-detect how many cells express each genes in the subset of cells
	early_prog <- detectGenes(early_prog, min_expr = 5)

	# Remove genes expressed in less than 10 cells
	early_prog <- early_prog[fData(early_prog)$num_cells_expressed >= 10, ]


	print("Perform DE analysis between XY and XX...")
	early_prog_DE_genes <- differentialGeneTest(
		early_prog, 
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	early_prog_sig_genes_0.05 <- subset(early_prog_DE_genes, qval < 0.05)
	early_prog_sig_genes_0.01 <- subset(early_prog_DE_genes, qval < 0.01)

	print(paste(nrow(early_prog_sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(early_prog_sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	print("Attribute sex to DE genes...")

	early_prog_cells_sex_0.05 <- sapply(strsplit(early_prog_cells, "_"), `[`, 2)
	names(early_prog_cells_sex_0.05) <- early_prog_cells
	early_prog_count <- count_matrix[,colnames(count_matrix) %in% early_prog_cells]

	early_prog_DE <- get_up_reg_clusters(
		count=early_prog_count, 
		clustering=early_prog_cells_sex_0.05, 
		DE_genes=early_prog_sig_genes_0.05
	)

	print("Done!")



	return(early_prog_DE)
}




early_prog_sex_DE <- DE_prog_lineage_early_prog_sex(
	rpkm_matrix=all,
	count_matrix=all_count,
	pseudotime=pseudotime,
	dyn_genes=gene_list,
	clusters=all_clustering
	)

write.csv(early_prog_sex_DE, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_DE_prog_lineage_early_prog_sex.csv"))




DE_prog_lineage_late_prog_sex <- function(
	rpkm_matrix=rpkm_matrix, 
	count_matrix=count_matrix, 
	pseudotime=pseudotime, 
	dyn_genes=dyn_genes,
	clusters=clusters
	){

	early_prog_lim <- 410

	late_prog_cells <- c(
		rownames(subset(pseudotime, pseudotime[,3]<=early_prog_lim))
	)


	print("Normalise data prior to DE analysis...")
	prep_de <- prepare_for_DE(
		count_matrix=count_matrix, 
		clustering=sapply(strsplit(colnames(count_matrix), "_"), `[`, 2), 
		stages=clusters
	)

	print("Select early progenitor cells...")
	#Select early prog cells
	late_prog <- prep_de[,colnames(prep_de) %in% late_prog_cells]

	# Select dynamic genes
	late_prog <- late_prog[rownames(late_prog) %in% dyn_genes,]

	print(dim(late_prog))

	# Re-detect how many cells express each genes in the subset of cells
	late_prog<- detectGenes(late_prog, min_expr = 5)

	# Remove genes expressed in less than 10 cells
	late_prog <- late_prog[fData(late_prog)$num_cells_expressed >= 10, ]


	print("Perform DE analysis between XY and XX...")
	late_prog_DE_genes <- differentialGeneTest(
		late_prog, 
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	late_prog_sig_genes_0.05 <- subset(late_prog_DE_genes, qval < 0.05)
	late_prog_sig_genes_0.01 <- subset(late_prog_DE_genes, qval < 0.01)

	print(paste(nrow(late_prog_sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(late_prog_sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	print("Attribute sex to DE genes...")

	late_prog_cells_sex_0.05 <- sapply(strsplit(late_prog_cells, "_"), `[`, 2)
	names(late_prog_cells_sex_0.05) <- late_prog_cells
	late_prog_count <- count_matrix[,colnames(count_matrix) %in% late_prog_cells]

	late_prog_DE <- get_up_reg_clusters(
		count=late_prog_count, 
		clustering=late_prog_cells_sex_0.05, 
		DE_genes=late_prog_sig_genes_0.05
	)

	print("Done!")



	return(late_prog_DE)
}




late_prog_sex_DE <- DE_prog_lineage_late_prog_sex(
	rpkm_matrix=all,
	count_matrix=all_count,
	pseudotime=pseudotime,
	dyn_genes=gene_list,
	clusters=all_clustering
	)

write.csv(late_prog_sex_DE, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_DE_prog_lineage_late_prog_sex.csv"))


#######
# Join DE result to gene clustering
#######

library(dplyr)


prog_genes <- read.csv(file="/media/windows/Users/zazooo/ownCloud/PhD/scRNAseq/Papers/Articles/Ovary_Dev/SupplementaryData/supplementary_data_5_progenitor_lineage_k-mean_25.csv")
prog_genes <- prog_genes[,c("genes", "clusters")]

early_prog <- early_prog_sex_DE[,c("genes", "cluster","qval")]
late_prog <- late_prog_sex_DE[,c("genes", "cluster","qval")]


join1 <- left_join(prog_genes,early_prog)
colnames(join1) <- c("genes", "clusters", "DE in Early prog", "qvalue1")
join2 <- left_join(join1,late_prog)
colnames(join2) <- c("genes", "clusters", "DE in Early prog", "qvalue1", "DE in Late prog.", "qvalue2")


write.csv(join2, file=paste0("../graph/all/", format(Sys.time(), "%Y%m%d"), "_test2.csv"))


########################################################################################



