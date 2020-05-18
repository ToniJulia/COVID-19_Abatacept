library(org.Hs.eg.db)
library(GO.db)
library(metap)
library(corrplot)
library(eulerr)
library(prabclus)

# Find BPs clusters to reduce redundancy =========================
# Load BP GO database
dictGO <- as.list(org.Hs.egGO2ALLEGS)
dictGO <- lapply(dictGO, function(x) {unique(x)})
dictGO <- dictGO[sapply(names(dictGO),function(x) length(dictGO[[x]])>9&length(dictGO[[x]])<301)]
dictTERMS <- as.list(GOTERM)
dictGO <- dictGO[sapply(names(dictGO), function(x){Ontology(dictTERMS[[x]])=='BP'})]

# Jaccard similarity matrix
genes <- unique(unlist(dictGO))
genes_in_bp_mat <- matrix(0,nrow=length(genes),ncol=length(dictGO),dimnames=list(genes,names(dictGO)))
for (path in names(dictGO)){
  path_genes <- dictGO[[path]]
  genes_in_bp_mat[path_genes,path] <- 1
}
Sim <- 1 - jaccard(genes_in_bp_mat)
rm(genes_in_bp_mat)

# Cluster BPs
tree <- hclust(d=as.dist(1-Sim))
cut <- cutree(tree, h=0.60)
table(table(cut))

bps_cluster <- data.frame(GOID = names(dictGO),
                  Description = sapply(names(dictGO), function(x) Term(dictTERMS[[x]])),
                  Cluster = cut[names(dictGO)],
                  Size = sapply(dictGO, function(x) length(x), USE.NAMES=F),
                  stringsAsFactors = FALSE)

# GSEA results antagonism =========================
# Load results
COVID <- readRDS("GSEA/COVID19_GSEA.RDS")
Abatacept <- readRDS("GSEA/RA_Abatacept_GranulocytesAdj_GSEA.RDS")

# Remove redundant BPs (keep the cluster-BP most strongly associated to COVID19 as representative)
COVID <- merge(bps_cluster[,2:4],COVID[,4:11],by=0)
Abatacept <- merge(bps_cluster[,2:4],Abatacept[,4:11],by=0)
names(COVID)[1] <- names(Abatacept)[1] <- "ID"
row.names(COVID) <- COVID$ID
row.names(Abatacept) <- Abatacept$ID
 
COVID <- COVID[order(COVID$pvalue), ]
Abatacept <- Abatacept[order(Abatacept$pvalue), ]
COVID <- COVID[COVID$ID %in% Abatacept$ID, ]
COVID <- COVID[!duplicated(COVID$Cluster), ]
Abatacept <- Abatacept[Abatacept$ID %in% COVID$ID, ]
n <- nrow(COVID)

# Significant and overlapping BPs
COVID <- COVID[COVID$p.adjust < 0.05, ]
Abatacept <- Abatacept[Abatacept$p.adjust < 0.05, ]
overlap <- merge(COVID[,c(1:4,6:8,10,11)],Abatacept[,c(6:8,10,11)],by=0,suffixes=c("_COVID","_Abatacept"))

write.csv(file = "Figures/SupplementaryTable1.csv", row.names = F,
          COVID[,c("ID","Description","Size","NES","pvalue","p.adjust","ratio")])

write.csv(file = "Figures/SupplementaryTable2.csv", row.names = F,
          Abatacept[,c("ID","Description","Size","NES","pvalue","p.adjust","ratio")])

write.csv(file = "Figures/SupplementaryTable3.csv", row.names = F,
          overlap[,c("ID","Description","Size","NES_COVID","pvalue_COVID","p.adjust_COVID","ratio_COVID",
                     "NES_Abatacept","pvalue_Abatacept","p.adjust_Abatacept","ratio_Abatacept")])

# Calculate antagonism
cov_pos_coherent <- sum(overlap$NES_COVID > 0 & overlap$NES_Abatacept < 0)
cov_neg_coherent <- sum(overlap$NES_COVID < 0 & overlap$NES_Abatacept > 0)
cov_pos_incoherent <- sum(overlap$NES_COVID > 0 & overlap$NES_Abatacept > 0)
cov_neg_incoherent <- sum(overlap$NES_COVID < 0 & overlap$NES_Abatacept < 0)

# Calculate antagonism probability
p0 <- (nrow(COVID[COVID$NES > 0, ])/n)*(nrow(Abatacept[Abatacept$NES < 0, ])/n) + (nrow(COVID[COVID$NES < 0, ])/n)*(nrow(Abatacept[Abatacept$NES > 0, ])/n)

p <- binom.test(x = cov_pos_coherent + cov_neg_coherent, n = n , p = p0, alternative = "greater")$p.value

print(paste0("n coherent = ",cov_pos_coherent + cov_neg_coherent, "; n incoherent = ", 
             cov_pos_incoherent + cov_neg_incoherent,"; analyzed BPs = ", n, "; binomial test probability = ", p))

# FIGURE: Euler diagram
cov_pos <- sum(COVID$NES > 0) - cov_pos_coherent - cov_pos_incoherent
cov_neg <- sum(COVID$NES < 0) - cov_neg_coherent - cov_neg_incoherent
abatacept_pos <- sum(Abatacept$NES > 0) - cov_neg_coherent - cov_pos_incoherent
abatacept_neg <- sum(Abatacept$NES < 0) - cov_pos_coherent - cov_neg_incoherent

fit1 <- euler(c("Up regulated in COVID-19" = cov_pos, "Down regulated in COVID-19" = cov_neg,
                "Up regulated by Abatacept" = abatacept_pos, "Down regulated by Abatacept" = abatacept_neg, 
                "Up regulated in COVID-19&Down regulated by Abatacept" = cov_pos_coherent,
                "Up regulated in COVID-19&Up regulated by Abatacept" = cov_pos_incoherent,
                "Down regulated in COVID-19&Up regulated by Abatacept" = cov_neg_coherent,
                "Down regulated in COVID-19&Down regulated by Abatacept" = cov_neg_incoherent)
)

jpeg("Figures/figure_4.jpeg", width=900, height=900)
  plot(fit1, shape="ellipse",
       labels = list(fontsize = 20, font = "plain"),
       edges ="azure4",
       quantities = list(fontsize = 20))
dev.off()



# FIGURE: Most antagonistic BPs
# Determine most different BPs by combining the P-value of both exposures (sum of the logs Fisher's method)
overlap$sum_pval <- apply(overlap[,c("pvalue_COVID","pvalue_Abatacept")], MARGIN = 1, function(x){y<-sumlog(x);y$p})
overlap <- overlap[order(overlap$sum_pval), ]
overlap <- overlap[1:20, ]

pvals <- overlap[,c("pvalue_COVID","pvalue_Abatacept")]
nes <- overlap[,c("NES_COVID","NES_Abatacept")]

row.names(pvals) <- row.names(nes) <- overlap$Description
colnames(nes) <- colnames(pvals) <- c("COVID-19 PBMCs","Abatacept in RA")

pdf("Figures/figure_5.pdf")
  corrplot::corrplot(as.matrix(nes), is.corr=F, p.mat=as.matrix(pvals), 
                     sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 1, pch.cex=0.7, 
                     insig="label_sig", tl.col = "black", cl.lim = c(-4,4), cl.pos = "n")
dev.off()
