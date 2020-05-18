library(corrplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GO.db)
library(gridExtra)
library(fgsea)

# Select curated BPs to analyze ==============================

bps1 <- c("GO:0002224","GO:0039528","GO:0060337","GO:0016197", # Virus sensing
          "GO:0035747","GO:0042267" # NKs
)

bps2 <- c("GO:0002456","GO:0002369","GO:0019882", # T cell and antigen processing and presentation
          "GO:0002573","GO:0042116","GO:0034341","GO:0071356", # Monocytes and macrophages activation and differentiation (& response to IFNG i TNF)
          "GO:0032612","GO:0032635","GO:0032640","GO:0032637",  # Monocyte and macrophage - produced cytokines 
          "GO:0098760","GO:0032613","GO:0006956", # Other associated cytokines + complement activation
          "GO:0019724", # B cell mediated immunity
          "GO:0030193" # Blood coagulation
          
)

# Determine genes in BP
dictGO <- as.list(org.Hs.egGO2ALLEGS)
dictTERMS <- as.list(GOTERM)
dictGO <- lapply(dictGO, function(x) {unique(x)})
dictGO <- dictGO[c(bps1,bps2)]
bps1 <- sapply(bps1, function(x) Term(dictTERMS[[x]]))
bps2 <- sapply(bps2, function(x) Term(dictTERMS[[x]]))
names(dictGO) <- c(bps1,bps2)
rm(dictTERMS)

genes_dict <- readRDS("Source_data/symbol2entrez.RDS")
genes_dict <- genes_dict[!is.na(genes_dict$ENTREZ_ID),]
rownames(genes_dict) <- genes_dict$ENTREZ_ID
dictGO <- sapply(dictGO, function(x) {genes_dict$SYMBOL[genes_dict$ENTREZ_ID %in% x]})

# Datasets 
dt <- gsea <- list(COVID=NULL,Abatacept=NULL)
dt[["COVID"]] <- readRDS("DEG/COVID19_DEG.RDS")
dt[["Abatacept"]] <- readRDS("DEG/RA_Abatacept_DEG.RDS")
dt[["Abatacept"]]$symbols <- row.names(dt[["Abatacept"]])

# GSEA of selected BPs
for (dt_name in c("COVID","Abatacept")){
  vals2test <- -log10(dt[[dt_name]]$PValue) * sign(dt[[dt_name]]$logFC) 
  names(vals2test) <- dt[[dt_name]]$symbols
  vals2test <- vals2test[order(vals2test,decreasing = T)]

  gsea[[dt_name]] <- fgsea(pathways=dictGO,stats=vals2test,nperm=1000000,minSize=1,maxSize=300,nproc=3)
  gsea[[dt_name]] <- gsea[[dt_name]][,c('pathway','NES','pval','padj')]
  names(gsea[[dt_name]]) <- c("Description","NES","pvalue","p.adjust")
  gsea[[dt_name]] <- as.data.frame(gsea[[dt_name]])
  row.names(gsea[[dt_name]]) <- gsea[[dt_name]]$Description
}


# FIGURE: Curated BPs modulation by Abatacept summary ====================================
COVID <- gsea[["COVID"]]
Abatacept <- gsea[["Abatacept"]]
nm1 <- "COVID-19 PBMCs"
nm2 <- "Abatacept in RA"
nm_ref <- "COVID-19 literature"

Descriptions1 <- COVID[bps1,"Description"]
Descriptions2 <- COVID[bps2,"Description"]
nes1_ref <- c(3.5,3.5,3.5,3.5,3.5,-3.5)
nes2_ref <- c(3.5,3.5,3.5,3.5,3.5,
              3.5,3.5,3.5,3.5,3.5,
              3.5,3.5,3.5,3.5,3.5,3.5)
pval1_ref <- rep(0.5,length(bps1))
pval2_ref <- rep(0.5,length(bps2))

summary_nes1 <- cbind(nes1_ref, COVID[bps1,"NES"], Abatacept[bps1,"NES"])
summary_pvals1 <- cbind(pval1_ref, COVID[bps1,"pvalue"], Abatacept[bps1,"pvalue"])
row.names(summary_nes1) <- row.names(summary_pvals1) <- Descriptions1

summary_nes2 <- cbind(nes2_ref, COVID[bps2,"NES"], Abatacept[bps2,"NES"])
summary_pvals2 <- cbind(pval2_ref, COVID[bps2,"pvalue"], Abatacept[bps2,"pvalue"])
row.names(summary_nes2) <- row.names(summary_pvals2) <- Descriptions2

colnames(summary_nes1) <- colnames(summary_pvals1) <- colnames(summary_nes2) <- colnames(summary_pvals2) <- c(nm_ref,nm1,nm2)

pdf("Figures/figure_2A.pdf")
  corrplot(as.matrix(summary_nes1), is.corr=F, p.mat=as.matrix(summary_pvals1), 
           sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 1, pch.cex=0.7, 
           insig="label_sig", tl.col = "black", cl.lim = c(-4,4), cl.pos = "n")
dev.off()

pdf("Figures/figure_2B.pdf")
  corrplot(as.matrix(summary_nes2), is.corr=F, p.mat=as.matrix(summary_pvals2), 
           sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 1, pch.cex=0.7, 
           insig="label_sig", tl.col = "black", cl.lim = c(-4,4), cl.pos = "n")
dev.off()

# Curated BPs modulation by Abatacept volcano plots =================================
# (base plot with a fraction of the points)
data <- dt[["Abatacept"]]
data$SIG <- ifelse(data$logFC < 0, 'NEG', 'POS')
data$PValue <- -log10(data$PValue)

n_parts <- round(nrow(data)/5)
del_smp <- c(sample(50:n_parts,(n_parts/2)),
             sample((n_parts+1):(n_parts*2),(n_parts/1.5)),
             sample(((n_parts*2)+1):(n_parts*3),(n_parts/1.3)),
             sample(((n_parts*3)+1):(n_parts*4),(n_parts/1.1)),
             sample(((n_parts*4)+1):(n_parts*5),(n_parts/1.05)))
base_data <- data[-del_smp, ]

# FIGURE: BP1 volcano plots ---------------
p <- vector(mode="list",length=3)
bps1[2] <- names(dictGO)[2] <- "cytoplasmic pattern recognition receptor \nsignaling pathway in response to virus"

for (line in 0:2){
  bp_line <- bps1[(1:2)+line*2]
  bp_line <- bp_line[!is.na(bp_line)]
  
  base_data_line <- NULL
  bp_data_line <- NULL
  for (bp in bp_line){
    base_data_line <- rbind(base_data_line,base_data)
    genes <- dictGO[[bp]]
    bp_data_line[[bp]] <- data[row.names(data) %in% genes, ]
    bp_data_line[[bp]]$BP <- bp
  }
  bp_data_line <- do.call(rbind,bp_data_line)
  bp_data_line$BP <- factor(bp_data_line$BP, levels=bp_line)
  base_data_line$BP <- factor(rep(bp_line,each=nrow(base_data)),levels=bp_line)
  
  
  p[[line+1]] <- ggplot() +
    geom_point(aes(x = logFC, y = PValue), data = base_data_line, color =  "darkgray", size = 3) +
    geom_point(aes(x = logFC, y = PValue, color = SIG), data = bp_data_line, size = 3, show.legend = F) +
    facet_grid(. ~ BP) +
    scale_color_manual(values=c("POS"="#0000cc","NEG"="#cc0000")) +
    scale_shape_identity() +
    geom_hline(yintercept=-log10(0.05),linetype="twodash",col="seagreen4") +
    xlab("Effect size") + 
    ylab("-log10 P Value") +
    theme(legend.position="none", 
          panel.border=element_rect(colour = "darkgray", fill=NA)) +
    theme_classic(base_size = 29, base_line_size = 1.25, base_rect_size = 1)
}

jpeg(file="Figures/figure_3A.jpeg",width = 1080, height = 1350)
  grid.arrange(p[[1]], p[[2]], p[[3]], ncol = 1, nrow = 3)
dev.off()

# FIGURE: BP2 volcano plots ---------------------
bps2[7] <- names(dictGO)[13] <- "cellular response to TNF"
p <- vector(mode="list",length=4)
for (line in 0:3){
  bp_line <- bps2[(1:4)+line*4]
  bp_line <- bp_line[!is.na(bp_line)]
  
  base_data_line <- NULL
  bp_data_line <- NULL
  for (bp in bp_line){
    base_data_line <- rbind(base_data_line,base_data)
    genes <- dictGO[[bp]]
    bp_data_line[[bp]] <- data[row.names(data) %in% genes, ]
    bp_data_line[[bp]]$BP <- bp
  }
  bp_data_line <- do.call(rbind,bp_data_line)
  bp_data_line$BP <- factor(bp_data_line$BP, levels=bp_line)
  base_data_line$BP <- factor(rep(bp_line,each=nrow(base_data)),levels=bp_line)
  
  
  p[[line+1]] <- ggplot() +
    geom_point(aes(x = logFC, y = PValue), data = base_data_line, color =  "darkgray", size = 3) +
    geom_point(aes(x = logFC, y = PValue, color = SIG), data = bp_data_line, size = 3, show.legend = F) +
    facet_grid(. ~ BP) +
    scale_color_manual(values=c("POS"="#0000cc","NEG"="#cc0000")) +
    scale_shape_identity() +
    geom_hline(yintercept=-log10(0.05),linetype="twodash",col="seagreen4") +
    xlab("Effect size") + 
    ylab("-log10 P Value") +
    theme(legend.position="none", 
          panel.border=element_rect(colour = "darkgray", fill=NA)) +
    theme_classic(base_size = 29, base_line_size = 1.25, base_rect_size = 1)
  
}

jpeg(file="Figures/figure_3B.jpeg",width = 1600, height = 1440)
  grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 1, nrow = 4)
dev.off()
