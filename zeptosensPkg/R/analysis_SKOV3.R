setwd("/Users/hwang29/Google Drive/targetscore_analysis")
x_2 <- read.csv( "data_TS/TS_anil/SKOV3.csv",row.names = 1)  #nProt=304

mab_to_genes <- read.table(system.file("targetscoreData", "antibodyMapFile.txt", package = "zeptosensPkg"),
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE
)

tcga_data <- read.csv("data_TS/TS_Anil/TCGA-OV-L4_1.csv",row.names = 1)

network <- zeptosensPkg::predict_bio_network(
    n_prot = dim(x_2)[2],
    proteomic_responses = x_2,
    max_dist = 1,
    mab_to_genes = mab_to_genes
)
prior <- network$wk  #sum(wk!=0)=834

#Index For proteins in Data/TCGA_Data
#index<-colnames(tcga_data)[which(colnames(tcga_data)%in%colnames(x_2))]
#tcga_data<-tcga_data[,index]
#write.csv(tcga_data,"data_TS/TS_Anil/TCGA-OV-L4_1.csv")
tcga_data<-read.csv("data_TS/TS_Anil/TCGA-OV-L4_1.csv",row.names = 1)
tcga_data[is.na(tcga_data)]<-0

network.Hyb<-zeptosensPkg::predict_hybrid_network(
    n_prot = 304 ,
    proteomic_responses = x_2,
    prior = prior,
    data =tcga_data 
)
write.csv(network.Hyb$wk,file="data_TS/TS_Anil/Hyb_network_ov.csv")
write.csv(network.Hyb$bic,file="data_TS/TS_Anil/Hyb_network_ov_bic.csv")

hybnet<-read.csv(file="data_TS/TS_Anil/Hyb_network.csv",row.names = 1)
hybnet_sif<-zeptosensPkg::create_sif_from_matrix(t_net = hybnet,col_genelist = colnames(hybnet),row_genelist = rownames(hybnet))
write.table(hybnet_sif,file="data_TS/TS_Anil/hybnet_sif.txt",quote = F,sep = "\t",row.names = F)

# Get functional score
fs_override_org <- readRDS(system.file("test_data_files", "fs_value_file.rds",
                                       package = "zeptosensPkg"
))

fs_value <- zeptosensPkg::get_fs_vals(
    n_prot = ncol(x_2), proteomic_responses = x_2,
    mab_to_genes = mab_to_genes, fs_override = fs_override_org
)

# SKOV3
proteomic_responses_2 <- x_2
write.table(proteomic_responses_2,file=paste0(Sample1,"_Zresp.txt"))
length(proteomic_responses_2)

ts <- array(0, dim = c(dim(proteomic_responses_2)[1], dim(proteomic_responses_2)[2]))
ts_p <- array(0, dim = c(dim(proteomic_responses_2)[1], dim(proteomic_responses_2)[2]))
ts_q <- array(0, dim = c(dim(proteomic_responses_2)[1], dim(proteomic_responses_2)[2]))

for (i in seq_len(dim(proteomic_responses_2)[1])) {
    results <- zeptosensPkg::get_target_score(
        wk = network.Hyb$wk,
        wks = network.Hyb$wks,
        dist_ind = network.Hyb$dist_ind,
        inter = network.Hyb$inter,
        n_dose = 1,
        n_prot = dim(proteomic_responses_2)[2],
        proteomic_responses = proteomic_responses_2[i, ],
        n_perm = 1000,
        verbose = FALSE,
        fs_dat = fs_value
    )
    ts[i, ] <- results$ts
    ts_p[i, ] <- results$pts
    ts_q[i, ] <- results$q
}

colnames(ts) <- colnames(proteomic_responses_2)
ts <- data.frame(rownames(proteomic_responses_2), ts)
write.csv(ts,file="data_TS/TS_anil/prglasso_Result_20191113/TS_SKOV3_20200519.csv",row.names = F)

colnames(ts_p) <- colnames(proteomic_responses_2)
ts_p <- data.frame(rownames(proteomic_responses_2), ts_p)
write.csv(ts_p,file="data_TS/TS_anil/prglasso_Result_20191113/TS_SKOV3_20200519_pvalue.csv",row.names = F)


colnames(ts_q) <- colnames(proteomic_responses_2)
ts_q <- data.frame(rownames(proteomic_responses_2), ts_q)
write.csv(ts_q,file="data_TS/TS_anil/prglasso_Result_20191113/TS_SKOV3_20200519_qvalue.csv",row.names = F)

#Plots
library(pheatmap)

TS_SKOV3=read.csv("data_TS/TS_anil/prglasso_Result_20191113/TS_SKOV3_20200519.csv",row.names = 1)
data <- as.matrix(TS_SKOV3)
data<-ifelse(data>=1.5,1.5,ifelse(data<=-1.5,-1.5,data))
bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,1.5,by=0.01))
p=pheatmap(data,
           scale = "none",
           color = c(colorRampPalette(colors = c("navy","white"))(length(seq(-1.5,-0.1,by=0.01))),colorRampPalette(colors = c("white","firebrick3"))(length(seq(0,1.5,by=0.01)))),
           legend_breaks=seq(-1.5,1.5,2),cellwidth = 1, cellheight =1, fontsize=1, fontsize_row=1,
           breaks=bk)
ggsave("data_TS/TS_anil/prglasso_Result_20191113/heatmap_SKOV3_20200519.pdf",p)
