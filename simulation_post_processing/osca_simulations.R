cd /scratch/local/daniel/daniel/methylSim/simulationsElior/forDaniel/simulations/osca


# create osca direct and download program
mkdir osca
cd osca
wget http://cnsgenomics.com/software/osca/download/osca_Linux.zip
unzip osca_Linux.zip



cd /scratch/local/monthly/daniel/forDaniel/simulations/osca
module add R/3.4.2
R
require(data.table)
files <- c("../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.01/simulated_datasim_data_n_2000_p_0.15_tau_0.01_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.03/simulated_datasim_data_n_2000_p_0.15_tau_0.03_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.05/simulated_datasim_data_n_2000_p_0.15_tau_0.05_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07/simulated_datasim_data_n_2000_p_0.15_tau_0.07_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.001/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.001_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.1/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.1_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_10/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_10_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.09/simulated_datasim_data_n_2000_p_0.15_tau_0.09_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.1_tau_0.07/simulated_datasim_data_n_2000_p_0.1_tau_0.07_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.3_tau_0.07/simulated_datasim_data_n_2000_p_0.3_tau_0.07_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.5_tau_0.07/simulated_datasim_data_n_2000_p_0.5_tau_0.07_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.7_tau_0.07/simulated_datasim_data_n_2000_p_0.7_tau_0.07_O",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.9_tau_0.07/simulated_datasim_data_n_2000_p_0.9_tau_0.07_O")
pheno <- c("../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.01/simulated_datasim_data_n_2000_p_0.15_tau_0.01_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.03/simulated_datasim_data_n_2000_p_0.15_tau_0.03_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.05/simulated_datasim_data_n_2000_p_0.15_tau_0.05_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07/simulated_datasim_data_n_2000_p_0.15_tau_0.07_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.001/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.001_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.1/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_0.1_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_10/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_10_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.09/simulated_datasim_data_n_2000_p_0.15_tau_0.09_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.1_tau_0.07/simulated_datasim_data_n_2000_p_0.1_tau_0.07_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.3_tau_0.07/simulated_datasim_data_n_2000_p_0.3_tau_0.07_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.5_tau_0.07/simulated_datasim_data_n_2000_p_0.5_tau_0.07_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.7_tau_0.07/simulated_datasim_data_n_2000_p_0.7_tau_0.07_y",
                "../simulated_data/simulated_datasim_data_n_2000_p_0.9_tau_0.07/simulated_datasim_data_n_2000_p_0.9_tau_0.07_y")




files <- c("../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_2/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_3/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_4/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_5/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_6/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_7/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_8/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_9/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_10/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_11/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_12/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_13/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_14/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
          "../simulated_data_15/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_O")
pheno <- c("../simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_O",
                    "../simulated_data_2/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_3/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_4/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_5/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_6/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_7/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_8/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_9/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_10/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_11/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_12/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_13/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_14/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y",
                    "../simulated_data_15/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_y")

for(i in 1:15){
    dat <- fread(paste(files[i],".dat",sep=''))
    dd <- data.frame("IID" = 1:ncol(dat), t(dat))
    names(dd)[2:ncol(dd)] <- paste("cg",1:(ncol(dd)-1),sep='')
    write.table(dd, paste(files[i],".txt",sep=''), row.names=FALSE, col.names=TRUE, quote=FALSE)
    system(paste("./OSCA_Linux --efile ",files[i],".txt --no-fid --methylation-beta --make-bod --out ",files[i], sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --make-orm --out ",files[i],sep=''))
    phen <- read.table(paste(pheno[i],".dat",sep=''))
    phen$FID <- 1:ncol(dat)
    phen$IID <- 1:ncol(dat)
    phen <- phen[,c(2,3,1)]
    write.table(phen, paste(pheno[i],".phen",sep=''), row.names=FALSE, col.names=FALSE, quote=FALSE)
    system(paste("./OSCA_Linux --orm ",files[i]," --reml --reml-pred-rand --pheno ",pheno[i],".phen  --out ",files[i],sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --blup-probe ",files[i],".indi.blp --out ",files[i],sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --mlma --pheno ",pheno[i],".phen --orm ",files[i]," --out ",files[i],sep=''))
}





cd /scratch/local/monthly/daniel/forDaniel/simulations/osca
module add R/3.4.2; R
require(data.table)
files <- c("../../simulations_no_sparse/simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_2/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_3/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_4/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_5/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_6/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_7/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_8/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_9/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_10/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_11/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_12/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_13/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_14/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O",
           "../../simulations_no_sparse/simulated_data_15/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_O")
pheno <- c("../../simulations_no_sparse/simulated_data/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_datasim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_2/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_2sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_3/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_3sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_4/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_4sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_5/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_5sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_6/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_6sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_7/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_7sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_8/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_8sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_9/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_9sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_10/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_10sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_11/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_11sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_12/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_12sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_13/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_13sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_14/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_14sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y",
          "../../simulations_no_sparse/simulated_data_15/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01/simulated_data_15sim_data_n_2000_p_0.15_tau_0.07_sp_1000_sigma_0.01_y")

for(i in 1:15){
    dat <- fread(paste(files[i],".dat",sep=''))
    dd <- data.frame("IID" = 1:nrow(dat), dat)
    names(dd)[2:ncol(dd)] <- paste("cg",1:(ncol(dd)-1),sep='')
    write.table(dd, paste(files[i],".txt",sep=''), row.names=FALSE, col.names=TRUE, quote=FALSE)
    system(paste("./OSCA_Linux --efile ",files[i],".txt --no-fid --methylation-beta --make-bod --out ",files[i], sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --make-orm --out ",files[i],sep=''))
    phen <- read.table(paste(pheno[i],".dat",sep=''))
    phen$FID <- 1:nrow(dat)
    phen$IID <- 1:nrow(dat)
    phen <- phen[,c(3,4,2)]
    write.table(phen, paste(pheno[i],".phen",sep=''), row.names=FALSE, col.names=FALSE, quote=FALSE)
    system(paste("./OSCA_Linux --orm ",files[i]," --reml --reml-pred-rand --pheno ",pheno[i],".phen  --out ",files[i],sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --blup-probe ",files[i],".indi.blp --out ",files[i],sep=''))
    system(paste("./OSCA_Linux --befile ",files[i]," --mlma --pheno ",pheno[i],".phen --orm ",files[i]," --out ",files[i],sep=''))
}
