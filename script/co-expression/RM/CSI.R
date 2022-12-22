#' Calculate CSI module activity over all cell types
#'
#' @param clusters_df -
#' @param regulonAUC -
#' @param metadata -
#' @keywords SCENIC, regulons, CSI activity
#' @import SCENIC
#' @import tidyverse
#' @import pheatmap
#' @import viridis
#' @export
#' @examples
#'

calc_csi_module_activity <- function(clusters_df,regulonAUC,metadata,cell_type_column){
  
  metadata$cell_type <- metadata[ , cell_type_column ]
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_sub <- regulonAUC#@assays@data@listData$AUC
  regulonAUC_sub <- regulonAUC_sub[regulons,]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_cluster" = c(),
                                     "mean_activity" = c(),
                                     "cell_type" = c())

  cell_type_counter <- 0
  regulon_counter <- 0
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1

      cell_type_aucs <- rowMeans(regulonAUC_sub[,rownames(subset(metadata,cell_type == ct))])
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }

  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_cluster)){
      csi_regulon <- subset(clusters_df,csi_cluster == cluster)

      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_cluster" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)

  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_cluster
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}


#' Calculates CSI values for all regulon pairs
#'
#' @param regulonAUC The AUC values for all regulons as calculated by SCENIC (content of file:3.4_regulonAUC.Rds).
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @import SCENIC
#' @import tidyverse
#' @import svMisc
#' @export
#' @examples
#' \donttest{
#' regulon_thresholds <- auc_thresh_kmeans(regulonAUC)
#' }

calculate_csi <- function(regulonAUC,
                          calc_extended = FALSE,
                          verbose = FALSE){

  compare_pcc <- function(vector_of_pcc,pcc,CONST=0.05){
        pcc_larger <- length(vector_of_pcc[vector_of_pcc > (pcc - CONST)])

        return(pcc_larger)
    
  }

  calc_csi <- function(reg,reg2,pearson_cor){
    test_cor <- pearson_cor[reg,reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)

    sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
    fraction_lower <- length(sums[sums == 0]) / total_n
    return(fraction_lower)
  }

  regulonAUC_sub <- regulonAUC#@assays@data@listData$AUC

  if(calc_extended == TRUE){
    regulonAUC_sub <- subset(regulonAUC_sub,grepl("extended",rownames(regulonAUC_sub)))
  } else if (calc_extended == FALSE){
    regulonAUC_sub <- subset(regulonAUC_sub,!grepl("extended",rownames(regulonAUC_sub)))
}

  regulonAUC_sub <- t(regulonAUC_sub)

  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    gather(regulon_2,pcc,-regulon_1) %>%
    mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))


  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)

  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))

  colnames(csi_regulons) <- c("regulon_1",
                              "regulon_2",
                              "CSI")

  num_regulons <- length(regulon_names)

  f <- 0
  for(reg in regulon_names){
    ## Check if user wants to print info
    if(verbose == TRUE){
      print(reg)
      }
    for(reg2 in regulon_names){
      f <- f + 1

      fraction_lower <- calc_csi(reg,reg2,pearson_cor)

      csi_regulons[f,] <- c(reg,reg2,fraction_lower)

    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}