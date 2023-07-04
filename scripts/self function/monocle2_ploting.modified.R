### plot_pseudotime_trajectory 220519
# changed the heatmap color theme
bks <- seq(-3.1,3.1, by = 0.1)
mycol<-colorRampPalette(colors = rev(c('#c83349','#e06377','#f1e3dd','#cfe0e8','#80ced6','#667292')),alpha=0.75)(length(bks) - 1)
#show_col(mycol)


# add parameter to define the manual annotation color


#The following code is swipped from colorRamps package which is used to make the pallette
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}


#' @importFrom grDevices rgb
rgb.tables <- function(n,
red = c(0.8, 0.25, 0.7),
green = c(0.5, 0.25, 0.7),
blue = c(0.2, 0.25, 0.7))
{
    rr <- do.call("table.ramp", as.list(c(n, red)))
    gr <- do.call("table.ramp", as.list(c(n, green)))
    br <- do.call("table.ramp", as.list(c(n, blue)))
    rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
rgb.tables(n,
red = c(0.8, 0.3, 0.7),
green = c(0.45, 0.3, 0.5),
blue = c(0.2, 0.3, 0.7))

blue2green2red <- matlab.like2


plot_pseudotime_heatmap <- function(cds_subset, 
                                    
                                    cluster_rows = TRUE,
                                    hclust_method = "ward.D2", 
                                    num_clusters = 6,
                                    
                                    hmcols = NULL, 
                                    
                                    add_annotation_row = NULL,
                                    add_annotation_col = NULL,
                                    show_rownames = FALSE, 
                                    use_gene_short_name = TRUE,
                                    
                                    annotation_color=NULL,
                                    
                                    norm_method = c("log", "vstExprs"), 
                                    scale_max=3, 
                                    scale_min=-3, 
                                    
                                    trend_formula = '~sm.ns(Pseudotime, df=3)',
                                    
                                    return_heatmap=FALSE,
                                    cores=1){
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime),length.out = 100)) 
  
  m <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                       relative_expr = T, new_data = newdata)
  

  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds_subset, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min

  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  #if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
   # hmcols <- blue2green2red(length(bks) - 1)
  #}
  #else {
 #   bks <- seq(-3.1,3.1, length.out = length(hmcols))
 # } 
  
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=cluster_rows, 
                 show_rownames=F, 
                 show_colnames=F, 
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 border_color = NA,
                 color=hmcols)

  if(cluster_rows) {
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  } else {
    annotation_row <- NULL
  }
  
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  if(!is.null(add_annotation_col)) {
    if(nrow(add_annotation_col) != 100) {
      stop('add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!')
    }
    annotation_col <- add_annotation_col
  } else {
    annotation_col <- NA
  }
 
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if(!is.null(annotation_row))
      row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  if(!is.null(annotation_row))
    row.names(annotation_row) <- row_ann_labels
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     annotation_color=annotation_color,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=mycol,#hmcols, 
                     border_color = NA,
                     silent=TRUE,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}



plot_genes_branched_heatmap <- function(cds_subset, 
                                        
                                        branch_point=1,
                                        branch_states=NULL,
                                        branch_labels = c("Cell fate 1", "Cell fate 2"), 
                                        cluster_rows = TRUE,
                                        hclust_method = "ward.D2", 
                                        num_clusters = 6,
                                        hmcols = NULL, 
                                        branch_colors = c('#979797', '#F05662', '#7990C8'), 
                                        cluster_colors=NULL,
                                        cluster_num=NULL,
                                        add_annotation_row = NULL,
                                        add_annotation_col = NULL,
                                        #annotation_colors=NULL,
                                        show_rownames = FALSE, 
                                        use_gene_short_name = TRUE,
                                        scale_max=3, 
                                        scale_min=-3, 
                                        norm_method = c("log", "vstExprs"), 
                                        
                                        trend_formula = '~sm.ns(Pseudotime, df=3) * Branch',
                                        
                                        return_heatmap=FALSE,
                                        cores = 1, ...) {
  
  cds <- NA
  new_cds <- buildBranchCellDataSet(cds_subset, 
                                    branch_states=branch_states, 
                                    branch_point=branch_point, 
                                    progenitor_method = 'duplicate',
                                    ...)
  
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  
  if(is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, 'State']
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  
  col_gap_ind <- 101
  # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  
  newdataA <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))   
  newdataB <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,  
                                    relative_expr = T, new_data = rbind(newdataA, newdataB))
  
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  
  #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == setdiff(pData(new_cds)$State, branch_states),])
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
  BranchB_num <- BranchA_num
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs') {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix=BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix=BranchB_exprs)
  }     
  else if(norm_method == 'log') {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  
  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
  
  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min
  
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  exp_rng <- range(heatmap_matrix) #bks is based on the expression range
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- colorRampPalette(colors = rev(c('#c83349','#e06377','#f1e3dd','#cfe0e8','#80ced6','#667292')),alpha=0.75)(length(bks) - 1)#blue2green2red(length(bks) - 1)
  }
  
  # prin  t(hmcols)
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=TRUE, 
                 show_rownames=F, 
                 show_colnames=F, 
                 #scale="row",
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols
                 #color=hmcols#,
                 # filename="expression_pseudotime_pheatmap.pdf",
  )
  #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, branch_labels, BranchA_num, BranchP_num, BranchB_num, file = 'heatmap_matrix')
  
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  
  if(!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(branch_labels[1], BranchA_num),
                                                                                      rep("Pre-branch",  2 * BranchP_num),
                                                                                      rep(branch_labels[2], BranchB_num)))
  
  colnames(annotation_col) <- "Cell Type"  
  
  if(!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
  }
  
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  
  annotation_colors=list("Cell Type"=branch_colors)
  if(!is.null(cluster_colors)){
      annotation_colors=list("Cell Type"=branch_colors,"Cluster"=cluster_colors)
      names(annotation_colors$Cluster)=cluster_num
    }
  names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols=FALSE, 
                     cluster_rows=TRUE, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     annotation_colors=annotation_colors,
                     gaps_col = col_gap_ind,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE)
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, heatmap_matrix = heatmap_matrix, 
                heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols, 
                annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col, 
                ph_res = ph_res))
  }
}




