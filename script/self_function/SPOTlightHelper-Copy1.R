SpatialScatterPie <- function(
  object,
  cell_types_all,
  cell.types = NULL,
  assay = "Spatial",
  images = NULL,
  image.alpha = 1,
  pie.alpha = 1,
  pie.scale = 1,
  cols = SelectColors(cell_types_all,palette="col50" )
) {
  
  ## If slice is not selected set, all slices will be used
  if ( is.null(images) ){
    images <- Seurat::Images(object, assay = assay)[1]
  }
  if (length(images) < 1) {
    stop("Could not find any spatial image information")
  }
  
  metadata_ds <- data.frame(object@meta.data)
  # metadata_ds <- Seurat::FetchData(object, vars = cell_types_all)
  
  # Change column names for consistency ##
  # [[:punct:]] - Any punctuation character: ! ' # S % & ' ( ) * + , - . / : ; < = > ? @ [ / ] ^ _ { | } ~
  colnames(metadata_ds) <- gsub("[[:punct:]]|[[:blank:]]", ".", colnames(metadata_ds), perl = TRUE)
  cell_types_all <- gsub("[[:punct:]]|[[:blank:]]", ".", cell_types_all, perl = TRUE)
  
  if (is.null(cell.types)) {
    cell.types <- cell_types_all
  } else {
    cell.types <- gsub( "[[:punct:]]|[[:blank:]]",  ".", cell.types, perl = TRUE)
  }
  
  # If not all cell types are in the cell types of interest we only want to keep those spots which have at least one of the cell types of interest
  slice <- images[1]
  metadata_ds <- metadata_ds %>% dplyr::filter( slice == slice )
  if (!all(cell_types_all %in% cell.types)) {
    metadata_ds <- metadata_ds %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::mutate(rsum = rowSums(.[, cell.types, drop = FALSE])) %>%
      dplyr::filter(rsum != 0) %>%
      dplyr::select("ID") %>%
      dplyr::left_join(metadata_ds %>% rownames_to_column("ID"),by = "ID") %>%
      tibble::column_to_rownames("ID")
  }
  
  # plots <- vector( mode = "list", length = length(images)  )
  ## Preprocess data
  scalefacor <- Seurat::ScaleFactors(object[[slice]])
  spatial_coord <- data.frame(object[[slice]]@coordinates) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::mutate(imagerow_scaled =  imagerow*scalefacor$lowres,#imagerow* scalefacor$lowres
                  imagecol_scaled =  imagecol*scalefacor$lowres) %>% #imagecol* scalefacor$lowres
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("ID"), by = "ID")
  
  ## Convert image to grob object
  #img <- Seurat::GetImage(object[[slice]], mode = "raw")
  #img_grob <- grid::rasterGrob(img, interpolate = FALSE,
  #                             width = grid::unit(1, "npc"),
  #                             height = grid::unit(1, "npc"))
  
  ## Plot spatial scatterpie plot
  scatterpie_plt <- suppressMessages(
    ggplot2::ggplot() +
      #ggplot2::annotation_custom(
      #  grob = img_grob,
      #  xmin = 0,
      #  xmax = ncol(img_grob$raster),
      #  ymin = 0,
      #  ymax = -nrow(img_grob$raster)
      #) +
      scatterpie::geom_scatterpie(data = spatial_coord,
                                  aes(x = imagecol_scaled,
                                      y = imagerow_scaled),
                                  cols = cell_types_all,
                                  color = NA,
                                  alpha = pie.alpha,
                                  pie_scale = pie.scale) +
      ggplot2::labs(title = ifelse(all(cell_types_all %in% cell.types),slice,paste(slice,cell.types,sep="_")), fill = "Celltype") +
      ggplot2::guides(fill = guide_legend(label.theme = element_text(size = 6), 
                                          override.aes = list(size = 1),
                                          ncol = 1)) +
      ggplot2::scale_y_reverse() +
      ggplot2::ylim(600, 0) +
      ggplot2::xlim(0, 600) +
      cowplot::theme_half_open(11, rel_small = 1) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed(ratio = 1,
                           xlim = NULL,
                           ylim = NULL,
                           expand = TRUE,
                           clip = "on") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), 
                     legend.key.size = unit(10, "pt"))  
  )
  plot <- scatterpie_plt + scale_fill_manual( values = cols)
  
  return(plot)
}

SelectColors <- function(
  object = NULL,
  palette = "col50",
  value = "celltype",
  n = NULL) {
  if (!is.null(object)) {
    if (class(object) == "data.frame") {
      colid <- ifelse(is.null(value), colnames(object)[1], value)
      names <- unique(object[[colid]])
    }
    if (is.factor(object)) {
      names <- levels(object)
    }
    if (is.vector(object)) {
      names <- unique(gsub("[[:punct:]]|[[:blank:]]", ".", object, perl = TRUE))
    }
    n <- length(names)
  } else if (!is.null(n)) {
    names <- NULL
  }
  
  colors2pick <- switch(palette,
                        seurat = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100),
                        ## ref: http://stackoverflow
                        # .com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
                        blindless = c(
                          "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
                          "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                          "#66A61E", "#E6AB02", "#A6761D", "#BC80BD", "#1F78B4",
                          "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                          "#8DD3C7", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#CAB2D6",
                          "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                          "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
                          "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#A6CEE3",
                          "#984EA3", "#FFFF33", "#A65628", "#F781BF", "#999999", "#FFED6F",
                          "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                          "#E5C494", "#B3B3B3", "#FFFFB3", "#BEBADA", "#FB8072",
                          "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#666666"
                        ),
                        col50 = c(
                          "#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                          "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                          "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                          "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                          "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                          "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                          "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                          "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                          "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                          "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"
                        ),
                        ditto = c(
                          "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                          "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
                          "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
                          "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                          "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"
                        ),
                        paired = brewer.pal(n = n, "Paired"),
                        colx22 = c(
                          "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4",
                          "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff",
                          "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
                          "#000075", "#808080", "#4f34ff", "#f340F0"
                        ),
                        jet = c(
                          "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                          "#FF7F00", "red", "#7F0000"
                        ),
                        tableau20 = c(
                          "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                          "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                          "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                          "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"
                        ),
                        tableau10medium = c(
                          "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                          "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                          "#CDCC5D", "#6DCCDA"
                        ),
                        colorblind10 = c(
                          "#006BA4", "#FF800E", "#ABABAB", "#595959",
                          "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                          "#FFBC79", "#CFCFCF"
                        ),
                        trafficlight = c(
                          "#B10318", "#DBA13A", "#309343", "#D82526",
                          "#FFC156", "#69B764", "#F26C64", "#FFDD71",
                          "#9FCD99"
                        ),
                        purplegray12 = c(
                          "#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
                          "#5F5A41", "#B4B19B", "#995688", "#D898BA",
                          "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"
                        ),
                        bluered12 = c(
                          "#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
                          "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
                          "#BD0A36", "#F4737A"
                        ),
                        greenorange12 = c(
                          "#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                          "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                          "#39737C", "#86B4A9", "#82853B", "#CCC94D"
                        ),
                        cyclic = c(
                          "#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
                          "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
                          "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
                          "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"
                        )
  )
  
  if (is.null(n)) {
    colors_use <- colors2pick
  } else {
    colors_use <- colors2pick[1:n]
  }
  if (!is.null(names)) {
    names(colors_use) <- names
  }
  return(colors_use)
}