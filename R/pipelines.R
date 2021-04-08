# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Defines general misty pipelines for 
#' their usage with Seurat objects 10x

#' Builds automatic MISTy pipelines for 10x Visium data in Seurat format
#' 
#' Warning: Feature IDs can't have "-" symbol
#' Warning: rownames must be the same in intra_df, para_df and geometry
#' Warning: assumes first element of each list to be the main view
#'
#' @param visium_slide: Seurat object containing spatial transcriptomics data
#' @param view_assays: named list containing the assays used to build the views
#' @param view_features: named list containing the features to be used, if NULL uses all features.
#' @param view_types: named list containing the type of view to be built: "intra", "juxta" or "para"
#' @param view_params: named list containing the params used for each view. Currently only accepts a numeric value representing 
#' the l parameter of the paraview or the number of neighbours used in the juxtaview
#' @param spot_ids: spot IDs to fit MISTy if null all spots are used
#' @param out_alias: folder name to be used in all MISTy outputs
run_misty_seurat <- function(visium.slide,
                         view.assays,
                         view.features,
                         view.types,
                         view.params,
                         spot.ids = NULL,
                         out.alias = "default"){
  
  # Geometry extraction --------------------------------  
  geometry <- visium.slide@images$slice1@coordinates[,c(2,3)]
  
  # Extracting data -------------------------------- 
  view.data <- map(view.assays, 
                   extract_seurat_data, 
                   geometry = geometry,
                   visium.slide = visium.slide)
  
  # Building pipeline --------------------------------
  build_misty_pipeline(view.data = view.data,
                       view.features = view.features,
                       view.types = view.types,
                       view.params = view.params,
                       geometry = geometry,
                       spot.ids = spot.ids,
                       out.alias = out.alias)
}

#' Builds automatic MISTy pipelines 
#' 
#' Warning: Feature IDs can't have "-" symbol
#' Warning: rownames must be the same in intra_df, para_df and geometry
#'
#' @param view_data: named list containing the data used to build the views
#' @param view_features: named list containing the features to be used in each view, if NULL uses all features.
#' @param view_types: named list containing the type of view to be built: "intra", "juxta" or "para"
#' @param view_params: named list containing the params used for each view. Currently only accepts a numeric value representing 
#' the l parameter of the paraview or the number of neighbours used in the juxtaview
#' @param geometry: a data frame with IDs as rows two coordinates as columns 
#' @param spot_ids: spot IDs to fit MISTy if null all rows are used
#' @param out_alias: folder name to be used in all MISTy outputs
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default"){
  
  # Adding all spots ids in case they are not defined ----------
  if(is.null(spot.ids) == TRUE){
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data ------------
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view ------------------------------------------
  views.main <- create_initial_view(view.data.filt[[1]][spot.ids, ], 
                                   unique.id = "main")
  
  # Create other views -------------------------------------------
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(view.data.filt[-1],
                         view.types[-1],
                         view.params[-1],
                         view.names[-1]), 
                    create_any_view, 
                    spot.ids = spot.ids,
                    geometry = geometry)
  
  pline.views <- add_views(views.main, 
                           unlist(all.views,recursive = F))
    
  # Run MISTy ----------------------------------
  run_misty(pline.views,out.alias)
}

################

# Helper functions

#'Extracts data from an specific assay in a Seurat object
#'and aligns the IDs to be identical to geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  data_df <- as.matrix(visium.slide@assays[[assay]]@data)
  data_df <- data_df %>%
    t %>% data.frame(check.names = F)
  data_df <- data_df[rownames(geometry), ]
  return(data_df)
}

#' Filters data to contain features of interest
filter_data_features <- function(data, 
                                 features) {
  if (is.null(features) == TRUE) {
    features <- colnames(data)
  }
  data_df <- data[, features]
  colnames(data_df) <- gsub("-", "_", colnames(data_df))
  return(data_df)
}

#' Builds views depending on the paramaters defined
#' @param data <- data.frame containing IDs in rows and features in columns
#' @param view_type <- type of view to be generated
#' @param view_param <- l or neighbor.thr parameter
#' @param spot_ids <- area to be considered in the analysis
#' @param geometry <- location data frame containing IDs in rows and coordinates in columns 
create_any_view <- function(data, 
                            view.type,
                            view.param,
                            view.name,
                            spot.ids,
                            geometry) {
  
  MISTy::clear_cache()
  
  if (view.type == "intra") {
    
    view.data.tmp <- create_initial_view(data, unique.id = view.name)
    
    # Spot specific view comes from the view above
    data.red <- view.data.tmp[["intraview"]]$data
    data.red <- data.red[spot.ids, ]
    
  } else if (view.type == "para") {
    
    view.data.tmp <- create_initial_view(data, unique.id = view.name) %>% 
      add_paraview(geometry,l = view.param)
    
    # Spot specific view comes from the view above -----------------------
    data.ix <- paste0("paraview.", view.param^2)
    data.red <- view.data.tmp[[data.ix]]$data
    rownames(data.red) <- rownames(view.data.tmp[["intraview"]]$data)
    data.red <- data.red[spot.ids, ]
    
  }else if(view.type == "juxta"){
    view.data.tmp <- create_initial_view(data, unique.id = view.name) %>% 
      add_juxtaview(positions = geometry,
                    neighbor.thr = view.param)
    
    # Spot specific view comes from the view above ----------------------
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data
    rownames(data.red) <- rownames(view.data.tmp[["intraview"]]$data)
    data.red <- data.red[spot.ids, ]
    
  }
  
  if (is.null(view.param) == TRUE) {
    MISTy.view <- create_view(paste0(view.name), 
                             data.red)
  } else {
    MISTy.view <- create_view(paste0(view.name,"_",view.param), 
                             data.red)
  }
  
  return(MISTy.view)
}






