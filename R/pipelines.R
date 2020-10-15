# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de
#
# Collection of wrapper functions to build automatic pipelines
# Suggest a Seurat wrapper for Visium data


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
MISTy_seurat = function(visium_slide,
                        view_assays,
                        view_features,
                        view_types,
                        view_params,
                        spot_ids = NULL,
                        out_alias = "default",
                        workers = 4){
  
  # Geometry extraction
  geometry = visium_slide@images$slice1@coordinates[,c(2,3)]
  
  # Extracting data
  view_data = map(view_assays, 
                  extract_seurat_data, 
                  geometry = geometry,
                  visium_slide = visium_slide)
  
  # Building pipeline
  
  build_MISTy_pipeline(view_data = view_data,
                       view_features = view_features,
                       view_types = view_types,
                       view_params = view_params,
                       geometry = geometry,
                       spot_ids = spot_ids,
                       out_alias = out_alias,
                       workers = workers)
  
}

#' Builds complete MISTy pipelines
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
build_MISTy_pipeline = function(view_data,
                                view_features,
                                view_types,
                                view_params,
                                geometry,
                                spot_ids = NULL,
                                out_alias = "default",
                                workers = 4){
  
  plan(multiprocess, workers = workers)
  
  clear_cache()
  
  #Checking areas to consider
  if(is.null(spot_ids)){
    spot_ids = rownames(view_data[[1]])
  }
  
  #First filter the features from the data
  view_data = map2(view_data,view_features,filter_data_features)
  
  #Create initial view
  #Defining useful data intra
  views_main = create_initial_view(view_data[[1]][spot_ids,], 
                                   unique.id = "main")
  
  #Create other views
  all_views = pmap(list(view_data[-1],
                        view_types[-1],
                        view_params[-1],
                        names(view_data[-1])), 
                   create_all_views, 
                   spot_ids = spot_ids,
                   geometry = geometry)
  
  
  # Define frankenstein view
  for(i in 1:length(all_views)){
    views_main = views_main %>%
      add_views(all_views[[i]])
  }
  
  # Run MISTy
  MISTy_run = run_misty(views_main,out_alias)
}

################

# Helper functions

#'Extracts data from an specific assay in a Seurat object
#'and aligns the IDs to be identical to geometry
extract_seurat_data = function(visium_slide,assay,geometry){
  data_df = as.matrix(visium_slide@assays[[assay]]@data)
  data_df = data_df %>%
    t %>% data.frame(check.names = F)
  data_df = data_df[rownames(geometry),]
  return(data_df)
}

#' Filters data to contain features of interested
filter_data_features = function(data, features){
  if(is.null(features)){
    features = colnames(data)
  }
  data_df = data[,features]
  colnames(data_df) = gsub("-","_", colnames(data_df))
  return(data_df)
}

#' Builds views depending on the paramaters defined
#' @param data = data.frame containing IDs in rows and features in columns
#' @param view_type = type of view to be generated
#' @param view_param = l or neighbor.thr parameter
#' @param spot_ids = area to be considered in the analysis
#' @param geometry = location data frame containing IDs in rows and coordinates in columns 
create_all_views = function(data, 
                            view_type,
                            view_param,
                            view_name,
                            spot_ids,
                            geometry){
  
  if(view_type == "intra"){
    
    view_data_tmp = create_initial_view(data, unique.id = view_name)
    
    # Spot specific view comes from the view above
    data_red = view_data_tmp[["intracellular"]]$data
    data_red = data_red[spot_ids,]
    
  }else if(view_type == "para"){
    
    view_data_tmp = create_initial_view(data, unique.id = view_name) %>% 
      add_paraview(geometry,l = view_param^2)
    
    # Spot specific view comes from the view above
    data_red = view_data_tmp[[3]]$data
    rownames(data_red) = rownames(view_data_tmp[["intracellular"]]$data) #we named rows just for easy access
    data_red = data_red[spot_ids,]
    
  }else if(view_type == "juxta"){
    view_data_tmp = create_initial_view(data, unique.id = view_name) %>% 
      add_juxtaview(positions = geometry,neighbor.thr = view_param)
    
    # Spot specific view comes from the view above
    data_red = view_data_tmp[[3]]$data
    rownames(data_red) = rownames(view_data_tmp[["intracellular"]]$data) #we named rows just for easy access
    data_red = data_red[spot_ids,]
    
  }
  
  return(create_view(paste0(view_name,"_",view_param),data_red))
  
}




























