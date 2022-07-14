#### Functions to extract sensitivity and specificity & overall AUC values for each fold from a blockCV::spatialBlock used inside a SDMtune::SDMmodelCV 

if (!require('pROC')) install.packages('pROC'); library('pROC')

extract_spec_sens_vals <- function(SDM_model_CV, Spatial_Blocks, SWD_data){

  list = list()
  
  k = length(SDM_model_CV@models)

  for(i in 1:k){
    # extract test data from each fold
    test <- Spatial_Blocks$folds[[i]][[2]]
    # extract P/A data and covariates for test points  
    test_df <- cbind(SWD_data@pa[test],SWD_data@data[test,])
    # extract each model
    sdm_model <- SDM_model_CV@models[[i]]
    # run prediction for each model
    prediction <- predict(sdm_model, data = test_df)
    # bind together the actual and predicted values
    df = as.data.frame(cbind(test_df[,1], prediction))
    list[[i]] <- df
    next
    
  }
  
  # run the pROC::roc() on each model
  roc_output <- lapply(1:length(list), function(x){roc(list[[x]]$V1, list[[x]]$prediction)})
  
  list2 <- list()
  
  for(i in 1:k){
  # extract the specificity and sensitivity values
  spec_sens <- as.data.frame(cbind(roc_output[[i]]$specificities, roc_output[[i]]$sensitivities))
  # rename columns
  names(spec_sens) <- c('specificities', 'sensitivities')
  # attach the model number
  spec_sens$model_no <- i
  list2[[i]] <- spec_sens
  next
  }
  
  # bind the rows together
  plot_df <- bind_rows(list2)
  return(plot_df)
}

extract_auc_vals <- function(SDM_model_CV, Spatial_Blocks, SWD_data){
  
  list = list()
  
  k = length(SDM_model_CV@models)
  
  for(i in 1:k){
    # extract test data from each fold
    test <- Spatial_Blocks$folds[[i]][[2]]
    # extract P/A data and covariates for test points  
    test_df <- cbind(SWD_data@pa[test],SWD_data@data[test,])
    # extract each model
    sdm_model <- SDM_model_CV@models[[i]]
    prediction <- predict(sdm_model, data = test_df)
    
    df = as.data.frame(cbind(test_df[,1], prediction))
    list[[i]] <- df
    next
    
  }
  
  roc_output <- lapply(1:length(list), function(x){roc(list[[x]]$V1, list[[x]]$prediction)})
  
  list2 <- list()
  
  for(i in 1:k){
    auc_vals <- as.data.frame(cbind(i, roc_output[[i]]$auc[1]))
    names(auc_vals) <- c('model_no', 'auc')
    list2[[i]] <- auc_vals
    next
  }
  return(bind_rows(list2))
}
