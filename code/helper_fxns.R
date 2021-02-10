## Helper functions

replace_parms <- function(new_vals, parms){
  val_names <- names(new_vals)
  for(i in 1:length(new_vals)){
    parms[val_names[i]] <- new_vals[i]  
  }
  return(parms)
}