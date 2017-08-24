#' Relabels the vector of labels to make sure it represents a valid partition
#' but it respects zero entries.
#' @return a cluster assignment vector that represents a valid partition
#' Code by Maria Lomeli

relabels <-function(c_block){
  if(sum(c_block==0)>0){
  minimum_posit <- 0 
  }
else{
   minimum_posit <-1
  }
  Max_Val <-max(c_block)
  aux_idx <-which(c_block>0)
  Min_Val <-min(c_block[aux_idx])
  c_block[aux_idx] <- c_block[aux_idx]-Min_Val+1
  cnt <- tabulate(c_block[aux_idx])
  idx_zeros <- which(cnt==0)
  
  minimum_posit <- 0 
  for (ell in 1:length(idx_zeros)){
    idx_help <-which(c_block>=idx_zeros[ell])
    c_block[idx_help] <- c_block[idx_help]-1
  }
  Max_Val <-max(c_block)
  cnt <- tabulate(c_block[aux_idx])
  idx_zeros <- which(cnt==0)
  while(length(idx_zeros)>0){
    for (ell in 1:length(idx_zeros)){
      idx_help <-which(c_block>=idx_zeros[ell])
      c_block[idx_help] <- c_block[idx_help]-1
    }
    Max_Val <-max(c_block)
    cnt <- tabulate(c_block[aux_idx])
    idx_zeros <- which(cnt==0)
  }
  return(c_block)
}



