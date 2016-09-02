get_ld =function(rsids,res_preparation){
  idx=match(rsids,res_preparation$data_set$rsid)
  ld_tmp=res_preparation$ld_matrix[idx,idx]
  return(ld_tmp)
}