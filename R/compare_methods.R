compare_logistic_results = function(freq_af, ld_matrix,dos_data){
  res.out= data.frame()
  for(i in seq(0.6,1,by=0.05)){
    correlations= c ()
    betas = c()
    ses = c()
    dos_data$RANDOM_BIN = dos_data$ENSG00000137601 > 6
    dos_data$RANDOM_BIN = ifelse(rbinom(n=nrow(dos_data),1,i), dos_data$RANDOM_BIN, 1- dos_data$RANDOM_BIN)
    m1 = summary(glm(dos_data$RANDOM_BIN ~ dos_data$rs10520157, family=binomial ))
    p_value = coef(m1)[2,4]
    var_y = var(dos_data$RANDOM_BIN)
    for(i in rsids){
      m1 = summary(glm(dos_data$RANDOM_BIN ~ dos_data$rs10520157 + dos_data[,i], family=binomial))
      if (colnames(dos_data)[i] != "rs10520157"){
        betas = c(betas,coef(m1)[3,1])
        ses = c(ses,coef(m1)[3,2])
      }else{
        betas = c(betas,NA)
        ses = c(ses,NA)
      }
    }
    beta_unadj =c()
    se_unadj = c()
    for(i in rsids){
      m1 = summary(glm(dos_data$RANDOM_BIN ~ dos_data[,i], family=binomial))
      beta_unadj = c(beta_unadj,coef(m1)[2,1])
      se_unadj = c(se_unadj,coef(m1)[2,2])
    }
    gg = table(dos_data$RANDOM_BIN)[1]/sum(table(dos_data$RANDOM_BIN))
    freq_af$beta = beta_unadj
    freq_af$se = se_unadj
    freq_af$Z =  beta_unadj/se_unadj
    res_preparation= prep_dataset_common(data_set = freq_af,ld_matrix= ld_matrix,hwe_variance = F) 
    g = cbind(df$res_step$beta_new, betas)
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    
    res_preparation = prep_dataset_common(data_set = freq_af,ld_matrix= ld_matrix,hwe_variance = T)
    #stepwise_results2 = stepwise_conditional_run(res_preparation = res_preparation,p_value_threshold = 9.186111e-06,colinear_threshold = 0.9)
   
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    
    # freq_af = merge(corr,frequencies,by=1)

    new_betas = logistic_to_linear(b = freq_af$beta, freq_af$se,phi = gg, ref_af=freq_af$FREQ1)
    freq_af$beta =  new_betas[,1]
    freq_af$se = new_betas[,2]
    freq_af$Z = freq_af$beta/freq_af$se
    
    res_preparation = prep_dataset_common(data_set = freq_af,ld_matrix= ld_matrix,hwe_variance = F)
    df = conditional_from_ids(rsids = "rs10520157", res_preparation)
    g = cbind(df$res_step$beta_new, betas)
    
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    g3 = linear_to_logistic(b = df$res_step$beta_new,df$res_step$se_new, phi=gg, theta=freq_af$FREQ1)
    g = cbind(betas/ses, g3[,1]/g3[,2])
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    
    res_preparation = prep_dataset_common(data_set = freq_af,ld_matrix= ld_matrix,hwe_variance = T)
    df = conditional_from_ids(rsids = "rs10520157", res_preparation)
    g = cbind(df$res_step$beta_new, betas)
    
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    g3 = linear_to_logistic(b = df$res_step$beta_new,df$res_step$se_new, phi=gg, theta=freq_af$FREQ1)
    g = cbind(betas/ses, g3[,1]/g3[,2])
    g = g[!is.na(g[,1]),]
    correlations = c(correlations, cor(g[,1],g[,2], method="spearman"))
    res.out = (rbind(res.out,c(var_y,p_value, correlations)))
  }
  colnames(res.out) = c("Var_bin", "rs10520157","no_step","no_step_hwe","step_no_convert","step_no_convert_hwe","step_log_to_lin","step_log_to_lin_hwe") 
  return(res.out) 
  
}