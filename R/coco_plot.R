library(ggplot2)
library(dplyr)
plot_stepwise = function(res_preparation, stepwise_results,one_per_step=T,z_plot=F){
  plot_list = c()
  num_plots = 1
  if(class(res_preparation) != "coco_data"){
    stop("res_preparation object must be of type coco_data")
  }
  if(class(stepwise_results) != "coco_stepwise"){
    stop("stepwise results object must be of type coco_stepwise")
  }
  if(is.null(stepwise_results$step_betas)){
    message("No stepwise betas found, generating them now.")
    cond_on = c()
    stepwise_results_list = list()
    for(snp in stepwise_results$stepwise_summary$rsid){
      cond_on = c(snp,cond_on)
      stepwise_results_list = c(stepwise_results_list, list(joint_from_ids(cond_on,res_preparation)))
    }
    stepwise_results$step_betas = stepwise_results_list
  }
  # Generate a plot per dataset.
  for(i in 1:nrow(stepwise_results$stepwise_summary)){
    # Must be the first locus
    if(i == 1){
      top_idx = which(max(abs(res_preparation$data_set$z)) == abs(res_preparation$data_set$z))[1]
      if(z_plot){
        p_logs= (abs(res_preparation$data_set$z))
      }else{
        p_logs= -log10(2*pnorm(abs(res_preparation$data_set$z), lower.tail=F))
        if(sum(is.infinite(p_logs))> 0){
          stop("Some P-values infinite cannot plot, convert to ZScore plotting with z_plot=T")
        }
      }
      max_p = max(p_logs,na.rm=T)
      fact_top = max(abs(res_preparation$data_set$z)) == abs(res_preparation$data_set$z)
    }else{
      top_idx = which(max(abs(stepwise_results$step_betas[[i-1]]$res_step$z_new), na.rm=T) == abs(stepwise_results$step_betas[[i-1]]$res_step$z_new))[1]
      if(z_plot){
        p_logs= abs(stepwise_results$step_betas[[i-1]]$res_step$z_new)
      }else{
        p_logs = -log10(2*pnorm(abs(stepwise_results$step_betas[[i-1]]$res_step$z_new), lower.tail=F))
      }
      fact_top = which(max(abs(stepwise_results$step_betas[[i-1]]$res_step$z_new),na.rm=T) == abs(stepwise_results$step_betas[[i-1]]$res_step$z_new))
    }
    ld_row = res_preparation$ld_matrix[top_idx,]^2
    ld_cut = cut(ld_row,breaks=c(0,.2,.4,.6,.8,1.0))
    plotting_data = data.frame(p_logs=p_logs,pos=res_preparation$data_set$pos,ld_cut=ld_cut, highlight_top=fact_top)
    plot_list[[num_plots]] = plotting_data %>% ggplot(aes(y=p_logs,x=pos,colour=ld_cut)) + geom_point() + theme_bw() + ggtitle(paste(res_preparation$data_set$rsid[top_idx]," MAX")) +
      ylab("-Log10 P-value") + xlab("Genomic Position") + scale_colour_discrete(name="Pairwise LD") +
        geom_point(data=plotting_data[top_idx,],aes(y=p_logs,x=pos), shape=18,size=4,colour="purple")
    num_plots = num_plots + 1
    if(!one_per_step){
      for(j in (i+1):nrow(stepwise_results$stepwise_summary)){
        if((i+1) >  nrow(stepwise_results$stepwise_summary)){
          next
        }
        snp_idx = which(res_preparation$data_set$rsid == stepwise_results$stepwise_summary$rsid[j])
        ld_row = res_preparation$ld_matrix[snp_idx,]^2 
        ld_cut = cut(ld_row,breaks=c(0,.2,.4,.6,.8,1.0))
        plotting_data$ld_cut = ld_cut
        plot_list[[num_plots]] = plotting_data %>% ggplot(aes(y=p_logs,x=pos,colour=ld_cut)) + geom_point()  + theme_bw() + ggtitle(paste(res_preparation$data_set$rsid[snp_idx], "Conditional")) +
               ylab("-Log10 P-value") + xlab("Genomic Position") + scale_colour_discrete(name="Pair-wise LD") + ylim(0,max_p + 1) + 
               geom_point(data=plotting_data[snp_idx,],aes(y=p_logs,x=pos), shape=18,size=4,colour="purple")
        num_plots = num_plots + 1
      }
      
    }
  }
  return(plot_list)
}



plot_all_but_one = function(res_preparation, all_but_one){
  if(class(res_preparation) != "coco_data"){
    stop("res_preparation object must be of type coco_data")
  }
  if(class(all_but_one) != "coco_all_but_one"){
    stop("stepwise results object must be of type coco_stepwise")
  }
  # Generate a plot per dataset.
  for(i in 1:length(all_but_one)){
    # Must be the first locus
    
    top_idx = which(max(abs(all_but_one[[i]]$res_step$z_new), na.rm=T) == abs(all_but_one[[i]]$res_step$z_new))[1]
    p_logs = -log10(2*pnorm(abs(all_but_one[[i]]$res_step$z_new), lower.tail=F))
    if(i == 1){
      max_p = max(p_logs,na.rm=T)
    }
    fact_top = which(max(abs(all_but_one[[i]]$res_step$z_new),na.rm=T) == abs(all_but_one[[i]]$res_step$z_new))
    print(fact_top)
    print(res_preparation$data_set$rsid[top_idx])
    ld_row = res_preparation$ld_matrix[top_idx,]^2
    ld_cut = cut(ld_row,breaks=c(0,.2,.4,.6,.8,1.0))
    plotting_data = data.frame(p_logs=p_logs,pos=res_preparation$data_set$pos,ld_cut=ld_cut, highlight_top=fact_top)
    print(plotting_data %>% ggplot(aes(y=p_logs,x=pos,colour=ld_cut)) + geom_point()  + theme_bw() + ggtitle(res_preparation$data_set$rsid[top_idx]) +
            ylab("-Log10 P-value") + xlab("Genomic Position") + scale_colour_discrete(name="Pair-wise LD") + ylim(0,max_p + 1))
  }
}