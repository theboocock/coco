### Conditional analysis using the matrix notation from GCTA-cojo paper.
### @author James Boocock
dos_data = setDF(fread("data/NEK1_Cauc_452_dos.txt", header=T))
# Beta estimates from a single SNP model.
# eq 4. in the GCTA-cojo paper
beta_hat = solve(diag(diag(XX))) %*% t(cbind(scale(dos_data$rs10520157,scale = F),scale(dos_data$rs4417927,scale = F))) %*% scale(dos_data$ENSG00000137601, scale=F)
# X'X in the GCTA-cojo paper. 
XX  = t(cbind(scale(dos_data$rs10520157,scale=F),scale(dos_data$rs4417927,scale = F))) %*% cbind(scale(dos_data$rs10520157,scale=F) , scale(dos_data$rs4417927,scale=F))
# eq 2. in the GCTA-cojo paper. For the estimate
beta_joint =solve(XX) %*% t(cbind(scale(dos_data$rs10520157,scale=F),scale(dos_data$rs4417927,scale = F))) %*% scale(dos_data$ENSG00000137601, scale=F)

# eq 2. in the GCTA-cojo paper. For the variance.
# calculated in TWO steps.
# 1. eq. 7. Calculate the joint variance from the new betas and the old betas.
# 2. eq. 12. Using the variance from eq 7. and the inverse of X'X calculate the standard error for each estimate.
# 
sigma_sq = (t(dos_data$ENSG00000137601)%*%(dos_data$ENSG00000137601) - t(beta_joint) %*% diag(diag(XX)) %*% beta_hat) / (452 - 3)
sigma_sq = sqrt(diag((sigma_sq[1]) *solve(XX)))
true_zs = beta_joint/sqrt(diag((sigma_sq[1]) *solve(XX)))
