### Conditional analysis using the matrix notation from GCTA-cojo paper.
### @author James Boocock
dos_data = setDF(fread("data/NEK1_Cauc_452_dos.txt", header=T))
beta_hat = solve(diag(diag(XX))) %*% t(cbind(scale(dos_data$rs10520157,scale = F),scale(dos_data$rs4417927,scale = F))) %*% scale(dos_data$ENSG00000137601, scale=F)
XX  = t(cbind(scale(dos_data$rs10520157,scale=F),scale(dos_data$rs4417927,scale = F))) %*% cbind(scale(dos_data$rs10520157,scale=F) , scale(dos_data$rs4417927,scale=F))
beta_joint =solve(XX) %*% t(cbind(scale(dos_data$rs10520157,scale=F),scale(dos_data$rs4417927,scale = F))) %*% scale(dos_data$ENSG00000137601, scale=F)
solve(XX) %*%diag(diag(XX)) %*% beta_hat

(sqrt(diag(XX))) %*% cor(cbind(scale(dos_data$rs10520157,scale = F),scale(dos_data$rs4417927,scale = F))) %*% sqrt(diag(XX))
cov_m = cov(cbind(scale(dos_data$rs10520157,scale = F),scale(dos_data$rs4417927,scale = F)))

sigma_sq = (t(dos_data$ENSG00000137601)%*%(dos_data$ENSG00000137601) - t(beta_joint) %*% diag(diag(XX)) %*% beta_hat) / (452 - 3)
sqrt(diag((sigma_sq[1]) *solve(XX)))
diag((sigma_sq[1]) *solve(XX))

true_zs = beta_joint/sqrt(diag((sigma_sq[1]) *solve(XX)))
