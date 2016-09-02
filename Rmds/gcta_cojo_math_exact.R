### Conditional analysis using the matrix notation from GCTA-cojo paper.
### First is conditional.
### @author James Boocock
summary(lm(scale(dos_data$ENSG00000137601, scale = F) ~dos_data$rs10520157+dos_data$rs4417927))

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
sigma_sq = (t(scale(dos_data$ENSG00000137601,scale = F))%*%(scale(dos_data$ENSG00000137601,scale = F)) - t(beta_joint) %*% diag(diag(XX)) %*% beta_hat) / (452 - 3)
sigma_sq = sqrt(diag((sigma_sq[1]) *solve(XX)))
true_zs = beta_joint/sigma_sq
true_zs

### 2. Conditional analysis of summary statistics.
### Conditional analysis is equivalent to doing a residual extraction.
summary(lm(residuals(lm(scale(dos_data$ENSG00000137601, scale = F) ~scale(dos_data$rs10520157 + dos_data$rs4565044,scale = F) )) ~ -1 + scale(dos_data$rs4417927,scale = F)))
### eq 4.
summary(lm(scale(dos_data$ENSG00000137601, scale = F) ~scale(dos_data$rs10520157 + dos_data$rs4565044,scale = F) ))
#beta_hat = solve(diag(diag(XX))) %*% t(cbind(scale(dos_data$rs10520157 + dos_data$rs4565044,scale = F),scale(dos_data$rs4417927,scale = F))) %*% scale(dos_data$ENSG00000137601, scale=F)

#### Matrices needed to estimate betas.
X1X= t(scale(dos_data$rs10520157, scale=F))%*%scale(dos_data$rs10520157,scale=F) 
X2X = t(scale(dos_data$rs4417927, scale=F)) %*% scale(dos_data$rs4417927, scale=F)
X2X1 = t(scale(dos_data$rs4417927, scale=F)) %*% scale(dos_data$rs10520157,scale=F) 
X1X2 = t(scale(dos_data$rs10520157, scale=F)) %*% scale(dos_data$rs4417927,scale=F) 

## eq.17.
beta_cond = solve(X2X)%*%  t(scale(dos_data$rs4417927, scale=F)) %*% dos_data$ENSG00000137601 - solve(X2X)%*%X2X1%*%solve(X1X) %*%  t(scale(dos_data$rs10520157,scale=F)) %*% dos_data$ENSG00000137601
beta_cond


## eq sigma_C^2

sigma_c = (((t(scale(dos_data$ENSG00000137601,scale=F))%*%(scale(dos_data$ENSG00000137601, scale=F)) - t(beta_hat[1]) %*% XX[1,1] %*% beta_hat[1])- (t(beta_cond[1]) %*% XX[2,2] %*% beta_hat[2])))/ 451

sigma_c
sigma_c2 =  sqrt(sigma_c * (XX[2,2] - X2X1%*%solve(X1X)%*%X1X2)/(XX[2,2]^2))
beta_cond/sigma_c2

0.007894774*5.82503

sigma_c = (((t(scale(dos_data$ENSG00000137601,scale=F))%*%(scale(dos_data$ENSG00000137601, scale=F)) - t(beta_hat[1]) %*% XX[1,1] %*% beta_hat[1])) /450)- (XX[2,2] %*% beta_cond^2) / (450)
sigma_c = (t(scale(dos_data$ENSG00000137601,scale=F))%*%(scale(dos_data$ENSG00000137601, scale=F)) - t(beta_hat[1]) %*% XX[1,1] %*% beta_hat[1] - (t(beta_cond) %*% XX[2,2] %*% beta_cond)) / 450


beta_cond/sqrt(sigma_c %*% solve(X2X))
coef(summary(lm(scale(residuals(lm(scale(dos_data$ENSG00000137601, scale = F) ~ dos_data$rs10520157)),scale=F) ~ dos_data$rs4417927)))

sqrt(sigma_c %*% solve(XX[2,2]) - sigma_c %*% solve(XX[2,2]) %*% X2X1 %*% solve(XX[1,1]) %*% t(X1X2) %*% solve(XX[2,2]))


dat_new = merge(dat$bed, dos_data,by.x="row.names",by.y="ID")


# Try again
t(residuals((lm(scale(dos_data$ENSG00000137601, scale = F) ~dos_data$rs10520157)))) %*% t(residuals((lm(scale(dos_data$ENSG00000137601, scale = F) ~dos_data$rs10520157))))


var(residuals(lm(residuals(lm(scale(dos_data$ENSG00000137601, scale = F) ~ dos_data$rs10520157)) ~ dos_data$rs4417927)))
summary(lm(scale(residuals(lm(scale(dos_data$ENSG00000137601, scale = F) ~ dos_data$rs10520157)),scale=F) ~ dos_data$rs4417927))
var(residuals((lm(residuals(lm(scale(dos_data$ENSG00000137601, scale = F) ~ -1 + scale(dos_data$rs10520157,scale = F))) ~ scale(dos_data$rs4417927, scale = F)))))

sigma_c = (((t(scale(dos_data$ENSG00000137601,scale=F)))%*%(scale(dos_data$ENSG00000137601, scale=F)) -XX[1,1] %*% beta_hat[1]^2))/(450)
sigma_c

sigma_c = (t(scale(dos_data$ENSG00000137601,scale=F))%*%(scale(dos_data$ENSG00000137601, scale=F)) - t(beta_hat[1]) %*% XX[1,1] %*% beta_hat[1] - (t(beta_cond) %*% XX[2,2] %*% beta_cond)) / 450
sqrt(sigma_c)
0.1097425

   #        - t(beta_cond) %*% XX[2,2] %*% beta_cond) /(452-)
#summary(lm(scale(dos_data$ENSG00000137601, scale = F) ~dos_data$rs10520157))

 # 
 # (452-1)
           
   #        /(452-3)



#sqrt(0.01199013573 %*% solve(X2X) - 0.01199013573*solve(X2X) %*% X2X1 %*% solve(X1X) %*% X2X1 %*% solve(X2X))



#%*%t(scale(dos_data$rs10520157,scale=F)) 

#iag((sigma_sq[1]) *solve(XX))
