# importing libraries -----------------------------------------------------
library(tidyverse)
library(MASS)
library(Rcpp)
# SNP selection (preprocessing) -------------------------------------------
data <- read.table('chromosom29.txt', sep = ',', header = FALSE)
data_cleaning <- function(data){
  data$V1 <- as.character(data$V1)
  n <- length(unique(data$V2))
  individuals <- unique(data$V2)
  n_occur <- data.frame(table(data$V1))
  n_occur$percentage <- ((n - n_occur$Freq)*100)/n
  
  df1 <- data[data$V1 %in% n_occur$Var1[n_occur$percentage == 0.0],]
  df2 <- data[data$V1 %in% n_occur$Var1[(n_occur$percentage < 1.0) & (n_occur$percentage > 0.0)],]
  most_frequent_allels <- data %>%
    group_by(V1) %>%
    count(V3,V4) %>%
    slice(which.max(n))
  rm(data)
  
  most_frequent_allels$V3 <- as.character(most_frequent_allels$V3)
  most_frequent_allels$V4 <- as.character(most_frequent_allels$V4)
  missing_individuals = as.character(individuals[!(individuals %in% df2$V2)])
  missing_markers <- length(unique(df2$V1))
  V2 <- c(replicate(missing_markers,missing_individuals))
  V11 <- c(replicate(length(missing_individuals),unique(df2$V1)))
  V3 <- c()
  V4 <- c()
  a <- c(unique(df2$V2))
  b <- c(unique(df2$V1))
  c <- expand.grid(b,a)
  df2$g <- paste(df2$V1,df2$V2)
  c$g <- paste(c$Var1,c$Var2)
  single_missing <- subset(c, !(g %in% df2$g),1:2)
  if (nrow(single_missing) != 0){
    colnames(single_missing) <- c('V1','V2')
    single_missing$V1 <- as.character(single_missing$V1)
    V2 <- c(V2,single_missing$V2)
    V11 <- c(V11,single_missing$V1)
  }
  for (i in 1:length(V11)){
    V3 <- c(V3,(most_frequent_allels[most_frequent_allels$V1 == V11[i],]$V3))
    V4 <- c(V4,(most_frequent_allels[most_frequent_allels$V1 == V11[i],]$V4))
  }
  single_missing_data <- data.frame('V1' = V11, 'V2'= V2, 'V3' = V3, 'V4' = V4)
  cleaned_data <- rbind(df1,df2[1:4],single_missing_data)
  adding_column <- function(data){
    tmp <- data[data$V3 == 'A' & data$V4 == 'B',]
    tmp$V5 <- 1
    tmp2 <- data[data$V4 == 'A',]
    tmp2$V5 <- 2
    tmp3 <- data[data$V3 == 'B',]
    tmp3$V5 <- 0
    tmp4 <- rbind(tmp,tmp2,tmp3)
    return(tmp4)
  }
  cleaned_data <- adding_column(cleaned_data)
  cleaned_data <- arrange(cleaned_data,V2,V1)
  return(cleaned_data)
}
data <- data_cleaning(data)
# Matrix G ----------------------------------------------------------------
G_Matrix <- function(data){
  m <- length(unique(data$V1))
  n <- length(unique(data$V2))
  M <- matrix(data$V5, ncol = m, nrow = n, byrow = TRUE)
  p = apply(M, 2, mean)/2
  p[p<0.5] = 1-p[p<0.5]
  P = 2 * (p-0.5)
  M = M - 1
  Z <- sweep(M,2,P)
  q <- 1 - p
  sum2pq <- 2*sum(p*q)
  Z2 <- Z%*%t(Z)
  G <- Z2/(sum2pq)
  return(G)
}

start_time <- Sys.time()
G <- G_Matrix(data)
end_time <- Sys.time()
time_of_G <- end_time - start_time
# importing phenotype data ------------------------------------------------
phenotype <- read.csv('DANESTU.csv', sep = ';', header = TRUE)
individuals <- read.table('individuals.txt')
phenotype$ID <- individuals
phenotype <- phenotype[order(phenotype$ID),1:8]
rm(individuals)
# setting up variables to BLUP --------------------------------------------
y1 <- matrix(phenotype$BV_MS)
y2 <- matrix(phenotype$BV_TEMP)
Z <- diag(1,7646,7646)
X <- matrix(1,7646,1)
sigma_y1 = var(y1)
sigma_y2 = var(y2)
sigma_k1 = 0.3*sigma_y1 #starting value for kegg effect
sigma_e1 = 0.7*sigma_y1 #starting value for error variance
sigma_k2 = 0.3*sigma_y2
sigma_e2 = 0.7*sigma_y2
# Cpp matrix multiplication -----------------------------------------------
sourceCpp("MatrixInverse.cpp")
# GBLUP -------------------------------------------------------------------
mme = function(y, X, Z, K, sigma_k, sigma_e) {
  G = ginv(K) * c(sigma_e)/c(sigma_k)
  C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
            cbind(t(Z)%*%X, t(Z)%*%Z+G))
  rhs = rbind(t(X)%*%y, t(Z)%*%y)
  mme = ginv(C)%*%rhs
  list(C = C, est = mme)
}
# finding sigma with EM algorithm -----------------------------------------------------------
EM = function(y, X, Z, K, sigma_k, sigma_e) {
  n = nrow(X)
  n_X = ncol(X)
  n_k = nrow(K)
  t = 1 #iteration number 1
  tmp = 0.1 #test for convergance
  print(t)
  while (tmp > 0.0001) {
    mme_new = mme(y, X, Z, K, sigma_k, sigma_e)
    C_new = ginv(mme_new$C)
    Ck = C_new[(n_X+1):(n_X+n_k), (n_X+1):(n_X+n_k)]
    mme2 = mme_new$est
    k = as.matrix(mme2[(n_X+1):(n_X+n_k)])
    sigma_k_new = (t(k)%*%ginv(K)%*%k + sum(diag(ginv(K)%*%Ck))*c(sigma_e))/n_k
    res = as.matrix(y-X%*%as.matrix(mme2[1:n_X]) - Z%*%as.matrix(mme2[(n_X+1):(n_X+n_k)]))
    X.tmp1 = cbind(X,Z) %*% C_new
    X.tmp2 = t(cbind(X,Z))
    sigma_e_new = (t(res)%*%res + sum(diag(eigenMatMult(X.tmp1,X.tmp2)))*c(sigma_e))/n
    tmp = max(abs(sigma_k - sigma_k_new), abs(sigma_e - sigma_e_new))
    sigma_k = sigma_k_new
    sigma_e = sigma_e_new
    t = t + 1
    write.table(c(sigma_k,sigma_e,t, tmp), "New_Results2.csv", col.names = FALSE, row.names = FALSE, append = TRUE, quote = FALSE)
  }
}
#EM(y1,X,Z,G,sigma_k1,sigma_e1)
#EM(y2,X,Z,G,sigma_k2,sigma_e2)
# sigma for my data -------------------------------------------------------
sigma_k1 <- 0.0000001301411442038
sigma_e1 <- 71.5344015633313
sigma_k2 <- 0.000000140127190931774
sigma_e2 <- 81.1969163466016
# calculation of heritability ---------------------------------------------
h2_1 = sigma_k1/(sigma_k1 + sigma_e1)
h2_2 = sigma_k2/(sigma_k2 + sigma_e2)
# estimating random effects -----------------------------------------------
u1 <- mme(y1,X,Z,G,sigma_k1,sigma_e1)
u2 <- mme(y2,X,Z,G,sigma_k2,sigma_e2)

u1 <- u1$est
u2 <- u2$est

a1 <- as.matrix(u1[2:7647])
b1 <- u1[1]

a2 <- as.matrix(u2[2:7647])
b2 <- u2[1]

write.table(a1, "MS_random_effects.csv", col.names = FALSE, row.names = FALSE)
write.table(a2, "TEMP_random_effects.csv", col.names = FALSE, row.names = FALSE)
# SNP effects -------------------------------------------------------------
k=0.2
A <- diag(1,7646,7646)
m <- length(unique(data$V1))
n <- length(unique(data$V2))

M <- matrix(data$V5, ncol = m, nrow = n, byrow = TRUE)
p = apply(M, 2, mean)/2
p[p<0.5] = 1-p[p<0.5]
P = 2 * (p-0.5)
M = M - 1
Z <- sweep(M,2,P)

f <- ((1-k)/(951))^(-1)
s <- 5 * t(Z) %*% solve(A)%*% Z
fs = f+s
invfs = ginv(fs)
snp1 <- (invfs * 5) %*% t(Z) %*% solve(A) %*% a1
snp2 <- (invfs * 5) %*% t(Z) %*% solve(A) %*% a2

write.table(snp1, "MS_SNP_effects.csv", col.names = FALSE, row.names = FALSE)
write.table(snp2, "TEMP_SNP_effects.csv", col.names = FALSE, row.names = FALSE)
# Finding efficient SNP ---------------------------------------------------
av1 <- mean(snp1)
sd1 <- sd(snp1)
Z1 = (snp1-av1)/sd1

p_value1 = 2*(1 - pnorm(abs(Z1),0,1))
bonferroni1 <- p.adjust(p_value1, "bonferroni")
efficientSNP1 <- which(bonferroni1 < 0.05)

av2 <- mean(snp2)
sd2 <- sd(snp2)
Z2 <- (snp2-av2)/sd2

p_value2 = 2*(1 - pnorm(abs(Z2),0,1))
bonferroni2 <- p.adjust(p_value2, "bonferroni")

efficientSNP2 <- which(bonferroni2 < 0.05)
efficientSNP <- data$V1[efficientSNP2]