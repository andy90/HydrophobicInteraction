# the dipole-dipole correlation function in k space
rm(list = ls())
library(foreach)
library(doParallel)
registerDoParallel(cores=2)
library(Rcpp)
sourceCpp("c_imp.cpp")
n_water <- 1000
df_corrs_dipoles <- readRDS(file = "df_corrs_dipoles")

structure_factor <- function(corrs){
  dis_corrs <- foreach(i = 1:n_water, .combine = 'rbind') %do% {
    sweep(corrs, 2, corrs[i, ])
  }
  
  k <- 1:50 *(2*pi) # the k space vectors we are interested in
  colSums(sapply(k, function(ki){
    colSums(cos(ki*dis_corrs))
  }))/3
}

pp_longitudinal <- function(corrs, dipoles){
  dis_corrs <- foreach(i = 1:n_water, .combine = 'rbind') %do% {
    sweep(corrs, 2, corrs[i, ])
  }
  
  mult_dipoles <- foreach(i = 1:n_water, .combine = 'rbind') %do% {
    sweep(dipoles, 2, dipoles[i, ], "*")
  }
  
  k <- 1:50 *(2*pi) # the k space vectors we are interested in
  colSums(sapply(k, function(ki){
    colSums(cos(ki*dis_corrs)*mult_dipoles)
  }))/3
}

# corrs <- as.matrix(df_corrs_dipoles[1001:2000, 2:4])
# dipoles <- as.matrix(df_corrs_dipoles[1001:2000, 5:7])
# a <- structure_factor(corrs)
# b <- structure_factorC(corrs)
# c <- pp_longitudinal(corrs, dipoles)
# d <- pp_longitudinalC(corrs, dipoles)
# 
# structure_factor(corrs)

s_t <- foreach(istep = 1:199, .combine = 'rbind') %dopar%{
  ibegin <- (istep-1)*n_water + 1
  iend <- istep*n_water 
  structure_factorC(as.matrix(df_corrs_dipoles[ibegin:iend, 2:4]))
}

colMeans(s_t)
saveRDS(s_t, file = "s_t")
plot(1:50 *(2*pi/31.1), colMeans(s_t))


pp_t <- foreach(istep = 1:199, .combine = 'rbind') %dopar%{
  ibegin <- (istep-1)*n_water + 1
  iend <- istep*n_water 
  pp_longitudinalC(as.matrix(df_corrs_dipoles[ibegin:iend, 2:4]), as.matrix(df_corrs_dipoles[ibegin:iend, 5:7]))
}

colMeans(pp_t)
saveRDS(pp_t, file = "pp_t")
plot(1:50 *(2*pi/31.1), colMeans(pp_t))
# df_corrs_dipoles %>%
#   groub_by(istep) %>%
#   
# #a <- apply(corrs, 1, function(Oxyz){
# #  sweep(corrs, 2, Oxyz)
# #}) 
# 
# dis_corrs <- foreach(i = 1:n_water, .combine = 'rbind') %do% {
#   sweep(corrs, 2, corrs[i, ])
# }
# 
# mult_dipoles <- foreach(i = 1:n_water, .combine = 'rbind') %do% {
#   dipoles*dipoles[i, ]
# }
# 
# k <- 1:10 *(2*pi) # the k space vectors we are interested in
# 
# colSums(cos(k[5]*dis_corrs)) # this is the structure factor
# colSums(cos(k[1]*dis_corrs)*mult_dipoles) # this is the structure factor
