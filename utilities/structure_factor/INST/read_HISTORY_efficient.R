# read the HISTORY file to get coordinates
rm(list = ls())
library(tidyverse)

to_vector <- function(string){
  string <- trimws(string)
  a <- strsplit(string, split = "\\s+")[[1]]
  as.numeric(a)
}

inputFile <- "HISTORY"
his_file  <- file(inputFile, open = "r")

readLines(con = his_file, n = 2, warn = FALSE)

nstep <- 199
n_water <- 1000 
qO <- -0.8476
qH <- 0.4238

df_corrs_dipoles <- matrix(0, nrow = (n_water*nstep), ncol = 7)
for (istep in 1:nstep){
  oneline <- readLines(con = his_file, n = 1, warn = FALSE)
  oneline <- readLines(con = his_file, n = 1, warn = FALSE)
  boxx <- to_vector(oneline)[1]
  oneline <- readLines(con = his_file, n = 1, warn = FALSE)
  boxy <- to_vector(oneline)[2]
  oneline <- readLines(con = his_file, n = 1, warn = FALSE)
  boxz <- to_vector(oneline)[3]
  
  pbc <- c(boxx, boxy, boxz)
  corrs <- readLines(con = his_file, n = 6*n_water, warn = FALSE)
  useless_ind <- grepl("OS|HS", corrs)
  corrs_new <- unname(t(sapply(corrs[!useless_ind], to_vector)))
  indOS <- (1:(3*n_water) %% 3 == 1)
  indHS1 <- (1:(3*n_water) %% 3 == 2)
  indHS2 <- (1:(3*n_water) %% 3 == 0)
  
  corrs_OS <- corrs_new[indOS,]
  corrs_HS1 <- corrs_new[indHS1,]
  corrs_HS2 <- corrs_new[indHS2,]
  
  dis_OH1 <- (corrs_HS1 - corrs_OS  ) - round(( corrs_HS1 - corrs_OS )/pbc)*pbc
  dis_OH2 <- (corrs_HS2 - corrs_OS ) - round((corrs_HS2 - corrs_OS )/pbc)*pbc
  
  dipole <- qH * (dis_OH1 + dis_OH2) # the dipole of a water is about 0.24
  
  ibegin <- (istep-1)*n_water + 1
  iend <- istep*n_water 
  df_corrs_dipoles[ibegin:iend,] <- cbind(rep(istep, n_water), corrs_OS/pbc, dipole)
}

colnames(df_corrs_dipoles) <- c("istep", "Ox", "Oy", "Oz", "px", "py", "pz")

saveRDS(as.data.frame(df_corrs_dipoles), file = "df_corrs_dipoles")
