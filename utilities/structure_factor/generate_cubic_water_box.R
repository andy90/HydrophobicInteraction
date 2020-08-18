rm(list = )
# generate water box
N = 1000 # number of water
sigma = 3.16 # sigma of water

n_water_x <- as.integer(N**(1/3)) + 1

O = c(0, 0, 0) # coordinates of reference oxygen
H1 = c(1, 0, 0)
H2 = c(cos(109.47/180*pi), sin(109.47/180*pi), 0)

corr_water <- rbind(O, H1, H2)

l = sigma # distance between water molecules

corrs <- numeric()
for (i in 1:n_water_x){
  for (j in 1:n_water_x){
    for (k in 1:n_water_x){
      dis_x <- i*l
      dis_y <- j*l
      dis_z <- k*l
      
      dis <- c(dis_x, dis_y, dis_z)
      corrs <- rbind(corrs, sweep(corr_water, 2 ,dis, "+"))
    }
  }
}

atoms <- rep(c("OS", "HS", "HS"), n_water_x**3)

df_corrs <- data.frame("atom" = atoms, corrs)

line <- sprintf("%-80s", "bulk spce water")
write(line,file="CONFIG")
line <- sprintf("%10d%10d", 0 , 1)
write(line,file="CONFIG", append = TRUE)
line <- sprintf("%20.12f%20.12f%20.12f", l*n_water_x, 0, 0)
write(line,file="CONFIG", append = TRUE)
line <- sprintf("%20.12f%20.12f%20.12f", 0, l*n_water_x, 0)
write(line,file="CONFIG", append = TRUE)
line <- sprintf("%20.12f%20.12f%20.12f",  0, 0, l*n_water_x)
write(line,file="CONFIG", append = TRUE)

for (i in 1:(3*N)){
  line <- sprintf("%-10s%8d",  df_corrs$atom[i], i)
  write(line,file="CONFIG", append = TRUE)
  line <- sprintf("%20.12f%20.12f%20.12f",  df_corrs[i,2], df_corrs[i,3], df_corrs[i,4])
  write(line,file="CONFIG", append = TRUE)
  
}
