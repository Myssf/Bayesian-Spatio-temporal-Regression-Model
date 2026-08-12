# setwd("D:/R/Group/Code/Simulation/")
#if running on my computer 
if(getwd()== "D:/R/Group/Code/Simulation"){
  path1 = "D:/R/Group/Code/Simulation/"
} else{ 
  path1 = "/work/chaom/ST_Sim/" #Running on cluster
}

library(MASS)
library(BayesLogit)
library(Matrix)
library(spdep)
library(invgamma)
library(usmap)
library(ggplot2)
library(sf)
library(coda)
library(dplyr)
library(units)
library(stringr)
library(ggmap)
library(RColorBrewer)
library(lwgeom)
library(terra)
library(Rcpp)
library(RcppArmadillo)
sourceCpp(paste0(path1, "step1.cpp"))

###########################
#ID
sim.id <- 6

###########################
#Iteration
n.iter <- 5000

###########################
# Load data
df.sim <- read.csv(paste0(path1,"sim.setup.csv"), fileEncoding="UTF-8-BOM")

sim.ind <- which(df.sim$sim.id==sim.id)

fips <- read.table(paste0(path1,'County FIPS Codes.txt'),sep = "\t",header = T)
fips.sc <- fips$FIPS[fips$State=="SC"]

#table(df.PCR$Pool.Size)
###########################
#1. Sample size

N <- df.sim$N[sim.ind] #automatically choose N
################################ 

#2. Generating surface from
gen.from <- df.sim$surface[sim.ind] 

################################ 

#3. Prevalence
p.seq <- c(0.09,0.18)

p.sim <- df.sim$p[sim.ind]

if(gen.from=="GPP"){
  if(p.sim==0.09){
    beta.true <- c(-2.5,-1.2,0.2,1) #p=0.09
  } else if(p.sim==0.18){
    beta.true <- c(-1.7,-1.0,0.2,1) #p=0.18
  }
}
#####################################

#4. Pool size

n.size.seq <- c(6,3)
n.size <- df.sim$n.size[sim.ind]

#####################################

#5. Number of knots m

m <- df.sim$m[sim.ind]
#####################################

#6. Sensitivity and specificity
#6.1 Sensitivity
sej.seq <- as.numeric(unlist(lapply(as.character(df.sim$sej),FUN=strsplit,",")[sim.ind]))
#6.2 Specificity
spj.seq <- as.numeric(unlist(lapply(as.character(df.sim$spj),FUN=strsplit,",")[sim.ind]))

#####################################

#7. Pool structure
#7.1 from same county
#7.2 randomly assigned
G.mat.seq <- c("sequential","random")

#G.mat.config <- "sequential"
G.mat.config <- "random"

#####################################

#8. Covariate estimation/Prevalence estimation

wd= getwd()
wd

#####################################

#9. Prediction
pred.from <- df.sim$pred[sim.ind] 


########################### Data Generation ################################

###### Individual Testing (IT) ###### 
#Using IT as Baseline testing

#Assume ""Sampling.Method" has 3 levels and 1 continuous 
#X <- model.matrix(~ Life.Stage + Sampling.Method)

n.b <- length(beta.true)
X <- matrix(NA, nrow = N, ncol = n.b)
X[,1] <- 1
Sampling.Method <- factor(sample(c("S1", "S2", "S3"), size = N, replace = TRUE, prob = c(1/3, 1/3, 1/3)))
# Sampling.Method (reference = "S1")
X[,2] <- ifelse(Sampling.Method == "S2", 1, 0)
X[,3] <- ifelse(Sampling.Method == "S3", 1, 0)
X[,4] <- rnorm(N, mean = 0, sd = 1)
beta <- beta.true

dim(X)

#Pool testing
#Split the individual testing data into batches
n.pools <- N/n.size  #200


n.assay <- df.sim$n.assay[sim.ind]
# df.sim$n.assay[sim.ind]

sens <- sej.seq[1:n.assay]
spec <- spj.seq[1:n.assay]

df.assay <- data.frame(id = c(1:n.assay), sens = sens, spec = spec) #sens.1 = 0.95; sens.2 = 0.98; spec.1 = 0.95; spec.2 = 0.98



######################################### Spatial random effects ############################################

# #Assign FIPS code (46 counties) so that each pool only contains test sample from individuals from same county
n.county <- length(fips.sc) #46
fips.sub <- fips.sc[1:n.county]
#fips.ind <- sort(rep(fips.sub,length.out=N))

# Assign FIPS code (46 counties) to a Time from 1:3
T <- 3
years.ind <- sample(1:T, size = N, replace = TRUE, prob = rep(1/T, T))

#Spatial random effects (from CAR model)
#Create W matrix and D matrix
if(!exists("shp")){
  shp <- st_read(paste0(path1,"cb_2018_us_county_500k.shp")[1]) #Need to also have .dbf and .shx in the folder w/ .shp
  #Subset shp with only the fips code appearing in fips.sub
  shp <- shp[shp$GEOID%in%fips.sub,]
}

neighbours <- poly2nb(shp)
W <- Matrix(as(nb2mat(neighbours,style = 'B', zero.policy = TRUE), "Matrix"), sparse = TRUE)
D <- Diagonal(n = ncol(W),colSums(W))

#inverse gamma
#Hyper Parameters for ell_s and ell_t
a.sigma <- 3
b.sigma <- 3

a.l.s <- 3
b.l.s <- 3

a.l.t <- 3
b.l.t <- 3



#Sample spatial locations  
shp.st <- shp
shp.sf <- st_as_sf(shp.st)

#Ground truth of ell_s and ell_t
ell_s.true <- ell_s <- 3
ell_s.sigma <- 0.3
ell_s.init <- 0.5
ell_s.mcmc <- rep(NA,n.iter)
ell_s.mcmc[1] <- ell_s.init

ell_t.true <- ell_t <- 2
ell_t.sigma <- 0.5
ell_t.init <- 0.5
ell_t.mcmc <- rep(NA,n.iter)
ell_t.mcmc[1] <- ell_t.init


sigma.sq_SE <- 0.5
sigma.sq.true <- sigma.sq_SE
sigma.sq.init <- 0.2

st_crs(shp.sf) <- 4326
loc <- st_sample(shp.sf, size = N, crs = 4326)

loc_geom <- st_sample(shp.sf, size = N)
attr_df <- data.frame(Year = years.ind[1:N])
loc <- st_sf(attr_df, geometry = loc_geom)
loc <- st_transform(loc, crs = 3086)
loc$Year


loc.coordinates <- st_coordinates(loc)

id <- sample(1:N, m/T)  
base_knot <- loc[id, ]
knot_data <- data.frame()
for (t in 1:T) {
  temp_knot <- base_knot
  temp_knot$Year <- t
  knot_data <- rbind(knot_data, temp_knot)
}
knot <- st_sf(knot_data)

scale <- 111000

dist.cs <- (st_distance(loc, knot)/ scale) ^2
dist.cs <- as.matrix(drop_units(dist.cs))
dim(dist.cs)

dist.ct <- outer(loc$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
dim(dist.ct)

hist(exp(-dist.cs/(2*ell_s)))
hist(exp(-dist.ct/(2*ell_t)))

c.mat <- sigma.sq_SE*exp(-dist.cs/(2*ell_s)-dist.ct/(2*ell_t)) 
hist(c.mat)

dim(dist.cs) # 2000 by 100
dim(dist.ct) # 2000 by 100
dim(c.mat) # 2000 by 100

dist.Cs <- (st_distance(knot,knot) / scale)^2
dist.Cs <- as.matrix(drop_units(dist.Cs))
dim(dist.Cs)

dist.Ct <- outer(knot$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
dim(dist.Ct)

C.prime <- exp(-dist.Cs/(2*ell_s)-dist.Ct/(2*ell_t)) 
C.mat <- sigma.sq_SE*C.prime  + diag(0.00001,m) #with a nugget
dim(dist.Cs) # 100 by 100
dim(dist.Ct) # 100 by 100
dim(C.mat) # 100 by 100



xi <- mvrnorm(1,rep(0,m),C.mat)   #m dim
xi <- xi - mean(xi)

if(gen.from=="GPP"){
  #(Option 2) Generate data from GPP
  xi <- mvrnorm(1,rep(0,m),C.mat)    #m dim
  xi.mean <- mean(xi)
  xi <- xi - xi.mean
  xi.true <- c.mat%*%solve(C.mat)%*%xi
  p <- as.vector(exp(X%*%beta + xi.true)/(1+exp(X%*%beta + xi.true)))
  p.true <- p
  mean(p.true)
}
hist(p.true)
mean(p.true)
#hist(as.vector(C.mat)/sigma.sq_SE)

#Make prevalence maps with grid 

# define dimensions, create grid

n.grid = 3000
cellsize <- 0.05*scale

if (is.na(st_crs(shp.sf))) {
  shp.sf <- st_set_crs(shp.sf, 4326)
}

total_years <- if (pred.from == "Yes") c(1, 2, 3, 4) else sort(unique(loc$Year))

# Transform to EPSG:3086
shp.sf <- st_transform(shp.sf, crs = 3086)
grid_list <- list()

for (yr in total_years) {
  grid <- shp.sf %>% 
    st_make_grid(cellsize = cellsize, what = "centers", n = c(n.grid/2, n.grid/2), square = FALSE) %>% 
    st_intersection(shp.sf)
  
  grid.sf <- st_as_sf(grid)
  grid.sf$Year <- yr
  grid_list[[as.character(yr)]] <- grid.sf
}


# Store grid coordinates once (same for all years/sims)
grid_coords <- do.call(rbind, lapply(grid_list, function(g) st_coordinates(g)))
lat <- st_coordinates(grid.sf)[,"Y"]
long <- st_coordinates(grid.sf)[,"X"]

write.table(data.frame(lat = lat, long = long, grid_id = 1:length(lat)), 
            file = paste0(path1, "grid_coordinates.txt"), row.names = FALSE, col.names = TRUE)


xi.grid.true_list <- list()
p.grid.true_list <- list()

# Create true p.grid
for (yr in total_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  # Spatial distances: grid points to knots
  dist.cs.grid <- (st_distance(grid.sf, knot) / scale)^2
  dist.cs.grid <- as.matrix(drop_units(dist.cs.grid))
  # Temporal distances: grid points to knot years
  dist.ct.grid <- outer(grid.sf$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
  # Cross-covariance matrix
  c.mat.grid <- sigma.sq.true * exp(-dist.cs.grid / (2 * ell_s.true) - dist.ct.grid / (2 * ell_t.true))
  # True spatial random effects
  xi.grid.true <- c.mat.grid %*% solve(C.mat) %*% xi
  xi.grid.true_list[[as.character(yr)]] <- xi.grid.true
  # True prevalence
  p.grid.true <- exp(beta[1] + xi.grid.true) / (1 + exp(beta[1] + xi.grid.true))
  p.grid.true_list[[as.character(yr)]] <- p.grid.true
}


max.lat = max(st_coordinates(shp.st)[,2])
min.lat = min(st_coordinates(shp.st)[,2])

max.long = max(st_coordinates(shp.st)[,1])
min.long = min(st_coordinates(shp.st)[,1])

# display the true prevalence for all years

plot.prev.true_list <- list()

for (yr in total_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  coords <- st_coordinates(grid.sf)
  grid.sf$long <- coords[, "X"]
  grid.sf$lat <- coords[, "Y"]
  grid.sf$p.grid.true <- p.grid.true_list[[as.character(yr)]]
  
  min.long <- min(st_coordinates(shp.sf)[, "X"]) - 1000
  max.long <- max(st_coordinates(shp.sf)[, "X"]) + 1000
  min.lat <- min(st_coordinates(shp.sf)[, "Y"]) - 1000
  max.lat <- max(st_coordinates(shp.sf)[, "Y"]) + 1000
  
  plot.prev.true <- ggplot() +
    geom_sf(data = grid.sf, aes(colour = p.grid.true), shape = 15, size = 2) +
    scale_colour_gradient(low = "white", high = "darkred", limits = c(0, 0.5 ), name = "True Prevalence") +
    coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
    geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.3) +
    theme_bw() +
    ggtitle(paste("True Prevalence in Year", yr))
  
  plot.prev.true_list[[as.character(yr)]] <- plot.prev.true
}

plot.prev.true_list

#Specify the true status Y_tilde
Y_tilde <- rbinom(N,1,p)
Y_tilde.true <- Y_tilde

#List indicating which individual tested by which assay
assay.ind.list <- sens.ind.list <- spec.ind.list <- rep(NA, n.assay)

#Sequential assignment of assay to individuals
for(i in 1:n.assay){
  assay.ind.list[((i-1)*(N/n.assay)+1) : (i*N/n.assay)] <- df.assay$id[df.assay$id[i]]
  sens.ind.list[((i-1)*(N/n.assay)+1) : (i*N/n.assay)] <- df.assay$sens[df.assay$id[i]]
  spec.ind.list[((i-1)*(N/n.assay)+1) : (i*N/n.assay)] <- df.assay$spec[df.assay$id[i]]
}



#Generate the individual test status
Y <- rep(NA,N)
for(i in 1:n.assay){
  Y[Y_tilde==1&sens.ind.list==df.assay$sens[i]] <- rbinom(length(Y_tilde[Y_tilde==1&sens.ind.list==df.assay$sens[i]]),1,df.assay$sens[i])
  Y[Y_tilde==0&spec.ind.list==df.assay$spec[i]] <- rbinom(length(Y_tilde[Y_tilde==0&spec.ind.list==df.assay$spec[i]]),1,1-df.assay$spec[i])
}

###### Pool Testing (PT) ######

#Create a list to save the indices
#Create a list for each array to record individual's id (i.e. array id to individual id's)  

pool_list <- list()
for (i in 1:(N/n.size)){
  pool_list[[i]] <- seq(i*n.size-(n.size-1),i*n.size)
}


#Create a list for each individual to record the pool being contributed to (i.e. individual id to pool id's) 
ind_list <- list()
for(i in 1:N){
  ind_list[[i]] <- which(lapply(lapply(lapply(pool_list,"%in%",i),as.numeric),sum)==1)
}

#True status for each pool
Z.tilde <- c() #n.pool by 1
for(i in 1:n.pools){
  Z.tilde[i] <- ifelse(sum(Y_tilde[unlist(pool_list[i])])>0,1,0)
}


#List indicating which pool tested by which assay
assay.pool.list <- sens.pool.list <- spec.pool.list <- rep(NA,n.pools)
assay.pool.list <- sample(1:n.assay, size = n.pools, replace=T, prob=rep(1/n.assay, n.assay))

for(i in 1:n.assay){
  sens.pool.list[assay.pool.list==i] <- df.assay$sens[i] 
  spec.pool.list[assay.pool.list==i] <- df.assay$spec[i] 
}

success_counts <- rep(0, n.iter)


Z <- rep(NA, n.pools)
for(i in 1:n.assay){
  #get(paste0("sens.",i))
  Z[assay.pool.list==i] <- rbinom(sum(assay.pool.list==i),1, df.assay$sens[i]*Z.tilde[assay.pool.list==i]+(1-df.assay$spec[i])*(1-Z.tilde[assay.pool.list==i]))
}

Z_0 <- Z
length(Z)
sum(Z)



n.tests <- length(Z)



ind_list <- list()
for(i in 1:N){
  ind_list[[i]] <- which(lapply(lapply(lapply(pool_list,"%in%",i),as.numeric),sum)==1)
}


#Gibbs sampler for beta and w
n.b <- length(beta)

#priors  
b.m <- rep(0,n.b)
b.cov <- diag(10,n.b)
b.cov.inv <- solve(b.cov)

#Initial values 
b <- matrix(NA,n.iter,n.b)
dim(b)
w <- matrix(NA,n.iter,N)
dim(w)

b[1,] <- rep(0.5,n.b)
w[1,] <- rep(0.5,N)
Y_tilde.mcmc <- matrix(0,n.iter,N)
Y_tilde <- rep(0,N)






xi.N <- c.mat%*%solve(C.mat)%*%xi
xi.mcmc  <- matrix(NA, n.iter, N)
xi.mcmc[1, ] <- xi.true

p.mcmc <- matrix(NA,n.iter,N)
p.mcmc[1,] <-  as.vector(exp(X%*%b[1,] + xi.N)/(1+exp(X%*%b[1,] + xi.N)))

# Initialize storage for p.grid.mcmc for each year
p.grid.mcmc_list <- lapply(total_years, function(yr) matrix(NA, nrow = n.iter, ncol = nrow(grid_list[[as.character(yr)]])))
names(p.grid.mcmc_list) <- total_years

xi.grid.mcmc_list <- lapply(total_years, function(yr) matrix(NA, nrow = n.iter, ncol = nrow(grid_list[[as.character(yr)]])))
names(xi.grid.mcmc_list) <- total_years


sigma.sq.mcmc <- rep(NA,n.iter)
sigma.sq.mcmc[1] <- sigma.sq <- sigma.sq.init

dist_cs_grid_list <- list()
dist_ct_grid_list <- list()

for (yr in total_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  dist.cs.grid <- (st_distance(grid.sf, knot) / scale)^2
  dist.cs.grid <- as.matrix(drop_units(dist.cs.grid))
  dist_cs_grid_list[[as.character(yr)]] <- dist.cs.grid
  dist.ct.grid <- outer(grid.sf$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
  dist_ct_grid_list[[as.character(yr)]] <- dist.ct.grid
}

omg <- Diagonal(x=w[1,])

xi.knots.mcmc <- matrix(NA, n.iter, m)

######## Use true sej and spj for debugging only (fixing sej, spj at true values) #######


sej <- sens.pool.list
spj <- spec.pool.list

###############################            RUN!!!!    #########################################   
total_years <- sort(unique(loc$Year))

# Initialize separate acceptance counters
ell_s.acc.count <- 0
ell_t.acc.count <- 0

run.start <- Sys.time()


for (t in 2:n.iter){
  #Step 1 Sample Y_tilde (Pool testing) 
  C.inv <- chol2inv(chol(C.mat))
  cC.inv <- c.mat %*% C.inv  
  cC_inv_xi <- cC.inv %*% xi
  
  #Call cpp function to sample Y_tilde
  result <- sample_Y_tilde(X, b, cC_inv_xi, Z, sej, spj, ind_list, pool_list, Y_tilde, t, N, success_counts)
  
  Y_tilde <- result$Y_tilde
  success_counts[t - 1] <- result$success_count  # Store success_count for this iteration
  Y_tilde.mcmc[t, ] <- Y_tilde
  
  #Step 2 Compute Z.tilde based on Y_tilde
  Z.tilde <- sapply(pool_list, function(p) {
    ifelse(sum(Y_tilde[unlist(p)]) > 0, 1, 0)
  })
  
  #Step 3 Sample sens and spec
  #Sensitivity & Specificity
  
  #sum_ZZt <- tapply(Z * Z.tilde, assay.pool.list, sum)
  #sum_oneZ_Zt <- tapply((1 - Z) * Z.tilde, assay.pool.list, sum)
  #sum_oneZ_oneZt <- tapply((1 - Z) * (1 - Z.tilde), assay.pool.list, sum)
  #sum_Z_oneZt <- tapply(Z * (1 - Z.tilde), assay.pool.list, sum)
  
  # Compute beta parameters per group
  #aem <- a.sens + sum_ZZt
  #bem <- b.sens + sum_oneZ_Zt
  #apm <- a.spec + sum_oneZ_oneZt
  #bpm <- b.spec + sum_Z_oneZt
  
  # Draw one beta random variable per group
  #se_per_group <- sapply(seq_along(aem), function(k) rbeta(1, aem[k], bem[k]))
  #sp_per_group <- sapply(seq_along(apm), function(k) rbeta(1, apm[k], bpm[k]))
  
  # Assign to sej and spj by mapping back to original positions via assay.pool.list
  
  #sej <- se_per_group[assay.pool.list]
  #spj <- sp_per_group[assay.pool.list]
  
  
  #sej.mcmc[t,] <- sej
  #spj.mcmc[t,] <- spj
  

  
  #Step 4 Sample beta
  
  M <- t(X)%*%omg%*%X+b.cov.inv
  k <- Y_tilde - 1/2
  h <- k/diag(omg)
  Q <- t(X)%*%omg%*%(h-cC_inv_xi) + b.cov.inv%*%b.m
  M.inv <- chol2inv(chol(M))   
  b[t,] <- as.vector(mvrnorm(n = 1,M.inv%*%Q,M.inv))
  #b[t,] <- beta #remove if update b
  
  ##Step 5 Sample w
  w[t,] <- apply(X = X%*%b[t,] + cC_inv_xi,1, FUN = rpg, num=1, h=1)
  
  #update omg and h
  omg <- Diagonal(x=w[t,])
  h <- k/diag(omg)
  
  
  #Step 6 Sample xi
  omg.cC.inv <- omg%*%cC.inv
  V <- t(cC.inv)%*%omg.cC.inv + C.inv
  d <- t(t(h-X%*%b[t,])%*%omg.cC.inv)
  #d <- as.vector((t(h)-t(b[t,])%*%t(X))%*%omg%*%cC.inv) #equivalent to above
  V.inv <- solve(V)
  #V.inv <- chol2inv(chol(V))
  
  
  #  if(gen.from=="trig"){
  #    xi <- mvrnorm(n = 1,mu= V.inv%*%d,Sigma = V.inv)
  #  } else if(gen.from=="GPP"){
  xi <- as.vector(mvrnorm(n = 1,mu= V.inv%*%d,Sigma = V.inv))
  xi <- xi - mean(xi)
  xi.knots.mcmc[t, ] <- xi
  #  } 
  
  cC_inv_xi <- cC.inv%*%xi
  xi.mcmc[t,] <- cC_inv_xi
  p.mcmc[t,] <- as.vector(exp(X%*%b[t,] + cC_inv_xi)/(1+exp(X%*%b[t,] + cC_inv_xi)))
  

  # Compute xi.grid and p.grid for each year
  for (yr in total_years) {
    dist.cs.grid <- dist_cs_grid_list[[as.character(yr)]]
    dist.ct.grid <- dist_ct_grid_list[[as.character(yr)]]
    c.mat.grid <- sigma.sq * exp(-dist.cs.grid / (2 * ell_s) - dist.ct.grid / (2 * ell_t))
    
    xi.grid <- c.mat.grid %*% C.inv %*% xi
    xi.grid.mcmc_list[[as.character(yr)]][t,] <- xi.grid
    
    eta.grid <- b[t,1] + xi.grid
    exp.eta.grid <- exp(eta.grid)
    p.grid.mcmc_list[[as.character(yr)]][t,] <- as.vector(exp.eta.grid / (1 + exp.eta.grid))
  }
  
  #Step 7 Sample sigma.sq
  #C.prime.inv <- solve(C.prime)
  C.prime.inv <- solve(C.prime + diag(0.00001,m))
  sigma.sq <- rinvgamma(1,
                        shape = a.sigma + m/2,
                        rate = b.sigma + 1/2*t(xi)%*%C.prime.inv%*%xi)
  sigma.sq.mcmc[t] <- sigma.sq
  
  
  # Step 8: Sample ell_s and ell_t separately (Metropolis-Hastings)
  
  # Step 8a: Sample ell_s
  u_s <- runif(1, 0, 1)
  ell_s.current <- ell_s.mcmc[t-1]
  ell_t.current <- ell_t.mcmc[t-1]  
  
  ell_s.log <- rnorm(1, log(ell_s.current), ell_s.sigma)
  ell_s.new <- exp(ell_s.log)
  
  C.prime.current <- exp(-dist.Cs / (2 * ell_s.current) - dist.Ct / (2 * ell_t.current))
  C.mat.current <- sigma.sq * C.prime.current + diag(0.00001, m)
  C.prime.new <- exp(-dist.Cs / (2 * ell_s.new) - dist.Ct / (2 * ell_t.current))
  C.mat.new <- sigma.sq * C.prime.new + diag(0.00001, m)
  c.mat.new <- sigma.sq * exp(-dist.cs / (2 * ell_s.new) - dist.ct / (2 * ell_t.current))
  
  eta <- X %*% b[t,] + xi.mcmc[t,]
  
  # Log posterior for current ell_s
  log.p.current <- -0.5 * (t(h) %*% omg %*% h - 2 * t(h) %*% omg %*% eta + t(eta) %*% omg %*% eta) -
    0.5 * determinant(C.mat.current, logarithm = TRUE)$modulus[1] -
    0.5 * t(xi) %*% chol2inv(chol(C.mat.current)) %*% xi -
    a.l.s * log(ell_s.current) - b.l.s / ell_s.current
  log.p.current <- as.numeric(log.p.current)
  
  # Log posterior for proposed ell_s
  log.p.new <- -0.5 * (t(h) %*% omg %*% h - 2 * t(h) %*% omg %*% eta + t(eta) %*% omg %*% eta) -
    0.5 * determinant(C.mat.new, logarithm = TRUE)$modulus[1] -
    0.5 * t(xi) %*% chol2inv(chol(C.mat.new)) %*% xi -
    a.l.s * log(ell_s.new) - b.l.s / ell_s.new
  log.p.new <- as.numeric(log.p.new)
  
  # Log proposal densities
  log.prop.current <- dnorm(log(ell_s.current), mean = log(ell_s.new), sd = ell_s.sigma, log = TRUE)
  log.prop.new <- dnorm(log(ell_s.new), mean = log(ell_s.current), sd = ell_s.sigma, log = TRUE)
  
  # Acceptance probability
  acc_s <- min(1, exp(log.p.new - log.p.current + log.prop.current - log.prop.new))
  
  if (u_s < acc_s) {
    ell_s.mcmc[t] <- ell_s <- ell_s.new
    C.prime <- C.prime.new
    C.mat <- C.mat.new
    c.mat <- c.mat.new
    ell_s.acc.count <- ell_s.acc.count + 1
  } else {
    ell_s.mcmc[t] <- ell_s <- ell_s.current
    C.prime <- C.prime.current
    C.mat <- C.mat.current
    c.mat <- sigma.sq * exp(-dist.cs / (2 * ell_s.current) - dist.ct / (2 * ell_t.current))
  }
  
  # Step 8b: Sample ell_t
  u_t <- runif(1, 0, 1)
  ell_t.current <- ell_t.mcmc[t-1]
  ell_t.log <- rnorm(1, log(ell_t.current), ell_t.sigma)
  ell_t.new <- exp(ell_t.log)
  
  # Compute covariance matrices for current and proposed ell_t
  C.prime.current <- exp(-dist.Cs / (2 * ell_s) - dist.Ct / (2 * ell_t.current))  # Use updated ell_s
  C.mat.current <- sigma.sq * C.prime.current + diag(0.00001, m)
  C.prime.new <- exp(-dist.Cs / (2 * ell_s) - dist.Ct / (2 * ell_t.new))
  C.mat.new <- sigma.sq * C.prime.new + diag(0.00001, m)
  c.mat.new <- sigma.sq * exp(-dist.cs / (2 * ell_s) - dist.ct / (2 * ell_t.new))
  
  # Log posterior for current ell_t
  log.p.current <- -0.5 * (t(h) %*% omg %*% h - 2 * t(h) %*% omg %*% eta + t(eta) %*% omg %*% eta) -
    0.5 * determinant(C.mat.current, logarithm = TRUE)$modulus[1] -
    0.5 * t(xi) %*% chol2inv(chol(C.mat.current)) %*% xi -
    a.l.t * log(ell_t.current) - b.l.t / ell_t.current
  log.p.current <- as.numeric(log.p.current)
  
  # Log posterior for proposed ell_t
  log.p.new <- -0.5 * (t(h) %*% omg %*% h - 2 * t(h) %*% omg %*% eta + t(eta) %*% omg %*% eta) -
    0.5 * determinant(C.mat.new, logarithm = TRUE)$modulus[1] -
    0.5 * t(xi) %*% chol2inv(chol(C.mat.new)) %*% xi -
    a.l.t * log(ell_t.new) - b.l.t / ell_t.new
  log.p.new <- as.numeric(log.p.new)
  
  # Log proposal density
  log.prop.current <- dnorm(log(ell_t.current), mean = log(ell_t.new), sd = ell_t.sigma, log = TRUE)
  log.prop.new <- dnorm(log(ell_t.new), mean = log(ell_t.current), sd = ell_t.sigma, log = TRUE)
  
  # Acceptance probability
  acc_t <- min(1, exp(log.p.new - log.p.current + log.prop.current - log.prop.new))
  
  if (u_t < acc_t) {
    ell_t.mcmc[t] <- ell_t <- ell_t.new
    C.prime <- C.prime.new
    C.mat <- C.mat.new
    c.mat <- c.mat.new
    ell_t.acc.count <- ell_t.acc.count + 1
  } else {
    ell_t.mcmc[t] <- ell_t <- ell_t.current
    C.prime <- C.prime.current
    C.mat <- C.mat.current
    c.mat <- sigma.sq * exp(-dist.cs / (2 * ell_s) - dist.ct / (2 * ell_t.current))
  }
  
  cat("| I:", t, 
      "| l_s :", ell_s.mcmc[t], sprintf("%.2f%%", 100 * ell_s.acc.count / t),
      "| l_t :", ell_t.mcmc[t], sprintf("%.2f%%", 100 * ell_t.acc.count / t), 
      "| sigma :", sigma.sq.mcmc[t],
      "| success :", success_counts[t - 1],
      "\n")
  
} 

#Remove burn-in and compute geweke statistics
burn.in <- 0.5 * n.iter
b.pred <- b[-(1:(burn.in)),]
b.geweke <- rep(NA,ncol(b.pred))
for(i in 1:ncol(b.pred)){
  b.geweke[i] <- unlist(geweke.diag(b.pred[,i]))["z.var1"]
}


ell_s.pred <- ell_s.mcmc[-(1:burn.in)]
ell_t.pred <- ell_t.mcmc[-(1:burn.in)]


sigma.sq.pred <- sigma.sq.mcmc[-(1:burn.in)]

xi.knots.pred <- xi.knots.mcmc[-(1:burn.in), ]

#---------------Prediction for yr 4----------------
total_years <- if (pred.from == "Yes") c(1, 2, 3, 4) else sort(unique(loc$Year))

p.grid.pred_list <- lapply(total_years, function(yr) p.grid.mcmc_list[[as.character(yr)]][-(1:burn.in),])
names(p.grid.pred_list) <- total_years
p.grid.hat_list <- lapply(p.grid.pred_list, function(x) apply(x, 2, mean))
names(p.grid.hat_list) <- total_years

xi.grid.pred_list <- lapply(total_years, function(yr) xi.grid.mcmc_list[[as.character(yr)]][-(1:burn.in),])
names(xi.grid.pred_list) <- total_years
xi.grid.hat_list <- lapply(xi.grid.pred_list, function(x) apply(x, 2, mean))
names(xi.grid.hat_list) <- total_years


# Prediction for year 4 (if pred.from == "Yes")
if (pred.from == "Yes") {
  
  # Loop over post-burn-in iterations
  for (i in 1:burn.in) {
    ell_s_i <- ell_s.pred[i]
    ell_t_i <- ell_t.pred[i]
    sig_i <- sigma.sq.pred[i]
    xi.knots_i <- xi.knots.pred[i, ]
    
    # Covariance at knots
    C.prime_i <- exp(-dist.Cs / (2 * ell_s_i) - dist.Ct / (2 * ell_t_i))
    C.mat_i <- sig_i * C.prime_i + diag(0.00001, m)
    C.inv_i <- chol2inv(chol(C.mat_i))  # Efficient inversion
    
    # Cross-covariance: grid4 to knots
    c.mat.grid_i <- sig_i * exp(- dist_cs_grid_list[[4]]/ (2 * ell_s_i) - dist_ct_grid_list[[4]] / (2 * ell_t_i))
    
    # Predicted xi at grid4
    xi.grid_i <- c.mat.grid_i %*% C.inv_i %*% xi.knots_i
    xi.grid.pred_list[[as.character(4)]][i,] <- xi.grid.mcmc_list[[as.character(4)]][i+burn.in,] <- xi.grid_i
    
    # Predicted prevalence
    eta.grid_i <- b.pred[i, 1] + xi.grid_i
    p.grid_i <- exp(eta.grid_i) / (1 + exp(eta.grid_i))
    p.grid.pred_list[[as.character(4)]][i,] <- p.grid.mcmc_list[[as.character(4)]][i+burn.in,] <- p.grid_i
  }
  
  # Posterior mean, quantiles, SD
  p.grid.hat_list[[as.character(4)]] <- colMeans(p.grid.pred_list[[as.character(4)]], na.rm = TRUE)
  p.grid.low <- apply(p.grid.pred_list[[as.character(4)]], 2, quantile, probs = 0.025, na.rm = TRUE)
  p.grid.upp <- apply(p.grid.pred_list[[as.character(4)]], 2, quantile, probs = 0.975, na.rm = TRUE)
  p.grid.cover <- as.vector(ifelse(p.grid.true_list[["4"]] >= p.grid.low & p.grid.true_list[["4"]]<= p.grid.upp, 1, 0))
  p.grid.sd <- apply(p.grid.pred_list[[as.character(4)]], 2, sd, na.rm = TRUE)
  
  # Compute bias for year 4
  p.grid.diff <- as.numeric(p.grid.true_list[["4"]]) - p.grid.hat_list[[as.character(4)]]
  p.grid.mse <- p.grid.diff^2 + p.grid.sd^2
  
  assign(paste0("p.grid.out_", 4), rbind(p.grid.diff, p.grid.cover, p.grid.hat_list[[as.character(4)]], p.grid.sd, p.grid.mse))
  
}


############################################################################
############################              SAVE!!!!    ###############################
run.end <- Sys.time()
run.time <- difftime(run.end, run.start, units='mins')
run.time 

#Trace-plots (full samples)
#Beta
for(i in 1:length(beta)){
  plot(b[,i],type = "l", main = paste0("Beta_",i))
}



Y_tilde.pred <- Y_tilde.mcmc[-(1:(0.5*n.iter)),]


ell_s.geweke <- unlist(geweke.diag(ell_s.pred))["z.var1"]
ell_t.geweke <- unlist(geweke.diag(ell_t.pred))["z.var1"]

sigma.sq.geweke <- unlist(geweke.diag(sigma.sq.pred))["z.var1"]  

xi.pred <- xi.mcmc[-(1:(0.5*n.iter)),]
xi.geweke <- rep(NA,ncol(xi.pred))

for(i in 1:ncol(xi.pred)){
  xi.geweke[i] <- unlist(geweke.diag(xi.pred[,i]))["z.var1"]
}

xi.grid.geweke <- vector("list", length(total_years))
names(xi.grid.geweke) <- total_years


for (yr in total_years) {
  xi.grid.geweke[[yr]] <- rep(NA, ncol(xi.grid.pred_list[[yr]])) 
  for (i in 1:ncol(xi.grid.pred_list[[yr]])) {
    xi.grid.geweke[[yr]][i] <- unlist(geweke.diag(xi.grid.pred_list[[yr]][,i]))["z.var1"]
  }
  xi.grid.geweke[[yr]] <- xi.grid.geweke[[yr]][!is.na(xi.grid.geweke[[yr]])]
}

#p.grid
p.grid.geweke <- vector("list", length(total_years))
names(p.grid.geweke) <- total_years


for (yr in total_years) {
  p.grid.geweke[[yr]] <- rep(NA,ncol(p.grid.pred_list[[yr]]))
  for(i in 1:ncol(p.grid.pred_list[[yr]])){
    p.grid.geweke[[yr]][i] <- unlist(geweke.diag(p.grid.pred_list[[yr]][,i]))["z.var1"]
  }
}


#b MCMC samples 
b.hat <- apply(b.pred,2,mean)
b.upp <- apply(b.pred,2,quantile,probs = 0.975) 
b.low <- apply(b.pred,2,quantile,probs = 0.025) 
b.sd <- apply(b.pred,2,sd)


b.diff <- beta - b.hat
b.cover <- ifelse(beta>=b.low&beta<=b.upp,1,0)
b.mse <- b.diff^2+b.sd^2
b.out <- round(rbind(beta,b.hat,b.diff,b.upp,b.low,b.sd,b.mse),4)
b.out

#Y_tilde MCMC samples
Y_tilde.hat <- apply(Y_tilde.pred,2,mean)
Y_tilde.upp <- apply(Y_tilde.pred,2,quantile,probs = 0.975)
Y_tilde.low <- apply(Y_tilde.pred,2,quantile,probs = 0.025)

Y_tilde.diff <- Y_tilde.true - Y_tilde.hat
Y_tilde.cover <- ifelse(Y_tilde.true>=Y_tilde.low&Y_tilde.true<=Y_tilde.upp,1,0)
Y_tilde.out <- rbind(Y_tilde.true,Y_tilde.hat, Y_tilde.diff,Y_tilde.cover,Y_tilde.true)
#Y_tilde.out

#plot(Y_tilde.true-colMeans(Y_tilde.pred),main="Residual plot_Y_tilde") 

#Sensitivity MCMC samples
#sem.hat <- apply(sem.pred,2,mean)
#sem.upp <- apply(sem.pred,2,quantile,probs=0.975)
#sem.low <- apply(sem.pred,2,quantile,probs=0.025)
#sem.sd <- apply(sem.pred,2,sd)

#sem.diff <- df.assay$sens - sem.hat
#sem.cover <- ifelse(df.assay$sens[unique(assay.pool.list)]>=sem.low&df.assay$sens[unique(assay.pool.list)]<=sem.upp,1,0)
#sem.cover <- ifelse(df.assay$sens>=sem.low&df.assay$sens<=sem.upp,1,0)
#sem.mse <- sem.diff^2 + sem.sd^2
#sem.out <- rbind(sem.hat,sem.diff,sem.cover,sem.sd,sem.mse)
#sem.out

#Specificity MCMC samples
#spm.hat <- apply(spm.pred,2,mean)
#spm.upp <- apply(spm.pred,2,quantile,probs=0.975)
#spm.low <- apply(spm.pred,2,quantile,probs=0.025)
#spm.sd <- apply(spm.pred,2,sd)

#spm.diff <- df.assay$spec - spm.hat
#spm.cover <- ifelse(df.assay$spec[unique(assay.pool.list)]>=spm.low&df.assay$spec[unique(assay.pool.list)]<=spm.upp,1,0)
#spm.cover <- ifelse(df.assay$spec>=spm.low&df.assay$spec<=spm.upp,1,0)
#spm.mse <- spm.diff^2 + spm.sd^2
#spm.out <- rbind(spm.hat,spm.diff,spm.cover,spm.sd,spm.mse)
#spm.out


#ell MCMC samples
ell_s.hat <- mean(ell_s.pred)
ell_s.upp <- quantile(ell_s.pred,0.975)
ell_s.low <- quantile(ell_s.pred,0.025)
ell_s.sd <- sd(ell_s.pred)
summary(ell_s.pred)


ell_t.hat <- mean(ell_t.pred)
ell_t.upp <- quantile(ell_t.pred,0.975)
ell_t.low <- quantile(ell_t.pred,0.025)
ell_t.sd <- sd(ell_t.pred)
summary(ell_t.pred)


ell_s.diff <- ell_s.true - ell_s.hat
ell_s.cover <- ifelse(ell_s.true>=ell_s.low&ell_s.true<=ell_s.upp,1,0)
ell_s.mse <- ell_s.diff^2 + ell_s.sd^2
ell_s.out <- rbind(ell_s.true, ell_s.hat, ell_s.diff, ell_s.cover, ell_s.sd, ell_s.mse)
ell_s.out

#ell_s plot
plot(ell_s.mcmc,main="ell_s.sq",type = "l") 


ell_t.diff <- ell_t.true - ell_t.hat
ell_t.cover <- ifelse(ell_t.true>=ell_t.low&ell_t.true<=ell_t.upp,1,0)
ell_t.mse <- ell_t.diff^2 + ell_t.sd^2
ell_t.out <- rbind(ell_t.true, ell_t.hat, ell_t.diff,ell_t.cover, ell_t.sd, ell_t.mse)
ell_t.out

#ell_t plot
plot(ell_t.mcmc,main="ell_t.sq",type = "l") 

#sigma.sq MCMC samples
sigma.sq.hat <- mean(sigma.sq.pred)
sigma.sq.upp <- quantile(sigma.sq.pred,0.975)
sigma.sq.low <- quantile(sigma.sq.pred,0.025)
sigma.sq.sd <- sd(sigma.sq.pred)

sigma.sq.diff <- sigma.sq.true - sigma.sq.hat
sigma.sq.cover <- ifelse(sigma.sq.true>=sigma.sq.low&sigma.sq.true<=sigma.sq.upp,1,0)
sigma.sq.mse <- sigma.sq.diff^2 + sigma.sq.sd^2
sigma.sq.out <- rbind(sigma.sq.hat, sigma.sq.true, sigma.sq.diff, sigma.sq.cover, sigma.sq.sd, sigma.sq.mse)
sigma.sq.out

#ell_sigma plot
plot(sigma.sq.mcmc,main="sigma.sq.sq",type = "l") 

#xi MCMC samples
dist.c <- (st_distance(loc,knot)/ scale)^2
dist.c <- as.matrix(drop_units(dist.c))

c.mat <- sigma.sq.hat*exp(-dist.cs/(2*ell_s)-dist.ct/(2*ell_t)) 

dist.C <- (st_distance(knot,knot) /scale)^2
dist.C <- as.matrix(drop_units(dist.C))

C.prime <- exp(-dist.Cs/(2*ell_s.hat)-dist.Ct/(2*ell_t.hat)) 
C.mat <- sigma.sq.hat*C.prime  + diag(0.00001,m)#with a nugget


# xi
xi.hat <- as.vector(apply(xi.pred,2,mean))
xi.upp <- as.vector(apply(xi.pred,2,quantile,probs=0.975))
xi.low <- as.vector(apply(xi.pred,2,quantile,probs=0.025))
xi.sd <- as.vector(apply(xi.pred,2,sd))

xi.diff <- as.vector(xi.true - xi.hat)
xi.cover <- as.vector(ifelse(xi.true>=xi.low&xi.true<=xi.upp,1,0))
xi.mse <- xi.diff^2 + xi.sd^2
xi.out <- rbind(xi.hat,xi.diff,xi.upp,xi.low,xi.cover,xi.sd,xi.mse)
#xi.out


# xi.grid
for (yr in total_years) {
  xi.grid.hat <- as.vector(apply(xi.grid.pred_list[[yr]], 2, mean))
  xi.grid.upp <- as.vector(apply(xi.grid.pred_list[[yr]], 2, quantile, probs = 0.975))
  xi.grid.low <- as.vector(apply(xi.grid.pred_list[[yr]], 2, quantile, probs = 0.025))
  xi.grid.sd <- as.vector(apply(xi.grid.pred_list[[yr]], 2, sd))
  
  # xi.grid.true
  if (is.list(xi.grid.true)) {
    xi.grid.true_yr <- xi.grid.true[[yr]]
  } else {
    xi.grid.true_yr <- xi.grid.true  
  }
  
  # Compute differences and coverage
  xi.grid.diff <- as.vector(xi.grid.true_yr - xi.grid.hat)
  xi.grid.mse <- xi.grid.diff^2 + xi.grid.sd^2
  xi.grid.cover <- as.vector(ifelse(xi.grid.true_yr >= xi.grid.low & xi.grid.true_yr <= xi.grid.upp, 1, 0))
  
  assign(paste0("xi.grid.out_", yr), rbind(xi.grid.diff, xi.grid.cover, xi.grid.hat, xi.grid.sd, xi.grid.mse))
}

#p
p.pred <- p.mcmc[-(1:(0.5*n.iter)),]
p.geweke <- rep(NA,N)
for(i in 1:N){
  p.geweke[i] <- unlist(geweke.diag(p.pred[,i]))["z.var1"]
}

p.hat <-  as.vector(apply(p.pred,2,mean, na.rm=TRUE))
p.upp <- as.vector(apply(p.pred,2,quantile,probs=0.975, na.rm=TRUE))
p.low <- as.vector(apply(p.pred,2,quantile,probs=0.025, na.rm=TRUE))
p.sd <- as.vector(apply(p.pred,2,sd))

p.diff <- as.vector(p.true - p.hat)
p.cover <- as.vector(ifelse(p.true>=p.low&p.true<=p.upp,1,0))
p.out <- rbind(p.diff,p.cover,p.hat, p.sd)


for (yr in total_years) {
  
  p.grid.upp <- as.vector(apply(p.grid.pred_list[[yr]], 2, quantile, probs = 0.975))
  p.grid.low <- as.vector(apply(p.grid.pred_list[[yr]], 2, quantile, probs = 0.025))
  p.grid.sd <- as.vector(apply(p.grid.pred_list[[yr]], 2, sd))
  
  
  if (is.list(p.grid.true)) {
    p.grid.true_yr <- p.grid.true[[yr]]
  } else {
    p.grid.true_yr <- p.grid.true
  }
  
  p.grid.diff <- as.vector(p.grid.true_yr - p.grid.hat_list[[yr]])
  p.grid.mse <- p.grid.diff^2 + p.grid.sd^2
  p.grid.cover <- as.vector(ifelse(p.grid.true_yr >= p.grid.low[[yr]] & p.grid.true_yr <= p.grid.upp[[yr]], 1, 0))
  
  assign(paste0("p.grid.out_", yr), rbind(p.grid.diff, p.grid.cover, p.grid.hat_list[[yr]], p.grid.sd[[yr]],p.grid.mse))
}


max.lat = max(st_coordinates(shp.st)[,2])
min.lat = min(st_coordinates(shp.st)[,2])

max.long = max(st_coordinates(shp.st)[,1])
min.long = min(st_coordinates(shp.st)[,1])

coords <- st_coordinates(grid.sf)

grid.sf$long <- coords[, "X"]  # Longitude in degrees
grid.sf$lat <- coords[, "Y"]   # Latitude in degrees

min.long <- min(st_coordinates(shp.sf)[, "X"]) - 1
max.long <- max(st_coordinates(shp.sf)[, "X"]) + 1

min.lat <- min(st_coordinates(shp.sf)[, "Y"]) - 1
max.lat <- max(st_coordinates(shp.sf)[, "Y"]) + 1

# choropleth map

plot_list <- list()
summary_list <- list()

for (yr in total_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  grid.sf$p.grid <- p.grid.hat_list[[as.character(yr)]]
  summary_list[[as.character(yr)]] <- summary(grid.sf$p.grid)
  coords <- st_coordinates(grid.sf)
  grid.sf$long <- coords[, "X"]
  grid.sf$lat <- coords[, "Y"]
  
  min.long <- min(st_coordinates(shp.sf)[, "X"]) - 1000
  max.long <- max(st_coordinates(shp.sf)[, "X"]) + 1000
  min.lat <- min(st_coordinates(shp.sf)[, "Y"]) - 1000
  max.lat <- max(st_coordinates(shp.sf)[, "Y"]) + 1000
  
  plot.prev <- ggplot() +
    geom_sf(data = grid.sf, aes(colour = p.grid), shape = 15, size = 3) +
    scale_colour_gradient(low = "white", high = "darkred", limits = c(0, 0.5), name = "Prevalence") +
    coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
    geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.2) +
    theme_bw() +
    ggtitle(paste("Predicted Prevalence in Year", yr))
  
  plot_list[[as.character(yr)]] <- plot.prev
}
plot.prev.true_list

plot_list


write.table(b.out,"b.out.txt",row.name=F,col.name=F)         
write.table(Y_tilde.out,"Y.tilde.out.txt",row.name=F,col.name=F)   

#write.table(sem.out,"sem.out.txt",row.name=F,col.name=F)
#write.table(spm.out,"spm.out.txt",row.name=F,col.name=F)


write.table(xi.out,"xi.out.txt",row.name=F,col.name=F) 
write.table(xi.grid.out_1,"xi.grid.out_1.txt",row.name=F,col.name=F) 
write.table(xi.grid.out_2,"xi.grid.out_2.txt",row.name=F,col.name=F) 
write.table(xi.grid.out_3,"xi.grid.out_3.txt",row.name=F,col.name=F) 
write.table(xi.grid.out_4,"xi.grid.out_4.txt",row.name=F,col.name=F)


write.table(p.grid.out_1,"p.grid.out_1.txt",row.name=F,col.name=F) 
write.table(p.grid.out_2,"p.grid.out_2.txt",row.name=F,col.name=F) 
write.table(p.grid.out_3,"p.grid.out_3.txt",row.name=F,col.name=F) 
write.table(p.grid.out_4,"p.grid.out_4.txt",row.name=F,col.name=F) 

write.table(ell_s.out,"ell_s.txt",row.name=F,col.name = F)
write.table(ell_t.out,"ell_t.txt",row.name=F,col.name = F)
write.table(sigma.sq.out,"sigma.sq.out.txt",row.name=F,col.name = F)

write.table(lat,"lat.txt",row.name=F,col.name=F)
write.table(long,"long.txt",row.name=F,col.name=F)
write.table(assay.pool.list,"assay.pool.list.txt",row.name=F,col.name=F)

write.table(b.geweke,"b.geweke.txt",row.name=F,col.name=F)
write.table(xi.geweke,"xi.geweke.txt",row.name=F,col.name=F) 

#for (i in 1:n.assay) {
#  write.table(sem.geweke[i],paste0("sem.geweke_", i, ".txt"), row.name=F,col.name=F)
#  write.table(spm.geweke[i],paste0("spm.geweke_", i, ".txt"), row.name=F,col.name=F)
#}


for (yr in total_years) {
  write.table(xi.grid.geweke[[as.character(yr)]],
              paste0("xi.grid.geweke_", yr, ".txt"), row.name=F,col.name=F)
  #write.table(p.grid.geweke[[as.character(yr)]],
  #paste0("p.grid.geweke_", yr, ".txt"), row.name=F,col.name=F)
  
}

write.table(p.out,"p.out.txt",row.name=F,col.name = F)
write.table(p.geweke,"p.geweke.txt",row.name=F,col.name = F)

write.table(ell_s.geweke,"ell_s.geweke.txt",row.name=F,col.name = F)
write.table(ell_t.geweke,"ell_t.geweke.txt",row.name=F,col.name = F)
write.table(sigma.sq.geweke,"sigma.sq.geweke.txt",row.name=F,col.name=F) 


#write.table(sens.pool.list, "sens.pool.list.txt", row.names = FALSE, col.names = FALSE)
#write.table(spec.pool.list, "spec.pool.list.txt", row.names = FALSE, col.names = FALSE)
# True Map
#plot.prev.true_list

# Estimated Map
#plot_list
