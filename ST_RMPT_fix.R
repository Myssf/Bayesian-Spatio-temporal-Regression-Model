
if (getwd() == "E:/R/Group/Code/Simulation") {
  path1 = "E:/R/Group/Code/Simulation/"
} else {
  path1 = "/work/chaom/ST_Sim/"    # cluster
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
# ID
sim.id <- 11
###########################
# Iteration
n.iter <- 5000
###########################
# Load data
df.sim  <- read.csv(paste0(path1, "sim.setup.csv"), fileEncoding = "UTF-8-BOM")
sim.ind <- which(df.sim$sim.id == sim.id)
fips    <- read.table(paste0(path1, "County FIPS Codes.txt"), sep = "\t", header = TRUE)
fips.sc <- fips$FIPS[fips$State == "SC"]

###########################
# 1. Sample size
N <- df.sim$N[sim.ind]
################################
# 2. Generating surface
gen.from <- df.sim$surface[sim.ind]
################################
# 3. Prevalence (target) + TRUE betas (8-length: real-data-like covariates)
p.sim <- df.sim$p[sim.ind]
# Order of covariate effects (intercept is auto-tuned below):
#   Intercept, TypeLarvae, TypeMale, TypeNymph,
#   NLCDForest, NLCDShrubland, NLCDWetlands, elevation
if (gen.from == "GPP") {
  if (p.sim == 0.09) {
    beta.true <- c(-2.5, -1.2, 0.2, 1)                       # legacy 4-cov
  } else if (p.sim == 0.18) {
    beta.true <- c(-1.7, -1.0, 0.2, 1)                       # legacy 4-cov
  } else {
    beta.true <- c(0, -0.4, 0.3, -0.5, 0.5, 0.3, 0.6, 0.5) 
  }
}
#####################################
# 4. Pool size / number of pools
#    Task: 700 pools, VARIABLE sizes summing to N (= 4500).
n.pools <- 700
generate_pool_sizes_table5 <- function(N, n.pools,
                                       outlier_sizes = c(60, 86, 553, 875, 1411),
                                       small_probs = c(288, 51, 198, 9, 120) /
                                         sum(c(288, 51, 198, 9, 120))) {
  n.outliers <- length(outlier_sizes)
  n.small.pools <- n.pools - n.outliers
  N.remaining <- N - sum(outlier_sizes)
  stopifnot(N.remaining >= n.small.pools)  # need at least 1 per small pool
  
  small_sizes <- seq_along(small_probs)  # 1:5
  sizes <- sample(small_sizes, size = n.small.pools, replace = TRUE, prob = small_probs)
  
  # nudge to hit exact sum
  diff <- N.remaining - sum(sizes)
  while (diff != 0) {
    idx <- sample(n.small.pools, 1)
    if (diff > 0 && sizes[idx] < max(small_sizes)) { sizes[idx] <- sizes[idx] + 1L; diff <- diff - 1L }
    else if (diff < 0 && sizes[idx] > 1)           { sizes[idx] <- sizes[idx] - 1L; diff <- diff + 1L }
  }
  
  n.size.vec <- sample(c(sizes, outlier_sizes))  # shuffle pool order
  stopifnot(sum(n.size.vec) == N, length(n.size.vec) == n.pools)
  n.size.vec
}

N <- 5000
n.pools <- 700
n.size.vec <- generate_pool_sizes_table5(N, n.pools)

stopifnot(sum(n.size.vec) == N, all(n.size.vec >= 1))
cat("Pool size summary:\n"); print(summary(n.size.vec))
cat("Pool size table:\n"); print(table(n.size.vec))
#####################################
# 5. Number of knots m  (make divisible by T = 4)
m.raw <- df.sim$m[sim.ind]
#####################################
# 6. Sensitivity and specificity
sej.seq <- as.numeric(unlist(lapply(as.character(df.sim$sej), FUN = strsplit, ",")[sim.ind]))
spj.seq <- as.numeric(unlist(lapply(as.character(df.sim$spj), FUN = strsplit, ",")[sim.ind]))
#####################################
# 7. Pool structure
G.mat.config <- "random"
#####################################
# 9. Prediction
pred.from <- df.sim$pred[sim.ind]

########################### Data Generation ################################
###### Covariate design matrix (real-data marginals, independent draws) ####
n.b <- length(beta.true)   # = 8 for the new scenarios


X <- matrix(NA, nrow = N, ncol = n.b)
X[, 1] <- 1
Type.f <- factor(sample(c("T1", "T2", "T3", "T4"), size = N, replace = TRUE,
                        prob = c(1/4, 1/4, 1/4, 1/4)))
X[, 2] <- ifelse(Type.f == "T2", 1, 0)
X[, 3] <- ifelse(Type.f == "T3", 1, 0)
X[, 4] <- ifelse(Type.f == "T4", 1, 0)
NLCD.f <- factor(sample(c("N1", "N2", "N3", "N4"), size = N, replace = TRUE,
                        prob = c(1/4, 1/4, 1/4, 1/4)))
X[, 5] <- ifelse(NLCD.f == "N2", 1, 0)
X[, 6] <- ifelse(NLCD.f == "N3", 1, 0)
X[, 7] <- ifelse(NLCD.f == "N4", 1, 0)
X[, 8] <- rnorm(N, mean = 0, sd = 1)
colnames(X) <- c("Intercept", "T2", "T3", "T4", "N2", "N3", "N4", "Cont")

write.table(X, "X.out.txt", row.name=F, col.name=F)

beta.true[1] <- qlogis(p.sim) - mean(X[, -1] %*% beta.true[-1])
beta <- beta.true

# Assays
n.assay <- df.sim$n.assay[sim.ind]
sens <- sej.seq[1:n.assay]
spec <- spj.seq[1:n.assay]
df.assay <- data.frame(id = 1:n.assay, sens = sens, spec = spec)

######################################### Spatial random effects ###########
n.county <- length(fips.sc)
fips.sub <- fips.sc[1:n.county]

# ---- T = 4 observed years; forecast year 5 ----
T <- 4
years.ind <- sample(1:T, size = N, replace = TRUE, prob = rep(1/T, T))

if (!exists("shp")) {
  shp <- st_read(paste0(path1, "cb_2018_us_county_500k.shp")[1])
  shp <- shp[shp$GEOID %in% fips.sub, ]
}
neighbours <- poly2nb(shp)
W <- Matrix(as(nb2mat(neighbours, style = "B", zero.policy = TRUE), "Matrix"), sparse = TRUE)
D <- Diagonal(n = ncol(W), colSums(W))

# Hyperparameters
a.sigma <- 3; b.sigma <- 3
a.l.s <- 3; b.l.s <- 3
a.l.t <- 3; b.l.t <- 3

shp.st <- shp
shp.sf <- st_as_sf(shp.st)

# True GP hyperparameters + MH inits
ell_s.true <- ell_s <- 3
ell_s.sigma <- 0.3; ell_s.init <- 0.5
ell_s.mcmc <- rep(NA, n.iter); ell_s.mcmc[1] <- ell_s.init
ell_t.true <- ell_t <- 2
ell_t.sigma <- 0.5; ell_t.init <- 0.5
ell_t.mcmc <- rep(NA, n.iter); ell_t.mcmc[1] <- ell_t.init
sigma.sq_SE <- 0.5; sigma.sq.true <- sigma.sq_SE; sigma.sq.init <- 0.2

st_crs(shp.sf) <- 4326
loc_geom <- st_sample(shp.sf, size = N)
attr_df  <- data.frame(Year = years.ind[1:N])
loc <- st_sf(attr_df, geometry = loc_geom)
loc <- st_transform(loc, crs = 3086)
loc.coordinates <- st_coordinates(loc)

# ---- Knots: m divisible by T; base knots replicated over the 4 years ----
n_base <- floor(m.raw / T)
m <- n_base * T
id <- sample(1:N, n_base)
base_knot <- loc[id, ]
knot_data <- data.frame()
for (t in 1:T) {
  temp_knot <- base_knot
  temp_knot$Year <- t
  knot_data <- rbind(knot_data, temp_knot)
}
knot <- st_sf(knot_data)
cat("Knots: base =", n_base, " total m =", m, "\n")

scale <- 111000
dist.cs <- (st_distance(loc, knot) / scale)^2
dist.cs <- as.matrix(drop_units(dist.cs))
dist.ct <- outer(loc$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
c.mat   <- sigma.sq_SE * exp(-dist.cs / (2 * ell_s) - dist.ct / (2 * ell_t))

dist.Cs <- (st_distance(knot, knot) / scale)^2
dist.Cs <- as.matrix(drop_units(dist.Cs))
dist.Ct <- outer(knot$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
C.prime <- exp(-dist.Cs / (2 * ell_s) - dist.Ct / (2 * ell_t))
C.mat   <- sigma.sq_SE * C.prime + diag(0.00001, m)

# True spatial effect + prevalence
xi <- mvrnorm(1, rep(0, m), C.mat)
xi <- xi - mean(xi)
if (gen.from == "GPP") {
  xi <- mvrnorm(1, rep(0, m), C.mat)
  xi <- xi - mean(xi)
  xi.true <- c.mat %*% solve(C.mat) %*% xi
  p <- as.vector(exp(X %*% beta + xi.true) / (1 + exp(X %*% beta + xi.true)))
  p.true <- p
}
cat("mean(p.true) =", round(mean(p.true), 4), " (target", p.sim, ")\n")

#####################################################################
# Prediction grid: observed years 1..T + forecast year 5
#####################################################################
n.grid = 3000
cellsize <- 0.05*scale

if (is.na(st_crs(shp.sf))) shp.sf <- st_set_crs(shp.sf, 4326)

pred_year   <- 5
total_years <- if (pred.from == "Yes") c(1:T, pred_year) else sort(unique(loc$Year))

shp.sf <- st_transform(shp.sf, crs = 3086)
grid_list <- list()
for (yr in total_years) {
  grid <- shp.sf %>%
    st_make_grid(cellsize = cellsize, what = "centers",
                 n = c(n.grid/2, n.grid/2), square = FALSE) %>%
    st_intersection(shp.sf)
  grid.sf <- st_as_sf(grid)
  grid.sf$Year <- yr
  grid_list[[as.character(yr)]] <- grid.sf
}

lat  <- st_coordinates(grid.sf)[, "Y"]
long <- st_coordinates(grid.sf)[, "X"]
write.table(data.frame(lat = lat, long = long, grid_id = 1:length(lat)),
            file = paste0(path1, "grid_coordinates.txt"),
            row.names = FALSE, col.names = TRUE)

# Grid<->knot distances (used for truth AND MCMC / prediction)
dist_cs_grid_list <- list(); dist_ct_grid_list <- list()
for (yr in total_years) {
  g <- grid_list[[as.character(yr)]]
  dcs <- (st_distance(g, knot) / scale)^2
  dist_cs_grid_list[[as.character(yr)]] <- as.matrix(drop_units(dcs))
  dist_ct_grid_list[[as.character(yr)]] <-
    outer(g$Year, knot$Year, FUN = function(x, y) abs(x - y)^2)
}

# TRUE prevalence surfaces (intercept + spatial only; covariates unavailable at grid)
xi.grid.true_list <- list(); p.grid.true_list <- list()
for (yr in total_years) {
  c.mat.grid <- sigma.sq.true *
    exp(-dist_cs_grid_list[[as.character(yr)]] / (2 * ell_s.true) -
          dist_ct_grid_list[[as.character(yr)]] / (2 * ell_t.true))
  xi.grid.true <- c.mat.grid %*% solve(C.mat) %*% xi
  xi.grid.true_list[[as.character(yr)]] <- xi.grid.true
  p.grid.true_list[[as.character(yr)]] <-
    exp(beta[1] + xi.grid.true) / (1 + exp(beta[1] + xi.grid.true))
}

# True-prevalence maps
plot.prev.true_list <- list()
for (yr in total_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  grid.sf$p.grid.true <- p.grid.true_list[[as.character(yr)]]
  min.long <- min(st_coordinates(shp.sf)[, "X"]) - 1000
  max.long <- max(st_coordinates(shp.sf)[, "X"]) + 1000
  min.lat  <- min(st_coordinates(shp.sf)[, "Y"]) - 1000
  max.lat  <- max(st_coordinates(shp.sf)[, "Y"]) + 1000
  plot.prev.true_list[[as.character(yr)]] <- ggplot() +
    geom_sf(data = grid.sf, aes(colour = p.grid.true), shape = 15, size = 2) +
    scale_colour_gradient(low = "white", high = "darkred",
                          limits = c(0, 1), name = "True Prevalence") +
    coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
    geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.3) +
    theme_bw() + ggtitle(paste("True Prevalence in Year", yr))
}

#####################################################################
# True individual status + pools + observed pool tests
#####################################################################
Y_tilde      <- rbinom(N, 1, p)
Y_tilde.true <- Y_tilde

# per-individual assay assignment (sequential)
assay.ind.list <- sens.ind.list <- spec.ind.list <- rep(NA, N)
for (i in 1:n.assay) {
  idxr <- ((i - 1) * (N / n.assay) + 1):(i * N / n.assay)
  assay.ind.list[idxr] <- df.assay$id[i]
  sens.ind.list[idxr]  <- df.assay$sens[i]
  spec.ind.list[idxr]  <- df.assay$spec[i]
}
Y <- rep(NA, N)
for (i in 1:n.assay) {
  Y[Y_tilde == 1 & sens.ind.list == df.assay$sens[i]] <-
    rbinom(sum(Y_tilde == 1 & sens.ind.list == df.assay$sens[i]), 1, df.assay$sens[i])
  Y[Y_tilde == 0 & spec.ind.list == df.assay$spec[i]] <-
    rbinom(sum(Y_tilde == 0 & spec.ind.list == df.assay$spec[i]), 1, 1 - df.assay$spec[i])
}

# ---- Build VARIABLE-size pools (contiguous blocks) ----
pool_list <- vector("list", n.pools)
idx <- 1L
for (i in 1:n.pools) {
  pool_list[[i]] <- idx:(idx + n.size.vec[i] - 1L)
  idx <- idx + n.size.vec[i]
}
stopifnot(idx - 1L == N)

ind_list <- vector("list", N)
for (i in 1:N) ind_list[[i]] <- which(vapply(pool_list, function(p) i %in% p, logical(1)))

# True pool status
Z.tilde <- vapply(pool_list, function(p) as.numeric(sum(Y_tilde[p]) > 0), numeric(1))

# assay per pool
assay.pool.list <- sample(1:n.assay, size = n.pools, replace = TRUE, prob = rep(1/n.assay, n.assay))
sens.pool.list  <- spec.pool.list <- rep(NA, n.pools)
for (i in 1:n.assay) {
  sens.pool.list[assay.pool.list == i] <- df.assay$sens[i]
  spec.pool.list[assay.pool.list == i] <- df.assay$spec[i]
}

# observed pool tests
Z <- rep(NA, n.pools)
for (i in 1:n.assay) {
  Z[assay.pool.list == i] <- rbinom(sum(assay.pool.list == i), 1,
                                    df.assay$sens[i] * Z.tilde[assay.pool.list == i] +
                                      (1 - df.assay$spec[i]) * (1 - Z.tilde[assay.pool.list == i]))
}
Z_0 <- Z
n.tests <- length(Z)
cat("n pools =", length(Z), " | # positive pools =", sum(Z), "\n")

#####################################################################
# Gibbs / MH setup
#####################################################################
n.b       <- length(beta)
b.m       <- rep(0, n.b)
b.cov     <- diag(10, n.b)
b.cov.inv <- solve(b.cov)

b <- matrix(NA, n.iter, n.b); b[1, ] <- rep(0.5, n.b)
w <- matrix(NA, n.iter, N);   w[1, ] <- rep(0.5, N)
Y_tilde.mcmc <- matrix(0, n.iter, N)
Y_tilde <- rep(0, N)

xi.N        <- c.mat %*% solve(C.mat) %*% xi
xi.mcmc     <- matrix(NA, n.iter, N); xi.mcmc[1, ] <- xi.true
p.mcmc      <- matrix(NA, n.iter, N)
p.mcmc[1, ] <- as.vector(exp(X %*% b[1, ] + xi.N) / (1 + exp(X %*% b[1, ] + xi.N)))

p.grid.mcmc_list  <- lapply(total_years, function(yr)
  matrix(NA, nrow = n.iter, ncol = nrow(grid_list[[as.character(yr)]])))
names(p.grid.mcmc_list) <- total_years
xi.grid.mcmc_list <- lapply(total_years, function(yr)
  matrix(NA, nrow = n.iter, ncol = nrow(grid_list[[as.character(yr)]])))
names(xi.grid.mcmc_list) <- total_years

sigma.sq.mcmc <- rep(NA, n.iter); sigma.sq.mcmc[1] <- sigma.sq <- sigma.sq.init
omg <- Diagonal(x = w[1, ])
xi.knots.mcmc <- matrix(NA, n.iter, m)

# Se/Sp FIXED at true values (single-/multi-assay debug)
sej <- sens.pool.list
spj <- spec.pool.list

# During the loop, fill grids for OBSERVED years only; year 5 is a
# dedicated post-burn-in forecast block.
years_obs <- 1:T

success_counts <- rep(0, n.iter)
ell_s.acc.count <- 0
ell_t.acc.count <- 0


###############################            RUN            ###################
run.start <- Sys.time()
for (t in 2:n.iter) {
  
    #Step 1 Sample Y_tilde (Pool testing) 
    C.inv <- chol2inv(chol(C.mat))
    cC.inv <- c.mat %*% C.inv  
    cC_inv_xi <- cC.inv %*% xi
    
    #Call cpp function to sample Y_tilde
    result <- sample_Y_tilde(X, b, cC_inv_xi, Z, sej, spj, ind_list, pool_list, Y_tilde, t, N, success_counts)
    
    Y_tilde <- result$Y_tilde
    success_counts[t - 1] <- result$success_count  # Store success_count for this iteration
    Y_tilde.mcmc[t, ] <- Y_tilde
  
  ## Step 2: Z.tilde
  Z.tilde <- sapply(pool_list, function(p) ifelse(sum(Y_tilde[unlist(p)]) > 0, 1, 0))
  
  ## Step 3: Se/Sp FIXED (no sampling)
  
  ## Step 4: sample beta
  M <- t(X) %*% omg %*% X + b.cov.inv
  k <- Y_tilde - 1/2
  h <- k / diag(omg)
  Q <- t(X) %*% omg %*% (h - cC_inv_xi) + b.cov.inv %*% b.m
  M.inv <- chol2inv(chol(M))
  b[t, ] <- as.vector(mvrnorm(1, M.inv %*% Q, M.inv))
  
  ## Step 5: sample w
  w[t, ] <- apply(X = X %*% b[t, ] + cC_inv_xi, 1, FUN = rpg, num = 1, h = 1)
  omg <- Diagonal(x = w[t, ])
  h   <- k / diag(omg)
  
  ## Step 6: sample xi
  omg.cC.inv <- omg %*% cC.inv
  V <- t(cC.inv) %*% omg.cC.inv + C.inv
  d <- t(t(h - X %*% b[t, ]) %*% omg.cC.inv)
  V.inv <- chol2inv(chol(V + diag(1e-8, m)))
  V.inv <- (V.inv + t(V.inv)) / 2
  xi <- as.vector(mvrnorm(1, mu = V.inv %*% d, Sigma = V.inv))
  xi <- xi - mean(xi)
  xi.knots.mcmc[t, ] <- xi
  
  cC_inv_xi <- cC.inv %*% xi
  xi.mcmc[t, ] <- cC_inv_xi
  p.mcmc[t, ]  <- as.vector(exp(X %*% b[t, ] + cC_inv_xi) /
                              (1 + exp(X %*% b[t, ] + cC_inv_xi)))
  
  ## grid xi + prevalence for OBSERVED years (intercept + spatial)
  for (yr in years_obs) {
    c.mat.grid <- sigma.sq *
      exp(-dist_cs_grid_list[[as.character(yr)]] / (2 * ell_s) -
            dist_ct_grid_list[[as.character(yr)]] / (2 * ell_t))
    xi.grid <- c.mat.grid %*% C.inv %*% xi
    xi.grid.mcmc_list[[as.character(yr)]][t, ] <- xi.grid
    eta.grid <- b[t, 1] + xi.grid
    p.grid.mcmc_list[[as.character(yr)]][t, ] <- as.vector(exp(eta.grid) / (1 + exp(eta.grid)))
  }
  
  ## Step 7: sigma^2
  C.prime.inv <- solve(C.prime + diag(0.00001, m))
  sigma.sq <- rinvgamma(1, shape = a.sigma + m/2,
                        rate = b.sigma + 0.5 * t(xi) %*% C.prime.inv %*% xi)
  sigma.sq.mcmc[t] <- sigma.sq
  
  eta <- X %*% b[t, ] + xi.mcmc[t, ]
  
  ## Step 8a: MH ell_s
  u_s <- runif(1)
  ell_s.current <- ell_s.mcmc[t-1]; ell_t.current <- ell_t.mcmc[t-1]
  ell_s.new <- exp(rnorm(1, log(ell_s.current), ell_s.sigma))
  C.prime.current <- exp(-dist.Cs/(2*ell_s.current) - dist.Ct/(2*ell_t.current))
  C.mat.current   <- sigma.sq * C.prime.current + diag(0.00001, m)
  C.prime.new     <- exp(-dist.Cs/(2*ell_s.new)     - dist.Ct/(2*ell_t.current))
  C.mat.new       <- sigma.sq * C.prime.new + diag(0.00001, m)
  c.mat.new       <- sigma.sq * exp(-dist.cs/(2*ell_s.new) - dist.ct/(2*ell_t.current))
  log.p.current <- -0.5*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) -
    0.5*determinant(C.mat.current, logarithm=TRUE)$modulus[1] -
    0.5*t(xi)%*%chol2inv(chol(C.mat.current))%*%xi -
    a.l.s*log(ell_s.current) - b.l.s/ell_s.current
  log.p.new <- -0.5*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) -
    0.5*determinant(C.mat.new, logarithm=TRUE)$modulus[1] -
    0.5*t(xi)%*%chol2inv(chol(C.mat.new))%*%xi -
    a.l.s*log(ell_s.new) - b.l.s/ell_s.new
  log.prop.current <- dnorm(log(ell_s.current), log(ell_s.new), ell_s.sigma, log=TRUE)
  log.prop.new     <- dnorm(log(ell_s.new), log(ell_s.current), ell_s.sigma, log=TRUE)
  acc_s <- min(1, exp(as.numeric(log.p.new) - as.numeric(log.p.current) + log.prop.current - log.prop.new))
  if (u_s < acc_s) {
    ell_s.mcmc[t] <- ell_s <- ell_s.new
    C.prime <- C.prime.new; C.mat <- C.mat.new; c.mat <- c.mat.new
    ell_s.acc.count <- ell_s.acc.count + 1
  } else {
    ell_s.mcmc[t] <- ell_s <- ell_s.current
    C.prime <- C.prime.current; C.mat <- C.mat.current
    c.mat <- sigma.sq * exp(-dist.cs/(2*ell_s.current) - dist.ct/(2*ell_t.current))
  }
  
  ## Step 8b: MH ell_t
  u_t <- runif(1)
  ell_t.current <- ell_t.mcmc[t-1]
  ell_t.new <- exp(rnorm(1, log(ell_t.current), ell_t.sigma))
  C.prime.current <- exp(-dist.Cs/(2*ell_s) - dist.Ct/(2*ell_t.current))
  C.mat.current   <- sigma.sq * C.prime.current + diag(0.00001, m)
  C.prime.new     <- exp(-dist.Cs/(2*ell_s) - dist.Ct/(2*ell_t.new))
  C.mat.new       <- sigma.sq * C.prime.new + diag(0.00001, m)
  c.mat.new       <- sigma.sq * exp(-dist.cs/(2*ell_s) - dist.ct/(2*ell_t.new))
  log.p.current <- -0.5*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) -
    0.5*determinant(C.mat.current, logarithm=TRUE)$modulus[1] -
    0.5*t(xi)%*%chol2inv(chol(C.mat.current))%*%xi -
    a.l.t*log(ell_t.current) - b.l.t/ell_t.current
  log.p.new <- -0.5*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) -
    0.5*determinant(C.mat.new, logarithm=TRUE)$modulus[1] -
    0.5*t(xi)%*%chol2inv(chol(C.mat.new))%*%xi -
    a.l.t*log(ell_t.new) - b.l.t/ell_t.new
  log.prop.current <- dnorm(log(ell_t.current), log(ell_t.new), ell_t.sigma, log=TRUE)
  log.prop.new     <- dnorm(log(ell_t.new), log(ell_t.current), ell_t.sigma, log=TRUE)
  acc_t <- min(1, exp(as.numeric(log.p.new) - as.numeric(log.p.current) + log.prop.current - log.prop.new))
  if (u_t < acc_t) {
    ell_t.mcmc[t] <- ell_t <- ell_t.new
    C.prime <- C.prime.new; C.mat <- C.mat.new; c.mat <- c.mat.new
    ell_t.acc.count <- ell_t.acc.count + 1
  } else {
    ell_t.mcmc[t] <- ell_t <- ell_t.current
    C.prime <- C.prime.current; C.mat <- C.mat.current
    c.mat <- sigma.sq * exp(-dist.cs/(2*ell_s) - dist.ct/(2*ell_t.current))
  }
  
  cat("| I:", t,
      "| l_s:", round(ell_s.mcmc[t],3), sprintf("%.1f%%",100*ell_s.acc.count/t),
      "| l_t:", round(ell_t.mcmc[t],3), sprintf("%.1f%%",100*ell_t.acc.count/t),
      "| sig:", round(sigma.sq.mcmc[t],3),
      "| succ:", success_counts[t-1], "\n")
}
run.end  <- Sys.time()
run.time <- difftime(run.end, run.start, units = "mins"); run.time

#####################################################################
# Burn-in + posterior summaries
#####################################################################
burn.in <- 0.5 * n.iter
b.pred  <- b[-(1:burn.in), ]
b.geweke <- apply(b.pred, 2, function(col) unlist(geweke.diag(col))["z.var1"])

ell_s.pred <- ell_s.mcmc[-(1:burn.in)]
ell_t.pred <- ell_t.mcmc[-(1:burn.in)]
sigma.sq.pred <- sigma.sq.mcmc[-(1:burn.in)]
xi.knots.pred <- xi.knots.mcmc[-(1:burn.in), ]

# observed-year grid posteriors
p.grid.pred_list  <- lapply(years_obs, function(yr) p.grid.mcmc_list[[as.character(yr)]][-(1:burn.in), ])
names(p.grid.pred_list) <- years_obs
p.grid.hat_list   <- lapply(p.grid.pred_list, function(x) apply(x, 2, mean))
names(p.grid.hat_list) <- years_obs
xi.grid.pred_list <- lapply(years_obs, function(yr) xi.grid.mcmc_list[[as.character(yr)]][-(1:burn.in), ])
names(xi.grid.pred_list) <- years_obs
xi.grid.hat_list  <- lapply(xi.grid.pred_list, function(x) apply(x, 2, mean))
names(xi.grid.hat_list) <- years_obs

#####################################################################
# YEAR-5 FORECAST + BIAS MAP
#####################################################################
if (pred.from == "Yes") {
  yc <- as.character(pred_year)
  n5 <- nrow(grid_list[[yc]])
  p.grid.pred_yr5  <- matrix(NA, burn.in, n5)
  xi.grid.pred_yr5 <- matrix(NA, burn.in, n5)
  
  for (i in 1:burn.in) {
    ell_s_i <- ell_s.pred[i]; ell_t_i <- ell_t.pred[i]; sig_i <- sigma.sq.pred[i]
    xik_i   <- xi.knots.pred[i, ]
    C.prime_i <- exp(-dist.Cs/(2*ell_s_i) - dist.Ct/(2*ell_t_i))
    C.mat_i   <- sig_i * C.prime_i + diag(0.00001, m)
    C.inv_i   <- chol2inv(chol(C.mat_i))
    c.mat.grid_i <- sig_i * exp(-dist_cs_grid_list[[yc]]/(2*ell_s_i) -
                                  dist_ct_grid_list[[yc]]/(2*ell_t_i))
    xi.grid_i <- c.mat.grid_i %*% C.inv_i %*% xik_i
    xi.grid.pred_yr5[i, ] <- xi.grid_i
    eta.grid_i <- b.pred[i, 1] + xi.grid_i
    p.grid.pred_yr5[i, ] <- as.vector(exp(eta.grid_i) / (1 + exp(eta.grid_i)))
  }
  # store into the master lists too
  xi.grid.pred_list[[yc]] <- xi.grid.pred_yr5
  p.grid.pred_list[[yc]]  <- p.grid.pred_yr5
  p.grid.hat_list[[yc]]   <- colMeans(p.grid.pred_yr5, na.rm = TRUE)
  xi.grid.hat_list[[yc]]  <- colMeans(xi.grid.pred_yr5, na.rm = TRUE)
  
  p.grid.low  <- apply(p.grid.pred_yr5, 2, quantile, probs = 0.025, na.rm = TRUE)
  p.grid.upp  <- apply(p.grid.pred_yr5, 2, quantile, probs = 0.975, na.rm = TRUE)
  p.grid.sd   <- apply(p.grid.pred_yr5, 2, sd, na.rm = TRUE)
  p.grid.true5 <- as.numeric(p.grid.true_list[[yc]])
  p.grid.diff <- p.grid.true5 - p.grid.hat_list[[yc]]        # BIAS
  p.grid.cover <- as.vector(ifelse(p.grid.true5 >= p.grid.low & p.grid.true5 <= p.grid.upp, 1, 0))
  p.grid.mse  <- p.grid.diff^2 + p.grid.sd^2
  assign(paste0("p.grid.out_", pred_year),
         rbind(p.grid.diff, p.grid.cover, p.grid.hat_list[[yc]], p.grid.sd, p.grid.mse))
  
  cat("\n--- Year-5 forecast ---\n")
  cat("mean |bias| :", round(mean(abs(p.grid.diff), na.rm=TRUE), 4), "\n")
  cat("coverage    :", round(mean(p.grid.cover, na.rm=TRUE), 3), "\n")
  cat("mean MSE    :", round(mean(p.grid.mse, na.rm=TRUE), 5), "\n")
}

# All years for downstream summaries
all_years <- total_years

#####################################################################
# Beta / hyperparameter recovery
#####################################################################
b.hat <- apply(b.pred, 2, mean)
b.upp <- apply(b.pred, 2, quantile, probs = 0.975)
b.low <- apply(b.pred, 2, quantile, probs = 0.025)
b.sd  <- apply(b.pred, 2, sd)
b.diff <- beta - b.hat
b.cover <- ifelse(beta >= b.low & beta <= b.upp, 1, 0)
b.mse <- b.diff^2 + b.sd^2
b.out <- round(rbind(beta, b.hat, b.diff, b.upp, b.low, b.sd, b.mse), 4)
if (n.b == 8) colnames(b.out) <- c("Intercept","Larvae","Male","Nymph",
                                   "Forest","Shrubland","Wetlands","Elevation")
cat("\n--- Beta recovery ---\n"); print(b.out)


#ell MCMC samples
ell_s.hat <- mean(ell_s.pred)
ell_s.upp <- quantile(ell_s.pred,0.975)
ell_s.low <- quantile(ell_s.pred,0.025)
ell_s.sd <- sd(ell_s.pred)


ell_t.hat <- mean(ell_t.pred)
ell_t.upp <- quantile(ell_t.pred,0.975)
ell_t.low <- quantile(ell_t.pred,0.025)
ell_t.sd <- sd(ell_t.pred)


ell_s.diff <- ell_s.true - ell_s.hat
ell_s.cover <- ifelse(ell_s.true>=ell_s.low&ell_s.true<=ell_s.upp,1,0)
ell_s.mse <- ell_s.diff^2 + ell_s.sd^2
ell_s.out <- rbind(ell_s.true, ell_s.hat, ell_s.diff, ell_s.cover, ell_s.sd, ell_s.mse)


ell_t.diff <- ell_t.true - ell_t.hat
ell_t.cover <- ifelse(ell_t.true>=ell_t.low&ell_t.true<=ell_t.upp,1,0)
ell_t.mse <- ell_t.diff^2 + ell_t.sd^2
ell_t.out <- rbind(ell_t.true, ell_t.hat, ell_t.diff,ell_t.cover, ell_t.sd, ell_t.mse)



#sigma.sq MCMC samples
sigma.sq.hat <- mean(sigma.sq.pred)
sigma.sq.upp <- quantile(sigma.sq.pred,0.975)
sigma.sq.low <- quantile(sigma.sq.pred,0.025)
sigma.sq.sd <- sd(sigma.sq.pred)

sigma.sq.diff <- sigma.sq.true - sigma.sq.hat
sigma.sq.cover <- ifelse(sigma.sq.true>=sigma.sq.low&sigma.sq.true<=sigma.sq.upp,1,0)
sigma.sq.mse <- sigma.sq.diff^2 + sigma.sq.sd^2
sigma.sq.out <- rbind(sigma.sq.hat, sigma.sq.true, sigma.sq.diff, sigma.sq.cover, sigma.sq.sd, sigma.sq.mse)

#####################################################################
# Bias / coverage / MSE tables for ALL years (obs + forecast)
#####################################################################

for (yr in all_years) {
  yc <- as.character(yr)
  pp <- p.grid.pred_list[[yc]]
  p.low <- apply(pp, 2, quantile, probs = 0.025, na.rm = TRUE)
  p.upp <- apply(pp, 2, quantile, probs = 0.975, na.rm = TRUE)
  p.sd  <- apply(pp, 2, sd, na.rm = TRUE)
  p.tru <- as.numeric(p.grid.true_list[[yc]])
  p.dif <- p.tru - p.grid.hat_list[[yc]]
  p.cov <- as.vector(ifelse(p.tru >= p.low & p.tru <= p.upp, 1, 0))
  p.mse <- p.dif^2 + p.sd^2
  assign(paste0("p.grid.out_", yr),
         rbind(p.dif, p.cov, p.grid.hat_list[[yc]], p.sd, p.mse))
  
  xp <- xi.grid.pred_list[[yc]]
  x.hat <- apply(xp, 2, mean); x.sd <- apply(xp, 2, sd)
  x.low <- apply(xp, 2, quantile, probs = 0.025)
  x.upp <- apply(xp, 2, quantile, probs = 0.975)
  x.tru <- as.numeric(xi.grid.true_list[[yc]])
  x.dif <- x.tru - x.hat
  x.cov <- as.vector(ifelse(x.tru >= x.low & x.tru <= x.upp, 1, 0))
  x.mse <- x.dif^2 + x.sd^2
  assign(paste0("xi.grid.out_", yr),
         rbind(x.dif, x.cov, x.hat, x.sd, x.mse))
}

#####################################################################
# Maps: estimated (obs years), forecast (yr5), bias (yr5)
#####################################################################
min.long <- min(st_coordinates(shp.sf)[, "X"]) - 1000
max.long <- max(st_coordinates(shp.sf)[, "X"]) + 1000
min.lat  <- min(st_coordinates(shp.sf)[, "Y"]) - 1000
max.lat  <- max(st_coordinates(shp.sf)[, "Y"]) + 1000

plot_list <- list()
for (yr in all_years) {
  grid.sf <- grid_list[[as.character(yr)]]
  grid.sf$p.grid <- p.grid.hat_list[[as.character(yr)]]
  ttl <- if (yr == pred_year) paste("Forecast Prevalence in Year", yr) else paste("Predicted Prevalence in Year", yr)
  plot_list[[as.character(yr)]] <- ggplot() +
    geom_sf(data = grid.sf, aes(colour = p.grid), shape = 15, size = 3) +
    scale_colour_gradient(low = "white", high = "darkred", limits = c(0, 1), name = "Prevalence") +
    coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
    geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.2) +
    theme_bw() + ggtitle(ttl)
}

# Year-5 bias map
g5 <- grid_list[[as.character(pred_year)]]
g5$p.grid.bias <- as.numeric(p.grid.true_list[[as.character(pred_year)]]) - p.grid.hat_list[[as.character(pred_year)]]
bias_range <- max(abs(range(g5$p.grid.bias, na.rm = TRUE)))
plot.bias_yr5 <- ggplot() +
  geom_sf(data = g5, aes(colour = p.grid.bias), shape = 15, size = 3) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-bias_range, bias_range), name = "Bias (True - Pred)") +
  coord_sf(xlim = c(min.long, max.long), ylim = c(min.lat, max.lat), expand = FALSE) +
  geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.3) +
  theme_bw() + ggtitle(paste("Bias in Forecast Prevalence for Year", pred_year))

#plot.prev.true_list
#plot_list
#plot.bias_yr5

#####################################################################
# SAVE
#####################################################################
write.table(b.out, "b.out.txt", row.name=F, col.name=F)
for (yr in all_years) {
  write.table(get(paste0("p.grid.out_", yr)),  paste0("p.grid.out_", yr, ".txt"),  row.name=F, col.name=F)
  write.table(get(paste0("xi.grid.out_", yr)), paste0("xi.grid.out_", yr, ".txt"), row.name=F, col.name=F)
}
write.table(ell_s.out,    "ell_s.txt",        row.name=F, col.name=F)
write.table(ell_t.out,    "ell_t.txt",        row.name=F, col.name=F)
write.table(sigma.sq.out, "sigma.sq.out.txt", row.name=F, col.name=F)
#write.table(b.geweke,     paste0(path1,"b.geweke.txt"),     row.name=F, col.name=F)
#write.table(lat,  paste0(path1,"lat.txt"),  row.name=F, col.name=F)
#write.table(long, paste0(path1,"long.txt"), row.name=F, col.name=F)
#ggsave(paste0(path1,"sim_yr5_bias.png"), plot.bias_yr5, width = 8, height = 6, dpi = 300)

