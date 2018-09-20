# This script simulates a dataset of 1 vessel and 3 bacterial families,
# with a few missing days.
library(compositions)
library(philr)
library(rstan)
library(tidyverse)
library(StanStateSpace)
library(ape)
library(abind)
library(ggtree)

set.seed(4)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


machine <- "local"
paths <- list()
paths[["local"]] <- list("results" = "~/Research/mdlm/results/2018-04-18_code_github/mdlm_toblerone/",
                        "mapping" = '~/Research/data/_data_raw/sequencing.2016.03.04/2016.03.25MappingFile.MergedPool.txt',
                        "phyloseq"= "~/Research/data/_data_derived/sequencing.2016.03.04/dada2_2016_10_18/hardac/phyloseq.rds")
paths[["hardac"]] <- list("results" = "/data/davidlab/users/jds/mdlm_publication_2/mdlm_toblerone/",
                          "mapping" = '/data/davidlab/users/jds/mdlm_publication_2/dada2/0_mapping/2016.03.25MappingFile.MergedPool.txt',
                          "phyloseq" = "/data/davidlab/users/jds/mdlm_publication_2/dada2/phyloseq.rds")

setwd(paths[[machine]]$results)
source("utils.R")

out.dir <- paste0(paths[[machine]]$results, "simulation/")
if (!dir.exists(out.dir)) dir.create(out.dir)
out <- function(filename){
  paste0(out.dir, filename)
}

# Simulation: Predefined Parameters ---------------------------------------

q <- 2 # number of ILR Coordinates (+1 = number of species)
r <- 1 # number of series
# T.length <- 200 # Length of time-series
# R.length <- 50 # Number of replicate samples 
T.length <- 200 # Length of time-series
R.length <- 25 # Number of replicate samples 
lambda <- 2000 # expected number of counts for each sample

# Define ILR Basis - random
tree <- named_rtree(q+1)
tax.names <- tree$tip.label[order(tree$tip.label)]
sbp <- phylo2sbp(tree)
contrast.matrix.sbp <- buildilrBasep(sbp, p = c(1,1,1)) # no taxa weights
contrast.matrix.sbp <- contrast.matrix.sbp[tax.names,]

# Define Starting Condition for longitudinal
theta0 <- c(1, -3)

# Define Covariance Structures, also cholesky of these as 
# they will be used to generate multivariate normal random variables

#tmp <- matrix(.5*(runif(q^2)-1), ncol=q) # Random matrix 
#W <- t(tmp)%*%tmp # Random Covariance Matrix
W <- rbind(c(.05, .01), c(.01, .05))
L_W <- t(chol(W))

# tmp <- matrix((runif(q^2)-1), ncol=q) # Random matrix
# V <- t(tmp)%*%tmp # Random Covariance Matrix
V <- rbind(c(.2, -.1), c(-.1, .2))
L_V <- t(chol(V))

# Simulation: Build Simulation - Longitudinal Data -----------------------------

# Draw Variation for Longitudinal Data
Z.W <- matrix(rnorm(T.length*q, 0, 1), nrow = T.length)
w <- t(apply(Z.W, 1, function(x) L_W%*%x)) # Draw MVN
Z.V <- matrix(rnorm(T.length*q, 0, 1), nrow = T.length)
v <- t(apply(Z.V, 1, function(x) L_V%*%x)) # Draw MVN

# Build time-series from variation and starting condition
theta <- apply(rbind(theta0, w), 2, cumsum)[-1,] # not including starting value
eta <- theta + v

# Now draw multinomial from eta
pi <- unclass(ilrInv(eta, contrast.matrix.sbp))
Y <- t(apply(pi, 1, function(x) rmultinom(1, rpois(1, lambda), prob = x)))

# Drop a few observations from Y
Y[c(15:16, 20),] <- NA

# Name Y
colnames(Y) <- tax.names

# Simulation: Build Simulation - Replicate Data --------------------------------

# Draw Variation for Longitudinal Data
Z.V <- matrix(rnorm(R.length*q, 0, 1), nrow = R.length)
v <- t(apply(Z.V, 1, function(x) L_V%*%x)) # Draw MVN

# Build time-series from variation and ending condition
eta.csme <- t(apply(v, 1, function(x) x+theta[T.length,]))

# Now draw multinomial from eta
pi.csme <- unclass(ilrInv(eta.csme, contrast.matrix.sbp))
Y.csme <- t(apply(pi.csme, 1, function(x) rmultinom(1, rpois(1, lambda), prob = x)))
colnames(Y.csme) <- tax.names


# Simulation: Write Results -----------------------------------------------

write.tree(tree, file=out("simulation.true.tree"))
write.table(theta, file=out("simulation.theta.tsv"), sep="\t")
write.table(eta, file=out("simulation.eta.tsv"), sep="\t")
write.table(Y, file=out("simulation.Y.tsv"), sep="\t")
write.table(W, file=out("simulation.W.tsv"), sep="\t")
write.table(V, file=out("simulation.V.tsv"), sep="\t")
write.table(eta.csme, file=out("simulation.eta.csme.tsv"), sep="\t")
write.table(Y.csme, file=out("simulation.Y.csme.tsv"), sep="\t")
write.table(contrast.matrix.sbp, file=out("simulation.contrast.matrix.sbp.true.tsv"), sep="\t")


theta.true <- theta
eta.true <- eta
W.true <- W
V.true <- V
eta.csme.true <- eta.csme
contrast.matrix.sbp.true <- contrast.matrix.sbp


# Build second random ILR basis for analysis ------------------------------

# Define ILR Basis - random
tree <- named_rtree(q+1)
tax.names <- tree$tip.label[order(tree$tip.label)]
sbp <- phylo2sbp(tree)
contrast.matrix.sbp <- buildilrBasep(sbp, p = c(1,1,1)) # no taxa weights
contrast.matrix.sbp <- contrast.matrix.sbp[tax.names,]

write.tree(tree, file=out("simulation.tree"))
write.table(contrast.matrix.sbp, file=out("simulation.contrast.matrix.sbp.tsv"), sep="\t")


# Analysis: Setup Simulation Data, Longitudinal ---------------------------

# Arrange CSME Samples
Tr <- rep(1, nrow(Y.csme))
num.replicates <- table(Tr)
max.replicates <- max(num.replicates)
d <- dim(Y)
Y.replicates <- array(NA, dim = c(max.replicates, r, q+1))
for (i in unique(Tr)) {
  Y.replicates[1:num.replicates[i],i,] <- Y.csme[Tr==i,]
}

# rearrange Y for stan code
Y <- array(Y, dim=c(T.length, 1, q+1)) 

# Add in CSME Samples
Y <- abind(Y, Y.replicates, along=1)

# Set up time-indexes for stan code
tt.total <- 1:(T.length + max.replicates)
tt.observed.longitudinal <- tt.total[c(array_any_cases(Y)[1:T.length], rep(FALSE, max.replicates))]
tt.replicate <- (1:max.replicates)+T.length
tt.sample <- 1:T.length

# Set up indicator variables for times
tt.observed.ind <- matrix(FALSE,length(tt.total), r)
for (i in 1:r) tt.observed.ind[,i] <- array_any_cases(Y[,i,])
tt.sample.ind <- rep(FALSE,length(tt.total)) # Expects this one is vector
tt.sample.ind[tt.sample] <- TRUE
tt.replicate.ind <- matrix(FALSE,length(tt.total), r)
tt.replicate.ind[tt.replicate, 1] <- TRUE

# Remove Unobserved Days in Longitudinal Samples
Y <- Y[c(tt.observed.longitudinal, tt.replicate),,,drop=F]

# Create Y.obs vector of complete cases
Y.obs <- matrix(0, dim(Y)[1], r)
for (i in 1:r) Y.obs[array_any_cases(Y[,i,]),i] <- 1

# Describe Features of Simulated Data -------------------------------------

# Proportion of counts < x
enframe(c(Y[,1,])) %>% 
  ggplot(aes(x = value)) +
  stat_ecdf(geom="step")+
  scale_x_log10() +
  ggtitle("ECDF of Counts")

ed <- ecdf(Y)

# Proportion = 0
ed(0)
# proportion <= 10
ed(10)

# Analysis: Setup Data for Rstan ------------------------------------------

br_dat <- within(list(), {
  N_timepoints_total <- length(tt.total)
  N_timepoints_observed <- dim(Y)[1]
  N_timepoints_sample <- length(tt.sample)
  
  r <- dim(Y)[2]
  N_species <- dim(Y)[3]
  p <- N_species-1
  
  TT_obs_ind <- tt.observed.ind # Time points with non-missing observations
  TT_obs_ind[,1] <- as.integer(tt.observed.ind[,1])
  TT_rep_ind <- tt.replicate.ind
  TT_rep_ind[,1] <- as.integer(tt.replicate.ind[,1])
  TT_sample_ind <- tt.sample.ind # A Vector (not a matrix )
  
  Y <- Y
  Y_obs <- Y.obs
  
  contrast_matrix <- contrast.matrix.sbp
  
  d_F <- 1
  d_G <- 1
  FF <- array(0, dim=c(1,p,p))
  FF[1,,] <- diag(p)
  
  GG <- array(0, dim=c(1,p,p))
  GG[1,,] <- diag(p)
  
  m0 <- matrix(0, r, p)
  C0 <- 100*diag(p)

  W_scale_mean_prior <-rep(-1.5, p)
  W_scale_var_prior <- rep(2, p)

  V_scale_mean_prior <- rep(-1,p)
  V_scale_var_prior <- rep(2,p)

  W_lkj_prior <- 1
  V_lkj_prior <- 1
})

inits <- list()
for (i in 1:4){
  inits[[i]] <- within(list(), {
    d <- dim(Y)
    p <- d[3]-1
    eta <- array(0, dim=c(d[1], d[2], p))
    r <- runif(1)
    for (j in 1:dim(eta)[1]){
      eta[j,,] <- unclass(ilr(Y[j,,]+r,contrast.matrix.sbp))
    }
  })
}


# Analysis: Run Stan ------------------------------------------------------

model <- stanc_builder("toblerone.stan", 
                       isystem=paths[[machine]]$results)
#write(model$model_code, file="tmp.stan")

m <- stan_model(model_code = model$model_code, verbose = FALSE)
fit <- vb(m, data=br_dat, init=inits[[1]])
#save(fit, file=out("tmp.fit.vb.RData"))
#load(out("tmp.fit.vb.RData"))

ex <- rstan::extract(fit)
samples <- sample(1:1000, 4)

inits <- list()
for (i in seq_along(samples)){
  inits[[i]] <- within(list(), {
    eta <- ex$eta[samples[i],,,]
    dim(eta) <- c(dim(eta)[1], 1, dim(eta)[2])
    theta <- ex$theta[samples[i],,,]
    dim(theta) <- c(dim(theta)[1], 1, dim(theta)[2])
    eta_csme <- ex$eta_csme[samples[i],,]
    W_scale <- ex$W_scale[samples[i],,]
    dim(W_scale) <- c(1, length(W_scale))
    V_scale <- ex$V_scale[samples[i],,]
    dim(V_scale) <- c(1, length(V_scale))
    V_corr <- ex$V_corr[samples[i],,,]
    dim(V_corr) <- c(1, dim(V_corr))
    W_corr <- ex$W_corr[samples[i],,,]
    dim(W_corr) <- c(1, dim(W_corr))
  })
}

rm(fit, ex, samples)

include.vars <- c("eta", "theta", "V", "W", "eta_csme")
fit <- stan(model_code = model$model_code, data=br_dat, chains=4, iter=2000,
            init = inits)
save(fit, file=out("tmp.fit.RData"))
#load(out("tmp.fit.RData"))
# Color Pallates ----------------------------------------------------------

# For Vessels and ... 
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Also Plot VAR(1) Covariance ---------------------------------------------

to.keep.index <- (1:length(tt.observed.longitudinal))[(tt.observed.longitudinal +1) %in% tt.observed.longitudinal]

ml.ilr <- unclass(ilr(Y[to.keep.index,,]+0.65, contrast.matrix.sbp))
ml.cov <- cov(ml.ilr[-1,], lag(ml.ilr)[-1,], na.rm=TRUE)  



# Analysis: Convergence Diag ----------------------------------------------

p <- stan_rhat(fit, bins=50)
ggsave(out("stan_rhat.pdf"), plot=p, height=4, width=4, units="in")



# Analysis: Plot Simulation Stuff -----------------------------------------

p <- Y[,,] %>% 
  `colnames<-`(c("t1", "t2", "t3")) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="Time") %>%
  gather(taxa, Count, -Time) %>% 
  mutate(Time = tt.observed.longitudinal[as.integer(Time)]) %>% 
  ggplot(aes(x=Time, y = Count)) +
  geom_line() +
  facet_grid(taxa~., scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=7), 
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=7), 
        axis.title.x = element_text(size=7),
        strip.text.y = element_text(size=7))
ggsave(out("simulation.counts.pdf"), plot=p, height=2, width=5, units="in")

# Analysis: Comparing Smoothing Estiamtes ---------------------------------

tidy_fit <- tidy_stan_summary(summary(fit))$all
labels <- setNames(colnames(sbp), seq_along(colnames(sbp)))

tidy_fit <- tidy_fit %>% 
  mutate(dim_3 = labels[dim_3])

# Caluclate Eta_ml - use 0.65 as pseudocount
eta_ml <- array(0, dim=dim(inits[[1]]$eta))
for (j in 1:dim(eta_ml)[1]){
  eta_ml[j,,] <- ilr(br_dat$Y[j,,]+0.65,contrast.matrix.sbp)
}
eta_ml <- tidy_array(eta_ml) %>% 
  rename(mean=var) %>% 
  mutate(dim_3 = labels[dim_3]) %>% 
  mutate(dim_1 = tt.observed.longitudinal[dim_1])

# Theta
tidy_theta_true <- theta.true %>% 
  ilrInv(contrast.matrix.sbp.true) %>% 
  ilr(contrast.matrix.sbp) %>% 
  unclass() %>% 
  tidy_array() %>% 
  mutate(dim_3 = labels[dim_2], 
         dim_2 = 1, 
         mean = var)

p <- tidy_fit %>% 
  filter(parameter =="theta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=tidy_theta_true, color=cbPalette[3], size=1.5) + 
  geom_line(data=eta_ml, color=cbPalette[1]) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill=cbPalette[4], alpha=0.5) +
  geom_line(color=cbPalette[2], size=1) +
  facet_grid(dim_3~., scales="free_y") +
  xlab("Time") + ylab("Balance Value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=7), 
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=7), 
        axis.title.x = element_text(size=7),
        strip.text.y = element_text(size=7))
ggsave(out("simulation.theta.pdf"), plot=p, height=2, width=7, units="in")

# Eta

# eta <- read.table(out("simulation.eta.tsv"), sep="\t", row.names = NULL) %>%
#   select(-row.names)

tidy_eta_true <- eta.true %>% 
  ilrInv(contrast.matrix.sbp.true) %>% 
  ilr(contrast.matrix.sbp) %>% 
  unclass() %>% 
  tidy_array() %>% 
  mutate(dim_3 = labels[dim_2], 
         dim_2 = 1,
         #dim_1 = tt.observed[dim_1],
         mean = var)

p <- tidy_fit %>% 
  filter(parameter =="eta") %>%
  mutate(dim_1 = tt.observed.longitudinal[dim_1]) %>% # Name dates accordingly
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=tidy_eta_true, color=cbPalette[3], size=1.5) + 
  geom_line(data=eta_ml, color=cbPalette[1]) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill=cbPalette[4], alpha=0.5) +
  geom_line(color=cbPalette[2], size=1) +
  facet_grid(dim_3~., scales="free_y") +
  xlab("Time") + ylab("Balance Value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=7), 
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=7), 
        axis.title.x = element_text(size=7),
        strip.text.y = element_text(size=7))
ggsave(out("simulation.eta.pdf"), plot=p, height=2, width=7, units="in")


# Analysis: Simplex View of Covariance ------------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))

W <- rstan::extract(fit, pars="W")$W
W.clr <- array(0, dim = c(dim(W)[1], dim(W)[2]+1, dim(W)[3]+1))
for (k in 1:dim(W)[1]){
  W.clr[k,,] <- contrast.matrix.sbp %*%  W[k,,] %*% t(contrast.matrix.sbp)
}
dimnames(W.clr) <- list(iterations=NULL, 
                        dim_2=labels_taxa, 
                        dim_3=labels_taxa)


V <- rstan::extract(fit, pars="V")$V
V.clr <- array(0, dim = c(dim(V)[1], dim(V)[2]+1, dim(V)[3]+1))
for (k in 1:dim(W)[1]){
  V.clr[k,,] <- contrast.matrix.sbp %*%  V[k,,] %*% t(contrast.matrix.sbp)
}
dimnames(V.clr) <- list(iterations=NULL, 
                        dim_2=labels_taxa, 
                        dim_3=labels_taxa)

produce_simplex_plot <- function(subset.families, true.value=NULL, 
                                 to.plot=c("W"), plot.isotropic=F){
  W.subset <- W.clr[,subset.families,subset.families]
  V.subset <- V.clr[,subset.families, subset.families]
  
  # based on page 83 of Compositional Data in R
  r <- sqrt(qchisq(p=0.95, df=2)) 
  
  neutral.element <- acomp(c(1,1,1))
  names(neutral.element) <- subset.families
  
  plot(neutral.element, pch=16, col="white", axes=TRUE)
  for (i in sample(1:dim(W)[1], size=400)){
    if (plot.isotropic){
      ellipses(mean=neutral.element,
               var= diag(rep(sum(diag(W.subset[i,,]))/3, 3)),
               r=r,
               col=scales::alpha(cbPalette[8], 0.3), 
               lwd=.6)
    }
    if ("W" %in% to.plot){
      ellipses(mean=neutral.element, 
               var=W.subset[i,,], 
               r=r, 
               col=scales::alpha(cbPalette[6], 0.3), 
               lwd=.6)      
    }
    if ("V" %in% to.plot){
      ellipses(mean=neutral.element,
               var=V.subset[i,,],
               r=r,
               col=scales::alpha(cbPalette[6], 0.3), 
               lwd=.6)
    }
    if ("both" %in% to.plot){
      ellipses(mean=neutral.element,
               var=V.subset[i,,] + W.subset[i,,],
               r=r,
               col=scales::alpha(cbPalette[8], 0.3), 
               lwd=.6)
    }
  } 
  if (!is.null(true.value)){
    true.value <- true.value[subset.families, subset.families]
    ellipses(mean=neutral.element, 
             var = true.value, 
             r=r, 
             col=cbPalette[3], 
             lwd=2)
  }
}

pdf(out("simplex.W.pdf"), height=5, width=5)
produce_simplex_plot(labels_taxa, to.plot = c("W"), 
                     true.value = ilrvar2clr(W.true, contrast.matrix.sbp.true))
ellipses(mean = acomp(c(1,1,1)), 
         var = ilrvar2clr(ml.cov, contrast.matrix.sbp), 
         r = sqrt(qchisq(p=0.95, df=2)), 
         col=cbPalette[1], 
         lwd=2)
plot(ilrInv(w, contrast.matrix.sbp), add=TRUE, pch = 16, cex=.4)
dev.off()

pdf(out("simplex.V.pdf"), height=5, width=5)
produce_simplex_plot(labels_taxa, to.plot = c("V"),
                     true.value = ilrvar2clr(V.true, contrast.matrix.sbp.true))
plot(ilrInv(rbind(eta-theta, v), contrast.matrix.sbp), add=TRUE, pch = 16, cex=.4)
ellipses(mean = acomp(c(1,1,1)),
         var = ilrvar2clr(ml.cov, contrast.matrix.sbp),
         r = sqrt(qchisq(p=0.95, df=2)),
         col=cbPalette[1], 
         lwd=2)
dev.off()


# Plot Tree ---------------------------------------------------------------

p.tree <- tree %>% 
  ggtree(size=2)+
  geom_label(aes(label=label), size=12)+
  xlim(c(-.1,1.5)) +
  ylim(c(.8,3.2))
p.tree <- rotate(p.tree, name.to.nn(tree, "n1"))
ggsave(out("simulation.tree.pdf"), plot=p.tree, height=4, width=4, units="in")
