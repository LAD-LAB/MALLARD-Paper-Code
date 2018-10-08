library(rstan)
library(phyloseq)
library(tidyverse)
library(stringr)
library(lubridate)
library(compositions)
library(philr)
# devtools::install_github("jrnold/ssmodels-in-stan/StanStateSpace")
library(StanStateSpace)
library(padr)
library(ggrepel)
library(gridExtra)
library(abind)
library(ape)
library(ggtree)
library(ggtern)
library(ks)
library(ggridges)
library(shapes)
# Requires driver (github.com/jsilve24/driver) be installed
library(driver)
library(treeio)

set.seed(4)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Set up pathnames used throughout script
# Files are stored in a created folder (indexed by slurm array.id)
# in the "results" directory
# 
# creating_better_family_tree/ is searched for in the "base" directory.
# Working path is set to "results".  
machine <- "local"
paths <- list()
paths[["local"]] <- list("base" = "~/Research/mdlm/results/2018-04-18_code_github/",
                        "results" = "~/Research/mdlm/results/2018-04-18_code_github/mdlm_toblerone/",
                        "mapping" = '~/Research/data/_data_raw/sequencing.2016.03.04/2016.03.25MappingFile.MergedPool.txt',
                        "dada2"= "~/Research/data/_data_derived/sequencing.2016.03.04/dada2_2016_10_18/hardac/",
                        "runnotes" = "~/Research/data/_data_raw/sequencing.2016.03.04/2017.10.26_ResearcherIdentity.csv")
# paths[["hardac"]] <- list("base"="/data/davidlab/users/jds/mdlm_publication_2/",
#                           "results" = "/data/davidlab/users/jds/mdlm_publication_2/mdlm_toblerone/",
#                           "mapping" = '/data/davidlab/users/jds/mdlm_publication_2/dada2/0_mapping/2016.03.25MappingFile.MergedPool.txt',
#                           "dada2" = "/data/davidlab/users/jds/mdlm_publication_2/dada2/", 
#                           "runnotes" = '/data/davidlab/users/jds/mdlm_publication_2/runnotes/2017.10.26_ResearcherIdentity.csv')

setwd(paths[[machine]]$results)

# A few utility functions used thorughout this analysis. 
source("utils.R")

# Read Params and get array ID
if (machine == "local") {
  array.id <- 1
} else {
  array.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}
params <- read.table('sensitivity_analysis_params.txt', stringsAsFactors = FALSE)
params <- params[array.id,]
##

# Create a function out() which returns propper file name for saved results
out.dir <- paste0(paths[[machine]]$results, array.id, "/")
if (!dir.exists(out.dir)) dir.create(out.dir)
out <- function(filename){
  paste0(out.dir, filename)
}


# Load Data ---------------------------------------------------------------

# Read in corrected mapping file 
mapping <- import_qiime_sample_data(paths[[machine]]$mapping)

# Read in Run Notes 
runnotes <- read_csv(paths[[machine]]$runnotes)
runnotes <- select(runnotes, X.SampleID, Researcher, Note, Comments) %>% 
  as.data.frame()
rownames(runnotes) <- runnotes$X.SampleID
mapping <- cbind(mapping, select(runnotes[as.character(mapping$X.SampleID), ], -X.SampleID))

# Load Sequencing batch Info and add to metadata
t1 <- readRDS(file.path(paths[[machine]]$dada2, "seqtab.s1.nochim.rds"))
t2 <- readRDS(file.path(paths[[machine]]$dada2, "seqtab.s2.nochim.rds"))
t1 <- rownames(t1)
t2 <- rownames(t2)
batch <- ifelse(rownames(mapping) %in% t1, 1, 2)
mapping$batch <- batch

# Read in Phyloseq Object (saved from dada2 scripts)
ps <- readRDS(file.path(paths[[machine]]$dada2, "phyloseq.rds"))
sample_data(ps) <- mapping

# Any Increase in B. Ovatus? -----------------------------------------------

tax_table(ps)[,"Species"][!is.na(tax_table(ps)[,"Species"])] %>% 
  as.data.frame() %>% 
  arrange(Species)

# From this result is is apparent that B. Ovatus was not found by dada2 
# or not assigned a taxonomy at the species level. 

# Preprocessing -----------------------------------------------------------

# Now parse Dates from Sample IDs 
# Note: Day 0 is 2016.11.24
# Daily samples were generally taken at 3pm
sample.ids <- as.character(sample_data(ps)$X.SampleID)
postinnoc <- str_detect(sample.ids, "PostInnoc")
times <- str_split(sample.ids, "V", simplify=TRUE)[,1] %>% 
  str_replace(regex("\\.(Days|Day)"), "") %>% 
  str_replace(regex("(?!(m|d))\\.$"), "d\\.15h\\.") %>% 
  str_replace("00md\\.15h\\.", "") %>% 
  str_replace("T\\.\\.", "2015\\.11\\.") %>% 
  str_replace_all(regex("d|h"), "") %>%
  ymd_h()+days(1) 
times[is.na(times)] <- ymd_h("2015.11.1.15")
# Now set time-series to begin 
#times <- times + days(23) # Removed so that day of month = experimental day for plots
sample_data(ps)$time <- times
sample_data(ps)$postinnoc <- postinnoc



# Investigate Retained by Bioreactor --------------------------------------

ps.tmp.total <- tax_glom(ps, "Family")

first5 <- sample_data(ps.tmp.total) %>% 
  as("data.frame") %>% 
  filter(time <= min(time) + days(5)) %>% 
  pull(X.SampleID) %>% 
  as.character()

last5 <- sample_data(ps.tmp.total) %>% 
  as("data.frame") %>% 
  filter(time >= max(time) - days(5)) %>% 
  pull(X.SampleID) %>% 
  as.character()

cs.first5 <- colSums(otu_table(ps.tmp.total)[first5,])
cs.last5 <- colSums(otu_table(ps.tmp.total)[last5,])

still.there <- names(cs.first5[cs.first5 > 0]) %in% names(cs.last5[cs.last5 > 0])
print("Percentage of Families in first 5 still in last 5 days:")
(sum(still.there)/length(still.there))*100

print("These families represent X percentage of total counts")
tmp.totals <- colSums(otu_table(ps.tmp.total))
sum(tmp.totals[names(cs.first5[!still.there])])/sum(tmp.totals)*100

# Clean up some variables
rm(ps.tmp.total, first5, last5, cs.first5, cs.last5, still.there, tmp.totals)


# Continue PreProcessing --------------------------------------------------

# Investigate Distribution of Counts per sample
hist(sample_sums(ps))
quantile(sample_sums(ps), probs=seq(.01, .1, by=0.01))

# Based on this I will remove all samples with fewer than 5000 reads
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>5000, ps)
sum(sample_sums(ps))/total.reads


# Now collapse to family level
ps <- tax_glom(ps, "Family")

# Now retain only the high abundance taxa
# hist(log(taxa_sums(ps)))

ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.90*length(x)), TRUE)
(remaining.reads <- sum(sample_sums(ps))/total.reads)

# Remove Duplicate Samples (ignoring the Noise estimation samples
# and postinnoculation samples)
# Heather and I are not sure why there are duplicates... but there are
duplicates.to.remove <- subset_samples(ps, Normal_Noise_Sample=="Normal" & postinnoc==FALSE) %>%
  sample_data() %>% 
  .[,c('time', 'Vessel')] %>% 
  cbind(., duplicate=duplicated(.)) %>% 
  rownames_to_column("SampleID") %>% 
  filter(duplicate==TRUE) 
ps <- prune_samples(!(sample_names(ps) %in% duplicates.to.remove$SampleID), ps)


# Replace with Manual Created Tree ----------------------------------------

tree <- ape::read.tree(file.path(paths[[machine]]$base,"creating_better_family_tree/manual_families.tree"))
ape::is.binary.tree(tree) # Tree must be binary for PhILR
ape::is.rooted(tree) # Tree must be rooted for PhILR
tree <- ape::makeNodeLabel(tree, method="number", prefix='n')
phy_tree(ps) <- tree


# Setup Time-Series Data --------------------------------------------------

# Keep all time-series datapoints
Y <- subset_samples(ps, Normal_Noise_Sample=="Normal" & 
                      postinnoc==FALSE) %>% 
  otu_table() %>% 
  as("matrix") %>%
  as.data.frame() %>%
  bind_cols(., as(sample_data(ps)[,c("time","Vessel")][rownames(.),], "data.frame")) %>%
  rename(t=time, series=Vessel) %>% 
  arrange(t, series)

# Index of Daily
tt.hourly <- subset_samples(ps, Normal_Noise_Sample=="Normal" &
                              SampleType=="Hourly" &
                              postinnoc==FALSE) %>%
  subset_samples(time<(ymd("2015-12-19")-days(23))) %>% 
  sample_data() %>% 
  .[["time"]] %>% 
  unique() %>% 
  as.data.frame() %>% 
  pad(interval="hour") %>% 
  .[[1]]

# Index of Hourly 
tt.daily <- Y[hour(Y$t)== 15, ] %>% 
  pad(interval="day", group="series") %>% 
  select(t) %>% 
  .[[1]] %>% 
  unique()

tt.total <- unique(Y$t) %>% 
  as.data.frame() %>% 
  pad(interval="hour") %>% 
  .[[1]]

# Pad Y to hourly
Y <- expand.grid(t = tt.total, series = factor(1:4)) %>% 
  left_join(Y, by=c("t", "series"))


# Set up Replicate Samples from End of Time-series Samples --------------------

# Create Vector of "Treatments" - this is just the vessel
Y_csme <- subset_samples(ps, Normal_Noise_Sample=="Noise_Estimation" &
                           SampleType=="Hourly" &
                           postinnoc==FALSE) %>%
  otu_table() %>%
  as("matrix")

Tr <-  subset_samples(ps, Normal_Noise_Sample=="Noise_Estimation" &
                        SampleType=="Hourly" &
                        postinnoc==FALSE) %>%
  sample_data() %>%
  .[["Vessel"]]
Tr <- as.numeric(Tr)

# Add in Time Info
replicate.datetime <- ymd_h("2015.12.22.15") - days(23)
Y_csme <- data.frame("t" = replicate.datetime, 
                    "series" = Tr, 
                     Y_csme)

# Expand for "missing" replicates
Y_csme <- Y_csme %>% 
  group_by(series) %>% 
  mutate(rep.n = 1:n()) %>% 
  ungroup() %>% 
  left_join(expand.grid("rep.n" = 1:max(.$rep.n), "series" = 1:4), ., 
             by=c("rep.n", "series"))
Y_csme$t[is.na(Y_csme$t)] <- replicate.datetime

n.replicates <- max(Y_csme$rep.n)
Y_csme <- dplyr::select(Y_csme, -rep.n) %>% 
  mutate(series = factor(series))


# Combine Longitudinal and Cross-Sectional --------------------------------
# Combining to one array indexed by sampling point rather than time point 
# (see methods). 

Y <- bind_rows(Y, Y_csme) # COMBINE!!!

tt.total <- filter(Y, series =="1")$t
tt.total.index <- seq_along(tt.total)

Y_indicators <- Y %>% 
  group_by(series) %>% 
  arrange(t) %>% 
  mutate(tt.total.index = 1:n(), 
         rep=duplicated(t)) %>%
  ungroup() %>% 
  mutate(cc = complete.cases(.)) %>% 
  select(series, tt.total.index, cc, rep)

tt.observed.ind <- Y_indicators %>% 
  select(-rep) %>% 
  spread(series, cc) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix()

tt.replicate.ind <- Y_indicators %>% 
  select(-cc) %>% 
  spread(series, rep) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix()

tt.observed <- tt.total[(rowSums(tt.observed.ind) > 0)]
tt.sample <- unique(as.POSIXct(format(c(unique(tt.observed), tt.hourly), tz="UTC", usetz = T), tz = "UTC"))

tt.sample.ind <- Y %>% 
  group_by(series) %>% 
  arrange(t) %>% 
  mutate(tt.total.index = 1:n(), 
         rep=duplicated(t)) %>%
  mutate(sample = !rep & (t %in% tt.sample)) %>% 
  ungroup() %>% 
  select(series, tt.total.index, sample) %>% 
  spread(series, sample) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix() %>% 
  rowSums() %>% 
  (function(x) x>0)

# Following pretty much just for plotting
tt.hourly.index <- match(tt.hourly, tt.total)
tt.daily.index <- match(tt.daily, tt.total)

# Convert to Array --------------------------------------------------------

# Now map Y to array
Y <- split(Y, Y$series) %>% 
  map(~dplyr::select(.x, -series, -t)) %>% 
  abind::abind(along=3) %>% 
  aperm(c(1,3,2)) 

# Remove time points/rows of Y that have no valid observations
Y <- Y[array_any_cases(Y),,]
Y.obs <- array_any_cases(Y, margin=c(1,2))

# Fill in NA with zeros
Y[is.na(Y)] <- 0

# Now set up basis
sbp <- phylo2sbp(phy_tree(ps))
contrast.matrix.sbp <- gsi.buildilrBase(sbp)

# Keep this in to warn in case dependent packages changes
# Currently relying on phyloseq to align these quantities. 
if (!all.equal(rownames(contrast.matrix.sbp),dimnames(Y)[[3]])) stop("Rownames of Y don't line up with colnames of sbp")


# Perturbation Series -----------------------------------------------------
# Create time-series of the B. Ovatus perturbation for plotting 

x <- data.frame(x=rep(0, length(tt.total)), t = tt.total) 
x[round_date(x$t, "hour") == (ymd_h("2015.12.16.16") - days(23)),'x'] <- 1
x.index <- x$x

t.perturb <- which(x==1) # timepoint where 
t.perturb.daily <- floor_date(tt.total[t.perturb], "day")
t.perturb.hourly <- floor_date(tt.total[t.perturb], "hour")


# Setup data for Rstan ----------------------------------------------------

br_dat <- within(list(),{
  N_timepoints_total <- length(tt.total)
  N_timepoints_observed <- dim(Y)[1]
  N_timepoints_sample <- length(tt.sample)
  
  r <- dim(Y)[2]
  N_species <- dim(Y)[3]
  p <- N_species-1
  
  # Time points with non-missing observations
  TT_obs_ind <- apply(tt.observed.ind, MARGIN = c(1,2), as.integer) 
  # Which are replicates
  TT_rep_ind <- apply(tt.replicate.ind, MARGIN = c(1,2), as.integer)
  # Which timepoints to sample
  TT_sample_ind <- as.integer(tt.sample.ind)
  
  Y <- Y
  # same as TT_obs_ind but indexes refer to location in Y matrix
  Y_obs <- apply(Y.obs, MARGIN = c(1,2), as.integer)

  contrast_matrix <- contrast.matrix.sbp
  
  d_F <- 1
  d_G <- 1
  FF <- array(0, dim=c(1,p,p))
  FF[1,,] <- diag(p)
  
  GG <- array(0, dim=c(1,p,p))
  GG[1,,] <- diag(p)
  
  ## Priors ## Note these are read in from the sensitivity analysis params file. 
  m0 <- matrix(0, r, p)
  C0 <- params$c0*diag(p)
  
  W_scale_mean_prior <-rep(params$W_scale_mean, p)
  W_scale_var_prior <- rep(params$W_scale_var, p)
  
  V_scale_mean_prior <- rep(params$V_scale_mean,p)
  V_scale_var_prior <- rep(params$V_scale_var,p)
  
  W_lkj_prior <- params$W_lkj
  V_lkj_prior <- params$V_lkj
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
    
    W_corr <- array(diag(p), dim=c(1,p,p))
    V_corr <- array(diag(p), dim=c(1,p,p))
  })
}


# Exploring Priors --------------------------------------------------------

# W.prior <- rlnorm(1500, br_dat$W_scale_mean_prior, br_dat$W_scale_var_prior)
# V.scale.prior <- rlnorm(1500, br_dat$V_scale_mean_prior, br_dat$V_scale_var_prior)
# signal.to.noise.ratio.prior <- (W.prior/v.prior)

# p <- data.frame("v"=v.prior,
#            "W"=W.prior,
#            "U_scale" = U.scale.prior,
#            "V_scale" = V.scale.prior,
#            "W_to_v" = signal.to.noise.ratio.prior) %>%
#   gather(Parameter, Value) %>%
#   ggplot(aes(x=Value)) +
#   geom_density(fill="darkgrey")+
#   facet_wrap(~Parameter) +
#   xlim(c(0,10)) +
#   scale_x_log10()+
#   theme_bw()
# ggsave("priors_log10.pdf", plot=p, height=5, width=6, units="in")

# Run Model ---------------------------------------------------------------

model <- stanc_builder("toblerone.stan", 
                       isystem=paths[[machine]]$results)
#write(model$model_code, file="tmp.stan") # for debugging stan code
m <- stan_model(model_code = model$model_code, verbose = FALSE)

fit <- vb(m, data=br_dat, init=inits[[1]])
#save(fit, file=out("tmp.fit.vb.RData")) # If variational model is to be saved
#load(out("tmp.fit.vb.RData"))           # "    " loaded

ex <- rstan::extract(fit)

# Randomly pick which of 4 variational samples to use as inits for MCMC
samples <- sample(1:1000, 4)

# Initial values for MCMC chains based on those 4 variational samples
inits <- list()
for (i in seq_along(samples)){
  inits[[i]] <- within(list(), {
   eta <- ex$eta[samples[i],,,]
   theta <- ex$theta[samples[i],,,]
   eta_csme <- ex$eta_csme[samples[i],,]
   W_scale <- ex$W_scale[samples[i],,]
   dim(W_scale) <- c(1, length(W_scale))
   V_scale <- ex$V_scale[samples[i],,]
   dim(V_scale) <- c(1, length(V_scale))
   V_corr <- ex$V_corr[samples[i],,,]
   dim(V_corr) <- c(1, dim(V_corr))
   W_corr <- ex$W_corr[samples[i],,,]
   dim(W_corr) <- c(1, dim(W_corr))
   mu_csme <- ex$mu_csme[samples[i],,]
  })
}

rm(fit, ex, samples) # To save memory remove unneded variational results

# Specify which parameters to sample in MCMC and run MCMC
include.vars <- c("eta", "theta", "V", "W", "eta_csme") 
fit <- stan(model_code = model$model_code, data=br_dat, chains=4, iter=2000,
            init = inits)
save(fit, file=out("tmp.fit.RData")) # Save stanfit object for reanalysis
#load(out("tmp.fit.RData")) # Load "    "


# Color Pallates ----------------------------------------------------------

# For Vessels and such 
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Look at Rhat for convergence support of Results -----------------------------

p <- stan_rhat(fit, bins=50)
ggsave(out("stan_rhat.pdf"), plot=p, height=4, width=4, units="in")

# More detailed traceplots and convergence diagnostics can be viewed from shinystan
# shinystan::launch_shinystan(fit)

# Look at Smoothing Estimates ---------------------------------------------
# This section is just some simple "ugly" plots to get a better sence of 
# model results. 

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
  mutate(dim_1 = tt.observed[dim_1])

## At Daily Scale ##
# Eta
p <- tidy_fit %>% 
  filter(parameter =="eta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  filter(dim_1 %in% tt.daily) %>%  # just take daily samples
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=filter(eta_ml, dim_1 %in% tt.daily), color="red") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.daily))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  ggtitle("eta") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("daily.eta.pdf"), plot=p, width=13.1, height=24, units="in")

# theta
p <- tidy_fit %>% 
  filter(parameter =="theta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  filter(dim_1 %in% tt.daily) %>%  # just take daily samples
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=filter(eta_ml, dim_1 %in% tt.daily), color="red") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.daily))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  ggtitle("theta") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("daily.theta.pdf"), plot=p, width=13.1, height=24, units="in")

## At Hourly Scale ##
# Eta
p <- tidy_fit %>% 
  filter(parameter =="eta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  filter(dim_1 %in% tt.hourly) %>%  # just take daily samples
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=filter(eta_ml, dim_1 %in% tt.hourly), color="red") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  ggtitle("eta") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("hourly.eta.pdf"), plot=p, width=13.1, height=24, units="in")


# theta
p <- tidy_fit %>% 
  filter(parameter =="theta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  filter(dim_1 %in% tt.hourly) %>%  # just take daily samples
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=filter(eta_ml, dim_1 %in% tt.hourly), color="red") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  ggtitle("theta") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("hourly.theta.pdf"), plot=p, width=13.1, height=24, units="in")

## Combined Hourly and Daily ##
# theta
p <- tidy_fit %>% 
  filter(parameter =="theta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=eta_ml, color="red") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5),fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  ggtitle("eta") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("combined.theta.pdf"), plot=p, width=13.1, height=24, units="in")

# Nice version for publication
eta_ml <- mutate(eta_ml, dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2])
p <- tidy_fit %>% 
  filter(parameter =="theta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2]) %>% 
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=eta_ml, color="black") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill="#BED7E6", alpha=.7) +
  geom_ribbon(aes(ymax=p75, ymin=p25), fill="#6DAFD4", alpha=.7) +
  geom_line(color="#FB0D1C") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  scale_x_datetime(date_labels="Day %d") +
  ylab("Balance Value (e.i.)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1), 
        axis.title.x = element_blank())
ggsave(out("combined.pub.theta.pdf"), plot=p, width=8, height=10.5, units="in")


# Nice version for publication
#eta_ml <- mutate(eta_ml, dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2])
p <- tidy_fit %>% 
  filter(parameter =="eta") %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2]) %>% 
  ggplot(aes(x=dim_1, y=mean)) +
  geom_line(data=eta_ml, color="black") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill="#BED7E6", alpha=.7) +
  geom_ribbon(aes(ymax=p75, ymin=p25), fill="#6DAFD4", alpha=.7) +
  geom_line(color="#FB0D1C") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_grid(dim_3~dim_2, scales="free_y") +
  scale_x_datetime(date_labels="Day %d") +
  ylab("Balance Value (e.i.)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1), 
        axis.title.x = element_blank())
ggsave(out("combined.pub.eta.pdf"), plot=p, width=8, height=10.5, units="in")

# 
# # # Nice version for publication
# # eta_ml <- tidy_array(eta_ml) %>% 
# #   rename(mean=var) %>% 
# mutate(dim_3 = labels[dim_3]) %>%
#   mutate(dim_1 = tt.observed[dim_1])

labs <- tax_table(ps)[,"Family"]
labs <- setNames(as.character(as.data.frame(labs)$Family), rownames(labs))

eta_ml_prop_log <- log(miniclo_array(br_dat$Y + 0.65, 3)) %>% 
  gather_array(mean, dim_1, dim_2, dim_3) %>% 
  mutate(dim_3 = labs[rownames(sbp)[dim_3]]) %>%
  mutate(dim_1 = tt.observed[dim_1]) %>% 
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2])
  

theta <- rstan::extract(fit, pars="theta")$theta
theta.prop <- array(0, dim=c(dim(theta)[1:3], dim(theta)[4]+1))
for(i in 1:dim(theta)[1]){
  for (j in 1:dim(theta)[2]){
    theta.prop[i,j,,] <- ilrInv(theta[i,j,,], contrast.matrix.sbp)  
  }
}
tidy_theta <- theta.prop %>% 
  gather_array(val, iter, dim_1, dim_2, dim_3) %>% 
  mutate(val = log(val)) %>% 
  group_by(dim_1, dim_2, dim_3) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(dim_3 = labs[rownames(sbp)[dim_3]]) %>% 
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2])
  
p <- ggplot(tidy_theta, aes(x=dim_1, y = mean)) +
  geom_line(data=eta_ml_prop_log, color="black") +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill="#BED7E6", alpha=.7) +
  geom_ribbon(aes(ymax=p75, ymin=p25), fill="#6DAFD4", alpha=.7) +
  geom_line(color="#FB0D1C") +
  facet_grid(dim_3~dim_2, scales="free_y") +
  scale_x_datetime(date_labels="Day %d") +
  ylab("Log Proportions") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1), 
        axis.title.x = element_blank(), 
        strip.text.y = element_text(size=10, angle=0))
ggsave(out("combined.pub.theta.prop.pdf"), plot=p, width=8, height=10.5, units="in")


# Plot tree ---------------------------------------------------------------

tax <- tax_table(ps)[,2:5] %>% 
  as.data.frame %>% 
  mutate(seq=rownames(.)) %>% 
  dplyr::select(seq, everything()) %>% 
  unite(d, everything() ,sep="\n") %>% 
  .[["d"]] %>% 
  setNames(rownames(tax_table(ps)))

o <- names(tax) %>% 
  substr(5, nchar(.)) %>% 
  as.integer() %>% 
  order

tax <- factor(tax, levels = tax[names(tax)[o]])

families <- as.data.frame(as(tax_table(ps), "matrix"))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))

tree <- phy_tree(ps)
tree$tip.label <- families[tree$tip.label]
tree$edge.length <- rep(NULL, nrow(tree$edge))
V.tmp <- contrast.matrix.sbp
rownames(V.tmp) <- families[rownames(V.tmp)]
p.tree <- ggtree(tree) + 
  geom_label2(aes(label=label, subset=!isTip), size = 3) + 
  geom_tiplab(size=3) +
  xlim(c(0, 9))
p.tree <- annotate_sbp(tree, V.tmp, p.tree, sep="\n")
#ggsave(out("family.tree.pdf"), plot=p.tree, width=3, height=6, units="in")



# Smoothing Estimates Figure ----------------------------------------------

p <- tidy_fit %>%
  filter(parameter == "theta") %>% 
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  mutate(dim_3=factor(dim_3, levels=c("n12", "n15", "n16", 
                                      "n2", "n13", "n14", 
                                      "n1", "n3", "n10"))) %>% # Order facets
  ggplot(aes(x=dim_1, y=mean)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5, fill = factor(dim_2)), alpha=0.5) +
  #geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_wrap(~dim_3) + 
  theme_bw() +
  guides(fill=guide_legend("Vessel")) +
  scale_fill_manual(values = cbPalette) +
  ylab("Balance Value") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90)) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("philr_combined.pdf"), plot=p, width=8, height=6, units="in")


p <- arrangeGrob(p.tree, p, ncol=2, widths=c(1,2))
ggsave(out("philr.basis.pdf"), plot=p, width=10, height=7, units="in")


# Just at hourly level
p <- tidy_fit %>%
  filter(parameter == "theta") %>% 
  mutate(dim_1 = tt.sample[dim_1]) %>% # Name dates accordingly
  filter(dim_1 %in% tt.hourly) %>% 
  mutate(dim_3=factor(dim_3, levels=c("n12", "n15", "n16", 
                                      "n2", "n13", "n14", 
                                      "n1", "n3", "n10"))) %>% # Order facets
  #filter(dim_3 == "n12") %>% 
  ggplot(aes(x=dim_1, y=mean)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5, fill = factor(dim_2)), alpha=0.5) +
  #geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_wrap(~dim_3, scales = "free_y") + 
  theme_bw() +
  guides(fill=guide_legend("Vessel")) +
  scale_fill_manual(values = cbPalette) +
  ylab("Balance Value") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90)) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("philr_hourly.pdf"), plot=p, width=8, height=6, units="in")

p <- arrangeGrob(p.tree, p, ncol=2, widths=c(1,2))
ggsave(out("philr.basis.hourly.pdf"), plot=p, width=10, height=7, units="in")


# 2:Plot n12 tree -----------------------------------------------------------
# Focusing on balance n12 for the moment
# Here is the tree 

tree <- phy_tree(ps)
tree.n12 <- extract.clade(tree, name.to.nn(tree, "n12"))
print(tax_table(ps)[tree.n12$tip.label,1:5])
tree.n12$tip.label <- paste0(tree.n12$tip.label, "  ")

p <- ggtree(tree.n12, size=2) +
  geom_tiplab() +
  xlim(0, 3.3)
p <- flip(p, name.to.nn(tree.n12, "n13"), name.to.nn(tree.n12, "n15"))
ggsave(out("tree.n12.pdf"), plot=p, height = 4, width = 1.5, units="in")


# 2:Fast Changes at Sub-Daily Intervals -------------------------------------
# Focusing on balance n12 for the moment
# Here are the dynamics 

t.overlap <- data.frame(`d` = tt.daily[tt.daily %in% tt.hourly])

# 23 days is subtracted to convert actual days of the month to "experimental days"
dates.changemedia <- c(ymd_hm("2015-12-11-12-34", tz="UTC") - days(23),
                       ymd_hm("2015-12-12-14-49", tz="UTC") - days(23),
                       ymd_hm("2015-12-14-15-20", tz="UTC") - days(23),
                       ymd_hm("2015-12-16-19-19", tz="UTC") - days(23),
                       ymd_hm("2015-12-18-15-30", tz="UTC") - days(23)) %>%
  enframe(value = "t")

# Figure out which taxa make up balance n12 (using PhILR function name.balance)
name.balance(phy_tree(ps), tax_table(ps), "n12", return.votes = c('up', 'down'))

# Look at balance n12 closer
p <- tidy_fit %>% 
  filter(parameter == "theta", 
         dim_3 == "n12") %>%
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2]) %>% 
  mutate(dim_1 = tt.sample[dim_1]) %>% 
  filter(dim_1 %in% tt.hourly) %>% 
  ggplot(aes(x=dim_1, y=mean, group=dim_2)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill="#BED7E6", alpha=1) +
  geom_ribbon(aes(ymax=p75, ymin=p25), fill="#6DAFD4", alpha=1) +
 # geom_vline(aes(xintercept=as.integer(d)), data = t.overlap, alpha=.7, size=1) +
  geom_line(color="#FB0D1C", size=.5) +
  #ggtitle("Phylum Bacteroidetes / Phylum Fusobacteria + Phylum Proteobacteria", 
  #        "Posterior 95% Credible Intervals for Theta with daily sampling indicated") +
  labs(fill = "Vessel")+
  xlab("") +
  ylab("Balance Value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_datetime(date_labels="Day %d") +
  facet_wrap(~dim_2, scales="free_y")
#ggsave(out("n12.hourly.pdf"), plot = p, height=3, width=7, units="in")
ggsave(out("n12.hourly.pdf"), plot = p, height=4, width=5.5, units="in")


p <- tidy_fit %>%
  filter(parameter == "theta",
         dim_3 == "n12",
         dim_2==4) %>%
  mutate(dim_1 = tt.sample[dim_1]) %>%
  filter(dim_1 %in% tt.hourly) %>%
  ggplot(aes(x=dim_1, y=mean, group=dim_2)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill=cbPalette[4], alpha=0.5) +
  geom_vline(aes(xintercept=as.integer(d)), data = t.overlap, alpha=.7, size=1) +
  geom_line(color="blue", size=.5) +
  #ggtitle("Phylum Bacteroidetes / Phylum Fusobacteria + Phylum Proteobacteria",
  #        "Posterior 95% Credible Intervals for Theta with daily sampling indicated") +
  labs(fill = "Vessel")+
  xlab("") +
  ylab("Balance Value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n12.v4only.hourly.pdf"), plot = p, height=2.8, width=4.5, units="in")


# Plot Researcher Sampling and Media Changes as well 
tmp <- select(sample_data(ps), time, Vessel, Researcher) %>% 
  mutate(Vessel = as.integer(Vessel)) %>% 
  mutate(Vessel = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[Vessel])
p <- tidy_fit %>% 
  filter(parameter == "theta", 
         dim_3 == "n12") %>%
  mutate(dim_2 = c("Vessel 1", "Vessel 2", "Vessel 3", "Vessel 4")[dim_2]) %>% 
  mutate(dim_1 = tt.sample[dim_1]) %>% 
  left_join(tmp, 
            by=c("dim_1" = "time", "dim_2" = "Vessel")) %>% 
  mutate(Researcher = factor(as.integer(factor(Researcher)))) %>% # Removing Initials
  filter(dim_1 %in% tt.hourly) %>% 
  ggplot(aes(x=dim_1, y=mean, group=dim_2)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill="#BED7E6", alpha=1) +
  geom_ribbon(aes(ymax=p75, ymin=p25), fill="#6DAFD4", alpha=1) +
  #geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill=cbPalette[4], alpha=0.5) +
  geom_vline(aes(xintercept=as.integer(d)), data = t.overlap, alpha=.7, size=1) +
  geom_line(aes(color = Researcher),  size=1) +
  geom_vline(aes(xintercept=as.integer(t)), data=dates.changemedia, color="red", size=2, linetype=2) +
  labs(fill = "Vessel")+
  xlab("") +
  ylab("Balance Value (e.i.)") +
  theme_bw()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_blank()) +
  scale_x_datetime(date_labels="Day %d") +
  facet_wrap(~dim_2)
#ggsave(out("n12.hourly.environmental_variables.pdf"), plot = p, height=2.8, width=4.5, units="in")
ggsave(out("n12.hourly.environmental_variables.pdf"), plot = p, height=4, width=6, units="in")


# Where there any run notes of interest in this section?
tmp <- select(sample_data(ps), time, Vessel, Note) %>% 
  mutate(Vessel = as.integer(Vessel))
# Plot Person Sampling 
p <- tidy_fit %>% 
  filter(parameter == "theta", 
         dim_3 == "n12", 
         dim_2==4) %>%
  mutate(dim_1 = tt.sample[dim_1]) %>% 
  left_join(tmp, 
            by=c("dim_1" = "time", "dim_2" = "Vessel")) %>% 
  filter(dim_1 %in% tt.hourly) %>% 
  ggplot(aes(x=dim_1, y=mean, group=dim_2)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5), fill=cbPalette[4], alpha=0.5) +
  geom_vline(aes(xintercept=as.integer(d)), data = t.overlap, alpha=.7, size=1) +
  geom_line(aes(color = Note),  size=1) +
  geom_vline(aes(xintercept=as.integer(t)), data=dates.changemedia, color="red", size=2, linetype=2) +
  labs(fill = "Vessel")+
  xlab("") +
  ylab("Balance Value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n12.hourly.runnotes.pdf"), plot = p, height=2.8, width=4.5, units="in")


# Scale of Balances 
contrast.matrix.sbp[,"n12"]

# Note that a 1 unit change in a balance is an equivalent amount of evidence
# to 1 taxa increasing 4 fold over another
exp(1/sqrt(2))^2
# Cite Evidence information paper.


# 2:Effective Time-Scale ----------------------------------------------------
# Figure 2B

V <- rstan::extract(fit, pars="V")$V
W <- rstan::extract(fit, pars="W")$W
trace.V <- rep(0, dim(V)[1])
trace.W <- rep(0, dim(V)[1])
for(i in 1:dim(V)[1]){
  trace.V[i] <- sum(diag(V[i,,]))
  trace.W[i] <- sum(diag(W[i,,]))
}

lags <- 1:24
trace.W.tidy <- as.list(lags) %>% 
  map(`*`, trace.W) %>% 
  do.call(rbind, .) %>% 
  t()
colnames(trace.W.tidy) <- lags
trace.W.tidy <- trace.W.tidy %>% 
  as.data.frame() %>% 
  gather(lag, trace) %>% 
  group_by(lag) %>% 
  summarize(p2.5 = quantile(trace, probs = 0.025), 
            p50 = quantile(trace, probs = 0.5), 
            p97.5 = quantile(trace, probs = 0.975))

trace.V.tidy <- trace.V %>% 
  enframe("lag", "trace") %>% 
  summarize(p2.5 = quantile(trace, probs = 0.025), 
            p50 = quantile(trace, probs = 0.5), 
            p97.5 = quantile(trace, probs = 0.975))
trace.V.tidy <- trace.V.tidy[rep(seq_len(nrow(trace.V.tidy)), each=length(lags)),] %>% 
  rownames_to_column(var="lag")

# Modified From here: 
# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_labels <- function(l) {
  # turn in to character string in non-scientific notation
  l <- format(l, scientific = FALSE)
  # Remove .0
  l <- gsub("\\.0", "", l)
  parse(text=l)
}

intersection <- 3.5
p <- ggplot(trace.W.tidy, aes(x = as.integer(lag), ymin = p2.5, y = p50, ymax = p97.5)) +
  geom_segment(x = intersection, xend = intersection, 
               y = -Inf, yend = trace.V.tidy$p50[1], linetype = "dotted") +
  geom_ribbon(data = trace.V.tidy, fill = cbPalette[4], alpha=0.5) +
  geom_line(data = trace.V.tidy, color = "blue", alpha=0.8) +
  geom_ribbon(fill=cbPalette[4], alpha=0.5) +
  geom_line(color = "blue", alpha=0.8) + 
  scale_x_continuous(breaks=c(intersection, seq(0, 25, by=5)), 
                     labels = fancy_labels) +
  theme_bw()+
  theme(axis.text.x = element_text(size=9, angle=90, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        panel.grid.minor= element_blank())
ggsave(out("bio.to.tech.noise.figure.pdf"), plot=p, height=2.25, width = 3, units="in")


# Visualize Abundance at Family Level ----------------------------
# A few plots to visualize raw counts at the family level. 

# col <- treepalette(as.data.frame(tax_table(ps)), c("Phylum", "Family"), method="HSV", 
#                    palette = cbPalette)
# col <- col[complete.cases(col),]

# [1] "Acidaminococcaceae"  "Bacteroidaceae"     
# [3] "Desulfovibrionaceae" "Enterobacteriaceae" 
# [5] "Fusobacteriaceae"    "Lachnospiraceae"    
# [7] "Porphyromonadaceae"  "Rikenellaceae"      
# [9] "Ruminococcaceae"     "Synergistaceae" 

col <- c("Synergistaceae" = "#E5CD34", "Fusobacteriaceae" = "#6C6011", 
  "Acidaminococcaceae" = "#2B5B51", "Ruminococcaceae" = "#1D9376", "Lachnospiraceae" = "#A3E4D5", 
  "Rikenellaceae" = "#4D2599", "Bacteroidaceae" = "#1E0748", "Porphyromonadaceae" = "#9C74E3", 
  "Desulfovibrionaceae" = "#6B2B0B", "Enterobacteriaceae" = "#E46628")


p <- driver::miniclo_array(Y, 3) %>% 
driver::gather_array(counts, t, vessel, taxa) %>% 
  mutate(t = tt.observed[t], 
         Family = as.character(as.data.frame(tax_table(ps))[['Family']])[taxa]) %>% 
  mutate(Family = factor(Family, levels = names(col))) %>% 
  mutate(vessel = paste("Vessel", vessel)) %>% 
  filter(t < ymd_hms("2015-11-29 15:00:00 UTC")) %>% 
  filter(t %in% tt.daily) %>%
  ggplot(aes(x = t, y = counts, fill = Family)) +
  geom_area() +
  facet_wrap(~vessel) +
  theme_minimal() +
  scale_fill_manual(values = col) +
  scale_x_datetime(date_labels="Day %d") +
  ylab("Proportions") +
  #guides(fill = element_text(size=16))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position="bottom", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14), 
        strip.text.x = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=16))
ggsave(out("daily.counts.pdf"), plot=p, height=7, width=12, units="in")



p <- driver::miniclo_array(Y, 3) %>% 
  driver::gather_array(counts, t, vessel, taxa) %>% 
  mutate(t = tt.observed[t], 
         Family = as.character(as.data.frame(tax_table(ps))[['Family']])[taxa]) %>% 
  mutate(Family = factor(Family, levels = names(col))) %>% 
  mutate(vessel = paste("Vessel", vessel)) %>% 
  filter(t < ymd_hms("2015-11-29 15:00:00 UTC")) %>% 
  filter(t %in% tt.hourly) %>%
  ggplot(aes(x = t, y = counts, fill = Family)) +
  geom_area() +
  geom_vline(xintercept=as.integer(t.perturb.hourly)) +
  facet_wrap(~vessel) +
  theme_minimal() +
  scale_fill_manual(values = col) +
  scale_x_datetime(date_labels="Day %d") +
  ylab("Proportions") +
  #guides(fill = element_text(size=16))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position="bottom", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14), 
        strip.text.x = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=16))
ggsave(out("hourly.counts.pdf"), plot=p, height=7, width=12, units="in")

# Y.na <- Y
# for (i in 1:ncol(Y.obs)){
#   Y.na[!Y.obs[,i],i,] <- NA
# }
# 
# for (i in 1:dim(Y)[1]){
#   for (j in 1:dim(Y)[2]){
#     for (k in 1:dim(Y)[3]){
#       if(is.na(Y.na[i,j,k])) Y.na[i,j,k] <- Y.na[i-1,j,k] - 1
#     }
#   }
# }
# 
# dat <- driver::gather_array(Y.na, counts, t, vessel, taxa) %>% 
#   mutate(t = tt.observed[t], 
#          Family = as.character(as.data.frame(tax_table(ps))[['Family']])[taxa]) #%>%
#   # filter(t < ymd_hms("2015-11-29 15:00:00 UTC")) %>% 
#   # filter(t %in% tt.hourly) %>% 
#   # ggplot(aes(x = t, y = counts, group = Family, fill = Family)) +
#   # stat_steamgraph() +
#   # facet_grid(vessel~.) +
#   # scale_x_datetime(date_labels="Day %d")
# 
# dat %>% 
#   filter(vessel == 2) %>% 
#   filter(t < ymd_hms("2015-11-29 15:00:00 UTC")) %>% 
#   filter(t %in% tt.daily) %>% 
#   arrange(t) %>% 
#   select(-taxa) %>%
#   #spread(Family, counts) %>% 
#   #gather(Family, counts, -t, -vessel) %>%
#   streamgraph(Family, counts, t, interactive=FALSE, offset="expand") %>% 
#   sg_axis_x(tick_units="day", tick_format="Day %d", tick_interval=5)
# 




# Daily
p <- subset_samples(ps, hour(time)==15 & 
                      Normal_Noise_Sample=="Normal" & 
                      postinnoc==FALSE) %>% 
  otu_table() %>% 
  as("matrix") %>%
  as.data.frame() %>%
  bind_cols(., as(sample_data(ps)[,c("time","Vessel")][rownames(.),], "data.frame")) %>%
  rename(t=time, series=Vessel) %>% 
  gather(seq, count, -t, -series) %>% 
  ggplot(aes(x=t, y=count)) +
  geom_line() +
  facet_grid(seq~series)+
  scale_y_log10() +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("counts.daily.pdf"), plot=p, width=13.1, height=24, units="in")

# Hourly
p <- subset_samples(ps, Normal_Noise_Sample=="Normal" &
                      SampleType=="Hourly" &
                      postinnoc==FALSE) %>%
  subset_samples(time<ymd("2015-12-19")) %>%
  otu_table() %>%
  as("matrix") %>%
  as.data.frame() %>%
  bind_cols(., as(sample_data(ps)[,c("time","Vessel")][rownames(.),], "data.frame")) %>%
  rename(t=time, series=Vessel) %>%
  gather(seq, count, -t, -series) %>%
  ggplot(aes(x=t, y=count)) +
  geom_line() +
  facet_grid(seq~series)+
  scale_y_log10() +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("counts.hourly.pdf"), plot=p, width=13.1, height=24, units="in")

# Combined Hourly and Daily
p <- subset_samples(ps, Normal_Noise_Sample=="Normal" & 
                      postinnoc==FALSE) %>% 
  otu_table() %>% 
  as("matrix") %>%
  as.data.frame() %>%
  bind_cols(., as(sample_data(ps)[,c("time","Vessel")][rownames(.),], "data.frame")) %>%
  rename(t=time, series=Vessel) %>% 
  gather(seq, count, -t, -series) %>% 
  ggplot(aes(x=t, y=count)) +
  geom_line() +
  facet_grid(seq~series)+
  scale_y_log10() +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("counts.combined.pdf"), plot=p, width=13.1, height=24, units="in")

# Investigate the Signal to Noise Ratio -----------------------------------

V <- rstan::extract(fit, pars="V")$V
W <- rstan::extract(fit, pars="W")$W
trace.V <- rep(0, dim(V)[1])
trace.W <- rep(0, dim(V)[1])
for(i in 1:dim(V)[1]){
  trace.V[i] <- sum(diag(V[i,,]))
  trace.W[i] <- sum(diag(W[i,,]))
}
signal.to.noise.ratio <- trace.W/trace.V
percent.bio.noise <- trace.W/(trace.V + trace.W)
quantile(signal.to.noise.ratio, probs = c(0.025, 0.25, .5, .75, .975)) %>% 
  write.table(file=out("signal_to_noise_ratio_posterior_summary.txt"))
quantile(percent.bio.noise, probs = c(0.025, 0.25, .5, .75, .975)) %>% 
  write.table(file=out("percent_bio_noise_summary.txt"))


# 2: Investigate Prior for Signal to Noise vs. Posterior ---------------------

rlnorm.square <- function(n, D, meanlog, sdlog){
  r <- rep(0, n)
  for (i in 1:n){
    x <- rlnorm(D, meanlog, sdlog)
    r[i] <- sum(x^2)
  }
  r
}

V.prior <- rlnorm.square(4000, dim(Y)[3], br_dat$V_scale_mean_prior[1], br_dat$V_scale_var_prior[1])
W.prior <- rlnorm.square(4000, dim(Y)[3], br_dat$W_scale_mean_prior[1], br_dat$W_scale_var_prior[1])
signal.to.noise.ratio.prior <- (W.prior/V.prior)


# Modified From here: 
# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # remove + values
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}

p <- data.frame(Posteior=signal.to.noise.ratio, Prior=signal.to.noise.ratio.prior) %>%
  gather(Distribution, Value) %>%
  ggplot(aes(x=Value, y=Distribution, fill = Distribution, height=..scaled..)) +
  geom_joy(scale=4) +
  scale_y_discrete(expand=c(0.01, 0)) +
  scale_x_log10(breaks=c(10^-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2), 
                labels=fancy_scientific) + 
  theme_joy(font_size = 6) +
  scale_fill_manual(values = cbPalette[6:7]) +
  guides(fill=FALSE) +
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_blank())
ggsave(out("signal_to_noise_ratio_prior.pdf"), plot=p, height=1.1, width=2.5, units="in")


# Variance Breakdown by Family --------------------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))

W <- rstan::extract(fit, pars="W")$W
V <- rstan::extract(fit, pars="V")$V

dW <- dim(W)
dV<- dim(V)
W.clr.diag <- array(0, dim=c(dW[1], dW[2]+1))
V.clr.diag <- array(0, dim=c(dV[1], dV[2]+1))
for (i in 1:dW[1]){
  W.clr.diag[i,] <- diag(contrast.matrix.sbp %*% W[i,,] %*% t(contrast.matrix.sbp))
  V.clr.diag[i,] <- diag(contrast.matrix.sbp %*% V[i,,] %*% t(contrast.matrix.sbp))
}


# Posterior distribution for CLR Signal to Noise Ratios
family.signal.to.noise.ratio <- W.clr.diag/V.clr.diag
colnames(family.signal.to.noise.ratio) <- families[labels_taxa]
p <- family.signal.to.noise.ratio %>% 
  as.data.frame() %>% 
  gather(CLR, BioVar.to.TechVar.Ratio) %>% 
  ggplot(aes(x=CLR, y=BioVar.to.TechVar.Ratio)) +
  geom_violin(scale="area", fill="grey") +
  #ggtitle("Biological Variation to Technical Variation Ratio by CLR Coordinates") +
  xlab("CLR Coordinate") +
  ylab("Tr(W)/Tr(V) Ratio") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(out("CLR.BioVar.To.TechVar.Ratio.pdf"), 
       plot=p, height=6, width=8, units="in")

# Biological Variation by Family CLR
family.ordering <- W.clr.diag %>% 
  colMeans() %>% 
  setNames(families[labels_taxa]) %>% 
  .[order(.)]

colnames(W.clr.diag) <- families[labels_taxa]
p <- W.clr.diag %>% 
  as.data.frame() %>% 
  gather(CLR, BioVar) %>%
  mutate(CLR = factor(CLR, levels = names(family.ordering))) %>% 
  ggplot(aes(x=CLR, y=BioVar)) +
  stat_summary_bin(fun.data = function(x) median_hilow(x, .95), geom=c("point")) +
  stat_summary_bin(fun.data = function(x) median_hilow(x, .95), geom=c("errorbar")) +
  xlab("CLR Coordinate") +
  ylab("Biological Variation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(out("CLR.BioVar.pdf"), plot=p, height=6, width=8, units="in")


# Technical Variation by Family CLR
family.ordering <- V.clr.diag %>% 
  colMeans() %>% 
  setNames(families[labels_taxa]) %>% 
  .[order(.)]

colnames(V.clr.diag) <- families[labels_taxa]
p <- V.clr.diag %>% 
  as.data.frame() %>% 
  gather(CLR, TechVar) %>%
  mutate(CLR = factor(CLR, levels = names(family.ordering))) %>% 
  ggplot(aes(x=CLR, y=TechVar)) +
  stat_summary_bin(fun.data = function(x) median_hilow(x, .95), geom=c("point")) +
  stat_summary_bin(fun.data = function(x) median_hilow(x, .95), geom=c("errorbar")) +
  xlab("CLR Coordinate") +
  ylab("Technical Variation [Tr(V)]") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(out("CLR.TechVar.pdf"), plot=p, height=6, width=8, units="in")


# Proportion of Bio Variation attributable to Top 4 ---------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))

W <- rstan::extract(fit, pars="W")$W

dW <- dim(W)
W.clr.diag <- array(0, dim=c(dW[1], dW[2]+1))
for (i in 1:dW[1]){
  W.clr.diag[i,] <- diag(contrast.matrix.sbp %*% W[i,,] %*% t(contrast.matrix.sbp))
}

focus <- c("seq_1", "seq_4", "seq_6", "seq_9")
colnames(W.clr.diag) <- labels_taxa
W.clr.diag %>% 
  miniclo() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="iter") %>% 
  gather(coord, val, -iter) %>% 
  mutate(top.4 = ifelse(coord %in% focus, TRUE, FALSE)) %>% 
  group_by(iter, top.4) %>% 
  summarize(sum.val = sum(val)) %>% 
  ungroup() %>% 
  group_by(top.4) %>% 
  summarize(mean = mean(sum.val), 
            p2.5 = quantile(sum.val, prob=0.025), 
            p97.5 = quantile(sum.val, prob=0.975)) %>% 
  write.table(file=out("percent.variance.top4.tsv"))


# 4:Investigate Innoculum Community -----------------------------------------

labels <- setNames(rownames(sbp), seq_along(rownames(sbp)))
theta.starting <- rstan::extract(fit, pars="theta")$theta[,1,,, drop=F] # Extract first datapoint

d <- dim(theta.starting)
theta_prop.starting <- array(0, dim=c(d[1:3], d[4]+1))
for (i in 1:d[1]){
  for (j in 1:d[3]){
    theta_prop.starting[i,,j,] <- unclass(clr(ilrInv(theta.starting[i,,j,], contrast.matrix.sbp)))
  }
}


day0.relabund <- drop(theta_prop.starting) %>% 
  tidy_array() %>% 
  mutate(dim_3 = labels[dim_3]) %>% 
  group_by(dim_3, dim_1) %>% 
  summarize(m.relabund = mean(var)) %>% # Mean rather than geometric mean (already in CLR)
  ungroup() %>% 
  rename(dim_2 = dim_3)

bio.var<- tidy_array(clr(miniclo(W.clr.diag))) %>% 
  mutate(dim_2 = labels[dim_2]) %>% 
  dplyr::select(-parameter) %>% 
  rename(bio.var = var)

d <- full_join(day0.relabund, bio.var, by=c("dim_1", "dim_2")) 

d.medians <- d %>% 
  group_by(dim_2) %>% 
  summarize(m.relabund = median(m.relabund), 
            bio.var = median(bio.var)) %>% 
  mutate(dim_2 = tax[dim_2]) %>% 
  mutate(dim_2 = str_split(dim_2, "\n", simplify=TRUE)[,5])

# Find Posterior Interval for the Regression
mods <- d %>% 
  group_by(dim_1) %>% 
  do(mod = lm(bio.var ~ m.relabund, data = .))

predict.m.relabund <- data.frame(m.relabund=seq(-4,4, by=.2))
preds <- mods$mod %>% 
  map(predict, predict.m.relabund) %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>% 
  gather(m.relabund, bio.var) %>% 
  mutate(m.relabund = as.integer(m.relabund)) %>% 
  mutate(m.relabund = predict.m.relabund$m.relabund[m.relabund]) %>% 
  group_by(m.relabund) %>% 
  summarize(p2.5 = quantile(bio.var, p=0.025),
            p50 = quantile(bio.var, p=0.5), 
            p97.5 = quantile(bio.var, p=0.975))

# Calculate Contours
contours <- d %>% 
  group_by(dim_2) %>% 
  do(kde = kde(cbind(.$m.relabund, .$bio.var), compute.cont=TRUE))

levels.to.plot <- c("5%", "25%", "50%", "75%", "95%")
contours <- as.list(levels.to.plot) %>% 
  map(function(x){
    contours %>% 
      do(contour = data.frame(contourLines(x = .$kde$eval.points[[1]],
                                           y = .$kde$eval.points[[2]],
                                           z = .$kde$estimate,
                                           levels=.$kde$cont[x])[[1]])) %>% 
      .$contour %>% 
      bind_rows(.id="dim_2") %>% 
      mutate(dim_2 = contours$dim_2[as.integer(dim_2)]) %>% 
      rename(m.relabund = x, 
             bio.var = y)
  }) %>% 
  bind_rows(.id = "pct") %>% 
  mutate(pct = levels.to.plot[as.integer(pct)])

p <- d %>% 
  ggplot(aes(x=m.relabund)) +
  geom_ribbon(data = preds, aes(ymin=p2.5, ymax=p97.5), fill="grey") +
  geom_line(data = preds, aes(y=p50), color="blue") + 
  geom_path(data = contours, aes(y = bio.var, fill=pct, color = dim_2)) +
  geom_label_repel(aes(y=bio.var, label=dim_2), data = d.medians) +
  xlab("Mean CLR Transformed Relative Abundance") +
  ylab("CLR Transformed Biological Variation") +
  ggtitle("Posterior Densities Day 1 Relative Abundance and Biological Variation", 
          paste("contours at:", paste(levels.to.plot, collapse = " "))) +
  xlim(c(-3, 4))+
  theme_bw() +
  theme(legend.position="none", 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14)) +
  scale_color_manual(values = rep("black", length(unique(d$dim_2))))
ggsave(out("RelativeAbund.vs.BioVar.pdf"), plot=p, height=8, width=8, units="in")


# Find Posterior Interval for the Regression
mods <- mods %>% 
  mutate(slope = summary(mod)$coeff[2]) %>% 
  dplyr::select(-mod) %>% 
  ungroup() 

mods %>% 
  summarise(p2.5 = quantile(slope, prob=0.025), 
            mean = mean(slope), 
            p97.5 = quantile(slope, prob=0.975)) %>% 
  write.table(file=out("posterior.summary.regression.coefficients.tsv"))


p <- mods %>% 
  ggplot(aes(x=slope)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  xlab("Regression Slope") +
  ylab("Count") +
  ggtitle("Posterior Distribution of Regression Coefficients", 
          "Regression between CLR(Relative Abundance) at Day 1 and CLR(Biological Variation)") +
  theme_bw()
ggsave(out("Posterior.Regression.Coefficients.pdf"), plot=p, height=5, width=6, units="in")


# Permutation Test on Relative Abundances ---------------------------------

# Create Permuted version of D where relative abundances have been permuted.
# Relies on bio.var already being defined above. 
permute.stat <- function(){
  day0.relabund <- drop(theta_prop.starting) %>% 
    .[,,sample(1:10)] %>% # Permute Relative Abundances
    tidy_array() %>% 
    mutate(dim_3 = labels[dim_3]) %>% 
    group_by(dim_3, dim_1) %>% 
    #summarize(gm.relabund = exp(mean(log(var)))) %>% 
    summarize(m.relabund = mean(var)) %>% 
    ungroup() %>% 
    rename(dim_2 = dim_3)
  
  d <- full_join(day0.relabund, bio.var, by=c("dim_1", "dim_2")) 
  
  stat <- d %>% 
    group_by(dim_1) %>% 
    do(mod = lm(bio.var ~ m.relabund, data = .)) %>% 
    mutate(slope = summary(mod)$coeff[2]) %>% 
    ungroup() %>% 
    summarise(p2.5 = quantile(slope, prob=0.025), 
              mean = mean(slope), 
              p97.5 = quantile(slope, prob=0.975))
  
  return(stat)
}

test.statistic <- mean(mods$slope)

# For each permutation, find posterior mean for regression
null.dist <- rep(0, 1000) %>% 
  as.list() %>% 
  map(~permute.stat())

null.dist %>% 
  bind_rows() %>% 
  write.table(out("null.dist.permutation.tsv"))

# Calculate "p-value" - not really a p-value
(perm.p.value <- sum(bind_rows(null.dist)$mean <= test.statistic)/length(null.dist)) 
write.table(perm.p.value, file = out("perm.p.value.txt"))

# 3: Variation Matrix Computation ---------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))
families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))

W <- rstan::extract(fit, pars="W")$W

variation.matrix.W <- array(0, dim=c(dim(W)[1], dim(W)[2]+1, dim(W)[3]+1))
for (k in 1:dim(W)[1]){
  clr.var <- contrast.matrix.sbp %*%  W[k,,] %*% t(contrast.matrix.sbp)
  for (i in 1:dim(variation.matrix.W)[2]){
    for (j in 1:dim(variation.matrix.W)[3]){
      variation.matrix.W[k,i,j] <- clr.var[i,i] + clr.var[j,j] - 2*clr.var[i,j]
    }
  }
}


# Perform biclustering (used to organize plots)
mean.variation.matrix.W <- apply(variation.matrix.W, c(2,3), mean)
rownames(mean.variation.matrix.W) <- colnames(mean.variation.matrix.W) <- tax[labels_taxa] 
pdf(out("heatmap.mean.W.pdf"), width = 15, height=15)
heatmap(mean.variation.matrix.W, 
        main="Heatmap from Mean of Biological Variation Matrix W",
        mar=c(10, 10), 
        symm=TRUE)
dev.off()

log.mean.variation.matrix.W <- log(mean.variation.matrix.W)
min.finite <- min(log.mean.variation.matrix.W[is.finite(log.mean.variation.matrix.W)])
log.mean.variation.matrix.W[is.infinite(log.mean.variation.matrix.W)] <- min.finite
pdf(out("heatmap.mean.W.log.pdf"), width = 15, height=15)
heatmap(log.mean.variation.matrix.W, 
        main="Heatmap from Mean of Biological Variation Matrix W, Log Scale",
        mar=c(10, 10), 
        symm=TRUE)
dev.off()


V <- rstan::extract(fit, pars="V")$V

variation.matrix.V <- array(0, dim=c(dim(V)[1], dim(V)[2]+1, dim(V)[3]+1))
for (k in 1:dim(V)[1]){
  clr.var <- contrast.matrix.sbp %*%  V[k,,] %*% t(contrast.matrix.sbp)
  for (i in 1:dim(variation.matrix.V)[2]){
    for (j in 1:dim(variation.matrix.V)[3]){
      variation.matrix.V[k,i,j] <- clr.var[i,i] + clr.var[j,j] - 2*clr.var[i,j]
    }
  }
}

p <- tidy_array(variation.matrix.V) %>% 
  mutate(dim_2 = tax[labels_taxa[dim_2]], 
         dim_3 = tax[labels_taxa[dim_3]]) %>% 
  ggplot(aes(x=var)) +
  geom_histogram(bins=100) +
  scale_x_log10() +
  facet_grid(dim_2 ~ dim_3) + 
  ggtitle("Technical Variation Array")
ggsave(out("variational_array_V.pdf"),plot=p, width=14, height=14, units="in")


# 3: Biological - Ward Clustering for Principle Balances - Biological ---------

families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))

trees <- list()
for (k in 1:dim(V)[1]){
  trees[[k]] <- as.phylo(hclust(as.dist(variation.matrix.W[k,,]), "ward.D")) 
  trees[[k]]$tip.label <- families[labels_taxa[1:(dim(V)[2]+1)][trees[[k]]$tip.label]]
}

boot.tree <- list()
consensus.tree <- consensus(trees, p=0.5)
consensus.tree <- multi2di(consensus.tree)
consensus.tree <- makeNodeLabel(consensus.tree, method="number", prefix='w')

### Now plot theta interms of consensus tree

ward.sbp <- phylo2sbp(consensus.tree)
contrast.matrix.ward <- gsi.buildilrBase(ward.sbp)
colnames(contrast.matrix.ward) <- paste0("w", 1:ncol(contrast.matrix.ward))

####

# Percent of variation explained by ward basis 

# convert W to ward basis
#Just look at traces for percent of total variation explained
W.ward <- array(0, dim=dim(W))
W.ward.percent <- matrix(0, dim(W.ward)[1], dim(W.ward)[2])
W.ward.percent.cum <- matrix(0, dim(W.ward)[1], dim(W.ward)[2])
for (i in 1:dim(W)[1]){
  W.ward[i,,] <- t(contrast.matrix.ward) %*% contrast.matrix.sbp %*% W[i,,] %*% t(contrast.matrix.sbp) %*% contrast.matrix.ward 
  tmp <- diag(W.ward[i,,]) 
  tmp <- tmp[order(tmp, decreasing = TRUE)]
  W.ward.percent[i,] <- tmp/sum(tmp)
  W.ward.percent.cum[i,] <- cumsum(W.ward.percent[i,])
}

W.ward.percent %>% 
  as.data.frame() %>% 
  gather(coord, percent) %>% 
  mutate(coord = substr(coord, 2, 2)) %>% 
  ggplot(aes(x = coord, y = percent)) +
  geom_jitter() +
  xlab("Number of Coordinates") +
  ylab("Percent Biological Variation Explained") +
  ggtitle("Ward Basis") +
  theme_bw()

W.ward.percent.cum %>% 
  as.data.frame() %>% 
  gather(coord, percent) %>% 
  mutate(coord = substr(coord, 2, 2)) %>% 
  ggplot(aes(x = coord, y = percent)) +
  geom_jitter() +
  xlab("Number of Coordinates") +
  ylab("Percent Biological Variation Explained") +
  ggtitle("Ward Basis") +
  theme_bw()

W.ward.percent.cum %>% 
  as.data.frame() %>% 
  gather(coord, percent) %>% 
  mutate(coord = substr(coord, 2, 2)) %>% 
  group_by(coord) %>% 
  summarise_posterior(percent)

####


#tr <- apeBoot(consensus.tree,  prop.clades(consensus.tree, trees)) # deprecated
tr <- treeio::as.treedata(consensus.tree, boot=prop.clades(consensus.tree, trees))
p.ward.tree <- ggtree(tr)
tmp <- p.ward.tree$data %>% 
  mutate(label = paste0(label, "\n(", signif(bootstrap/length(trees), 3), ")")) %>% 
  filter(!isTip) %>% 
  .[["label"]] 
p.ward.tree$data[,"label"] <- as.character(p.ward.tree$data[,"label"])
p.ward.tree$data[!p.ward.tree$data$isTip, "label"] <-  tmp
l.tmp <- p.ward.tree$data %>% 
  filter(!isTip)
l.tmp <- setNames(l.tmp$label, paste0("w", 1:length(l.tmp$label)))
V.tmp <- contrast.matrix.ward
colnames(V.tmp) <- l.tmp[colnames(V.tmp)]
consensus.tree.tmp <- consensus.tree
consensus.tree.tmp$node.label <- l.tmp[consensus.tree.tmp$node.label]
p.ward.tree <- annotate_sbp(consensus.tree.tmp, V.tmp, p.ward.tree, sep="\n")

p.ward.tree <- p.ward.tree + 
  geom_label2(aes(label=label, subset=!isTip), size = 2) + 
  geom_tiplab(size=3) +
  xlim(c(0, 9))
ggsave(out("ward.consensus.tree.pdf"), plot=p.ward.tree, width=3, height=6, units="in")

###
theta <- rstan::extract(fit, pars="theta")$theta

d <- dim(theta)
theta_ward <- array(0, dim=c(d[1:3], ncol(contrast.matrix.ward)))
for (i in 1:d[1]){
  for (j in 1:d[3]){
    theta_ward[i,,j,] <- unclass(ilr(ilrInv(theta[i,,j,], contrast.matrix.sbp), contrast.matrix.ward))
  }
}

# Posterior Samples
tidy_theta_ward <- tidy_array(theta_ward) %>% 
  group_by(dim_2, dim_3, dim_4) %>% 
  summarize(mean = mean(var), 
            p2.5=quantile(var, prob=c(.025)), 
            p97.5=quantile(var, prob=c(.975))) %>% 
  ungroup() %>% 
  mutate(dim_2 = tt.sample[dim_2]) %>% 
  mutate(dim_4 = colnames(contrast.matrix.ward)[dim_4])

p <- tidy_theta_ward %>%
  ggplot(aes(x=dim_2, y=mean)) +
  geom_ribbon(aes(ymax=p97.5, ymin=p2.5, fill = factor(dim_3)), alpha=0.5) +
  #geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly))+
  facet_wrap(~dim_4) + 
  theme_bw() +
  #ggtitle("Ward Clustering Consensus Basis")  +
  guides(fill=guide_legend("Vessel")) +
  ylab("Balance Value") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90)) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = cbPalette) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("WardConsensus_combined.pdf"), plot=p, width=8, height=6, units="in")


p <- arrangeGrob(p.ward.tree, p, ncol=2, widths=c(1,2))
ggsave(out("ward.basis.pdf"), plot=p, width=10, height=7, units="in")


tmp.w1 <- tidy_theta_ward %>% 
  filter(dim_4 %in% c("w1"))
tmp.w2 <- tidy_theta_ward %>% 
  filter(dim_4 %in% c("w2"))
full_join(tmp.w1, tmp.w2, by=c("dim_2", "dim_3"), suffix=c(".w1", ".w2")) %>% 
  ggplot(aes(x = mean.w1, y = mean.w2)) +
  geom_path() +
  facet_wrap(~dim_3) +
  theme_bw()


# 3: Biological Make Heatmap with Text -----------------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))
families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))


od <- order_tips(consensus.tree)

tidy_variation.matrix.W <- tidy_array(variation.matrix.W) %>% 
  mutate(dim_2 = families[labels_taxa[dim_2]], 
         dim_3 = families[labels_taxa[dim_3]]) %>% 
  mutate(dim_2 = factor(dim_2, levels=od), 
         dim_3 = factor(dim_3, levels=od))

mean_tidy_variation.matrix.W <- tidy_variation.matrix.W%>% 
  group_by(dim_2, dim_3) %>% 
  summarize(p2.5 = quantile(var, prob = 0.025), 
            p50 = quantile(var, prob=0.5), 
            p97.5 = quantile(var, prob=0.975)) %>% 
  ungroup()

lt <- lower_triangle_factors(mean_tidy_variation.matrix.W$dim_2, c("dim_2", "dim_3")) %>% 
  transmute(d23 = paste0(dim_2,"_", dim_3)) %>% 
  .$d23

tidy_variation.matrix.W <- tidy_variation.matrix.W %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt)

mean_tidy_variation.matrix.W <- mean_tidy_variation.matrix.W %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt) %>% 
  mutate(dim_2 = factor(dim_2, levels = rev(levels(dim_2))))

p <- mean_tidy_variation.matrix.W %>% 
  ggplot(aes(x = dim_3, y = dim_2, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)), 
            color="white", size= 7, nudge_y = .1) +
  geom_text(aes(label = paste0("(", signif(p2.5, 1), "-",signif(p97.5, 1), ")")), 
            color="white", size=5.3, nudge_y = -.1)+
  theme_minimal() +
  # scale_fill_gradient(low = "black", high = cbPalette[7], trans="log", 
  #                     breaks = c(0.006, 0.018, 0.05)) +
  scale_fill_gradient(low = "black", high = cbPalette[2]) +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank())

# Separate Legend from Plot and Save Seperately
parts <- g_separate_legend(p)

# Save Legend
pdf(out("heatmap.text.legend.pdf"))
grid.arrange(parts$legend)
dev.off()

p <- parts$main
ggsave(out("heatmap_variational_array_text_W.pdf"),plot=p, width=14, height=14, units="in")


 # Create Zoom-in
p <-mean_tidy_variation.matrix.W %>% 
  filter(dim_2 == "Lachnospiraceae", 
         dim_3 == "Bacteroidaceae") %>% 
  ggplot(aes(x = dim_3, y = dim_2, fill = p50)) +
  geom_tile(show.legend=FALSE) +
  geom_text(aes(label = signif(p50, 2)), 
            color="white", size= 18, nudge_y = .1) +
  geom_text(aes(label = paste0("(", signif(p2.5, 2), "-",signif(p97.5, 2), ")")), 
            color="white", size=15, nudge_y = -.1)+
  theme_minimal() +
  # scale_fill_gradient(low = "black", high = cbPalette[7], trans="log", 
  #                     breaks = c(0.006, 0.018, 0.05)) +
  scale_fill_gradient(low = "black", high = cbPalette[2], 
                      limits = c(min(mean_tidy_variation.matrix.W$p50), 
                                 max(mean_tidy_variation.matrix.W$p50))) +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_blank()) +
  guides(fill=NULL)
ggsave(out("heatmap_variational_array_text_W_zoomin.pdf"),plot=p, width=6, height=6, units="in")



# 3: Biological - Make Heatmap with posteriors---------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))
families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))


od <- order_tips(consensus.tree)

tidy_variation.matrix.W <- tidy_array(variation.matrix.W) %>% 
  mutate(dim_2 = families[labels_taxa[dim_2]], 
         dim_3 = families[labels_taxa[dim_3]]) %>% 
  mutate(dim_2 = factor(dim_2, levels=od), 
         dim_3 = factor(dim_3, levels=od))

mean_tidy_variation.matrix.W <- tidy_variation.matrix.W%>% 
  group_by(dim_2, dim_3) %>% 
  summarize(p50 = quantile(var, prob=0.5)) %>% 
  ungroup()

lt <- lower_triangle_factors(mean_tidy_variation.matrix.W$dim_2, c("dim_2", "dim_3")) %>% 
  transmute(d23 = paste0(dim_2,"_", dim_3)) %>% 
  .$d23

tidy_variation.matrix.W <- tidy_variation.matrix.W %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt)

mean_tidy_variation.matrix.W <- mean_tidy_variation.matrix.W %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt) 

p <-tidy_variation.matrix.W %>% 
  ggplot() +
  geom_rect(data=mean_tidy_variation.matrix.W, 
            aes(fill = p50), 
            xmin=-Inf, xmax = +Inf, ymin=-Inf, ymax=+Inf)+
  geom_density(aes(x=var), fill = "grey") +
  scale_x_log10() +
  # scale_fill_gradient2(low=scales::muted("blue"), mid="white", high=scales::muted("red"), 
  #                      midpoint = median(mean_tidy_variation.matrix.W$p50))+
  scale_fill_gradient2(low="darkblue",mid="white", high="red", 
                       midpoint = 0.04, limits = c(min(mean_tidy_variation.matrix.W$p50), 
                                                   max(mean_tidy_variation.matrix.W$p50)))+
  facet_grid(dim_2 ~ dim_3)+
  ggtitle("Biological Variation Matrix") +
  xlab("Var(Log(X/Y))")+
  #theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

# Separate Legend from Plot and Save Seperately
parts <- g_separate_legend(p)

# Save Legend
pdf(out("heatmap.legend.pdf"))
grid.arrange(parts$legend)
dev.off()

p <- parts$main
ggsave(out("heatmap_variational_array_W.pdf"),plot=p, width=14, height=14, units="in")


# Create Zoom-in
xlimit <- c(0.001, 0.1)
ylimit <- c(0, 9)
p <-tidy_variation.matrix.W %>% 
  filter(dim_2 == "Lachnospiraceae", 
         dim_3 == "Bacteroidaceae") %>% 
  ggplot() +
  geom_rect(data=filter(mean_tidy_variation.matrix.W, 
                        dim_2 == "Lachnospiraceae", 
                        dim_3 == "Bacteroidaceae"), 
            aes(fill = p50), 
            xmin=-Inf, xmax = +Inf, ymin=-Inf, ymax=+Inf)+
  geom_density(aes(x=var), fill = "grey") +
  ylim(ylimit) +
  scale_x_log10(limits = xlimit) +
  # scale_fill_gradient2(low=scales::muted("blue"), mid="white", high=scales::muted("red"), 
  #                      midpoint = median(mean_tidy_variation.matrix.W$p50))+
  scale_fill_gradient2(low="darkblue",mid="white", high="red", 
                       midpoint = 0.04, limits = c(min(mean_tidy_variation.matrix.W$p50), 
                                                   max(mean_tidy_variation.matrix.W$p50)))+
  ggtitle("Biological Variation Matrix") +
  xlab("Var(Log(X/Y))")+
  #theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())
ggsave(out("heatmap_variational_array_W_zoomin.pdf"),plot=p, width=6, height=6, units="in")


# 3: Technical - Ward Clustering for Principle Balances -----------------------

families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))

trees <- list()
for (k in 1:dim(V)[1]){
  trees[[k]] <- as.phylo(hclust(as.dist(variation.matrix.V[k,,]), "ward.D")) 
  trees[[k]]$tip.label <- families[labels_taxa[1:(dim(V)[2]+1)][trees[[k]]$tip.label]]
}

consensus.tree <- consensus(trees, p=0.5)
consensus.tree <- multi2di(consensus.tree)
consensus.tree <- makeNodeLabel(consensus.tree, method="number", prefix='w')

### Now plot theta interms of consensus tree

ward.sbp <- phylo2sbp(consensus.tree)
contrast.matrix.ward <- gsi.buildilrBase(ward.sbp)
colnames(contrast.matrix.ward) <- paste0("w", 1:ncol(contrast.matrix.ward))

####


tr <- apeBoot(consensus.tree,  prop.clades(consensus.tree, trees))
p.ward.tree <- ggtree(tr)
p.ward.tree$data[!p.ward.tree$data$isTip, "label"] <-  p.ward.tree$data %>% 
  mutate(label = paste0(label, "\n(", signif(bootstrap/length(trees), 3), ")")) %>% 
  filter(!isTip) %>% 
  .[["label"]]
l.tmp <- p.ward.tree$data %>% 
  filter(!isTip)
l.tmp <- setNames(l.tmp$label, paste0("w", 1:length(l.tmp$label)))
V.tmp <- contrast.matrix.ward
colnames(V.tmp) <- l.tmp[colnames(V.tmp)]
consensus.tree.tmp <- consensus.tree
consensus.tree.tmp$node.label <- l.tmp[consensus.tree.tmp$node.label]
p.ward.tree <- annotate_sbp(consensus.tree.tmp, V.tmp, p.ward.tree, sep="\n")

p.ward.tree <- p.ward.tree + 
  geom_label2(aes(label=label, subset=!isTip), size = 2) + 
  geom_tiplab(size=2) +
  xlim(c(0, 12))
ggsave(out("ward.consensus.tree.V.pdf"), plot=p.ward.tree, width=3, height=6, units="in")


# 3: Technical - Make Heatmap with Text ---------------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))
families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))


od <- order_tips(consensus.tree)

tidy_variation.matrix.V <- tidy_array(variation.matrix.V) %>% 
  mutate(dim_2 = families[labels_taxa[dim_2]], 
         dim_3 = families[labels_taxa[dim_3]]) %>% 
  mutate(dim_2 = factor(dim_2, levels=od), 
         dim_3 = factor(dim_3, levels=od))

mean_tidy_variation.matrix.V <- tidy_variation.matrix.V%>% 
  group_by(dim_2, dim_3) %>% 
  summarize(p2.5 = quantile(var, prob = 0.025), 
            p50 = quantile(var, prob=0.5), 
            p97.5 = quantile(var, prob=0.975)) %>% 
  ungroup()

lt <- lower_triangle_factors(mean_tidy_variation.matrix.V$dim_2, c("dim_2", "dim_3")) %>% 
  transmute(d23 = paste0(dim_2,"_", dim_3)) %>% 
  .$d23

tidy_variation.matrix.V <- tidy_variation.matrix.V %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt)

mean_tidy_variation.matrix.V <- mean_tidy_variation.matrix.V %>% 
  filter(paste0(dim_2,"_", dim_3) %in% lt) %>% 
  mutate(dim_2 = factor(dim_2, levels = rev(levels(dim_2))))

p <- mean_tidy_variation.matrix.V %>% 
  ggplot(aes(x = dim_3, y = dim_2, fill = p50)) +
  geom_tile() +
  geom_text(aes(label = signif(p50, 2)), 
            color="white", size= 7, nudge_y = .1) +
  geom_text(aes(label = paste0("(", signif(p2.5, 1), "-",signif(p97.5, 1), ")")), 
            color="white", size=5.3, nudge_y = -.1)+
  theme_minimal() +
  # scale_fill_gradient(low = "black", high = cbPalette[7], trans="log", 
  #                     breaks = c(0.006, 0.018, 0.05)) +
  scale_fill_gradient(low = "black", high = cbPalette[2], name="Median") +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title=element_text(size=12) , 
        legend.text=element_text(size=11))

# Separate Legend from Plot and Save Seperately
parts <- g_separate_legend(p)

# Save Legend
pdf(out("heatmap.text.legend.V.pdf"))
grid.arrange(parts$legend)
dev.off()

p <- parts$main
ggsave(out("heatmap_variational_array_text_V.pdf"),plot=p, width=14, height=14, units="in")


# 3:Simplex View of Variation Array -----------------------------------------

labels_taxa <- rownames(contrast.matrix.sbp) %>% setNames(seq_along(.))
families <- as.data.frame(tax_table(ps))
families <- as.character(families$Family) %>% 
  setNames(rownames(families))

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


produce_simplex_plot <- function(subset.families, to.plot=c("W"), plot.isotropic=F){
  W.subset <- W.clr[,subset.families,subset.families]
  V.subset <- V.clr[,subset.families, subset.families]
  
  # based on page 83 of Compositional Data in R
  r <- sqrt(qchisq(p=0.95, df=2)) 
  
  neutral.element <- acomp(c(1,1,1))
  names(neutral.element) <- subset.families
  
  plot(neutral.element, pch=16, col="white", axes=TRUE)
  for (i in sample(1:dim(W)[1], size=100)){
    if (plot.isotropic){
      ellipses(mean=neutral.element,
               var= diag(rep(sum(diag(W.subset[i,,]))/3, 3)),
               r=r,
               col=scales::alpha(cbPalette[8], 0.3), 
               lwd=.3)
    }
    if ("W" %in% to.plot){
      ellipses(mean=neutral.element, 
               var=W.subset[i,,], 
               r=r, 
               col=scales::alpha(cbPalette[6], 0.3), 
               lwd=.3)      
    }
    if ("V" %in% to.plot){
      ellipses(mean=neutral.element,
               var=V.subset[i,,],
               r=r,
               col=scales::alpha(cbPalette[7], 0.3), 
               lwd=.3)
    }
    if ("both" %in% to.plot){
      ellipses(mean=neutral.element,
               var=V.subset[i,,] + W.subset[i,,],
               r=r,
               col=scales::alpha(cbPalette[8], 0.3), 
               lwd=.3)
    }
  } 
}


pdf(out("simplex_2_10_9.pdf"), width = 5, height=5)
subset.families <- c("seq_2", "seq_10", "seq_9")
produce_simplex_plot(subset.families, to.plot="W", plot.isotropic=T)
dev.off()

pdf(out("simplex_2_10_9.both.pdf"), width=5, height=5)
subset.families <- c("seq_2", "seq_10", "seq_9")
produce_simplex_plot(subset.families, to.plot=c("W", "V","both"))
dev.off()

pdf(out("simplex_2_10_4.both.pdf"), width=5, height=5)
subset.families <- c("seq_2", "seq_10", "seq_4")
produce_simplex_plot(subset.families, to.plot=c("W", "V","both"))
dev.off()


# 5:Special Balances to Investigate -----------------------------------------
# Figure 5B-D

n.rikenellacae <- rep(-1, nrow(sbp))
names(n.rikenellacae) <- rownames(sbp)
n.oral <- n.rikenellacae
n.bacteroidaceae <- n.rikenellacae
n.enterobacteriaceae <- n.rikenellacae

n.rikenellacae["seq_1"] <- 1
n.oral[c("seq_4", "seq_6")] <- 1
n.bacteroidaceae["seq_2"] <- 1
n.enterobacteriaceae["seq_9"] <- 1
n.lachnospiraceae.bacteroidaceae <- rep(0, nrow(sbp))
names(n.lachnospiraceae.bacteroidaceae) <- rownames(sbp)
n.lachnospiraceae.bacteroidaceae["seq_2"] <- 1
n.lachnospiraceae.bacteroidaceae["seq_10"] <- -1

sbp.special <- cbind(n.rikenellacae, n.oral, n.bacteroidaceae, n.lachnospiraceae.bacteroidaceae, 
                     n.enterobacteriaceae)
contrast.matrix.special <- gsi.buildilrBase(sbp.special)


theta <- rstan::extract(fit, pars="theta")$theta

d <- dim(theta)
theta_special <- array(0, dim=c(d[1:3], ncol(contrast.matrix.special)))
for (i in 1:d[1]){
  for (j in 1:d[3]){
    theta_special[i,,j,] <- unclass(ilr(ilrInv(theta[i,,j,], contrast.matrix.sbp), contrast.matrix.special))
  }
}

# Posterior Samples
tidy_theta_special <- tidy_array(theta_special) %>% 
  group_by(dim_2, dim_3, dim_4) %>% 
  summarize(mean = mean(var), 
            p2.5=quantile(var, prob=c(.025)), 
            p97.5=quantile(var, prob=c(.975))) %>% 
  ungroup() %>% 
  mutate(dim_2 = tt.sample[dim_2]) %>% 
  mutate(dim_4 = colnames(contrast.matrix.special)[dim_4])

# Which days was nadir of n.oral?
tidy_theta_special %>% 
  filter(dim_4 == "n.oral") %>% 
  filter(dim_3 == 1) # Just look at first vessel to make it easy to see.  
  
  
tt.starvation <- c(ymd_h("2015-12-5-15", tz = "UTC") - days(23), ymd_h("2015-12-7-15", tz = "UTC")-days(23))
p <- tidy_theta_special %>% 
  filter(dim_4 == "n.rikenellacae") %>% 
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_rect(xmin=as.integer(tt.starvation[1]), xmax = as.integer(tt.starvation[2]), 
            ymin=-Inf, ymax=Inf, fill="lightgrey", alpha=1) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_blank()) +
  ggtitle("Rikenellacae / Remaining Taxa", "Posterior 95% credible interval") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Vessel"), fill=guide_legend(title="Vessel")) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.rikenellacae.combined.pdf"), plot=p, width=8, height=4, units="in")



## combined plot
tmp.rect <- data.frame(xmin=as.integer(tt.starvation[1]), xmax = as.integer(tt.starvation[2]), 
                       ymin=-Inf, ymax=Inf, dim_4 = "n.rikenellacae", 
                       mean=1, dim_3= 1 )
p <- tidy_theta_special %>% 
  filter(dim_4 %in% c("n.rikenellacae", "n.oral", "n.enterobacteriaceae")) %>% 
  mutate(dim_4 = factor(dim_4, levels = c("n.rikenellacae", "n.oral", "n.enterobacteriaceae"))) %>% 
  ggplot() +
  geom_rect(data = tmp.rect, xmin=as.integer(tt.starvation[1]), xmax = as.integer(tt.starvation[2]), 
            ymin=-Inf, ymax=Inf, fill="lightgrey", alpha=1) +
  # geom_rect(data = tmp.rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            # fill = "lightgrey", alpha=1) +
  geom_ribbon(aes(x=dim_2, group = dim_3, ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.6) +
  facet_wrap(~dim_4, nrow=3, strip.position = "top", scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_blank()) +
  ggtitle("Posterior 95% credible interval") +
  ylab("Balance Value") +
  guides(color=guide_legend(title="Vessel"), fill=guide_legend(title="Vessel")) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = cbPalette) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("special.balances.combined.pdf"), plot=p, width=6, height=8, units="in")

##

p <- tidy_theta_special %>% 
  filter(dim_4 == "n.oral") %>% 
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  #geom_line(aes(color=factor(dim_3))) +
  #facet_grid(dim_3~., scales="free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_blank()) +
  ggtitle("Fusobacteriaceae + Synergistaceae / Remaining Taxa", "Posterior 95% credible interval") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Vessel"), fill=guide_legend(title="Vessel")) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.oral.combined.pdf"), plot=p, width=8, height=4, units="in")


p <- tidy_theta_special %>% 
  filter(dim_4 == "n.enterobacteriaceae") %>% 
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_blank()) +
  ggtitle("Enterobacteriaceae / Remaining Taxa", "Posterior 95% credible interval") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Vessel"), fill=guide_legend(title="Vessel")) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.enterobacteriaceae.combined.pdf"), plot=p, width=8, height=4, units="in")

# Media Changes and person sampling dates defined above - dates.changemedia
tmp <- select(sample_data(ps), time, Vessel, Researcher) %>% 
  mutate(Vessel = as.integer(Vessel))
p <- tidy_theta_special %>%
  left_join(tmp, 
            by=c("dim_2" = "time", "dim_3" = "Vessel")) %>% 
  mutate(Researcher = factor(as.integer(factor(Researcher)))) %>% # Removing Initials
  filter(dim_4 == "n.enterobacteriaceae") %>%
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  geom_line(aes(color = Researcher), size = 1) +
  geom_vline(aes(xintercept=as.integer(t)), data=dates.changemedia) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ggtitle("Enterobacteriaceae / Remaining Taxa", "Posterior 95% credible interval With Media Changes Indicated") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Researcher"), fill=guide_legend(title="Vessel")) +
  scale_fill_manual(values = cbPalette) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.enterobacteriaceae.environmental_variables.combined.pdf"), plot=p, width=8, height=5, units="in")

# Sequencing batch?
tmp <- select(sample_data(ps), time, Vessel, batch) %>% 
  mutate(Vessel = as.integer(Vessel))
p <- tidy_theta_special %>%
  left_join(tmp, 
            by=c("dim_2" = "time", "dim_3" = "Vessel")) %>% 
  mutate(batch = factor(as.integer(factor(batch)))) %>% # Removing Initials
  filter(dim_4 == "n.enterobacteriaceae") %>%
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  geom_line(aes(color = batch), size = 1) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ggtitle("Enterobacteriaceae / Remaining Taxa", "Posterior 95% credible interval With Media Changes Indicated") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Sequencing Batch"), fill=guide_legend(title="Vessel")) +
  scale_fill_manual(values = cbPalette) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.enterobacteriaceae.sequencing_batch.combined.pdf"), plot=p, width=8, height=5, units="in")


# Check for pertinent runnotes
tmp <- select(sample_data(ps), time, Vessel, Note) %>% 
  mutate(Vessel = as.integer(Vessel))
p <- tidy_theta_special %>%
  left_join(tmp, 
            by=c("dim_2" = "time", "dim_3" = "Vessel")) %>% 
  filter(dim_4 == "n.enterobacteriaceae") %>%
  ggplot(aes(x=dim_2, y=mean, group=dim_3)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5, fill=factor(dim_3)), alpha=0.5) +
  geom_line(aes(color = Note), size = 1) +
  geom_vline(aes(xintercept=as.integer(t)), data=dates.changemedia) +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ggtitle("Enterobacteriaceae / Remaining Taxa", "Posterior 95% credible interval With Media Changes Indicated") +
  xlab("") + ylab("Balance Value") +
  guides(color=guide_legend(title="Vessel"), fill=guide_legend(title="Vessel")) +
  scale_fill_manual(values = cbPalette) +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.enterobacteriaceae.runnotes.combined.pdf"), plot=p, width=8, height=5, units="in")


p <- tidy_theta_special %>% 
  filter(dim_2 %in% tt.hourly) %>% 
  filter(dim_4 == "n.bacteroidaceae") %>% 
  ggplot(aes(x=dim_2, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly)) +
  facet_grid(dim_3~., scales="free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ggtitle("Bacteroidaceae / Remaining Taxa", "Posterior mean and 95% credible interval") +
  xlab("") + ylab("Balance Value") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.bacteroidaceae.hourly.pdf"), plot=p, width=8, height=5, units="in")


p <- tidy_theta_special %>% 
  #filter(dim_2 %in% tt.hourly) %>% 
  filter(dim_4 == "n.lachnospiraceae.bacteroidaceae") %>% 
  ggplot(aes(x=dim_2, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.7) +
  geom_line(color="blue") +
  geom_vline(xintercept=as.integer(t.perturb.hourly)) +
  facet_grid(dim_3~., scales="free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0)) +
  ggtitle(" Bacteroidaceae / Lachnospiraceae", "Posterior mean and 95% credible interval") +
  xlab("") + ylab("Balance Value") +
  scale_x_datetime(date_labels="Day %d")
ggsave(out("n.bacteroidaceae.lachnospiraceae.combined.pdf"), plot=p, width=8, height=5, units="in")



# Percent of variation explained by "special balances" (as defined above) -----

#specifically the following above created balances 
contrast.special3 <- contrast.matrix.special[,c("n.rikenellacae", 
                                                "n.oral", 
                                                "n.enterobacteriaceae")]

# Orthogonalize with QR to compare total variation of these three
# to total variation of philr basis (note total variation is 
# a matrix invariant to this comparison is appropriate)
contrast.special3 <- qr.Q(qr(contrast.special3))
W <- rstan::extract(fit, pars="W")$W
totvar.all <- totvar.special3 <- rep(0, dim(W)[1])
for (i in 1:dim(W)[1]){
 totvar.all[i] <- sum(diag(W[i,,]))
 W.special3 <- ilrvar2ilrvar(W[i,,], contrast.matrix.sbp, contrast.special3)
 totvar.special3[i] <- sum(diag(W.special3))
}

enframe(totvar.special3/totvar.all) %>% 
  summarise_posterior(value)

# Comparing Overlap between V and W by Distance ---------------------------

split.along.dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

W <- rstan::extract(fit, pars="W")$W
V <- rstan::extract(fit, pars="V")$V

# Just use 500 samples from the posterior distribution of W and V. 
samsize <- 500
sam <- sample(1:dim(W)[1], size = samsize)

WV <- c(split.along.dim(W, 1)[sam], split.along.dim(V,1)[sam])
WV <- map(WV, cov2cor) # convert to correlation matricies to compare shape only
labs <- c(rep("W", samsize), rep("V", samsize))

# Calcualte distance matrix between matricies
dist.wv <- matrix(0, length(WV), length(WV))
ntmp <- length(WV)
for (i in 1:ntmp){
  for (j in 1:i){
    dist.wv[i,j] <- distcov(WV[[i]],WV[[j]], method = "Riemannian")
    if (i %% 100 == 0 && j == 1) print(i)
  }
}
dist.wv <- dist.wv + t(dist.wv)
write.table(dist.wv, file = out("dist.wv.csv"))

# Calculate blocking totals given labeling
# d is a distance matrix 
# l is a vector of categorical labels (2 groups)
blocks <- function(d, l){
  o <- matrix(0, 2,2 ) # Matrix of block sums 
  lw <- l=="W"
  lv <- l=="V"
  o[1,1] <- sum(dist.wv[lw, lw])  ## WW block
  o[2,1] <- o[1,2] <- sum(dist.wv[lw, lv])  ## WV and VW blocks
  o[2,2] <- sum(dist.wv[lv,lv])
  o
}

# Create Function to take block output and make statistic out of it
# Statitic will be the ratio between the sum of the diagonal blocks
# and the sum of the off-diagonal blocks
# Low stat means that most distance is between the groups
# High stat means most distance is within the groups
blockstat <- function(b){
  (b[1,1]+b[2,2])/(2*b[1,2])
}

# Calculate test statistic
teststat <- blocks(dist.wv, labs)
teststat <- blockstat(teststat)

# Create Permutation scheme 
nperm <- 1000 
null.dist <- rep(0, nperm)
for (i in 1:nperm){
  b <- blocks(dist.wv, sample(labs)) # build blocks off of permuted labels
  null.dist[i] <- blockstat(b)
}


# Compare Test statistic to null distribution
# result represents the posterior probability 
# that the two distributions are the same under this 
# statistic and permutation scheme
ecdf(null.dist)(teststat)

# Create plot
logit <- function(x){ # helper function x is a vector
  log(x/(1-x))
}

p <- ggplot(enframe(logit(null.dist)), aes(x=value)) + 
  geom_density(fill="grey") +
  geom_vline(xintercept=logit(teststat), color="red") +
  xlim(c(-1, 10)) +
  xlab("Logit(Sum of Diagonal Blocks / Sum of Off Diagonal Blocks)") +
  theme_bw() +
  ggtitle("Comparing Correlation Matricies W and V", "Using Riemmanian Distance")
ggsave(out("comparing_corr_riemmanian_distance.pdf"), plot=p, height=3, width=5, units="in")


# PhILR PCoA  -------------------------------------------------------------

# Create a version of Y that is easier to work with for PCoA
Yperm <- aperm(Y, c(3, 2, 1))
d <- dim(Yperm)
meta <- data.frame(R=rep(1:d[2], times=d[3]), 
                   tt=rep(1:d[3], each=d[2]))
Yperm <- t(matrix(Yperm, d[1], prod(d[-1])))
colnames(Yperm) <- taxa_names(ps)

# Convert to philr basis with pseudocount
Ypp <- philr(Yperm+0.65, phy_tree(ps))

# Calculate distance (aitchison distance)
daitch <- as.matrix(dist(Ypp))

# Calculate pcoa decomposition
decomp <- ape::pcoa(daitch)

# Get Batch info
tmp <- select(sample_data(ps), time, Vessel, batch) %>% 
  mutate(R = as.integer(Vessel)) %>% 
  rename(Time=time) %>% 
  mutate( R = paste("Vessel", R))

# Combine data together prior for plotting
meta <- meta %>% 
  cbind(decomp$vectors[,1:2]) %>% 
  mutate(replicate=duplicated(tt.observed)[tt], 
         observed=Y.obs[cbind(tt, R)], 
         Time=tt.observed[tt], 
         R = paste("Vessel", R)) %>% 
  filter(observed==TRUE) %>%
  left_join(tmp, by=c("Time", "R")) %>% 
  mutate(Batch=factor(batch))

meta.norep <- filter(meta, replicate==FALSE)
meta.rep <- filter(meta, replicate==TRUE)
  
scale_colour_datetime <- function (..., low = "#132B43", high = "#56B1F7", mid="blue", space = "Lab", 
                                   na.value = "grey50", guide = "colourbar") 
{
  ggplot2:::datetime_scale("colour", "time", 
                           palette = scales::div_gradient_pal(low,mid, 
                                                              high, space), 
                           na.value = na.value, guide = guide, ...)
}

p <- meta.norep %>% 
  ggplot(aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Time)) +
  geom_point(data=meta.rep) +
  facet_wrap(~R) +
  theme_bw() +
  xlab(paste0("Axis 1 [", 100*signif(decomp$values$Relative_eig[1], 3), "%]")) +
  ylab(paste0("Axis 2 [", 100*signif(decomp$values$Relative_eig[2], 3), "%]")) +
  scale_colour_datetime(date_labels="Day %d", low="green", high="black", mid="red", 
                        rescaler=ggplot2:::mid_rescaler(mid=as.integer(ymd_hms("2015-11-16 15:00:00"))))
ggsave(out("pcoa_aitchison_time.pdf"), plot=p, width=7.5, height=5, units="in")


p <- meta.norep %>% 
  mutate(Batch=factor(batch)) %>% 
  ggplot(aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Batch)) +
  geom_point(data=meta.rep, aes(color=Batch)) +
  facet_wrap(~R) +
  theme_bw() +
  xlab(paste0("Axis 1 [", 100*signif(decomp$values$Relative_eig[1], 3), "%]")) +
  ylab(paste0("Axis 2 [", 100*signif(decomp$values$Relative_eig[2], 3), "%]")) +
  scale_color_brewer(palette = "Set1")
ggsave(out("pcoa_aitchison_batch.pdf"), plot=p, width=7, height=5, units="in")


p <- meta.norep %>% 
  mutate(Batch=factor(batch)) %>% 
  mutate(R = str_extract(R, "[:digit:]")) %>% 
  mutate(Vessel=R) %>% 
  ggplot(aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Vessel)) +
  xlab(paste0("Axis 1 [", 100*signif(decomp$values$Relative_eig[1], 3), "%]")) +
  ylab(paste0("Axis 2 [", 100*signif(decomp$values$Relative_eig[2], 3), "%]")) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") 
ggsave(out("pcoa_aitchison_vessels.pdf"), plot=p, width=7, height=3, units="in")


# Original Computing Environment ------------------------------------------

devtools::session_info()

