library(tidyverse)

machine <- "hardac"
paths <- list()
paths[["local"]] <- list("results" = "~/Research/mdlm/results/2018-03-21_publication_2/",
                        "mapping" = '~/Research/data/_data_raw/sequencing.2016.03.04/2016.03.25MappingFile.MergedPool.txt',
                        "phyloseq"= "~/Research/data/_data_derived/sequencing.2016.03.04/dada2_2016_10_18/hardac/phyloseq.rds")
paths[["hardac"]] <- list("base"="/data/davidlab/users/jds/mdlm_publication/",
                          "results" = "/data/davidlab/users/jds/mdlm_publication_2/mdlm_toblerone/",
                          "mapping" = '/data/davidlab/users/jds/mdlm_publication_2/dada2/0_mapping/2016.03.25MappingFile.MergedPool.txt',
                          "phyloseq" = "/data/davidlab/users/jds/mdlm_publication_2/dada2/phyloseq.rds", 
                          "runnotes" = '/data/davidlab/users/jds/mdlm_publication_2/runnotes/2017.10.26_ResearcherIdentity.csv')

setwd(paths[[machine]]$results)
source("utils.R")

# Read Params
params.all <- read.table('sensitivity_analysis_params.txt', stringsAsFactors = FALSE)

# Helper Function for File paths
out <- function(filename,i){
  out.dir <- paste0(paths[[machine]]$results, i, "/")
  paste0(out.dir, filename)
}



# Colnames for Parameters -------------------------------------------------
params <- params.all

param.names <- c("Base")
for (i in 2:nrow(params)){
  p <- colnames(params)[!(params[1,] == params[i,])]
  v <- params[i,][!(params[1,] == params[i,])]
  p <- paste0(p, ": ",v)
  param.names <- c(param.names, p)
}

param.names.expr <- c("Base")
for (i in 2:nrow(params)){
  p <- colnames(params)[!(params[1,] == params[i,])]
  #p <- params.exprs[p]
  v <- params[i,][!(params[1,] == params[i,])]
  v <- paste0(": ", unname(v))
  if (p == "V_lkj"){
    p <- bquote(zeta^V~.(v))
  } else if (p == "W_lkj") {
    p <- bquote(zeta^W~.(v))
  } else if (p == "V_scale_var") {
    p <- bquote(xi^V[i]~.(v)) 
  } else if (p == "W_scale_var"){
    p <- bquote(xi^W[i]~.(v))
  } else if (p == "V_scale_mean"){
    p <- bquote(tau^V[i]~.(v))
  } else if (p == "W_scale_mean"){
    p <- bquote(tau^V[i]~.(v))
  } else if (p == "c0") {
    p <- bquote(C[0]~.(v))
  }
  param.names.expr <- c(param.names.expr, p)
}
names(param.names.expr) <- param.names


# Percent Bio Noise Figure ------------------------------------------------

bio.noise <- matrix(0, nrow(params), 5)
for (i in 1:nrow(params)){
  #for (i in 1:1){
  tmp <- t(read.table(out("percent_bio_noise_summary.txt", i)))
  bio.noise[i,] <- tmp
}
colnames(bio.noise) <- colnames(tmp)

p <- bio.noise %>% 
  as.data.frame() %>% 
  rownames_to_column(var="Parameter") %>% 
  mutate(Parameter = factor(Parameter, levels=as.character(1:nrow(params)))) %>% 
  mutate(Parameter = factor(param.names[Parameter], level=param.names)) %>% 
  ggplot(aes(y = Parameter, x = `50%`, xmin = `2.5%`, xmax = `97.5%`)) +
  geom_errorbarh() +
  geom_point() +
  xlab("Percentage of Total Variation Attributable to Biological Sources")+
  theme_bw() +
  scale_y_discrete(labels=as.expression(param.names.expr))
ggsave("bio.noise.sensitivity.pdf", plot=p, height=6, width=7, units="in")


# Linear Regression Figure ------------------------------------------------


lin.reg <- matrix(0, nrow(params), 4)
for (i in 1:nrow(params)){
  tmp <- as.matrix(read.table(out("posterior.summary.regression.coefficients.tsv", i)))
  p.val <- c(read.table(out("perm.p.value.txt", i)))$x
  lin.reg[i,1:3] <- tmp
  lin.reg[i,4] <- p.val # note not realy a p-value (see main text)
}
colnames(lin.reg) <- c(colnames(tmp), "p.val")

p <- lin.reg %>% 
  as.data.frame() %>% 
  rownames_to_column(var="Parameter") %>% 
  mutate(Parameter = factor(Parameter, levels=as.character(1:nrow(params)))) %>% 
  mutate(Parameter = factor(param.names[Parameter], level=param.names)) %>% 
  ggplot(aes(y = Parameter, x = mean, xmin = p2.5, xmax = p97.5)) +
  geom_errorbarh() +
  geom_point() +
  #geom_text(nudge_y = .35, size = 4) +
  xlab("Regression Slope") +
  theme_bw()+
  scale_y_discrete(labels=as.expression(param.names.expr))
ggsave("lin.reg.sensitivity.pdf", plot=p, height=6, width=7, units="in")
