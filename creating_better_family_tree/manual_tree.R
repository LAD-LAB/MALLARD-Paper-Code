library(ape)
library(ggtree)

setwd("~/Research/mdlm/results/2017-03-20_csme_dlm_combo/creating_better_family_tree/")


# Based on First 1000 cultured Gastrointestinal microbes, FEMS microbiology review 2014

tree_txt<- c("(6,(((25, 15)Negativicutes,((34 ,(73, 30)Bacilli), (117, (35, (5, (10, 54))))Clostridia))Firmicutes,((1, (2,14))Bacteroides, ((23, (47, 9))Proteo, 4))));")
tree <- read.tree(text=tree_txt)
tree$tip.label <- paste0("seq_", tree$tip.label)


ggtree(tree) +
  geom_label2(aes(label=label))

write.tree(tree, file="manual_families.tree")
