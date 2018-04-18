# Assumes tree's internal nodes and tips are named
# can pass smaller contrast matrix to subset which are annotated
# tip names on tree and V matrix must align if V pre-specified
annotate_sbp <- function(tr, V, p=NULL, sep="\n", orientation="vertical"){
  if (!setequal(tr$tip.label, rownames(V))) stop("mismatch between tip.label of tree and rownames of V")
  if (!setequal(tr$node.label, colnames(V))) stop("mismatch between node.label of tree and colnames of V")
  if (is.null(p)) p <- ggtree(tr)
  if (orientation =="vertical"){
    dec <- TRUE
  } else if (orientation=="horizontal") {
    dec <- FALSE
  }
  d <- p$data
  n.tip <- ape::Ntip(tr)
  n.node <- ape::Nnode(tr)
  n.numbers <- (n.tip+1):(n.node+n.tip)
  children <- phangorn::Children(tr, n.numbers)
  children <- lapply(children, function(x) x[order(c(d[x,"y"]),decreasing=dec)])
  names(children) <- tr$node.label
  V.sign <- sign(V)
  tips <- phangorn::Descendants(tr, 1:nrow(d), type="tips")
  tips <- lapply(tips, function(x) tr$tip.label[x])

  l <- list()
  for (n in names(children)){
    signs <- sapply(children[[n]], function(x) unique(sign(V[tips[[x]], n])))
    signs <- ifelse(signs==1, "+", "-")
    l[[n]] <- paste(signs[1], n, signs[2], sep=sep)
  }
  l <- unlist(l)
  d.order <- d$label[d$label %in% names(l)]
  d$label[d$label %in% names(l)] <- l[d.order]
  p$data <- d
  return(p)
}


tidy_array <- function(a){
  d <- dim(a)
  l <- list()
  for (i in 1:length(d)){
    l[[paste0('dim_',i)]] <- 1:d[i]
  }
  tidy <- expand.grid(l)
  tidy[["var"]] <- a[as.matrix(tidy)]
  tidy[["parameter"]] <- lazyeval::expr_text(a)
  
  return(tidy)
}

array_any_cases <- function(a, margin=1){
  return(apply(a, margin, function(x) any(!is.na(x))))
}

# a is vector of factors
lower_triangle_factors <- function(a, VarNames=c("Var1", "Var2")){
  if(!is.factor(a)) a <- factor(a)
  l <- levels(a)
  z <- matrix(0, length(l), length(l))
  r <- row(z)[lower.tri(z)]
  c <- col(z)[lower.tri(z)]
  
  idx <- cbind(l[r],l[c])
  colnames(idx) <- VarNames
  idx <- as.data.frame(idx)
  if (is.factor(a)){
    idx[[VarNames[1]]] <- factor(idx[[VarNames[1]]], levels=l)
    idx[[VarNames[2]]] <- factor(idx[[VarNames[2]]], levels=l)
  }
  idx
}

# Order of Phylogenetic Tips on a tree
order_tips <- function(tree){
  d=fortify(tree)
  dd = subset(d, isTip)
  dd$label[order(dd$y, decreasing=TRUE)]
}

# Extract legend from a ggplot plot
g_separate_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  main <- a.gplot + guides(color='none')
  list(legend=legend,
       main=main)
}

bool_ordered <- function(x){
  all.equal(x, x[order(x)])
} 

# coda functions ----------------------------------------------------------

ilr_array <- function(x, V){ # Always works on first dimention as margin
  y <- array(0, dim = c(dim(x)[1], dim(x)[2], dim(x)[3]-1))
  for (i in 1:dim(x)[1]){
    y[i,,] <- unclass(ilr(x[i,,], V))
  }
  return(y)
}
ilrInv_array <- function(y, V){ # Always works on first dimention as margin
  x <- array(0, dim = c(dim(y)[1], dim(y)[2], dim(y)[3]+1))
  for (i in 1:dim(y)[1]){
    x[i,,] <- unclass(ilrInv(y[i,,], V))
  }
  return(x)
}
clr_array <- function(x, V){ # Always works on first dimention as margin
  y <- array(0, dim = c(dim(x)[1], dim(x)[2], dim(x)[3]))
  for (i in 1:dim(x)[1]){
    y[i,,] <- unclass(clr(x[i,,], V))
  }
  return(y)
}
