#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ape); library(nlme); library(phytools)})

argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(key,default=NULL){i<-match(paste0("--",key),argv); if(is.na(i)) default else argv[[i+1L]]}
envfile <- getarg("environment"); treefile <- getarg("tree"); outdir <- getarg("out-dir")
min_n <- as.integer(getarg("min-species-n","5"))
if(is.null(envfile)||is.null(outdir)) stop("Required: --environment and --out-dir")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

traits <- c(
  "orientation_angle_degrees_median","corolla_lab_lightness_median",
  "corolla_lab_chroma_median","corolla_hue_sin_median","corolla_hue_cos_median",
  "shape_aspect_ratio_median","shape_circularity_median","shape_solidity_median",
  "shape_width_cv_median"
)
env <- c("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15")

if(is.null(treefile) || !file.exists(treefile)){
  write.csv(data.frame(status="skipped",reason="No Newick tree supplied"),file.path(outdir,"phylogenetic_status.csv"),row.names=FALSE)
  quit(save="no",status=0)
}

df <- read.csv(envfile,check.names=FALSE)
need <- c("taxon_name",traits,env); miss <- setdiff(need,names(df)); if(length(miss)) stop("Missing: ",paste(miss,collapse=", "))
for(v in c(traits,env)) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
keep <- names(which(table(df$taxon_name)>=min_n)); df <- df[df$taxon_name%in%keep,,drop=FALSE]
sp <- aggregate(df[c(traits,env)],list(taxon_name=df$taxon_name),median,na.rm=TRUE)
rownames(sp) <- sp$taxon_name

tree <- read.tree(treefile)
normalize_name <- function(x) gsub("[^A-Za-z0-9_]","_",gsub(" ","_",trimws(x)))
tree$tip.label <- normalize_name(tree$tip.label)
sp$tree_name <- normalize_name(sp$taxon_name)
common <- intersect(tree$tip.label,sp$tree_name)
if(length(common)<20) stop("Fewer than 20 species overlap between tree and data: ",length(common))
tree <- drop.tip(tree,setdiff(tree$tip.label,common))
sp <- sp[match(tree$tip.label,sp$tree_name),,drop=FALSE]

lambda_rows <- list(); pgls_rows <- list()
for(trait in traits){
  x <- sp[[trait]]; names(x) <- tree$tip.label
  ok <- is.finite(x)
  if(sum(ok)<20) next
  tr <- drop.tip(tree,tree$tip.label[!ok]); xx <- x[ok]
  lam <- tryCatch(phylosig(tr,xx,method="lambda",test=TRUE),error=function(e)e)
  if(!inherits(lam,"error")) lambda_rows[[trait]] <- data.frame(trait=trait,n_species=length(xx),lambda=unname(lam$lambda),logL=unname(lam$logL),p_value=unname(lam$P))
  for(pred in env){
    dat <- data.frame(y=sp[[trait]],x=sp[[pred]],species=tree$tip.label)
    dat <- dat[is.finite(dat$y)&is.finite(dat$x),,drop=FALSE]
    if(nrow(dat)<20) next
    tr2 <- drop.tip(tree,setdiff(tree$tip.label,dat$species))
    dat <- dat[match(tr2$tip.label,dat$species),,drop=FALSE]
    dat$y <- as.numeric(scale(dat$y)); dat$x <- as.numeric(scale(dat$x)); rownames(dat)<-dat$species
    fit <- tryCatch(gls(y~x,data=dat,correlation=corPagel(0.5,phy=tr2,fixed=FALSE),method="ML"),error=function(e)e)
    if(inherits(fit,"error")) next
    tt <- summary(fit)$tTable
    lamfit <- coef(fit$modelStruct$corStruct,unconstrained=FALSE)
    pgls_rows[[paste(trait,pred)]] <- data.frame(trait=trait,predictor=pred,n_species=nrow(dat),beta=tt["x","Value"],se=tt["x","Std.Error"],t=tt["x","t-value"],p_value=tt["x","p-value"],lambda=unname(lamfit))
  }
}
if(length(lambda_rows)) write.csv(do.call(rbind,lambda_rows),file.path(outdir,"phylogenetic_signal_pagel_lambda.csv"),row.names=FALSE)
if(length(pgls_rows)){
  z <- do.call(rbind,pgls_rows); z$q_bh <- p.adjust(z$p_value,method="BH")
  write.csv(z,file.path(outdir,"phylogenetic_gls_trait_environment.csv"),row.names=FALSE)
}
write.csv(data.frame(status="ok",species_overlap=length(common),tree_tips=length(tree$tip.label)),file.path(outdir,"phylogenetic_status.csv"),row.names=FALSE)
