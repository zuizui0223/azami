#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(INLA); library(Matrix); library(sf)})

traits <- c(
  "orientation_angle_degrees_median","corolla_lab_lightness_median",
  "corolla_lab_chroma_median","corolla_hue_sin_median",
  "corolla_hue_cos_median","shape_aspect_ratio_median",
  "shape_circularity_median","shape_solidity_median",
  "shape_width_cv_median"
)
climate <- c("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15")
topography <- c("topo_elevation","topo_slope","topo_roughness")
soil <- paste0("soil_",c("bdod","cec","cfvo","clay","sand","silt","nitrogen","phh2o","soc","ocd"),"_0_30cm")
requested <- c(climate,topography,soil)
model_groups <- list(
  climate=climate,
  climate_topography=c(climate,topography),
  climate_soil=c(climate,soil),
  full=c(climate,topography,soil)
)

argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(key,default=NULL){i<-match(paste0("--",key),argv); if(is.na(i)) default else argv[[i+1L]]}
infile <- getarg("environment"); outdir <- getarg("out-dir")
if(is.null(infile)||is.null(outdir)) stop("Required: --environment and --out-dir")
min_species_n <- as.integer(getarg("min-species-n","5"))
cor_cutoff <- as.numeric(getarg("cor-cutoff","0.85"))
max_edge <- as.numeric(getarg("mesh-max-edge-km","1000"))
cutoff <- as.numeric(getarg("mesh-cutoff-km","100"))
offset <- as.numeric(getarg("mesh-offset-km","1500"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

df <- read.csv(infile,check.names=FALSE)
need <- c("taxon_name","latitude","longitude",traits)
miss <- setdiff(need,names(df)); if(length(miss)) stop("Missing required columns: ",paste(miss,collapse=", "))
for(v in intersect(c("latitude","longitude",traits,requested),names(df))) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
df <- df[is.finite(df$latitude)&is.finite(df$longitude),,drop=FALSE]
keep_species <- names(which(table(df$taxon_name)>=min_species_n)); df <- df[df$taxon_name%in%keep_species,,drop=FALSE]

available <- requested[requested%in%names(df)]
available <- available[vapply(available,function(p) sum(is.finite(df[[p]]))>=max(100,0.8*nrow(df)) && sd(df[[p]],na.rm=TRUE)>0,logical(1))]
if(length(intersect(climate,available))<2) stop("Fewer than two covered climate predictors")

Xall <- matrix(NA_real_,nrow(df),length(available),dimnames=list(NULL,available))
for(j in seq_along(available)){
  v <- df[[available[j]]]
  centered <- v-ave(v,df$taxon_name,FUN=function(z) mean(z,na.rm=TRUE))
  s <- sd(centered,na.rm=TRUE)
  if(is.finite(s)&&s>0) Xall[,j] <- centered/s
}
usable <- colnames(Xall)[vapply(seq_len(ncol(Xall)),function(j) sum(is.finite(Xall[,j]))>=max(100,0.8*nrow(Xall)),logical(1))]
Xall <- Xall[,usable,drop=FALSE]

prune_group <- function(candidates){
  candidates <- intersect(candidates,colnames(Xall))
  if(!length(candidates)) return(character())
  xx <- Xall[,candidates,drop=FALSE]
  C <- suppressWarnings(cor(xx,use="pairwise.complete.obs"))
  selected <- character()
  for(p in candidates){
    if(!length(selected) || all(abs(C[p,selected,drop=TRUE])<cor_cutoff,na.rm=TRUE)) selected <- c(selected,p)
  }
  selected
}
selected_groups <- lapply(model_groups,prune_group)
selection_rows <- do.call(rbind,lapply(names(model_groups),function(g){
  data.frame(model_group=g,predictor=model_groups[[g]],status=ifelse(model_groups[[g]]%in%selected_groups[[g]],"used",ifelse(model_groups[[g]]%in%colnames(Xall),"dropped_correlation","missing_or_low_coverage")),stringsAsFactors=FALSE)
}))
write.csv(selection_rows,file.path(outdir,"spde_predictor_selection_by_group.csv"),row.names=FALSE)

pts <- st_as_sf(df,coords=c("longitude","latitude"),crs=4326,remove=FALSE)
loc <- st_coordinates(st_transform(pts,"+proj=eqearth +datum=WGS84 +units=m +no_defs"))/1000
mesh <- inla.mesh.2d(loc=loc,max.edge=c(max_edge,max_edge*2),cutoff=cutoff,offset=c(offset,offset*2))
spde <- inla.spde2.pcmatern(mesh,alpha=2,prior.range=c(500,0.5),prior.sigma=c(1,0.01))
A <- inla.spde.make.A(mesh,loc=loc)
spatial_index <- inla.spde.make.index("spatial",mesh$n)

fixed_rows <- list(); hyper_rows <- list(); model_rows <- list()
for(trait in traits){
  union_preds <- unique(unlist(selected_groups,use.names=FALSE))
  base_work <- is.finite(df[[trait]]) & complete.cases(Xall[,union_preds,drop=FALSE])
  counts <- table(df$taxon_name[base_work]); ok_species <- names(counts[counts>=min_species_n])
  base_work <- base_work & df$taxon_name%in%ok_species
  if(sum(base_work)<100 || length(ok_species)<3){
    for(g in names(selected_groups)) model_rows[[paste(trait,g)]] <- data.frame(trait=trait,model_group=g,status="insufficient",n_observations=sum(base_work),n_species=length(ok_species))
    next
  }
  y <- df[[trait]][base_work]; y <- (y-mean(y))/sd(y)
  species <- as.integer(factor(df$taxon_name[base_work]))
  Awork <- A[base_work,,drop=FALSE]
  for(g in names(selected_groups)){
    preds <- selected_groups[[g]]
    if(!length(preds)){
      model_rows[[paste(trait,g)]] <- data.frame(trait=trait,model_group=g,status="no_predictors",n_observations=length(y),n_species=length(ok_species)); next
    }
    X <- as.data.frame(Xall[base_work,preds,drop=FALSE]); names(X) <- make.names(names(X),unique=TRUE)
    stk <- inla.stack(data=list(y=y),A=list(1,Awork),effects=list(c(list(Intercept=rep(1,length(y)),species=species),X),spatial_index),tag="est")
    form <- as.formula(paste("y ~ 0 + Intercept +",paste(names(X),collapse=" + "),"+ f(species,model='iid') + f(spatial,model=spde)"))
    fit <- tryCatch(inla(form,data=inla.stack.data(stk),family="gaussian",control.predictor=list(A=inla.stack.A(stk),compute=FALSE),control.compute=list(waic=TRUE,dic=TRUE,cpo=TRUE,config=FALSE),verbose=FALSE),error=function(e)e)
    key <- paste(trait,g)
    if(inherits(fit,"error")){
      model_rows[[key]] <- data.frame(trait=trait,model_group=g,status=paste0("error: ",conditionMessage(fit)),n_observations=length(y),n_species=length(ok_species)); next
    }
    fx <- fit$summary.fixed; fx$term <- rownames(fx); fx$trait <- trait; fx$model_group <- g; rownames(fx)<-NULL
    fx$posterior_tail_area_two_sided <- 2*pnorm(-abs(fx$mean/fx$sd))
    fixed_rows[[key]] <- fx
    hp <- fit$summary.hyperpar; hp$parameter <- rownames(hp); hp$trait <- trait; hp$model_group <- g; rownames(hp)<-NULL
    hyper_rows[[key]] <- hp
    model_rows[[key]] <- data.frame(trait=trait,model_group=g,status="ok",n_observations=length(y),n_species=length(ok_species),n_predictors=ncol(X),n_mesh_vertices=mesh$n,waic=fit$waic$waic,dic=fit$dic$dic,cpo_failures=sum(fit$cpo$failure>0,na.rm=TRUE))
    rm(fit,stk); gc()
  }
}

models <- do.call(rbind,model_rows)
models$delta_waic_within_trait <- ave(models$waic,models$trait,FUN=function(x) x-min(x,na.rm=TRUE))
models$delta_dic_within_trait <- ave(models$dic,models$trait,FUN=function(x) x-min(x,na.rm=TRUE))
write.csv(models,file.path(outdir,"spde_model_group_summary.csv"),row.names=FALSE)

if(length(fixed_rows)){
  fixed <- do.call(rbind,fixed_rows)
  non_intercept <- fixed$term!="Intercept"
  fixed$q_bh_posterior_tail_global <- NA_real_
  fixed$q_bh_posterior_tail_within_group <- NA_real_
  fixed$q_bh_posterior_tail_global[non_intercept] <- p.adjust(fixed$posterior_tail_area_two_sided[non_intercept],method="BH")
  for(g in unique(fixed$model_group)){
    idx <- non_intercept & fixed$model_group==g
    fixed$q_bh_posterior_tail_within_group[idx] <- p.adjust(fixed$posterior_tail_area_two_sided[idx],method="BH")
  }
  fixed$credible_95_excludes_zero <- fixed$`0.025quant`>0 | fixed$`0.975quant`<0
  write.csv(fixed,file.path(outdir,"spde_fixed_effects_by_group.csv"),row.names=FALSE)
  envfx <- fixed[non_intercept,,drop=FALSE]
  stability <- do.call(rbind,lapply(split(envfx,paste(envfx$trait,envfx$term,sep="||")),function(z){
    data.frame(trait=z$trait[1],term=z$term[1],n_groups=nrow(z),n_positive=sum(z$mean>0),n_negative=sum(z$mean<0),n_credible_95=sum(z$credible_95_excludes_zero),sign_consistent=length(unique(sign(z$mean[z$mean!=0])))<=1,min_abs_mean=min(abs(z$mean)),max_abs_mean=max(abs(z$mean)),stringsAsFactors=FALSE)
  }))
  write.csv(stability,file.path(outdir,"spde_effect_stability_across_groups.csv"),row.names=FALSE)
}
if(length(hyper_rows)) write.csv(do.call(rbind,hyper_rows),file.path(outdir,"spde_hyperparameters_by_group.csv"),row.names=FALSE)
write.csv(data.frame(n_observations=nrow(df),n_species=length(unique(df$taxon_name)),n_mesh_vertices=mesh$n,max_edge_km=max_edge,cutoff_km=cutoff,offset_km=offset,correlation_cutoff=cor_cutoff,n_models_requested=length(traits)*length(model_groups)),file.path(outdir,"spde_run_metadata.csv"),row.names=FALSE)
