#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(TMB); library(INLA); library(Matrix); library(sf)})

traits <- c("orientation_angle_degrees_median","corolla_lab_lightness_median","corolla_lab_chroma_median","corolla_hue_sin_median","corolla_hue_cos_median","shape_aspect_ratio_median","shape_circularity_median","shape_solidity_median","shape_width_cv_median")
climate <- c("chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15")
topo <- c("topo_elevation","topo_slope","topo_roughness")
soil <- paste0("soil_",c("bdod","cec","cfvo","clay","sand","silt","nitrogen","phh2o","soc","ocd"),"_0_30cm")
requested_predictors <- c(climate,topo,soil)
argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(key, default=NULL){i<-match(paste0("--",key),argv); if(is.na(i)) default else argv[[i+1L]]}
log_stage <- function(x){cat(sprintf("[%s] %s\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),x)); flush.console()}

env_file <- getarg("environment"); out_dir <- getarg("out-dir")
if(is.null(env_file)||is.null(out_dir)) stop("Required: --environment PATH --out-dir DIR")
n_factors <- as.integer(getarg("n-factors","1")); min_species_n <- as.integer(getarg("min-species-n","5"))
max_edge <- as.numeric(getarg("mesh-max-edge-km","1500")); cutoff <- as.numeric(getarg("mesh-cutoff-km","100")); offset <- as.numeric(getarg("mesh-offset-km","1500"))
max_observations <- as.integer(getarg("max-observations","0")); max_iterations <- as.integer(getarg("max-iterations","50"))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

log_stage("read environment table")
df <- read.csv(env_file,check.names=FALSE)
need <- c("taxon_name","latitude","longitude",traits); miss <- setdiff(need,names(df)); if(length(miss)) stop("Missing required columns: ",paste(miss,collapse=", "))
predictors <- intersect(requested_predictors,names(df)); if(length(predictors)<3) stop("Fewer than three usable environmental predictors")
for(v in c("latitude","longitude",traits,predictors)) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
df <- df[is.finite(df$latitude)&is.finite(df$longitude),,drop=FALSE]
keep <- names(which(table(df$taxon_name)>=min_species_n)); df <- df[df$taxon_name%in%keep,,drop=FALSE]
usable <- predictors[vapply(predictors,function(p){v<-df[[p]]; sum(is.finite(v))>0 && sd(v,na.rm=TRUE)>0},logical(1))]
complete_n <- vapply(usable,function(p)sum(is.finite(df[[p]])),numeric(1)); predictors <- usable[complete_n>=max(100,0.5*nrow(df))]
dropped <- setdiff(requested_predictors,predictors); if(length(predictors)<3) stop("Fewer than three sufficiently covered predictors")
df <- df[complete.cases(df[,predictors,drop=FALSE]) & rowSums(is.finite(as.matrix(df[,traits,drop=FALSE])))>0,,drop=FALSE]
if(max_observations>0 && nrow(df)>max_observations){set.seed(20260714); idx<-unlist(lapply(split(seq_len(nrow(df)),df$taxon_name),function(z){n<-max(1,round(max_observations*length(z)/nrow(df))); sample(z,min(length(z),n))})); if(length(idx)>max_observations) idx<-sample(idx,max_observations); df<-df[sort(unique(idx)),,drop=FALSE]}
if(nrow(df)<100||length(unique(df$taxon_name))<3) stop("Insufficient data")
log_stage(sprintf("prepared rows=%d species=%d predictors=%d",nrow(df),length(unique(df$taxon_name)),length(predictors)))

X <- matrix(NA_real_,nrow(df),length(predictors),dimnames=list(NULL,predictors)); final_predictors<-character()
for(j in seq_along(predictors)){v<-df[[predictors[j]]]; c0<-v-ave(v,df$taxon_name,FUN=mean); s<-sd(c0); if(is.finite(s)&&s>0){X[,j]<-c0/s; final_predictors<-c(final_predictors,predictors[j])}}
X <- X[,seq_along(final_predictors),drop=FALSE]; predictors<-final_predictors
Y <- matrix(NA_real_,nrow(df),length(traits),dimnames=list(NULL,traits)); scale_meta<-data.frame(trait=traits,mean=NA_real_,sd=NA_real_)
for(t in seq_along(traits)){v<-df[[traits[t]]]; m<-mean(v,na.rm=TRUE); s<-sd(v,na.rm=TRUE); if(!is.finite(s)||s<=0) stop("No trait variation: ",traits[t]); Y[,t]<-(v-m)/s; scale_meta$mean[t]<-m; scale_meta$sd[t]<-s}
species_levels<-sort(unique(df$taxon_name)); species_index<-match(df$taxon_name,species_levels)-1L

log_stage("build projected mesh")
pts<-st_as_sf(df,coords=c("longitude","latitude"),crs=4326,remove=FALSE); loc<-st_coordinates(st_transform(pts,"+proj=eqearth +datum=WGS84 +units=m +no_defs"))/1000
mesh<-inla.mesh.2d(loc=loc,max.edge=c(max_edge,max_edge*2),cutoff=cutoff,offset=c(offset,offset*2)); spde<-inla.spde2.matern(mesh,alpha=2); A<-inla.spde.make.A(mesh,loc=loc)
long<-which(is.finite(Y),arr.ind=TRUE); y<-Y[long]; obs_index<-as.integer(long[,1]-1L); trait_index<-as.integer(long[,2]-1L)
n_traits<-length(traits); n_factors<-min(n_factors,n_traits); n_species<-length(species_levels)
log_stage(sprintf("mesh vertices=%d trait values=%d factors=%d",mesh$n,length(y),n_factors))

cpp<-file.path("ch1_global","v2","tmb_multivariate_spatial.cpp"); log_stage("compile TMB template"); TMB::compile(cpp,flags="-O1"); dyn.load(TMB::dynlib(sub("\\.cpp$","",cpp)))
pars<-list(beta=matrix(0,ncol(X),n_traits),species_re=matrix(0,n_species,n_traits),log_sigma_species=rep(log(.5),n_traits),log_sigma_obs=rep(log(.7),n_traits),omega=matrix(0,mesh$n,n_factors),loadings=matrix(0,n_traits,n_factors),log_kappa=rep(log(1/1000),n_factors),log_tau=rep(0,n_factors))
for(k in seq_len(n_factors)) pars$loadings[k,k]<-.5
lm<-matrix(seq_len(n_traits*n_factors),n_traits,n_factors); for(t in seq_len(n_traits)) for(k in seq_len(n_factors)) if(k>t) lm[t,k]<-NA_integer_; lm<-factor(lm)
log_stage("construct MakeADFun")
obj<-MakeADFun(data=list(y=as.numeric(y),obs_index=obs_index,trait_index=trait_index,species_index=as.integer(species_index),X=X,A=as(A,"TsparseMatrix"),M0=as(spde$param.inla$M0,"TsparseMatrix"),M1=as(spde$param.inla$M1,"TsparseMatrix"),M2=as(spde$param.inla$M2,"TsparseMatrix"),n_traits=n_traits,n_factors=n_factors),parameters=pars,random=c("species_re","omega"),map=list(loadings=lm),DLL="tmb_multivariate_spatial",silent=TRUE)
log_stage(sprintf("MakeADFun complete fixed_parameters=%d",length(obj$par)))
log_stage("evaluate initial objective"); initial_objective<-obj$fn(obj$par); log_stage(sprintf("initial objective=%g",initial_objective))
log_stage("evaluate initial gradient"); initial_gradient<-obj$gr(obj$par); log_stage(sprintf("initial max_abs_gradient=%g",max(abs(initial_gradient))))
log_stage(sprintf("start nlminb max_iterations=%d",max_iterations)); opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=max_iterations,eval.max=max_iterations*3)); log_stage(sprintf("nlminb complete convergence=%d objective=%g",opt$convergence,opt$objective))
fit<-sdreport(obj,par.fixed=opt$par,getJointPrecision=FALSE); ph<-obj$env$parList(opt$par)
beta_out<-expand.grid(predictor=predictors,trait=traits,stringsAsFactors=FALSE); beta_out$estimate_standardized_trait<-as.vector(ph$beta); beta_orig<-ph$beta; for(t in seq_len(n_traits)) beta_orig[,t]<-beta_orig[,t]*scale_meta$sd[t]; beta_out$estimate_original_trait_units_per_1sd_within_species_predictor<-as.vector(beta_orig)
write.csv(beta_out,file.path(out_dir,"tmb_multivariate_environment_effects.csv"),row.names=FALSE)
write.csv(data.frame(requested_predictor=requested_predictors,status=ifelse(requested_predictors%in%predictors,"used","dropped")),file.path(out_dir,"tmb_predictor_availability.csv"),row.names=FALSE)
load_out<-expand.grid(trait=traits,factor=seq_len(n_factors),stringsAsFactors=FALSE); load_out$loading<-as.vector(ph$loadings); write.csv(load_out,file.path(out_dir,"tmb_spatial_factor_loadings.csv"),row.names=FALSE)
write.csv(scale_meta,file.path(out_dir,"tmb_trait_scaling.csv"),row.names=FALSE)
summary<-data.frame(n_observations=nrow(df),n_trait_values=length(y),n_species=n_species,n_traits=n_traits,n_predictors=ncol(X),n_spatial_factors=n_factors,n_mesh_vertices=mesh$n,initial_objective=initial_objective,initial_max_abs_gradient=max(abs(initial_gradient)),optimizer_convergence=opt$convergence,objective=opt$objective,max_abs_gradient=max(abs(obj$gr(opt$par))))
write.csv(summary,file.path(out_dir,"tmb_multivariate_model_summary.csv"),row.names=FALSE)
saveRDS(list(optimizer=opt,sdreport=fit,parameters=ph,mesh=mesh,traits=traits,predictors=predictors,dropped_predictors=dropped,trait_scaling=scale_meta,species_levels=species_levels),file.path(out_dir,"tmb_multivariate_model_fit.rds")); print(summary)
