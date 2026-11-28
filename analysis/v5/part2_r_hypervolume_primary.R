# Standard stochastic Gaussian hypervolumes, not fixed-grid HDR approximation.
# Pilot comparison; no unregistered significance search. No raw observations exported.
suppressPackageStartupMessages(library(hypervolume))
args <- commandArgs(trailingOnly=TRUE)
out <- "work/submission/part2_r_hypervolume"
input <- "work/submission/private/r_hypervolume_inputs.csv"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fit <- function(x, name, h=.5, q=.95, target=60000L) {
  stopifnot(is.matrix(x), nrow(x)>=20, all(is.finite(x)))
  b <- estimate_bandwidth(x,method="fixed",value=rep(h,ncol(x)))
  v <- hypervolume_gaussian(x,name=name,
    samples.per.point=max(20L,ceiling(target/nrow(x))),kde.bandwidth=b,
    sd.count=3,quantile.requested=q,quantile.requested.type="probability",
    verbose=FALSE)
  stopifnot(is(v,"Hypervolume"),is.finite(v@Volume),v@Volume>0,
    ncol(v@RandomPoints)==ncol(x),nrow(v@RandomPoints)>500)
  v
}
compare <- function(a,b,h=.5,q=.95,target=60000L,cap=12000L,seed=20261128L,
                    export=FALSE) {
  set.seed(seed)
  v1 <- fit(a,"Reference",h,q,target)
  v2 <- fit(b,"Comparison",h,q,target)
  z <- hypervolume_set(v1,v2,num.points.max=cap,check.memory=FALSE,verbose=FALSE)
  stopifnot(is(z,"HypervolumeList"))
  s <- z@HVList
  w <- hypervolume_overlap_statistics(z)
  va <- v1@Volume;vb <- v2@Volume;i <- s$Intersection@Volume
  ua <- s$Unique_1@Volume;ub <- s$Unique_2@Volume;u <- s$Union@Volume
  stopifnot(all(is.finite(c(va,vb,i,ua,ub,u))),
    i>=0,ua>=0,ub>=0,abs(va-i-ua)<1e-7*max(1,va),
    abs(vb-i-ub)<1e-7*max(1,vb),abs(u-va-vb+i)<1e-7*max(1,u))
  stat <- data.frame(reference_n=nrow(a),comparison_n=nrow(b),dimensions=ncol(a),
    bandwidth=h,coverage=q,seed=seed,target=target,set_cap=cap,
    reference_points=nrow(v1@RandomPoints),comparison_points=nrow(v2@RandomPoints),
    reference_volume=va,comparison_volume=vb,volume_ratio=vb/va,
    intersection=i,union=u,reference_unique=ua,comparison_unique=ub,
    jaccard=i/u,sorensen=2*i/(va+vb),comparison_unique_fraction=ub/vb,
    turnover=2*min(ua,ub)/(va+vb),net_volume_difference=abs(ua-ub)/(va+vb))
  stopifnot(abs(stat$jaccard-w[["jaccard"]])<1e-8,
            abs(stat$sorensen-w[["sorensen"]])<1e-8)
  pts <- NULL
  if(export){
    # Export only stochastic volume points, never @Data or original observations.
    set.seed(seed+101L)
    pts <- lapply(c("Unique_1","Intersection","Unique_2"),function(k){
      p <- s[[k]]@RandomPoints
      if(nrow(p)==0) return(NULL)
      ix <- sample.int(nrow(p),min(nrow(p),4000L))
      d <- as.data.frame(p[ix,,drop=FALSE]);names(d) <- paste0("x",seq_len(ncol(p)))
      d$component <- k
      d
    })
    pts <- do.call(rbind,pts)
  }
  list(stat=stat,points=pts)
}
if("--self-test" %in% args) {
  set.seed(119L)
  a <- matrix(rnorm(180),60,3);b <- a+matrix(rep(c(.2,.5,0),each=60),60,3)
  x <- compare(a,b,target=6000L,cap=1500L,export=TRUE)
  stopifnot(nrow(x$stat)==1,nrow(x$points)>0, x$stat$jaccard>=0,x$stat$jaccard<=1)
  cat("Actual R hypervolume Gaussian/set synthetic smoke test passed; version ",
      as.character(packageVersion("hypervolume")),"\n")
  quit(save="no",status=0)
}
d <- read.csv(input,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(setequal(unique(d$species),c("Cirsium arvense","Cirsium vulgare")))
meta <- unique(d[c("species","space","design")])
# Fixed priority: exact original primary fits, without changing h/q/seed/target/cap.
meta <- meta[(meta$space=="trait" & meta$design=="native_novel") |
             (meta$space=="climate_raw" & meta$design=="native_all"),,drop=FALSE]
stopifnot(nrow(meta)==4L)
stats <- list();points <- list()
add <- function(row,a,b,h=.5,q=.95,dim=ncol(a),seed=20261128L,
                scenario="primary",target=60000L,cap=12000L,export=FALSE) {
  aa <- a[,seq_len(dim),drop=FALSE];bb <- b[,seq_len(dim),drop=FALSE]
  ans <- compare(aa,bb,h,q,target,cap,seed,export)
  z <- cbind(row,scenario=scenario,ans$stat)
  stats[[length(stats)+1L]] <<- z
  if(export && !is.null(ans$points)) {
    pp <- ans$points
    # Uniform columns across dimensions, for straightforward downstream CSV.
    for(k in setdiff(paste0("x",1:5),names(pp))) pp[[k]] <- NA_real_
    pp <- pp[c(paste0("x",1:5),"component")]
    points[[length(points)+1L]] <<- cbind(row[rep(1,nrow(pp)),,drop=FALSE],
                                          dimensions=dim,pp)
  }
  cat("completed",row$species,row$space,row$design,scenario,dim,
      "ratio",signif(z$volume_ratio,4),"J",signif(z$jaccard,4),"\n")
  flush.console()
  write.csv(do.call(rbind,stats),file.path(out,"comparisons.csv"),row.names=FALSE)
  if(length(points)) write.csv(do.call(rbind,points),file.path(out,"modeled_set_points.csv"),row.names=FALSE)
}
for(k in seq_len(nrow(meta))) {
  row <- meta[k,,drop=FALSE]
  g <- d[d$species==row$species & d$space==row$space & d$design==row$design,,drop=FALSE]
  labs <- if(startsWith(row$design,"matched_")) c("analogous","novel") else
          if(row$design=="native_novel") c("native","novel") else c("native","introduced")
  dim <- if(row$space=="trait") 3L else 5L
  a <- as.matrix(g[g$group==labs[1],paste0("x",seq_len(dim)),drop=FALSE])
  b <- as.matrix(g[g$group==labs[2],paste0("x",seq_len(dim)),drop=FALSE])
  stopifnot(nrow(a)>=20,nrow(b)>=20)
  if(startsWith(row$design,"matched_")) stopifnot(nrow(a)==nrow(b))
  target <- if(dim==5L) 100000L else 60000L
  add(row,a,b,dim=dim,target=target,export=TRUE)
}
z <- do.call(rbind,stats)
stopifnot(nrow(z)==4L,all(is.finite(z$volume_ratio)))
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
cat("Four original primary comparisons complete; no new inferential test.\n")
