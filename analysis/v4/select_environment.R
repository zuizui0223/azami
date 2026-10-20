# Select once on the full environment-only native table, before any trait fit.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- args[1]
output <- args[2]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
metadata <- jsonlite::read_json(file.path(input,"selection_input.json"), simplifyVector=TRUE)
x <- read.csv(file.path(input,"environment_only.csv"), check.names=FALSE)
stopifnot(identical(names(x),metadata$candidates), nrow(x)==metadata$selection_n)
stopifnot(all(is.finite(as.matrix(x))), all(vapply(x,sd,numeric(1))>0))
selection <- usdm::vifstep(x, th=metadata$threshold, keep=NULL, size=nrow(x)+1)
selected <- as.character(selection@results$Variables)
excluded <- as.character(selection@excluded)
# Reconstruct each exclusion's VIF table through the public package API.
trace <- list()
remaining <- names(x)
for (v in excluded) {
  values <- usdm::vif(x[,remaining,drop=FALSE],size=nrow(x)+1)
  stopifnot(v %in% remaining)
  stopifnot(values$VIF[match(v,values$Variables)] >= max(values$VIF)-1e-7)
  trace[[length(trace)+1]] <- list(removed=v, vif=values)
  remaining <- setdiff(remaining,v)
}
stopifnot(setequal(remaining,selected),all(selection@results$VIF<=metadata$threshold))
result <- list(status="selected_before_trait_models", input=metadata,
               selected=selected, excluded_in_order=excluded, removal_trace=trace,
               final_vif=selection@results,
               correlation=list(columns=names(x),values=unname(cor(x))),
               usdm_version=as.character(packageVersion("usdm")),
               r_version=R.version.string,
               git_sha=Sys.getenv("GITHUB_SHA"),run_id=Sys.getenv("GITHUB_RUN_ID"),
               interpretation="VIF selection reduces predictor redundancy; not causal confounding control. Groups need not all survive.")
jsonlite::write_json(result,file.path(output,"environment_selection.json"),
                     pretty=TRUE,auto_unbox=TRUE,digits=16)
write.csv(selection@results,file.path(output,"selected_environment_vif.csv"),row.names=FALSE)
writeLines(c("# 環境変数の選定",
             "6群・9候補の全自生域完全行にusdm::vifstepを適用。閾値10、強制保持なし、抽出なし。",
             paste("採択：",paste(selected,collapse=", ")),
             paste("除外順：",paste(excluded,collapse=", ")),
             "形質の値・測定可否は選定に使用していない。群ごとに1変数残る保証はない。",
             "選定後の種内／種間モデル行列の診断は別途必要。VIF低下は因果的交絡除去の証明ではない。"),
           file.path(output,"environment_selection_ja.md"),useBytes=TRUE)
