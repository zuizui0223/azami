library(httr)
library(jsonlite)

# GloBI API エンドポイント
globi_url <- "https://api.globalbioticinteractions.org/interaction"

# Cirsium に visitedBy の関係がある種
response <- GET(globi_url,
                query = list(
                  sourceTaxon = "Cirsium",
                  interactionType = "visitedBy",
                  limit = 1000
                )
)

# JSON をパース
json_txt <- content(response, as = "text")
json_list <- fromJSON(json_txt)

# DataFrame に変換
df_visitors <- as.data.frame(json_list$data, stringsAsFactors = FALSE)
colnames(df_visitors) <- json_list$columns

# 訪花昆虫の学名
visitors <- unique(df_visitors$target_taxon_name)
visitors <- visitors[!is.na(visitors)]

print(visitors)
genera <- sub(" .*", "", visitors)
unique_genera <- unique(genera)
print(unique_genera)

# 2) 沖縄県ポリゴン（WKT）
okinawa_wkt <- 
  "POLYGON((122.5 24.0,122.5 28.5,131.0 28.5,131.0 24.0,122.5 24.0))"

# 3) 種名 → taxonKey 取得関数
get_key <- function(species_name) {
  kb <- name_backbone(name = species_name, rank = "species")
  kb$usageKey
}

# 4) 1種分のオカレンス取得関数
get_occ <- function(species_name) {
  key <- get_key(species_name)
  if (is.null(key)) return(NULL)
  
  # ページング：5000件ずつ
  # まず総数を確認
  total <- occ_search(
    taxonKey      = key,
    geometry      = okinawa_wkt,
    hasCoordinate = TRUE,
    limit         = 0
  )$meta$count %||% 0L
  
  if (total == 0) return(NULL)
  
  offsets <- seq(0, total-1, by = 5000)
  pages <- map(offsets, function(off) {
    occ_search(
      taxonKey      = key,
      geometry      = okinawa_wkt,
      hasCoordinate = TRUE,
      limit         = 5000,
      start         = off
    )$data
  })
  
  df <- bind_rows(pages)
  df$queried_species <- species_name
  df
}

# 5) 全種分を一気に取得＆結合
okinawa_visitors <- visitors %>%
  set_names() %>%    # 名前付きベクトルにして…
  map_dfr(~ {
    tryCatch(
      get_occ(.x),
      error = function(e) NULL
    )
  })

# 6) 確認
nrow(okinawa_visitors)
head(okinawa_visitors)

library(dplyr)
library(readr)

okinawa_visitors_grp <- okinawa_visitors %>%
  mutate(
    group = case_when(
      order == "Hymenoptera" ~ "trichromatic",                      # ハチ目
      family == "Pieridae"   ~ "trichromatic",                      # シロチョウ科
      genus == "Macroglossum" ~ "tetrachromatic",                  # 昼行性スズメガ
      genus == "Hyles" ~ "simplified",                             # 夜行性スズメガ
      order == "Lepidoptera" ~ "tetrachromatic",                   # その他チョウ目
      order == "Diptera"     ~ "simplified",                       # ハエ目
      TRUE                   ~ NA_character_
    )
  )

# 件数確認
okinawa_visitors_grp %>% count(group)

# 各グループごとに CSV 保存
okinawa_visitors_grp %>% 
  filter(group == "trichromatic") %>% 
  write_csv("~/Desktop/okivis_trichromatic.csv")

okinawa_visitors_grp %>% 
  filter(group == "tetrachromatic") %>% 
  write_csv("~/Desktop/okivis_tetrachromatic.csv")

okinawa_visitors_grp %>% 
  filter(group == "simplified") %>% 
  write_csv("~/Desktop/okivis_simplified.csv")




install.packages(c("terra", "sf", "usdm", "raster"))  # 初回のみ
library(terra)
library(sf)
library(usdm)
library(raster)

# 1) 沖縄県全域ポリゴン（WKT → sf → terra vect）
okinawa_wkt <- "POLYGON((122.5 24.0,122.5 28.5,131.0 28.5,131.0 24.0,122.5 24.0))"
okinawa_sf  <- st_as_sfc(okinawa_wkt, crs = 4326)
okinawa_vec <- vect(okinawa_sf)

# 2) フォルダ内の .tif を全読み込み＆沖縄でマスク
env_dir   <- "~/Desktop/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

# 基準ラスタとして1枚目を読み込み
r0 <- rast(env_files[1]) |>
  crop(okinawa_vec) |>
  mask(okinawa_vec)

# 残りのラスタを同じ範囲＆解像度に合わせて読み込み
rasters <- list(r0)
for (f in env_files[-1]) {
  r <- rast(f) |>
    crop(okinawa_vec) |>
    mask(okinawa_vec) |>
    resample(r0, method = "bilinear")
  rasters <- c(rasters, list(r))
}

# 3) SpatRaster のリストを一気にスタック化
env_stack <- rast(rasters)
print(env_stack) 

# ── 前提：env_stack が raster::RasterStack で読み込まれていること ──
# （以下は raster::stack() で得られたオブジェクトです）

library(usdm)

# 1) RasterStack → data.frame に変換（NA を除く、フルデータを使う場合）
df_all <- as.data.frame(env_stack, na.rm = TRUE)

# 2) VIF 分析（閾値 th = 10）
vif_res <- vifstep(df_all, th = 10)



# 4) 結果確認
print(vif_res)          # 除外・残存変数一覧
print(vif_res@results)  # 各変数の VIF 値

# 5) env_stack_r から残したいレイヤーを抽出
keep_vars    <- vif_res@results$Variables
env_selected <- subset(env_stack_r, keep_vars)

# 6) 確認と書き出し
print(env_selected)
plot(env_selected)
writeRaster(
  env_selected,
  filename  = "~/Desktop/env/okinawa_env_selected.tif",
  format    = "GTiff",
  overwrite = TRUE
)
message("VIF選択後の環境スタックを出力しました: ~/Desktop/env/okinawa_env_selected.tif")


library(SSDM)
library(readr)
library(raster)

#── 色嗜好性グループのファイルリスト ──
groups <- list(
  trichromatic   = "~/Desktop/okivis_trichromatic.csv",
  tetrachromatic = "~/Desktop/okivis_tetrachromatic.csv",
  simplified     = "~/Desktop/okivis_simplified.csv"
)

#── 各グループについてループ処理 ──
for (grp in names(groups)) {
  occ_path <- groups[[grp]]
  
  # (1) CSV読み込み
  occ_df <- read_csv(occ_path, show_col_types = FALSE)
  
  # (2) SSDMモデリング
  ssdm_mod <- stack_modelling(
    algorithms  = c("GLM", "RF"),
    Occurrences = as.data.frame(occ_df),
    Env         = env_stack_r,  # ← 事前に読み込まれているraster stack
    Xcol        = "decimalLongitude",
    Ycol        = "decimalLatitude",
    Spcol       = "species",
    rep         = 1,
    cores       = 8
  )
  
  # (3) GeoTIFF保存（ファイル名＝グループ名）
  tif_fn <- paste0("~/Desktop/", grp, "_diversity.tif")
  writeRaster(
    ssdm_mod@diversity.map,
    filename = tif_fn,
    overwrite = TRUE
  )
  
  # (4) PNGプロット保存（ファイル名＝グループ名）
  png_fn <- paste0("~/Desktop/", grp, "_diversity.png")
  png(png_fn, width = 1000, height = 800)
  plot(ssdm_mod@diversity.map,
       main = paste0("Diversity Map - ", grp))
  dev.off()
  
  message("出力完了: ", grp, " → ", tif_fn, ", ", png_fn)
}


#── 必要パッケージ ──
library(rgbif)
library(dplyr)
library(purrr)
library(readr)


okinawa_wkt <- "POLYGON((122.5 24.0,122.5 28.5,131.0 28.5,131.0 24.0,122.5 24.0))"

cirsium_key <- name_backbone(name = "Cirsium", rank = "genus")$usageKey

message("Downloading up to 200,000 records for genusKey = ", cirsium_key)
occ_data <- occ_search(
  taxonKey      = cirsium_key,
  country       = "JP",
  hasCoordinate = TRUE,
  geometry      = okinawa_wkt,
  limit         = 200000
)$data


message("Downloaded records: ", ifelse(is.null(occ_data), 0, nrow(occ_data)))
head(occ_data)


if (!is.null(occ_data) && nrow(occ_data)>0) {
  write_csv(occ_data, "~/Desktop/okinawa_cirsium.csv")
  saveRDS(occ_data, "~/Desktop/okinawa_cirsium.rds")
  message("Saved to ~/Desktop/okinawa_cirsium.csv and .rds")
} else {
  message("No Cirsium occurrences found for Okinawa.")
}


library(SSDM)
library(raster)
library(dplyr)
library(readr)

okinawa_cirsium <- read_csv("~/Desktop/okinawa_cirsium.csv")

occ_irum <- okinawa_cirsium %>%
  filter(specificEpithet == "irumtiense") %>%
  transmute(decimalLongitude,
            decimalLatitude,
            species = specificEpithet) %>%
  distinct() %>%
  as.data.frame()


env_selected <- stack("~/Desktop/env/okinawa_env_selected.tif")

mod_glm_irum <- modelling(
  algorithm   = "GLM",
  Occurrences = occ_irum,
  Env         = env_selected,
  Xcol        = "decimalLongitude",
  Ycol        = "decimalLatitude",
  Spcol       = "species",
  rep         = 1,
  PA.strategy = "random",
  PA.nb.absences = 1000,
  cores       = 8
)

mod_rf_irum  <- modelling(
  algorithm   = "RF",
  Occurrences = occ_irum,
  Env         = env_selected,
  Xcol        = "decimalLongitude",
  Ycol        = "decimalLatitude",
  Spcol       = "species",
  rep         = 1,
  PA.strategy = "random",
  PA.nb.absences = 1000,
  cores       = 8
)


ens_mod_irum <- ensemble_modelling(
  models      = list(GLM = mod_glm_irum, RF = mod_rf_irum),
  Occurrences = occ_irum,
  Env         = env_selected,
  Xcol        = "decimalLongitude",
  Ycol        = "decimalLatitude",
  Spcol       = "species",
  algorithms  = c("GLM", "RF"),
  stat        = "mean"
)


pred_irum <- SSDM::project(
  ens_mod_irum,
  env_selected,
  output.format = "rasters",
  uncertainty   = TRUE
)


raster::writeRaster(
  pred_irum$projection,
  filename  = "~/Desktop/cirsium_irumtiense_suitability.tif",
  format    = "GTiff",
  overwrite = TRUE
)
raster::writeRaster(
  pred_irum$uncertainty,
  filename  = "~/Desktop/cirsium_irumtiense_uncertainty.tif",
  format    = "GTiff",
  overwrite = TRUE
)


png("~/Desktop/cirsium_irumtiense_suitability.png", width=800, height=600)
raster::plot(
  pred_irum$projection,
  main="Ensemble SDM: Cirsium irumtiense",
  xlab="Longitude", ylab="Latitude"
)
dev.off()

library(raster)
library(stats)


brev <- raster("~/Desktop/cirsium_brevicule_suitability.tif")  # シマアザミ
irum <- raster("~/Desktop/cirsium_irumtiense_suitability.tif") # イリオモテアザミ

tri   <- raster("~/Desktop/trichromatic_diversity.tif")
tetra <- raster("~/Desktop/tetrachromatic_diversity.tif")
simp  <- raster("~/Desktop/simplified_diversity.tif")


stk_brev <- stack(brev, tri, tetra, simp)
names(stk_brev) <- c("brev", "tri", "tetra", "simp")

df_brev <- na.omit(as.data.frame(getValues(stk_brev)))

mod_brev <- glm(brev ~ tri + tetra + simp, data = df_brev, family = quasibinomial())
summary(mod_brev)


stk_irum <- stack(irum, tri, tetra, simp)
names(stk_irum) <- c("irum", "tri", "tetra", "simp")

df_irum <- na.omit(as.data.frame(getValues(stk_irum)))

mod_irum <- glm(irum ~ tri + tetra + simp, data = df_irum, family = quasibinomial())
summary(mod_irum)

ggplot(df_brev_long, aes(x = visitor_diversity, y = suitability)) +
  geom_point(alpha = 0.3, color = "gray20") +
  geom_smooth(
    method = "glm", 
    method.args = list(family = "quasibinomial"), 
    se = TRUE, 
    color = "black"
  ) +
  facet_wrap(
    ~ visitor_group,
    scales = "free_x",
    labeller = as_labeller(c(
      tri  = "Trichromatic visitors",
      tetra = "Tetrachromatic visitors",
      simp = "Simplified vision visitors"
    ))
  ) +
  labs(
    x = "Visitor diversity",
    y = expression("Suitability of " * italic("C. brevicaule"))
  ) +
  theme_minimal()

ggplot(df_irum_long, aes(x = visitor_diversity, y = suitability)) +
  geom_point(alpha = 0.3, color = "gray20") +
  geom_smooth(
    method = "glm", 
    method.args = list(family = "quasibinomial"), 
    se = TRUE, 
    color = "black"
  ) +
  facet_wrap(
    ~ visitor_group,
    scales = "free_x",
    labeller = as_labeller(c(
      tri  = "Trichromatic visitors",
      tetra = "Tetrachromatic visitors",
      simp = "Simplified vision visitors"
    ))
  ) +
  labs(
    x = "Visitor diversity",
    y = expression("Suitability of " * italic("C. irumtiense"))
  ) +
  theme_minimal()

library(dplyr)
library(ggplot2)

df_plot <- df_irum_long %>%
  filter(visitor_group %in% c("tetra","simp")) %>%
  mutate(
    visitor_taxon = recode(visitor_group,
                           "tetra" = "Lepidoptera (Tetrachromatic)",
                           "simp"  = "Diptera (Simplified)")
  )

ggplot(df_plot, aes(x = visitor_diversity, y = suitability)) +
  geom_point(alpha = 0.35, color = "gray20", size = 1.6) +
  geom_smooth(
    method = "glm",
    method.args = list(family = quasibinomial()),
    se = TRUE, color = "black", linewidth = 0.7
  ) +
  facet_wrap(~ visitor_taxon, scales = "free_x") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(
    title = "Relationship between pollinator diversity and habitat suitability",
    subtitle = "Diptera (simplified) vs.Lepidoptera (tetrachromatic)",
    x = "Visitor diversity (Shannon index, per site)",
    y = expression("Predicted habitat suitability for " * italic("C. irumtiense") * " (0–1)"),
    caption = "Each dot = one site. 
x-axis = pollinator diversity (Shannon index from observations). 
y-axis = habitat suitability predicted by MaxEnt models.
Facets = main pollinator group with associated visual system."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid.minor = element_blank()
  )


library(terra)
library(dplyr)


irum2  <- resample(irum,  brev)

tri2   <- resample(tri,   brev)
tetra2 <- resample(tetra, brev)
simp2  <- resample(simp,  brev)


stk <- c(brev, irum2, tri2, tetra2, simp2)
names(stk) <- c("brev", "irum", "tri", "tetra", "simp")


xy   <- as.data.frame(terra::xyFromCell(brev, 1:ncell(brev)))
vals <- as.data.frame(terra::values(stk))

dat <- bind_cols(xy, vals)
colnames(dat)[1:2] <- c("x", "y")


dat <- dat %>%
  filter(!is.na(brev), !is.na(irum),
         !is.na(tri), !is.na(tetra), !is.na(simp)) %>%
  mutate(
    total_suit = brev + irum
  ) %>%
  filter(total_suit > quantile(total_suit, 0.6, na.rm = TRUE)) %>% 
  mutate(
    species_bin = ifelse(irum > brev, 1L, 0L)  
  )

table(dat$species_bin)
dat_glm <- dat %>%
  mutate(
    z_tri   = as.numeric(scale(tri)),
    z_tetra = as.numeric(scale(tetra)),
    z_simp  = as.numeric(scale(simp))
  )

m_glm <- glm(
  species_bin ~ z_tri + z_tetra + z_simp,
  data   = dat_glm,
  family = binomial(link = "logit")
)

summary(m_glm)


crs_brev <- crs(brev, proj = TRUE)

pts_sf <- st_as_sf(dat, coords = c("x", "y"), crs = crs_brev)
pts_sf <- st_transform(pts_sf, 3857)   

coords_utm <- st_coordinates(pts_sf)
dat$x_utm <- coords_utm[,1]
dat$y_utm <- coords_utm[,2]


set.seed(123)

N_TARGET <- 25000   

if (nrow(dat) > N_TARGET) {
  dat_sub <- dat %>% slice_sample(n = N_TARGET)
} else {
  dat_sub <- dat
}


dat_sub <- dat_sub %>%
  mutate(
    z_tri   = as.numeric(scale(tri)),
    z_tetra = as.numeric(scale(tetra)),
    z_simp  = as.numeric(scale(simp))
  )

coords <- as.matrix(dat_sub[, c("x_utm", "y_utm")])

mesh <- inla.mesh.2d(
  loc      = coords,
  max.edge = c(50e3, 200e3),  
  cutoff   = 20e3,           
  offset   = c(50e3, 200e3)  
)


plot(mesh); points(coords[sample(1:nrow(coords), 2000),], pch=".", col="red")

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(200e3, 0.5),  
  prior.sigma = c(1, 0.01)      
)

spde_index <- inla.spde.make.index(
  name   = "spatial",
  n.spde = spde$n.spde
)


A_obs <- inla.spde.make.A(
  mesh = mesh,
  loc  = coords
)



stack_est <- inla.stack(
  data = list(
    species_bin = dat_sub$species_bin
  ),
  A = list(A_obs, 1),
  effects = list(
    spatial = spde_index,
    data.frame(
      Intercept = 1,
      z_tri     = dat_sub$z_tri,
      z_tetra   = dat_sub$z_tetra,
      z_simp    = dat_sub$z_simp
    )
  ),
  tag = "est"
)


formula_inla <- species_bin ~ 
  1 + z_tri + z_tetra + z_simp +
  f(spatial, model = spde)


fit_inla_vi <- inla(
  formula_inla,
  family = "binomial",
  data   = inla.stack.data(stack_est),
  control.predictor = list(
    A       = inla.stack.A(stack_est),
    compute = TRUE
  ),
  control.compute = list(
    dic  = TRUE,
    waic = TRUE,
    cpo  = TRUE
  ),
  verbose = TRUE
)

summary(fit_inla_vi)

fit_inla_vi$summary.fixed

fit_inla_vi$summary.hyperpar

sp_field <- fit_inla_vi$summary.random$spatial

mesh_df <- data.frame(
  x = mesh$loc[,1],
  y = mesh$loc[,2],
  mean = sp_field$mean
)

ggplot(mesh_df, aes(x, y, color = mean)) +
  geom_point(size = 1) +
  scale_color_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "Spatial effect"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Spatial random field (INLA–SPDE)")


idx <- inla.stack.index(stack_est, "est")$data
pred <- fit_inla_vi$summary.fitted.values[idx, "mean"]

dat_sub$pred_pink <- pred

library(ggplot2)

ggplot(dat_sub, aes(x_utm, y_utm, color = pred_pink)) +
  geom_point(size = 0.6) +
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
    limits = c(0, 1),
    name = "Predicted\nP(pink)"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Predicted pink-flower probability (INLA–SPDE)")



library(INLA)
library(dplyr)
library(ggplot2)


idx_est <- inla.stack.index(stack_est, "est")$data
y_obs   <- dat_sub$species_bin  # 1=イリオモテ側, 0=シマ側


tjur_R2 <- function(p, y) {
  mean(p[y == 1]) - mean(p[y == 0])
}


p_full <- fit_inla_vi$summary.fitted.values[idx_est, "mean"]
R2_full <- tjur_R2(p_full, y_obs)


stack_fix <- inla.stack(
  data = list(species_bin = dat_sub$species_bin),
  A    = list(1),
  effects = list(
    data.frame(
      Intercept = 1,
      z_tri   = dat_sub$z_tri,
      z_tetra = dat_sub$z_tetra,
      z_simp  = dat_sub$z_simp
    )
  ),
  tag = "fix"
)

form_fix <- species_bin ~ 1 + z_tri + z_tetra + z_simp

fit_fix <- inla(
  form_fix,
  family = "binomial",
  data   = inla.stack.data(stack_fix),
  control.predictor = list(
    A       = inla.stack.A(stack_fix),
    compute = TRUE
  ),
  control.compute = list(dic = TRUE, waic = TRUE)
)

idx_fix <- inla.stack.index(stack_fix, "fix")$data
p_fix   <- fit_fix$summary.fitted.values[idx_fix, "mean"]
R2_fix  <- tjur_R2(p_fix, y_obs)



stack_sp <- inla.stack(
  data = list(
    species_bin = dat_sub$species_bin
  ),
  A = list(
    A_obs,   
    1       
  ),
  effects = list(
    spatial = spde_index,
    data.frame(
      Intercept = rep(1, nrow(dat_sub))   
    )
  ),
  tag = "sp"
)

form_sp <- species_bin ~ 1 + f(spatial, model = spde)

fit_sp <- inla(
  form_sp,
  family = "binomial",
  data   = inla.stack.data(stack_sp),
  control.predictor = list(
    A       = inla.stack.A(stack_sp),
    compute = TRUE
  ),
  control.compute = list(
    dic  = TRUE,
    waic = TRUE
  )
)

idx_sp <- inla.stack.index(stack_sp, "sp")$data
p_sp   <- fit_sp$summary.fitted.values[idx_sp, "mean"]
R2_sp  <- tjur_R2(p_sp, y_obs)

df_r2 <- tibble(
  component = c("Full model", "Pollinator vision only", "Spatial structure only"),
  R2        = c(R2_full, R2_fix, R2_sp)
)

print(df_r2)


p_r2 <- ggplot(df_r2, aes(x = component, y = R2)) +
  geom_col(fill = "grey70", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", R2)),
            vjust = -0.5, size = 4) +
  ylim(0, max(df_r2$R2) * 1.15) +
  labs(
    title = "Variance partition: pollinator vision vs spatial structure",
    x     = NULL,
    y     = "Tjur's R² (explained variation in pink vs white)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 10, hjust = 1)
  )

p_r2

ggsave("~/Desktop/inla_variance_partition.png",
       p_r2, width = 6.5, height = 4.2, dpi = 600)

