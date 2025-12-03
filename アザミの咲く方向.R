library(terra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(readr)
library(SSDM)
library(raster)
library(usdm)

# 行政範囲でマスクする準備
japan_sf <- ne_countries(country = "Japan", returnclass = "sf")
japan_vect <- vect(japan_sf)

# フォルダ内の .tif を全部読む
env_files <- list.files("~/Desktop/env", pattern = "\\.tif$", full.names = TRUE)

# 1枚目をベースにする
r1 <- rast(env_files[1])
r1 <- crop(r1, japan_vect) |> mask(japan_vect)

# 同じ範囲・解像度に合わせて全部読み込む
raster_list <- list(r1)
for (i in 2:length(env_files)) {
  r <- rast(env_files[i]) |> crop(japan_vect) |> mask(japan_vect)
  r <- resample(r, r1, method = "bilinear")
  raster_list[[i]] <- r
}

# スタックにする
env_stack <- rast(raster_list)
print(env_stack)

stacked_rasters <- rast(raster_list) 
vif_result <- vifstep(stacked_rasters)
print(vif_result)
keep_vars <- vif_result@results$Variables
keep_vars <- vif_result@results$Variables

env_stack_vif <- env_stack[[keep_vars]]
Env_r <- raster::stack(env_stack_vif)



clean_occ <- function(df, g){
  df %>%
    dplyr::filter(group == g) %>%
    dplyr::group_by(species) %>%
    dplyr::filter(n() > 20) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lon_r = round(longitude, 4),
      lat_r = round(latitude, 4)
    ) %>%
    dplyr::distinct(species, lon_r, lat_r, .keep_all = TRUE) %>%
    dplyr::select(species, longitude, latitude)
}



run_ssdm <- function(occ_df, model_name){
  message("=== Running SSDM for ", model_name, " ===")
  
  ssdm_mod <- stack_modelling(
    algorithms   = c("GLM", "RF"),
    Occurrences  = as.data.frame(occ_df),
    Env          = Env_r,
    Xcol         = "longitude",
    Ycol         = "latitude",
    Spcol        = "species",
    rep          = 1,
    cores        = 8,
    verbose      = TRUE
  )
  
  return(ssdm_mod)
}

# ============================
# 4. Cirsium の SSDM（group0 / group1）
# ============================

cirsium_df <- read_csv("~/Desktop/cirsium2.csv")

occ_g0 <- clean_occ(cirsium_df, 0)
occ_g1 <- clean_occ(cirsium_df, 1)

ssdm_g0 <- run_ssdm(occ_g0, "Cirsium_Group0")
ssdm_g1 <- run_ssdm(occ_g1, "Cirsium_Group1")

plot(ssdm_g0@diversity.map, main="SSDM Stack (Group0)")
plot(ssdm_g1@diversity.map, main="SSDM Stack (Group1)")

terra::writeRaster(ssdm_g0@diversity.map, "~/Desktop/ssdm_group0_stack.tif", overwrite=TRUE)
terra::writeRaster(ssdm_g1@diversity.map, "~/Desktop/ssdm_group1_stack.tif", overwrite=TRUE)

# ============================
# 5. Visitor の SSDM（group1〜5）
# ============================

visitor_df <- read_csv("~/Desktop/visitor2.csv")

library(dplyr)
library(SSDM)

### ============================
### 1. genus からギルド分類
### ============================
visitor_df <- visitor_df %>%
  mutate(genus = sub(" .*", "", species))

# ---- Guild definitions ----
bee_genus <- c("Bombus","Xylocopa","Apis","Megachile","Eucera","Halictus","Lasioglossum")
butter_genus <- c(
  "Papilio","Pieris","Vanessa","Argynnis","Aglais","Fabriciana","Erynnis",
  "Polygonia","Nymphalis","Cynthia","Speyeria","Aporia","Erebia","Araschnia"
)
hover_genus <- c("Episyrphus","Eristalis","Melanostoma","Syrphus","Sphaerophoria","Eupeodes")
macro_genus <- c("Macroglossum")

bee_df   <- visitor_df %>% filter(genus %in% bee_genus)
butter_df <- visitor_df %>% filter(genus %in% butter_genus)
hover_df <- visitor_df %>% filter(genus %in% hover_genus)
macro_df <- visitor_df %>% filter(genus %in% macro_genus)

### ============================
### 2. 重複除去＋20点未満の種を除外
### ============================
clean_occ <- function(df){
  df %>%
    dplyr::mutate(
      lon_r = round(longitude, 4),
      lat_r = round(latitude, 4)
    ) %>%
    dplyr::distinct(species, lon_r, lat_r, .keep_all = TRUE) %>%
    dplyr::group_by(species) %>%
    dplyr::filter(n() >= 20) %>%
    dplyr::ungroup() %>%
    dplyr::select(species, longitude, latitude)
}

occ_bee   <- clean_occ(bee_df)
occ_butter <- clean_occ(butter_df)
occ_hover <- clean_occ(hover_df)
occ_macro <- clean_occ(macro_df)

run_ssdm <- function(occ, label){
  message("=== Running SSDM for: ", label, " ===")
  
  ssdm_mod <- SSDM::stack_modelling(
    algorithms   = c("GLM", "RF"),
    Occurrences  = as.data.frame(occ),
    Env          = Env_r,
    Xcol         = "longitude",
    Ycol         = "latitude",
    Spcol        = "species",
    rep          = 1,
    cores        = 8,
    verbose      = TRUE
  )
  
  # Plot
  plot(ssdm_mod@diversity.map, main = paste("SSDM:", label))
  
  # Save raster
  out_tif <- paste0("~/Desktop/ssdm_", label, "_richness.tif")
  terra::writeRaster(ssdm_mod@diversity.map, out_tif, overwrite = TRUE)
  
  message("Saved: ", out_tif)
  
  return(ssdm_mod)
}

### ============================
### 4. 実行
### ============================
ssdm_bee    <- run_ssdm(occ_bee, "bee")
ssdm_butter <- run_ssdm(occ_butter, "butterfly")
ssdm_hover  <- run_ssdm(occ_hover, "hoverfly")
ssdm_macro  <- run_ssdm(occ_macro, "macroglossum")

library(terra)
library(dplyr)
library(sf)
library(INLA)  # すでにインストール済み前提

## ---- 1) ラスタ読み込み ----
r_down   <- rast("~/Desktop/ssdm_group0_stack.tif")          # 下向き
r_up     <- rast("~/Desktop/ssdm_group1_stack.tif")          # 上向き

r_bee    <- rast("~/Desktop/ssdm_bee_richness.tif")          # ハチ
r_butter <- rast("~/Desktop/ssdm_butterfly_richness.tif")    # チョウ
r_hover  <- rast("~/Desktop/ssdm_hoverfly_richness.tif")     # ハナアブ

## ---- 2) 解像度・範囲を down に合わせる ----
r_bee2    <- resample(r_bee,    r_down)
r_butter2 <- resample(r_butter, r_down)
r_hover2  <- resample(r_hover,  r_down)

## ---- 3) stack にまとめる ----
stk <- c(r_down, r_up, r_bee2, r_butter2, r_hover2)
names(stk) <- c("down", "up", "bee", "butter", "hover")

## ---- 4) xy + 値を data.frame に ----
xy   <- as.data.frame(xyFromCell(r_down, 1:ncell(r_down)))
vals <- as.data.frame(values(stk))
dat  <- bind_cols(xy, vals)
colnames(dat)[1:2] <- c("x", "y")   # ← ここが座標名

## ---- 5) 欠損除去 & 弱いピクセルをフィルタ ----
dat <- dat %>%
  filter(!is.na(down), !is.na(up),
         !is.na(bee), !is.na(butter), !is.na(hover)) %>%
  mutate(total_flower = down + up) %>%
  # 花の豊富度がある程度ある場所だけ（上位 60%）
  filter(total_flower > quantile(total_flower, 0.6, na.rm = TRUE)) %>%
  mutate(
    orient_bin = ifelse(down > up, 1L, 0L),           # 下向き優勢=1, 上向き優勢=0
    orient_idx = (down - up) / (down + up + 1e-6)     # -1〜+1の連続指標
  )

table(dat$orient_bin)

## =========================
## 3. 座標を sf で投影
## =========================

# ★ポイント１：coords は c("x","y")
# ★ポイント２：CRS は r_down に合わせる方が安全
crs_down <- crs(r_down, proj=TRUE)  # 文字列として取得

pts_sf <- st_as_sf(dat, coords = c("x", "y"), crs = crs_down)
pts_sf <- st_transform(pts_sf, 3857)   # 平面座標（単位: m）

coords_utm <- st_coordinates(pts_sf)
dat$x_utm <- coords_utm[,1]
dat$y_utm <- coords_utm[,2]

set.seed(123)

N_TARGET <- 25000  # 好きに調整（PC きつければ 10000 くらいでもOK）

if (nrow(dat) > N_TARGET) {
  dat_sub <- dat %>% slice_sample(n = N_TARGET)
} else {
  dat_sub <- dat
}

## 共変量をスケール
dat_sub <- dat_sub %>%
  mutate(
    z_bee    = as.numeric(scale(bee)),
    z_butter = as.numeric(scale(butter)),
    z_hover  = as.numeric(scale(hover))
  )

## =========================
## 4. メッシュ作成
## =========================
coords <- as.matrix(dat_sub[, c("x_utm", "y_utm")])

mesh <- inla.mesh.2d(
  loc       = coords,
  max.edge  = c(50e3, 200e3),  # メッシュの細かさ（50km, 200km）
  cutoff    = 20e3,            # 近接点をまとめる距離（20km）
  offset    = c(50e3, 200e3)   # 領域から外に余白
)

# 必要なら確認
# plot(mesh); points(coords[sample(1:nrow(coords), 2000),], pch=".", col="red")

## =========================
## 5. SPDE モデル定義（PC prior）
## =========================
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(200e3, 0.5),  # P(range < 200km) = 0.5
  prior.sigma = c(1, 0.01)      # P(sigma > 1) = 0.01
)

spde_index <- inla.spde.make.index(
  name   = "spatial",
  n.spde = spde$n.spde
)

## 観測点用 A 行列
A_obs <- inla.spde.make.A(
  mesh = mesh,
  loc  = coords
)

## =========================
## 6. INLA stack 作成
## =========================
stack_est <- inla.stack(
  data = list(
    orient_bin = dat_sub$orient_bin
  ),
  A = list(A_obs, 1),
  effects = list(
    spatial = spde_index,
    data.frame(
      Intercept = 1,
      z_bee     = dat_sub$z_bee,
      z_butter  = dat_sub$z_butter,
      z_hover   = dat_sub$z_hover
    )
  ),
  tag = "est"
)

## =========================
## 7. モデル式 & 当てる
## =========================
formula_inla <- orient_bin ~ 
  1 + z_bee + z_butter + z_hover +
  f(spatial, model = spde)

fit_inla <- inla(
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

summary(fit_inla)

## 固定効果（ギルド効果）
fit_inla$summary.fixed

## 空間場の分散・レンジ
fit_inla$summary.hyperpar







