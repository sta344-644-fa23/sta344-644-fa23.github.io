library(sf)
library(sp)
library(dplyr)

load(url("http://www.stat.osu.edu/~pfc/teaching/5012_spatial_statistics/datasets/Scottish_shapefile.RData"))
lips <- read.table("http://www.stat.osu.edu/~pfc/teaching/5012_spatial_statistics/datasets/Scottish_lip_cancer.txt", header=T)

lip_cancer = st_as_sf(Scottish.shapefile) %>% 
  rename(Name = NAME, District = ID) %>% 
  full_join(lips) %>%
  rename(id = District, District = Name)

rm(lips, Scottish.shapefile)

st_crs(lip_cancer) = st_crs(lip_cancer)

saveRDS(lip_cancer, file=here::here("data/lip_cancer.rds"))

