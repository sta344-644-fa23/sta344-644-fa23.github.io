library(forecast)
library(dplyr)

aus_wine = data.frame(
  date  = time(wineind) %>% as.double(),
  sales = c(wineind)
) %>% tbl_df()

save(aus_wine, file="aus_wine.Rdata")