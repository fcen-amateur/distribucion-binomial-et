library(tidyverse)
setwd("~/maestria/teoria-estadistica/distribucion-binomial-et")

tipos_columnas <- cols(
  param = col_character(),
  metodo = col_character(),
  khat = col_double(),
  phat = col_double(),
  k = col_double(),
  p = col_double(),
  n = col_double()
)

dir_raiz <- "sims"
archivos <- dir(dir_raiz, full.names = T)
df <- map_df(archivos[sample(1:60000, 1000)], read_csv, col_types = tipos_columnas, .id = "sim_id")


resumen_tita <- df %>%
#  slice(1:100) %>%
  filter(param == 'tita') %>%
  mutate(ecmhat = (khat - k)^2,
         ecmhat_rel = sqrt(ecmhat)/k)

resumen_tita %>%
  ggplot(aes(k, ecmhat_rel, color = metodo)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = F) +
  coord_cartesian(ylim = c(0,2))

resumen_tita %>%
  group_by((khat == Inf | khat < 0), metodo) %>%
  summarise(mu_ecmhat_rel = mean(ecmhat_rel))

resumen_tita %>%
  ggplot(aes(k, ecmhat_rel, color = metodo)) +
  geom_point() +
  coord_cartesian(ylim = c(0,2))
  
resumen_tita %>%
  select(-c(sim_id, param, ecmhat)) %>%
  arrange(desc(ecmhat_rel)) 

resumen_tita %>%
  filter(khat != Inf, khat > 0) %>%
  group_by(metodo) %>%
  nest() %>%
  mutate(tendencia = map(data, ~lm(ecmhat_rel ~ k, .)))


