# =============================================================================
# Análisis de las comunidades ecológicas (apartado 3.1 del manuscrito)
#
# Entrada : results_1000_with_and_without.csv
#           1.000 redes de 10 especies simuladas sin y con competencia por
#           recursos (fusión de los recursos de dos especies).
# Salida  : cifras del apartado 3.1 en consola y en resultados_ecologia.txt
#           Tabla 1 y Tabla suplementaria S2 en formato CSV
#
# Ejecutar desde la carpeta raíz del repositorio:  Rscript analisis_ecologia.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(lme4)      # modelos lineales mixtos
  library(geepack)   # ecuaciones de estimación generalizadas
  library(survival)  # regresión logística condicional
})

# ── Rutas ────────────────────────────────────────────────────────────────────
DATOS      <- "datos/results_1000_with_and_without.csv"
RESULTADOS <- "resultados"
dir.create(RESULTADOS, showWarnings = FALSE)

# ── Utilidades para informar de los resultados ───────────────────────────────
salida <- file.path(RESULTADOS, "resultados_ecologia.txt")
cat("", file = salida)
seccion <- function(titulo) {
  txt <- paste0("\n", strrep("=", 78), "\n", titulo, "\n", strrep("=", 78), "\n")
  cat(txt); cat(txt, file = salida, append = TRUE)
}
dato <- function(etiqueta, valor) {
  txt <- paste0("  ", format(etiqueta, width = 64), " ", valor, "\n")
  cat(txt); cat(txt, file = salida, append = TRUE)
}
pct <- function(x, d = 2) sprintf(paste0("%.", d, "f%%"), 100 * x)
fp  <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))
ic  <- function(b, se, f = exp) sprintf("%.2f (%.2f–%.2f)", f(b), f(b - 1.96 * se), f(b + 1.96 * se))


# =============================================================================
# 1. Preparación de los datos
# =============================================================================
raw <- read.csv(DATOS, stringsAsFactors = FALSE)

# Los rasgos de cada especie se miden en el escenario sin competencia, porque la
# fusión de recursos modifica las interacciones entrantes de las competidoras.
# Los enlaces sintróficos se cuentan a la vez entre los productores y entre los
# consumidores de una especie, así que se descomponen en enlaces no recíprocos y
# sintróficos.
basal <- raw %>%
  filter(network_competition == 0) %>%
  transmute(network_id, species,
            focal   = sp_competition,              # especie del par que compite
            ext_res = incoming_external_res,       # recibe el aporte del quimiostato
            rep_rate,                              # tasa de replicación
            incoming_resources,                    # suministro total
            outgoing_resources,                    # subproductos liberados (requerimiento)
            incoming_degree,                       # número de productores
            in_uni  = incoming_degree - syntrophy_links,   # productores no recíprocos
            out_uni = outgoing_degree - syntrophy_links,   # consumidores no recíprocos
            syn     = syntrophy_links,                     # enlaces sintróficos
            in_rate  = ifelse(incoming_degree > 0, incoming_resources / incoming_degree, 0),  # suministro medio por productor
            out_rate = ifelse(outgoing_degree > 0, outgoing_resources / outgoing_degree, 0))  # liberación media por consumidor

COV <- c("rep_rate", "incoming_resources", "outgoing_resources",
         "in_uni", "out_uni", "syn", "in_rate", "out_rate")

d <- raw %>%
  select(network_id, species, comp = network_competition, abundance) %>%
  inner_join(basal, by = c("network_id", "species")) %>%
  mutate(ext = as.integer(abundance == 0)) %>%
  arrange(network_id, comp, species)
for (v in COV) d[[paste0("z_", v)]] <- as.numeric(scale(d[[v]]))   # unidades de desviación típica

seccion("Diseño")
dato("Redes", n_distinct(d$network_id))
dato("Observaciones (especie x red x escenario)", nrow(d))


# =============================================================================
# 3.1.1  La abundancia de cada especie aumenta con el suministro disponible
# =============================================================================
seccion("3.1.1 Abundancia, suministro y requerimiento")

# Elasticidades en el escenario sin competencia. Se excluyen la especie del
# quimiostato y las que no liberan subproductos, cuyo suministro o requerimiento
# no recoge la matriz de interacciones.
el <- d %>% filter(comp == 0, abundance > 0, ext_res == 0,
                   outgoing_resources > 0, incoming_resources > 0)
m_el <- lmer(log(abundance) ~ log(incoming_resources) + log(outgoing_resources) + log(rep_rate) +
               (1 | network_id), data = el)
ce <- summary(m_el)$coefficients
dato("Especies analizadas", nrow(el))
dato("Elasticidad del suministro (estimación ± EE)", sprintf("%.2f ± %.2f", ce[2, 1], ce[2, 2]))
dato("Elasticidad de la liberación de subproductos", sprintf("%.2f ± %.2f", ce[3, 1], ce[3, 2]))
dato("Elasticidad de la tasa de replicación", sprintf("%.2f ± %.2f", ce[4, 1], ce[4, 2]))

# Modelo de extinción (Tabla 1): GEE binomial con redes como clústeres
FIX <- paste("comp * focal +", paste(paste0("z_", COV), collapse = " + "), "+ ext_res")
m1 <- geeglm(as.formula(paste("ext ~", FIX)), id = network_id, data = d,
             family = binomial, corstr = "exchangeable")
s1 <- summary(m1)$coefficients

# Modelo de abundancia de las supervivientes (Tabla 1): modelo lineal mixto
s2 <- d %>% filter(abundance > 0)
m2 <- lmer(as.formula(paste("log(abundance) ~", FIX, "+ (1 | network_id)")), data = s2)
c2 <- summary(m2)$coefficients

# Simetría de los efectos del suministro y de la liberación sobre la extinción
b <- coef(m1); V <- vcov(m1)
k <- c("z_incoming_resources", "z_outgoing_resources")
suma <- sum(b[k]); ee_suma <- sqrt(V[k[1], k[1]] + V[k[2], k[2]] + 2 * V[k[1], k[2]])
dato("Cambio de abundancia por +1 DE de suministro", sprintf("%+.1f%%", 100 * (exp(c2["z_incoming_resources", 1]) - 1)))
dato("Cambio de abundancia por +1 DE de liberación", sprintf("%+.1f%%", 100 * (exp(c2["z_outgoing_resources", 1]) - 1)))
dato("Suma de coeficientes de extinción (suministro + liberación)", sprintf("%.3f ± %.3f", suma, ee_suma))
dato("Test de Wald de simetría, P", sprintf("%.2f", pchisq((suma / ee_suma)^2, 1, lower.tail = FALSE)))

# Tabla 1
nombres_t1 <- c(comp = "Scenario with competition", focal = "Competing species",
                `comp:focal` = "Scenario with competition × competing species",
                z_rep_rate = "Replication rate", z_incoming_resources = "Total resource supply",
                z_outgoing_resources = "Total by-product release", z_in_uni = "Non-reciprocal producers",
                z_out_uni = "Non-reciprocal consumers", z_syn = "Syntrophic links",
                z_in_rate = "Mean supply per producer", z_out_rate = "Mean release per consumer",
                ext_res = "Chemostat inflow")
tabla1 <- data.frame(
  Predictor = nombres_t1,
  `Odds ratio of extinction` = ic(s1[names(nombres_t1), "Estimate"], s1[names(nombres_t1), "Std.err"]),
  `P value (extinction)` = fp(s1[names(nombres_t1), "Pr(>|W|)"]),
  `Change in abundance (%)` = ic(c2[names(nombres_t1), 1], c2[names(nombres_t1), 2], function(x) 100 * (exp(x) - 1)),
  `P value (abundance)` = fp(2 * pnorm(-abs(c2[names(nombres_t1), 3]))),
  check.names = FALSE, row.names = NULL)
write.csv(tabla1, file.path(RESULTADOS, "Tabla1_extincion_abundancia.csv"), row.names = FALSE)
dato("Tabla 1 guardada en", file.path(RESULTADOS, "Tabla1_extincion_abundancia.csv"))


# =============================================================================
# 3.1.2  Compartir recursos lleva a la exclusión competitiva
# =============================================================================
seccion("3.1.2 Exclusión competitiva")

tasas <- d %>% group_by(focal, comp) %>% summarise(tasa = mean(ext), .groups = "drop")
tasa <- function(f, c) tasas$tasa[tasas$focal == f & tasas$comp == c]
dato("Extinción sin competencia: competidoras / resto", sprintf("%s / %s", pct(tasa(1, 0)), pct(tasa(0, 0))))
dato("Odds ratio de las competidoras sin competencia (P)",
     sprintf("%s (P = %.2f)", ic(s1["focal", 1], s1["focal", 2]), s1["focal", 4]))
dato("Extinción de las competidoras con competencia", pct(tasa(1, 1)))
dato("Odds ratio de la interacción escenario × competidora", ic(s1["comp:focal", 1], s1["comp:focal", 2]))

destino <- d %>% filter(focal == 1) %>% group_by(comp, network_id) %>%
  summarise(n_ext = sum(ext), .groups = "drop")
tab_dest <- destino %>% count(comp, n_ext) %>% pivot_wider(names_from = comp, values_from = n, names_prefix = "comp_")
dato("Redes: ninguna / una / ambas competidoras extintas (con competencia)",
     paste(tab_dest$comp_1, collapse = " / "))
dato("Redes en que ambas persisten sin competencia", tab_dest$comp_0[tab_dest$n_ext == 0])

# Regresión logística condicional (Tabla S2) en las redes con exclusión
excl <- destino %>% filter(comp == 1, n_ext == 1) %>% pull(network_id)
ab0 <- d %>% filter(comp == 0) %>% select(network_id, species, ab_basal = abundance)
pares <- d %>% filter(comp == 1, focal == 1, network_id %in% excl) %>%
  inner_join(ab0, by = c("network_id", "species")) %>%
  mutate(sobrevive = 1L - ext, z_ab_basal = as.numeric(scale(ab_basal)))
m3 <- clogit(sobrevive ~ z_rep_rate + z_outgoing_resources + z_out_uni + z_syn + z_out_rate +
               z_ab_basal + strata(network_id), data = pares)
s3 <- summary(m3)
acierto <- pares %>% mutate(lp = as.numeric(predict(m3, type = "lp"))) %>% group_by(network_id) %>%
  summarise(ok = lp[sobrevive == 1] > lp[sobrevive == 0], .groups = "drop") %>% pull(ok) %>% mean()
regla_rep <- pares %>% group_by(network_id) %>%
  summarise(ok = rep_rate[sobrevive == 1] > rep_rate[sobrevive == 0], .groups = "drop") %>% pull(ok) %>% mean()
dato("Redes con exclusión", length(excl))
dato("Acierto del modelo al identificar a la superviviente", pct(acierto, 1))
dato("Odds ratio de sobrevivir: tasa de replicación (+1 DE)",
     ic(s3$coefficients["z_rep_rate", 1], s3$coefficients["z_rep_rate", 3]))
dato("Odds ratio de sobrevivir: liberación de subproductos (+1 DE)",
     ic(s3$coefficients["z_outgoing_resources", 1], s3$coefficients["z_outgoing_resources", 3]))
dato("Odds ratio de sobrevivir: enlaces sintróficos (P)",
     sprintf("%.2f (P = %.2f)", s3$coefficients["z_syn", 2], s3$coefficients["z_syn", 5]))
dato("Sobrevive la competidora con mayor tasa de replicación", pct(regla_rep, 1))

nombres_s2 <- c(z_rep_rate = "Replication rate", z_outgoing_resources = "Total by-product release",
                z_out_uni = "Non-reciprocal consumers", z_syn = "Syntrophic links",
                z_out_rate = "Mean release per consumer", z_ab_basal = "Abundance without competition")
tablaS2 <- data.frame(Trait = nombres_s2,
                      `Odds ratio of survival` = sprintf("%.2f (%.2f–%.2f)", s3$conf.int[names(nombres_s2), 1],
                                                         s3$conf.int[names(nombres_s2), 3], s3$conf.int[names(nombres_s2), 4]),
                      `P value` = fp(s3$coefficients[names(nombres_s2), 5]), check.names = FALSE, row.names = NULL)
write.csv(tablaS2, file.path(RESULTADOS, "TablaS2_logit_condicional.csv"), row.names = FALSE)

ganadoras <- pares %>% filter(sobrevive == 1, ab_basal > 0) %>% mutate(ratio = abundance / ab_basal)
dato("Aumento de abundancia de la superviviente (mediana del cociente)", sprintf("%.2f", median(ganadoras$ratio)))
dato("Aumento de abundancia estimado por el modelo (interacción)",
     sprintf("%+.1f%%", 100 * (exp(c2["comp:focal", 1]) - 1)))

# Coexistencia: asimetría en el índice de capacidad competitiva
cf3 <- coef(m3)
asim <- d %>% filter(comp == 0, focal == 1) %>%
  mutate(h = cf3["z_rep_rate"] * z_rep_rate + cf3["z_outgoing_resources"] * z_outgoing_resources +
           cf3["z_out_uni"] * z_out_uni) %>%
  group_by(network_id) %>% summarise(asim = abs(diff(h)), .groups = "drop") %>%
  inner_join(destino %>% filter(comp == 1), by = "network_id") %>% filter(n_ext <= 1)
mw <- wilcox.test(asim ~ n_ext, data = asim)
write.csv(asim %>% select(network_id, asim, n_ext), file.path(RESULTADOS, "asimetria_competidoras.csv"), row.names = FALSE)
dato("Redes en que ambas competidoras persisten", sum(asim$n_ext == 0))
dato("Asimetría mediana: coexisten / exclusión",
     sprintf("%.2f / %.2f", median(asim$asim[asim$n_ext == 0]), median(asim$asim[asim$n_ext == 1])))
dato("Test de Mann-Whitney, P", sprintf("%.1e", mw$p.value))
coex_tramos <- asim %>%
  mutate(tramo = cut(asim, c(-Inf, .5, 1, 2, 3, 5, Inf),
                     labels = c("< 0.5", "0.5-1", "1-2", "2-3", "3-5", "> 5"))) %>%
  group_by(tramo) %>% summarise(coexisten = mean(n_ext == 0), redes = n(), .groups = "drop")
dato("Coexistencia en el tramo más simétrico / más asimétrico",
     sprintf("%s / %s", pct(coex_tramos$coexisten[1], 1), pct(tail(coex_tramos$coexisten, 1), 1)))


# =============================================================================
# 3.1.3  Las especies con un único productor son las más vulnerables
# =============================================================================
seccion("3.1.3 Extinciones secundarias")

dato("Extinción del resto de especies: sin / con competencia", sprintf("%s / %s", pct(tasa(0, 0)), pct(tasa(0, 1))))
dato("Odds ratio del escenario con competencia", ic(s1["comp", 1], s1["comp", 2]))

sec <- d %>% filter(focal == 0) %>%
  select(network_id, species, incoming_degree, comp, ext) %>%
  pivot_wider(names_from = comp, values_from = ext, names_prefix = "ext_c") %>%
  mutate(secundaria = ext_c0 == 0 & ext_c1 == 1) %>%
  left_join(destino %>% filter(comp == 1) %>% select(network_id, n_ext), by = "network_id")
dato("Extinciones secundarias (total)", sum(sec$secundaria))
por_destino <- sec %>% group_by(n_ext) %>%
  summarise(ext = sum(secundaria), redes = n_distinct(network_id), .groups = "drop") %>%
  mutate(por_red = ext / redes)
dato("En redes con exclusión (total; por red)",
     sprintf("%d; %.3f", por_destino$ext[por_destino$n_ext == 1], por_destino$por_red[por_destino$n_ext == 1]))
dato("En redes con coexistencia (total; por red)",
     sprintf("%d; %.3f", por_destino$ext[por_destino$n_ext == 0], por_destino$por_red[por_destino$n_ext == 0]))
por_prod <- sec %>% filter(ext_c0 == 0, incoming_degree > 0) %>%
  mutate(productores = pmin(incoming_degree, 4)) %>% group_by(productores) %>%
  summarise(extintas = sum(secundaria), especies = n(), tasa = mean(secundaria), .groups = "drop")
for (i in seq_len(nrow(por_prod)))
  dato(sprintf("Extinción secundaria con %s productor(es)", ifelse(por_prod$productores[i] == 4, "4 o más", por_prod$productores[i])),
       sprintf("%s (%d de %d)", pct(por_prod$tasa[i], 1), por_prod$extintas[i], por_prod$especies[i]))
dato("Odds ratio de extinción por +1 DE de suministro medio por productor",
     ic(s1["z_in_rate", 1], s1["z_in_rate", 2]))

com <- d %>% group_by(network_id, comp) %>%
  summarise(riqueza = sum(abundance > 0), total = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = comp, values_from = c(riqueza, total), names_sep = "_c")
dato("Riqueza media: sin / con competencia", sprintf("%.2f / %.2f", mean(com$riqueza_c0), mean(com$riqueza_c1)))
dato("Abundancia total media: sin / con competencia (cambio)",
     sprintf("%.0f / %.0f (%+.1f%%)", mean(com$total_c0), mean(com$total_c1),
             100 * (mean(com$total_c1) / mean(com$total_c0) - 1)))
dato("Test de Wilcoxon pareado de la abundancia total, P",
     fp(wilcox.test(com$total_c0, com$total_c1, paired = TRUE)$p.value))

cat("\nResultados guardados en", salida, "\n")
