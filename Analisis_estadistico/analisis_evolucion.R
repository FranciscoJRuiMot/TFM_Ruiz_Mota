# =============================================================================
# Análisis de las comunidades evolutivas (apartado 3.2 del manuscrito)
#
# Entradas: results_final.csv  abundancia media, riesgo de extinción y tiempo de
#                              gestación de cada especie en las 10 comunidades
#                              ancestrales (evo_id = 0) y en las 10 réplicas
#                              evolutivas de cada una (evo_id = 1-10)
#           network_50.csv     universo de interacciones (72 x 72)
# Salidas : cifras del apartado 3.2 en consola y en resultados_evolucion.txt
#           Tablas 2 y 3 y Tabla suplementaria S1 en formato CSV
#           topologia_redes.csv y modelos_nulos_*.csv (los usa figuras_manuscrito.R)
#
# Ejecutar desde la carpeta raíz del repositorio:  Rscript analisis_evolucion.R
# Los modelos nulos (300 redes por red observada) tardan unos 10 minutos.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(lme4)      # modelos mixtos
  library(geepack)   # ecuaciones de estimación generalizadas
  library(igraph)    # propiedades de red
})

# ── Rutas y parámetros ───────────────────────────────────────────────────────
DATOS      <- "datos/results_final.csv"
UNIVERSO   <- "datos/network_50.csv"
RESULTADOS <- "resultados"
N_NULO     <- 300      # redes simuladas por red observada y modelo nulo
dir.create(RESULTADOS, showWarnings = FALSE)

# Especie que recibe el aporte del quimiostato en cada comunidad ancestral
QUIMIOSTATO <- c(`1` = "EQU_NOT", `2` = "OR_NOT", `3` = "ORN_NAND", `4` = "ANDN_NOT", `5` = "OR_NAND",
                 `6` = "NOT_NAND", `7` = "ANDN_NOT", `8` = "NAND_NOT", `9` = "NAND_NOT", `10` = "EQU_NOT")

# ── Utilidades para informar de los resultados ───────────────────────────────
salida <- file.path(RESULTADOS, "resultados_evolucion.txt")
cat("", file = salida)
seccion <- function(titulo) {
  txt <- paste0("\n", strrep("=", 78), "\n", titulo, "\n", strrep("=", 78), "\n")
  cat(txt); cat(txt, file = salida, append = TRUE)
}
dato <- function(etiqueta, valor) {
  txt <- paste0("  ", format(etiqueta, width = 64), " ", valor, "\n")
  cat(txt); cat(txt, file = salida, append = TRUE)
}
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
fp  <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))
ic  <- function(b, se) sprintf("%.2f (%.2f–%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))


# =============================================================================
# 1. Preparación de los datos
# =============================================================================
U    <- as.matrix(read.csv(UNIVERSO, row.names = 1, check.names = FALSE))
SP   <- rownames(U)
PROD <- rowSums(U)                                    # subproductos liberados por cada especie

d <- read.csv(DATOS, stringsAsFactors = FALSE) %>%
  mutate(sp      = sub("+", "_", species, fixed = TRUE),   # X+Y en los datos = X_Y en el universo
         rep     = 30 / mean_gest_time,                   # tasa de replicación (eficacia en Avida)
         present = mean_abundance > 0,
         stable  = extinction_risk == 0,                  # presente en todos los registros finales
         fase    = ifelse(evo_id == 0, "ancestral", "evolucionada"))

chem   <- data.frame(eco_id = as.integer(names(QUIMIOSTATO)), chem = unname(QUIMIOSTATO))
S9     <- d %>% filter(evo_id == 0, present)               # 9 especies de cada comunidad ancestral
S18    <- d %>% filter(evo_id == 0)                        # 18 especies submuestreadas
anc_sp <- S9 %>% select(eco_id, sp)

# Suministro que reciben las especies x desde las especies presentes S con abundancias w
suministro <- function(S, w, x = S) as.numeric(t(U[S, x, drop = FALSE]) %*% w)

# Especies de cada comunidad evolucionada clasificadas por su origen
ev <- d %>% filter(evo_id > 0, present) %>%
  inner_join(chem, by = "eco_id") %>%
  left_join(anc_sp %>% mutate(anc = TRUE), by = c("eco_id", "sp")) %>%
  mutate(origen = factor(case_when(sp == chem ~ "quimiostato", !is.na(anc) ~ "residente", TRUE ~ "nueva"),
                         c("quimiostato", "residente", "nueva")))

# Propiedades de cada comunidad
metricas_comunidad <- function(S, w, ch) {
  p <- w / sum(w)
  data.frame(S = length(S),
             cap = sum(rowSums(U[S, S, drop = FALSE]) * w) / sum(PROD[S] * w),  # captura ponderada por abundancia
             tot = sum(w),                                                      # abundancia total
             share_chem = w[S == ch] / sum(w),                                  # fracción de la especie del quimiostato
             J = -sum(p * log(p)) / log(length(S)))                             # equitatividad de Pielou
}
M <- d %>% filter(present) %>% inner_join(chem, by = "eco_id") %>%
  group_by(eco_id, evo_id, fase) %>%
  group_modify(~ metricas_comunidad(.x$sp, .x$mean_abundance, .x$chem[1])) %>% ungroup()
Me <- M %>% filter(fase == "evolucionada")

# Test de Wilcoxon pareado: cada comunidad ancestral frente a la media de sus réplicas
media_eco <- function(tab, vars) tab %>% group_by(eco_id, fase) %>% summarise(across(all_of(vars), mean), .groups = "drop")
p_pareado <- function(tab_eco, v) {
  x <- tab_eco %>% select(eco_id, fase, all_of(v)) %>% pivot_wider(names_from = fase, values_from = all_of(v))
  if (sd(x$evolucionada - x$ancestral) == 0) return(NA)
  wilcox.test(x$ancestral, x$evolucionada, paired = TRUE, exact = FALSE)$p.value
}
Mf <- media_eco(M, c("S", "cap", "tot", "share_chem", "J"))


# =============================================================================
# 3.2.1  Las especies nuevas aumentan la riqueza y la abundancia total
# =============================================================================
seccion("3.2.1 Riqueza y abundancia total")

retencion <- ev %>% filter(origen != "nueva") %>% count(eco_id, evo_id, name = "retenidas")
persist   <- d %>% filter(evo_id > 0, present) %>% group_by(eco_id, evo_id) %>% summarise(n = sum(stable), .groups = "drop")
dominante <- ev %>% group_by(eco_id, evo_id) %>% slice_max(mean_abundance, n = 1, with_ties = FALSE) %>% ungroup()
dato("Riqueza tras la evolución: media (rango)", sprintf("%.1f (%d-%d)", mean(Me$S), min(Me$S), max(Me$S)))
dato("Especies presentes en todos los registros finales (media)", sprintf("%.1f", mean(persist$n)))
dato("Especies ancestrales retenidas: media (mínimo)", sprintf("%.2f (%d)", mean(retencion$retenidas), min(retencion$retenidas)))
dato("Experimentos en que domina la especie del quimiostato", sprintf("%d de %d", sum(dominante$sp == dominante$chem), nrow(dominante)))
dato("Fracción de la abundancia de esa especie: ancestral / evolucionada",
     sprintf("%.2f / %.2f", mean(Mf$share_chem[Mf$fase == "ancestral"]), mean(Mf$share_chem[Mf$fase == "evolucionada"])))
dato("Equitatividad de Pielou: ancestral / evolucionada",
     sprintf("%.2f / %.2f", mean(Mf$J[Mf$fase == "ancestral"]), mean(Mf$J[Mf$fase == "evolucionada"])))
t0 <- mean(Mf$tot[Mf$fase == "ancestral"]); t1 <- mean(Mf$tot[Mf$fase == "evolucionada"])
dato("Abundancia total: ancestral / evolucionada (cambio)", sprintf("%.0f / %.0f (%+.0f%%)", t0, t1, 100 * (t1 / t0 - 1)))
dato("Test de Wilcoxon pareado, P", sprintf("%.4f", p_pareado(Mf, "tot")))
dato("Captura ponderada de subproductos: ancestral / evolucionada",
     sprintf("%.3f / %.3f", mean(Mf$cap[Mf$fase == "ancestral"]), mean(Mf$cap[Mf$fase == "evolucionada"])))
dato("Spearman entre captura y abundancia total (evolucionadas)", sprintf("%.2f", cor(Me$cap, Me$tot, method = "spearman")))


# =============================================================================
# 3.2.2  Evolución paralela
# =============================================================================
seccion("3.2.2 Repetibilidad y evolución paralela")

m_S <- lmer(S ~ 1 + (1 | eco_id), data = Me)
v_S <- as.data.frame(VarCorr(m_S))$vcov
dato("Varianza de la riqueza explicada por la comunidad ancestral", pct(v_S[1] / sum(v_S)))

# Similitud de Jaccard entre pares de experimentos
sets <- d %>% filter(evo_id > 0, present) %>% group_by(eco_id, evo_id) %>%
  summarise(s = list(sp), .groups = "drop") %>% arrange(eco_id, evo_id)
anc_set <- split(anc_sp$sp, anc_sp$eco_id)
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
n <- nrow(sets); Jall <- Jnew <- matrix(NA, n, n)
for (i in 1:n) for (j in 1:n) {
  Jall[i, j] <- jac(sets$s[[i]], sets$s[[j]])
  a <- setdiff(sets$s[[i]], anc_set[[as.character(sets$eco_id[i])]])
  b <- setdiff(sets$s[[j]], anc_set[[as.character(sets$eco_id[j])]])
  Jnew[i, j] <- if (length(union(a, b)) > 0) jac(a, b) else NA
}
lab <- sets$eco_id; iu <- upper.tri(Jall); misma <- outer(lab, lab, "==")[iu]
dif_jac <- function(J, l) { s <- outer(l, l, "==")[iu]; mean(J[iu][s], na.rm = TRUE) - mean(J[iu][!s], na.rm = TRUE) }
set.seed(1)
nulo_jac <- replicate(2000, dif_jac(Jnew, sample(lab)))
dato("Jaccard, todas las especies: misma comunidad / distinta",
     sprintf("%.2f / %.2f", mean(Jall[iu][misma]), mean(Jall[iu][!misma])))
dato("Jaccard, especies nuevas: misma comunidad / distinta",
     sprintf("%.2f / %.2f", mean(Jnew[iu][misma], na.rm = TRUE), mean(Jnew[iu][!misma], na.rm = TRUE)))
dato("Test de permutación (2.000 permutaciones), P",
     sprintf("%.4f", (sum(nulo_jac >= dif_jac(Jnew, lab)) + 1) / 2001))

# Evolución paralela: especie nueva presente en al menos 5 de las 10 réplicas
paralelas <- ev %>% filter(origen == "nueva") %>% count(eco_id, sp, name = "replicas")
frac_par <- function(x) sum(x[x >= 5]) / sum(x)
por_eco <- paralelas %>% group_by(eco_id) %>% summarise(f = frac_par(replicas), .groups = "drop")
entre <- paralelas %>% filter(replicas >= 5) %>% count(sp, name = "comunidades") %>% arrange(desc(comunidades))
dato("Apariciones de especies nuevas en paralelo (total)", pct(frac_par(paralelas$replicas)))
dato("Apariciones en paralelo: rango entre comunidades", sprintf("%s-%s", pct(min(por_eco$f), 0), pct(max(por_eco$f), 0)))
dato("Especies en paralelo en al menos 5 comunidades", sum(entre$comunidades >= 5))
dato("Especie en paralelo en más comunidades", sprintf("%s (%d de 10)", entre$sp[1], entre$comunidades[1]))

# Índice de repetibilidad de Venkataram y Kryazhimskiy (2023)
indice_rep <- function(dd) {
  X <- dd %>% group_by(eco_id, evo_id) %>% mutate(p = mean_abundance / sum(mean_abundance)) %>% ungroup() %>%
    select(eco_id, evo_id, sp, p) %>% pivot_wider(names_from = sp, values_from = p, values_fill = 0) %>%
    arrange(eco_id, evo_id)
  l <- X$eco_id; S <- 1 - as.matrix(dist(as.matrix(X[, -(1:2)]))) / sqrt(2)
  u <- upper.tri(S); s <- outer(l, l, "==")[u]
  c(dentro = mean(S[u][s]), entre = mean(S[u][!s]))
}
r_con <- indice_rep(ev); r_sin <- indice_rep(ev %>% filter(origen != "quimiostato"))
dato("Índice de repetibilidad: misma comunidad / distinta", sprintf("%.2f / %.2f", r_con["dentro"], r_con["entre"]))
dato("Ídem sin la especie del quimiostato", sprintf("%.2f / %.2f", r_sin["dentro"], r_sin["entre"]))


# =============================================================================
# 3.2.3  Las especies nuevas son funcionalmente simples y ocupan nichos vacíos
# =============================================================================
seccion("3.2.3 Complejidad booleana y emergencia")

# Complejidad: mínimo número de operaciones NAND para calcular cada función (Lenski et al., 2003)
CX <- c(NOT = 1, NAND = 1, AND = 2, ORN = 2, OR = 3, ANDN = 3, NOR = 4, XOR = 4, EQU = 5)
fun <- do.call(rbind, strsplit(SP, "_"))
spinfo <- data.frame(sp = SP, cx = CX[fun[, 1]] + CX[fun[, 2]],
                     indeg = colSums(U > 0), outdeg = rowSums(U > 0))

occ <- d %>% filter(evo_id > 0, present) %>% count(sp, name = "runs") %>%
  right_join(spinfo, by = "sp") %>% mutate(runs = coalesce(runs, 0L))
rho <- cor.test(occ$runs, occ$cx, method = "spearman", exact = FALSE)
dato("Especies del universo presentes en algún experimento", sum(occ$runs > 0))
dato("Spearman entre frecuencia y complejidad (P)", sprintf("%.2f (%.0e)", rho$estimate, rho$p.value))
dato("Presencia máxima de una especie (experimentos)", max(occ$runs))

# Candidatas: las 63 especies ausentes de cada comunidad ancestral, en cada experimento
cand <- d %>% filter(evo_id > 0) %>% distinct(eco_id, evo_id) %>% group_by(eco_id, evo_id) %>%
  group_modify(~ {
    a9 <- S9 %>% filter(eco_id == .y$eco_id); x <- setdiff(SP, a9$sp)
    pres <- d %>% filter(eco_id == .y$eco_id, evo_id == .y$evo_id, present)
    data.frame(sp = x,
               emerge    = as.integer(x %in% pres$sp),
               establece = as.integer(x %in% pres$sp[pres$stable]),
               sup0      = suministro(a9$sp, a9$mean_abundance, x),          # suministro potencial
               excluida  = as.integer(x %in% S18$sp[S18$eco_id == .y$eco_id]))
  }) %>% ungroup() %>% inner_join(spinfo, by = "sp") %>%
  mutate(z_cx = as.numeric(scale(cx)), z_sup = as.numeric(scale(log(sup0 + 0.1))),
         z_in = as.numeric(scale(indeg)), z_out = as.numeric(scale(outdeg)))
em_cx <- cand %>% group_by(cx) %>% summarise(p = mean(emerge))
dato("Candidatas (especie x experimento)", nrow(cand))
dato("Emergen: complejidad mínima / máxima (8-9)",
     sprintf("%s / %s", pct(em_cx$p[em_cx$cx == min(em_cx$cx)], 0), pct(mean(cand$emerge[cand$cx >= 8]), 0)))

m_em <- glmer(emerge ~ z_cx + z_sup + excluida + z_in + z_out + (1 | eco_id / evo_id),
              data = cand, family = binomial, control = glmerControl(optimizer = "bobyqa"))
m_es <- update(m_em, establece ~ .)
ce <- summary(m_em)$coefficients; cs <- summary(m_es)$coefficients
nombres_t2 <- c(z_cx = "Boolean complexity", z_sup = "Potential resource supply from the ancestral community",
                excluida = "Excluded during ecological assembly",
                z_in = "Potential producers in the interaction universe",
                z_out = "Potential consumers in the interaction universe")
tabla2 <- data.frame(Predictor = nombres_t2,
                     `Odds ratio of emergence` = ic(ce[names(nombres_t2), 1], ce[names(nombres_t2), 2]),
                     `P value (emergence)` = fp(ce[names(nombres_t2), 4]),
                     `Odds ratio of establishment` = ic(cs[names(nombres_t2), 1], cs[names(nombres_t2), 2]),
                     `P value (establishment)` = fp(cs[names(nombres_t2), 4]), check.names = FALSE, row.names = NULL)
write.csv(tabla2, file.path(RESULTADOS, "Tabla2_emergencia.csv"), row.names = FALSE)
for (k in names(nombres_t2)[c(1, 2, 4)])
  dato(paste("OR emergencia / establecimiento:", nombres_t2[k]),
       sprintf("%s / %s", ic(ce[k, 1], ce[k, 2]), ic(cs[k, 1], cs[k, 2])))

r_anc <- with(S18 %>% filter(!is.na(rep)) %>% inner_join(spinfo, by = "sp"), cor.test(rep, cx, method = "spearman", exact = FALSE))
r_evo <- with(ev %>% filter(!is.na(rep)) %>% inner_join(spinfo, by = "sp"), cor.test(rep, cx, method = "spearman", exact = FALSE))
dato("Spearman tasa de replicación-complejidad: genomas ancestrales (P)", sprintf("%.2f (%.2f)", r_anc$estimate, r_anc$p.value))
dato("Spearman tasa de replicación-complejidad: comunidades evolucionadas", sprintf("%.2f", r_evo$estimate))


# =============================================================================
# 3.2.4  Las especies nuevas son raras pese a replicarse más rápido
# =============================================================================
seccion("3.2.4 Propiedades de las especies nuevas")

evs <- ev %>% filter(origen != "quimiostato") %>% mutate(nueva = as.integer(origen == "nueva"))
# (los avisos de ajuste singular indican que la varianza entre réplicas es prácticamente nula)
m_ab  <- suppressMessages(lmer(log(mean_abundance) ~ nueva + (1 | eco_id / evo_id), data = evs))
m_rep <- suppressMessages(lmer(log(rep) ~ nueva + (1 | eco_id / evo_id), data = evs %>% filter(!is.na(rep))))
m_st  <- glmer(stable ~ nueva + (1 | eco_id / evo_id), data = evs, family = binomial,
               control = glmerControl(optimizer = "bobyqa"))
efecto <- function(m) { s <- summary(m)$coefficients["nueva", ]; c(b = s[[1]], p = 2 * pnorm(-abs(s[[1]] / s[[2]]))) }
e_ab <- efecto(m_ab); e_rep <- efecto(m_rep); e_st <- efecto(m_st)
dato("Abundancia de las nuevas frente a residentes (P)", sprintf("%+.0f%% (%s)", 100 * (exp(e_ab["b"]) - 1), fp(e_ab["p"])))
dato("Tasa de replicación de las nuevas frente a residentes (P)", sprintf("%+.0f%% (%s)", 100 * (exp(e_rep["b"]) - 1), fp(e_rep["p"])))
dato("Persistentes: nuevas / residentes (P)",
     sprintf("%s / %s (%s)", pct(mean(evs$stable[evs$nueva == 1])), pct(mean(evs$stable[evs$nueva == 0])), fp(e_st["p"])))

sup_evo <- ev %>% group_by(eco_id, evo_id) %>% mutate(sup = suministro(sp, mean_abundance)) %>% ungroup() %>%
  filter(origen != "quimiostato", sup > 0, !is.na(rep))
m_sup <- lmer(log(mean_abundance) ~ log(sup) + log(rep) + (1 | eco_id / evo_id), data = sup_evo)
cf <- summary(m_sup)$coefficients
dato("Especies analizadas", nrow(sup_evo))
dato("Elasticidad del suministro (estimación ± EE)", sprintf("%.2f ± %.3f", cf[2, 1], cf[2, 2]))
dato("Elasticidad de la tasa de replicación", sprintf("%.2f ± %.2f", cf[3, 1], cf[3, 2]))
dato("R2 (ajustados frente a observados)", sprintf("%.2f", cor(fitted(m_sup), log(sup_evo$mean_abundance))^2))

directa <- function(ch, sp) U[cbind(ch, sp)] > 0     # ¿consume subproductos de la especie del quimiostato?
pos <- d %>% filter(present) %>% inner_join(chem, by = "eco_id") %>% filter(sp != chem) %>%
  left_join(anc_sp %>% mutate(anc = TRUE), by = c("eco_id", "sp")) %>%
  mutate(grupo = case_when(evo_id == 0 ~ "residente", !is.na(anc) ~ "residente_evo", TRUE ~ "nueva"),
         dir = directa(chem, sp))
azar <- mean(sapply(unique(chem$chem), function(c) mean(U[c, setdiff(SP, c)] > 0)))
dato("Consumidoras directas del quimiostato: residentes / nuevas / azar",
     sprintf("%s / %s / %s", pct(mean(pos$dir[pos$grupo == "residente"]), 0),
             pct(mean(pos$dir[pos$grupo == "nueva"]), 0), pct(azar, 0)))


# =============================================================================
# 3.2.5  Cambios en la tasa de replicación de las especies presentes
# =============================================================================
seccion("3.2.5 Tasa de replicación de las especies residentes")

nuevas_exp <- ev %>% group_by(eco_id, evo_id) %>% summarise(n_new = sum(origen == "nueva"), .groups = "drop")
tasa <- S9 %>% select(eco_id, sp, rep_anc = rep, ab_anc = mean_abundance) %>%
  inner_join(ev %>% select(eco_id, evo_id, sp, rep_evo = rep, origen, chem), by = c("eco_id", "sp")) %>%
  filter(!is.na(rep_anc), !is.na(rep_evo)) %>%
  inner_join(nuevas_exp, by = c("eco_id", "evo_id")) %>%
  mutate(ratio = rep_evo / rep_anc, sube = as.integer(ratio > 1.05))

# Fracción del suministro de cada residente que procede de especies nuevas
frac_nuevas <- ev %>% group_by(eco_id, evo_id) %>% group_modify(~ {
  tot <- suministro(.x$sp, .x$mean_abundance)
  nv <- .x$origen == "nueva"
  de_nuevas <- if (any(nv)) suministro(.x$sp[nv], .x$mean_abundance[nv], .x$sp) else 0
  data.frame(sp = .x$sp, frac_new = ifelse(tot > 0, de_nuevas / tot, NA))
}) %>% ungroup()
tasa <- tasa %>% left_join(frac_nuevas, by = c("eco_id", "evo_id", "sp"))

q <- tasa %>% filter(origen == "quimiostato")
dato("Especie del quimiostato: aumenta / aumenta más de un 5%", sprintf("%s / %s", pct(mean(q$ratio > 1), 0), pct(mean(q$sube), 0)))
dato("Especie del quimiostato: cociente mediano", sprintf("%.2f", median(q$ratio)))
dato("Especie del quimiostato: Spearman cociente-nº de especies nuevas",
     sprintf("%.2f", cor(q$ratio, q$n_new, method = "spearman")))

r <- tasa %>% filter(origen == "residente") %>%
  mutate(directa = as.integer(directa(chem, sp)), eco = factor(eco_id),
         experimento = interaction(eco_id, evo_id, drop = TRUE),
         z_ab = as.numeric(scale(log(ab_anc))), z_new = as.numeric(scale(n_new))) %>%
  arrange(experimento)
dato("Residentes analizadas (especie x experimento)", nrow(r))
dato("Aumentan más de un 5%: consumidoras directas / resto",
     sprintf("%s / %s", pct(mean(r$sube[r$directa == 1])), pct(mean(r$sube[r$directa == 0]))))
dato("Cociente mediano cuando aumentan más de un 5%", sprintf("%.2f", median(r$ratio[r$sube == 1])))
dato("Conservan exactamente su tasa ancestral", pct(mean(abs(r$ratio - 1) < 1e-9), 0))

# GEE con la comunidad ancestral como efecto fijo y cada experimento como clúster
g <- geeglm(sube ~ directa + z_ab + z_new + eco, id = experimento, data = r,
            family = binomial, corstr = "independence")
sg <- summary(g)$coefficients
dato("OR consumidora directa, dentro de cada comunidad (P)", sprintf("%.1f (%s)", exp(sg["directa", 1]), fp(sg["directa", 4])))
dato("OR por +1 DE del nº de especies nuevas (P)", sprintf("%.2f (%s)", exp(sg["z_new", 1]), fp(sg["z_new", 4])))
for (grupo in c(1, 0)) {
  sub <- r %>% filter(directa == grupo) %>% mutate(eco = droplevels(eco))
  gg <- summary(geeglm(sube ~ z_new + eco, id = experimento, data = sub, family = binomial, corstr = "independence"))$coefficients
  dato(sprintf("OR nº de especies nuevas: %s (P)", ifelse(grupo == 1, "consumidoras directas", "resto de residentes")),
       sprintf("%.2f (%s)", exp(gg["z_new", 1]), fp(gg["z_new", 4])))
}
dato("Spearman cociente-fracción del suministro procedente de nuevas",
     sprintf("%.2f", cor(r$ratio, r$frac_new, method = "spearman", use = "complete.obs")))


# =============================================================================
# 3.2.6  La emergencia de especies nuevas altera la topología de la red
# =============================================================================
seccion("3.2.6 Topología de las redes")

nodf <- function(M) {                    # anidamiento (Almeida-Neto et al., 2008), filas y columnas
  pares <- function(X) {
    k <- rowSums(X); O <- X %*% t(X); ii <- which(upper.tri(O), arr.ind = TRUE)
    ki <- k[ii[, 1]]; kj <- k[ii[, 2]]; kmin <- pmin(ki, kj)
    v <- ifelse(ki != kj & kmin > 0, O[ii] / kmin, 0); c(sum(v), length(v))
  }
  r <- pares(M); cc <- pares(t(M)); 100 * (r[1] + cc[1]) / (r[2] + cc[2])
}
incoherencia <- function(W) {           # incoherencia trófica (MacKay et al., 2020)
  s_in <- colSums(W); s_out <- rowSums(W)
  h <- as.numeric(MASS::ginv(diag(s_in + s_out) - W - t(W)) %*% (s_in - s_out))
  sum(W * (outer(h, h, function(a, b) b - a) - 1)^2) / sum(W)
}
metricas_topo <- function(S) {
  W <- U[S, S, drop = FALSE]; B <- (W > 0) * 1; n <- length(S); L <- sum(B)
  G <- graph_from_adjacency_matrix(B, mode = "directed"); Gu <- as.undirected(G, mode = "collapse")
  kin <- colSums(B); kout <- rowSums(B); tc <- triad_census(G); conn <- sum(tc[-1])
  c(C = L / (n * (n - 1)), rec = sum(B * t(B)) / L,
    cv_in = sd(kin) / mean(kin), cv_out = sd(kout) / mean(kout),
    assort = assortativity_degree(G, directed = TRUE), trans = transitivity(Gu, type = "global"),
    path = mean_distance(G, directed = TRUE), nodf = nodf(B), modul = modularity(cluster_louvain(Gu)),
    scc = max(components(G, mode = "strong")$csize) / n,
    ffl = tc[9] / conn, cyc = tc[10] / conn, full = tc[16] / conn,   # motivos de tres especies
    capu = sum(W) / sum(PROD[S]), F0 = incoherencia(W))
}

set.seed(2026)
TOPO <- d %>% filter(present) %>% group_by(eco_id, evo_id, fase) %>%
  group_modify(~ {
    S <- .x$sp; G <- graph_from_adjacency_matrix((U[S, S] > 0) * 1, mode = "directed")
    as.data.frame(t(c(metricas_topo(S), LS = sum(U[S, S] > 0) / length(S), vconn = vertex_connectivity(G))))
  }) %>% ungroup()
write.csv(TOPO, file.path(RESULTADOS, "topologia_redes.csv"), row.names = FALSE)

# Tabla 3: propiedades estructurales y funcionales
filas_t3 <- list(
  S = "Number of species", C = "Connectance", LS = "Links per species", vconn = "Vertex connectivity",
  rec = "Reciprocity", cv_in = "In-degree heterogeneity", cv_out = "Out-degree heterogeneity",
  assort = "Degree assortativity", trans = "Transitivity", path = "Mean path length", nodf = "Nestedness",
  modul = "Modularity", scc = "Largest strongly connected component (fraction of species)",
  ffl = "Feed-forward loops (fraction of connected triads)", cyc = "Three-species cycles (fraction of connected triads)",
  full = "Fully reciprocal triads (fraction of connected triads)", F0 = "Trophic incoherence",
  capu = "By-product capture", cap = "By-product capture weighted by producer abundance",
  tot = "Total abundance", share_chem = "Share of the species receiving the chemostat inflow",
  J = "Evenness of abundances (Pielou index)")
Tm <- media_eco(TOPO, setdiff(names(TOPO), c("eco_id", "evo_id", "fase"))) %>%
  inner_join(Mf, by = c("eco_id", "fase"))
tabla3 <- data.frame(Property = unlist(filas_t3),
                     Ancestral = sapply(names(filas_t3), function(v) mean(Tm[[v]][Tm$fase == "ancestral"])),
                     Evolved = sapply(names(filas_t3), function(v) mean(Tm[[v]][Tm$fase == "evolucionada"])),
                     `P value` = sapply(names(filas_t3), function(v) p_pareado(Tm, v)),
                     check.names = FALSE, row.names = NULL)
write.csv(tabla3, file.path(RESULTADOS, "Tabla3_propiedades_redes.csv"), row.names = FALSE)
for (v in c("LS", "capu", "vconn", "C")) dato(paste(filas_t3[[v]], ": ancestral / evolucionada"),
  sprintf("%.3f / %.3f (P = %.4f)", tabla3$Ancestral[names(filas_t3) == v], tabla3$Evolved[names(filas_t3) == v],
          tabla3$`P value`[names(filas_t3) == v]))
dato("Conectancia del universo", sprintf("%.3f", sum(U > 0) / (72 * 71)))

# Coherencia del cambio entre comunidades y convergencia de las propiedades intensivas
signif <- tabla3 %>% mutate(v = names(filas_t3)) %>%
  filter(`P value` < 0.05, v %in% names(TOPO), v != "S") %>% pull(v)
anc_m <- Tm %>% filter(fase == "ancestral") %>% arrange(eco_id)
evo_m <- Tm %>% filter(fase == "evolucionada") %>% arrange(eco_id)
coherencia <- sapply(signif, function(v) max(sum(evo_m[[v]] > anc_m[[v]]), sum(evo_m[[v]] < anc_m[[v]])))
dato("Propiedades con cambio significativo", paste(signif, collapse = ", "))
dato("Mínimo de comunidades con la misma dirección de cambio", sprintf("%d de 10", min(coherencia)))
for (v in c("C", "trans")) dato(paste("Desviación típica entre comunidades:", filas_t3[[v]]),
  sprintf("%.3f -> %.3f", sd(anc_m[[v]]), sd(evo_m[[v]])))
for (v in c("trans", "ffl")) dato(paste("Comunidades en que aumenta / disminuye:", filas_t3[[v]]),
  sprintf("%d / %d", sum(evo_m[[v]] > anc_m[[v]]), sum(evo_m[[v]] < anc_m[[v]])))

# Tabla suplementaria S1: definiciones de las propiedades estructurales
tablaS1 <- data.frame(Property = c("Connectance", "Links per species", "Vertex connectivity", "Reciprocity",
  "In-degree heterogeneity", "Out-degree heterogeneity", "Degree assortativity", "Transitivity", "Mean path length",
  "Nestedness", "Modularity", "Largest strongly connected component", "Feed-forward loops", "Three-species cycles",
  "Fully reciprocal triads", "By-product capture", "Trophic incoherence"),
  `Code variable` = c("C", "LS", "vconn", "rec", "cv_in", "cv_out", "assort", "trans", "path", "nodf", "modul",
                      "scc", "ffl", "cyc", "full", "capu", "F0"), check.names = FALSE)
write.csv(tablaS1, file.path(RESULTADOS, "TablaS1_propiedades_red_variables.csv"), row.names = FALSE)

# Modelos nulos: A = especies al azar del universo (siempre con la del quimiostato);
# B = residentes + nuevas al azar; C = residentes + nuevas con probabilidad
# proporcional al suministro potencial desde la comunidad ancestral
MET_NULO <- setdiff(names(metricas_topo(SP[1:10])), "scc")
S9ab <- S9 %>% select(eco_id, sp, mean_abundance)
set.seed(11)
Z <- d %>% filter(present) %>% group_by(eco_id, evo_id) %>% group_modify(~ {
  e <- .y$eco_id; k <- .y$evo_id; S <- .x$sp; n <- length(S); ch <- QUIMIOSTATO[as.character(e)]
  obs <- metricas_topo(S)[MET_NULO]
  zf <- function(sim) (obs - colMeans(sim)) / apply(sim, 2, sd)
  simA <- t(replicate(N_NULO, metricas_topo(c(ch, sample(setdiff(SP, ch), n - 1)))[MET_NULO]))
  out <- data.frame(metrica = MET_NULO, nulo = "A", z = zf(simA))
  if (k > 0) {
    a <- S9ab %>% filter(eco_id == e); ret <- intersect(S, a$sp); pool <- setdiff(SP, a$sp)
    nn <- n - length(ret); w <- suministro(a$sp, a$mean_abundance, pool) + 1e-6
    simB <- t(replicate(N_NULO, metricas_topo(c(ret, sample(pool, nn)))[MET_NULO]))
    simC <- t(replicate(N_NULO, metricas_topo(c(ret, sample(pool, nn, prob = w)))[MET_NULO]))
    out <- rbind(out, data.frame(metrica = MET_NULO, nulo = "B", z = zf(simB)),
                      data.frame(metrica = MET_NULO, nulo = "C", z = zf(simC)))
  }
  out
}) %>% ungroup() %>% mutate(fase = ifelse(evo_id == 0, "ancestral", "evolucionada"))
write.csv(Z, file.path(RESULTADOS, "modelos_nulos_por_red.csv"), row.names = FALSE)

# Desviación media por comunidad, test t frente a cero y corrección de Benjamini-Hochberg
Zs <- Z %>% filter(is.finite(z)) %>%
  group_by(fase, nulo, metrica, eco_id) %>% summarise(z = mean(z), .groups = "drop") %>%
  group_by(fase, nulo, metrica) %>% summarise(z_media = mean(z), p = t.test(z)$p.value, .groups = "drop") %>%
  group_by(fase, nulo) %>% mutate(p_ajustada = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(comparacion = ifelse(fase == "ancestral", "Ancestral vs null A", paste("Evolved vs null", nulo)))
write.csv(Zs, file.path(RESULTADOS, "modelos_nulos_resumen.csv"), row.names = FALSE)

tabla_z <- Zs %>% mutate(celda = sprintf("%+.2f%s", z_media, ifelse(p_ajustada < 0.05, "*", ""))) %>%
  select(metrica, comparacion, celda) %>% pivot_wider(names_from = comparacion, values_from = celda)
txt <- paste(capture.output(print(as.data.frame(tabla_z), row.names = FALSE)), collapse = "\n")
cat("\n  Desviación estandarizada media frente a cada modelo nulo (* P ajustada < 0.05)\n", txt, "\n")
cat("\n  Desviación estandarizada media frente a cada modelo nulo (* P ajustada < 0.05)\n", txt, "\n",
    file = salida, append = TRUE)

cat("\nResultados guardados en", salida, "\n")
