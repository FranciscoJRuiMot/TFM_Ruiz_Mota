# =============================================================================
# Figuras del manuscrito
#
# Figura 1 y Figura suplementaria S1: comunidades ecológicas (apartado 3.1)
# Figuras 2-5 y Figura suplementaria S2: comunidades evolutivas (apartado 3.2)
#
# Necesita los datos de simulación y algunos resultados intermedios que generan
# los scripts de análisis, así que hay que ejecutar antes:
#   Rscript analisis_ecologia.R
#   Rscript analisis_evolucion.R
# Ejecutar desde la carpeta raíz del repositorio:  Rscript figuras_manuscrito.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

# ── Rutas ────────────────────────────────────────────────────────────────────
DATOS_ECO  <- "datos/results_1000_with_and_without.csv"
DATOS_EVO  <- "datos/results_final.csv"
UNIVERSO   <- "datos/network_50.csv"
RESULTADOS <- "resultados"
FIGURAS    <- "figuras"
dir.create(FIGURAS, showWarnings = FALSE)

# Especie que recibe el aporte del quimiostato en cada comunidad ancestral
QUIMIOSTATO <- c(`1` = "EQU_NOT", `2` = "OR_NOT", `3` = "ORN_NAND", `4` = "ANDN_NOT", `5` = "OR_NAND",
                 `6` = "NOT_NAND", `7` = "ANDN_NOT", `8` = "NAND_NOT", `9` = "NAND_NOT", `10` = "EQU_NOT")

# ── Estilo común ─────────────────────────────────────────────────────────────
theme_set(theme_classic(base_size = 11) + theme(plot.tag = element_text(face = "bold", size = 13),
          legend.position = "top", legend.title = element_blank(), axis.title = element_text(size = 9.5)))
C0 <- "#4C72B0"; C1 <- "#C44E52"; GR <- "grey45"                                   # sin / con competencia
C_ANC <- "#4C72B0"; C_NEW <- "#55A868"; C_CHEM <- "#DD8452"; C_EVO <- "#C44E52"    # evolución
LN <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE) # etiquetas de ejes logarítmicos

# =============================================================================
# Figura 1 y Figura S1: comunidades ecológicas
# =============================================================================
raw <- read.csv(DATOS_ECO)
b <- raw %>% filter(network_competition == 0) %>%
  transmute(network_id, species, focal = sp_competition, ext_res = incoming_external_res, rep_rate,
            incoming_resources, outgoing_resources, incoming_degree, out_uni = outgoing_degree - syntrophy_links)
d <- raw %>% select(network_id, species, comp = network_competition, abundance) %>%
  inner_join(b, by = c("network_id", "species")) %>% mutate(ext = as.integer(abundance == 0))
x <- d %>% filter(comp == 0, abundance > 0, incoming_resources > 0, outgoing_resources > 0, ext_res == 0) %>%
  mutate(ratio = incoming_resources / outgoing_resources)
e1a <- ggplot(x, aes(ratio, abundance)) + geom_point(alpha = .1, size = .5, colour = C0) +
  geom_smooth(method = "lm", colour = "black", linewidth = .7, se = FALSE) + scale_x_log10(labels = LN) + scale_y_log10(labels = LN) +
  labs(x = "Supply relative to by-product release", y = "Abundance at steady state")
e1b <- ggplot(x, aes(rep_rate, abundance)) + geom_point(alpha = .1, size = .5, colour = C0) +
  geom_smooth(method = "lm", colour = "black", linewidth = .7, se = FALSE) + scale_y_log10(labels = LN) +
  labs(x = "Replication rate", y = "Abundance at steady state")
t22 <- d %>% group_by(focal, comp) %>% summarise(r = mean(ext), .groups = "drop") %>%
  mutate(sp = factor(ifelse(focal == 1, "Competing\nspecies", "Non-competing\nspecies"), c("Non-competing\nspecies", "Competing\nspecies")),
         sc = factor(ifelse(comp == 1, "With competition", "Without competition"), c("Without competition", "With competition")))
e1c <- ggplot(t22, aes(sp, r, fill = sc)) + geom_col(position = position_dodge(.75), width = .65) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * r)), position = position_dodge(.75), vjust = -.4, size = 2.7) +
  scale_fill_manual(values = c(C0, C1)) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .12))) +
  labs(x = NULL, y = "Proportion of species extinct")
fo <- d %>% filter(comp == 1, focal == 1)
k1 <- fo %>% group_by(network_id) %>% summarise(k = sum(ext))
pr <- fo %>% filter(network_id %in% k1$network_id[k1$k == 1]) %>% mutate(fate = factor(ifelse(ext == 1, "Excluded", "Survivor"), c("Survivor", "Excluded")))
vio <- function(df, v, lab) ggplot(df, aes(fate, .data[[v]], fill = fate)) + geom_violin(alpha = .4, colour = NA) +
  geom_boxplot(width = .15, outlier.size = .3) + scale_fill_manual(values = c(C0, C1), guide = "none") + labs(x = NULL, y = lab)
e1d <- vio(pr, "rep_rate", "Replication rate"); e1e <- vio(pr, "outgoing_resources", "Total by-product release")
w <- raw %>% select(network_id, species, network_competition, abundance) %>%
  pivot_wider(names_from = network_competition, values_from = abundance, names_prefix = "a") %>%
  inner_join(b %>% select(network_id, species, focal, incoming_degree), by = c("network_id", "species"))
nf <- w %>% filter(focal == 0, a0 > 0, incoming_degree > 0) %>% mutate(vic = a1 == 0)
e1f <- nf %>% mutate(g = factor(pmin(incoming_degree, 4), 1:4, c("1", "2", "3", "4 or more"))) %>% group_by(g) %>%
  summarise(p = mean(vic), n = n(), v = sum(vic), .groups = "drop") %>% ggplot(aes(g, p)) + geom_col(fill = C1, width = .6) +
  geom_text(aes(label = paste0(v, " / ", n)), vjust = -.4, size = 2.6) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .12))) +
  labs(x = "Number of producers", y = "Proportion of secondary extinctions")
F1 <- (e1a | e1b | e1c) / (e1d | e1e | e1f) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIGURAS, "Figure1.jpg"), F1, width = 11, height = 7.4, dpi = 300)

## Suplementaria S1 (ecología)
fate <- d %>% filter(focal == 1) %>% group_by(comp, network_id) %>% summarise(k = sum(ext), .groups = "drop") %>% count(comp, k) %>%
  mutate(sc = factor(ifelse(comp == 1, "With competition", "Without competition"), c("Without competition", "With competition")),
         lab = factor(k, 0:2, c("Both persist", "One extinct", "Both extinct")))
s1a <- ggplot(fate, aes(lab, n, fill = sc)) + geom_col(position = position_dodge(.75), width = .65) +
  geom_text(aes(label = n), position = position_dodge(.75), vjust = -.4, size = 2.7) + scale_fill_manual(values = c(C0, C1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) + labs(x = "Fate of the competing pair", y = "Number of networks")
ab0 <- d %>% filter(comp == 0) %>% select(network_id, species, ab0 = abundance)
win <- pr %>% filter(ext == 0) %>% inner_join(ab0, by = c("network_id", "species")) %>% filter(ab0 > 0)
s1b <- ggplot(win, aes(ab0, abundance)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = GR) +
  geom_point(alpha = .3, size = .7, colour = C0) + scale_x_log10(labels = LN) + scale_y_log10(labels = LN) +
  labs(x = "Abundance without competition", y = "Abundance with competition")
# Asimetría en capacidad competitiva de cada par (calculada en analisis_ecologia.R)
co <- read.csv(file.path(RESULTADOS, "asimetria_competidoras.csv")) %>%
  mutate(g = cut(asim, c(-Inf, .5, 1, 2, 3, 5, Inf), labels = c("< 0.5", "0.5 to 1", "1 to 2", "2 to 3", "3 to 5", "> 5"))) %>%
  group_by(g) %>% summarise(p = mean(n_ext == 0), n = n(), .groups = "drop")
s1c <- ggplot(co, aes(g, p)) + geom_col(fill = C0, width = .6) + geom_text(aes(label = paste0("n = ", n)), vjust = -.4, size = 2.6) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .12))) +
  labs(x = "Asymmetry in competitive ability", y = "Proportion of pairs that coexist")
com <- d %>% group_by(network_id, comp) %>% summarise(R = sum(abundance > 0), B = sum(abundance), .groups = "drop") %>%
  mutate(sc = factor(ifelse(comp == 1, "With competition", "Without competition"), c("Without competition", "With competition")))
s1d <- com %>% count(sc, R) %>% group_by(sc) %>% mutate(p = n / sum(n)) %>% ungroup() %>% filter(R >= 5) %>%
  ggplot(aes(factor(R), p, fill = sc)) + geom_col(position = position_dodge(.8), width = .75) + scale_fill_manual(values = c(C0, C1)) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) + labs(x = "Number of surviving species", y = "Proportion of networks")
s1e <- ggplot(com, aes(sc, B, fill = sc)) + geom_hline(yintercept = 12100, linetype = 2, colour = GR) + geom_violin(alpha = .4, colour = NA) +
  geom_boxplot(width = .13, outlier.size = .3) + scale_fill_manual(values = c(C0, C1), guide = "none") + labs(x = NULL, y = "Total community abundance")
S1 <- (s1a | s1b | s1c) / (s1d | s1e | plot_spacer()) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "top")
ggsave(file.path(FIGURAS, "FigureS1.jpg"), S1, width = 11, height = 7.4, dpi = 300)

# =============================================================================
# Figuras 2-4: comunidades evolutivas
# =============================================================================
U <- as.matrix(read.csv(UNIVERSO, row.names = 1, check.names = FALSE)); SP <- rownames(U); PROD <- rowSums(U)
d <- read.csv(DATOS_EVO, stringsAsFactors = FALSE) %>%
  mutate(sp = sub("+", "_", species, fixed = TRUE), rep = 30 / mean_gest_time, present = mean_abundance > 0, stable = extinction_risk == 0)
sumin <- function(S, w, x = S) as.numeric(t(U[S, x, drop = FALSE]) %*% w)
S9 <- d %>% filter(evo_id == 0, present)
chem <- data.frame(eco_id = as.integer(names(QUIMIOSTATO)), chem = unname(QUIMIOSTATO))
anc_sp <- S9 %>% select(eco_id, sp)
CX <- c(NOT = 1, NAND = 1, AND = 2, ORN = 2, OR = 3, ANDN = 3, NOR = 4, XOR = 4, EQU = 5)
fun <- do.call(rbind, strsplit(SP, "_")); spinfo <- data.frame(sp = SP, cx = CX[fun[, 1]] + CX[fun[, 2]])
ORIG <- c("Chemostat inflow" = C_CHEM, "Resident" = C_ANC, "Novel" = C_NEW)

## ---- Figura 2: riqueza, redes y abundancia total
Me <- d %>% filter(evo_id > 0, present) %>% count(eco_id, evo_id, name = "S")
f2a <- ggplot(Me, aes(factor(eco_id), S)) + geom_hline(yintercept = 9, linetype = 2, colour = C_ANC) +
  geom_boxplot(fill = C_EVO, alpha = .25, outlier.shape = NA, width = .6) + geom_jitter(width = .12, size = 1.2, alpha = .7, colour = C_EVO) +
  labs(x = "Ancestral community", y = "Number of species after evolution")
circ <- function(k, r) { a <- seq(0, 2 * pi, length.out = k + 1)[-1]; cbind(r * cos(a), r * sin(a)) }
dibuja <- function(e, k) {
  g <- d %>% filter(eco_id == e, evo_id == k, present) %>% left_join(anc_sp %>% filter(eco_id == e) %>% mutate(anc = TRUE), by = c("eco_id", "sp"))
  ch <- chem$chem[chem$eco_id == e]; ea <- !is.na(g$anc); xy <- matrix(0, nrow(g), 2)
  xy[ea, ] <- circ(sum(ea), .45); if (any(!ea)) xy[!ea, ] <- circ(sum(!ea), 1)
  nodes <- data.frame(x = xy[, 1], y = xy[, 2], ab = g$mean_abundance,
    o = factor(ifelse(g$sp == ch, "Chemostat inflow", ifelse(ea, "Resident", "Novel")), names(ORIG)))
  B <- which(U[g$sp, g$sp] > 0, arr.ind = TRUE)
  ed <- data.frame(x = xy[B[, 1], 1], y = xy[B[, 1], 2], xend = xy[B[, 2], 1], yend = xy[B[, 2], 2])
  ggplot() + geom_segment(data = ed, aes(x, y, xend = xend, yend = yend), colour = "grey40", alpha = .12, linewidth = .3) +
    geom_point(data = nodes, aes(x, y, colour = o, size = sqrt(ab))) + scale_colour_manual(values = ORIG, drop = FALSE) +
    scale_size(range = c(1.3, 8), guide = "none") + coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) + theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(), plot.tag = element_text(face = "bold", size = 13)) }
M <- d %>% filter(present) %>% group_by(eco_id, evo_id) %>% group_modify(~ {
  S <- .x$sp; w <- .x$mean_abundance; data.frame(cap = sum(rowSums(U[S, S, drop = FALSE]) * w) / sum(PROD[S] * w), tot = sum(w)) }) %>%
  ungroup() %>% mutate(fase = factor(ifelse(evo_id == 0, "Ancestral", "Evolved"), c("Ancestral", "Evolved")))
f2d <- ggplot(M, aes(cap, tot, colour = fase, size = fase)) + geom_point(alpha = .75) +
  scale_colour_manual(values = c(Ancestral = C_ANC, Evolved = C_EVO)) + scale_size_manual(values = c(Ancestral = 2.4, Evolved = 1.4), guide = "none") +
  labs(x = "Fraction of released by-products with a consumer present", y = "Total community abundance") + theme(legend.position = "bottom")
F2 <- (f2a + dibuja(3, 0) + dibuja(3, 2) + f2d + plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")) & theme(legend.position = "bottom")
ggsave(file.path(FIGURAS, "Figure2.jpg"), F2, width = 9.5, height = 8.8, dpi = 300)

## ---- Figura 3: paralelismo y emergencia
sets <- d %>% filter(evo_id > 0, present) %>% anti_join(anc_sp, by = c("eco_id", "sp")) %>%
  group_by(eco_id, evo_id) %>% summarise(s = list(sp), .groups = "drop") %>% arrange(eco_id, evo_id)
n <- nrow(sets); J <- matrix(NA, n, n)
for (i in 1:n) for (j in 1:n) { a <- sets$s[[i]]; bq <- sets$s[[j]]; J[i, j] <- length(intersect(a, bq)) / length(union(a, bq)) }
f3a <- expand.grid(i = 1:n, j = 1:n) %>% mutate(J = J[cbind(i, j)]) %>% ggplot(aes(i, j, fill = J)) + geom_raster() +
  scale_fill_gradient(low = "white", high = C_EVO, name = "Similarity") +
  geom_vline(xintercept = seq(10.5, 90.5, 10), colour = "grey60", linewidth = .3) + geom_hline(yintercept = seq(10.5, 90.5, 10), colour = "grey60", linewidth = .3) +
  coord_equal() + scale_x_continuous(expand = c(0, 0), breaks = seq(5.5, 95.5, 10), labels = 1:10) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(5.5, 95.5, 10), labels = 1:10) +
  labs(x = "Ancestral community", y = "Ancestral community") + theme(legend.position = "right", legend.title = element_text(size = 9))
occ <- d %>% filter(evo_id > 0, present) %>% count(sp, name = "runs") %>% right_join(spinfo, by = "sp") %>% mutate(runs = coalesce(runs, 0L))
f3b <- ggplot(occ, aes(cx, runs)) + geom_jitter(width = .12, height = 0, size = 1.6, alpha = .7, colour = C_EVO) +
  stat_summary(fun = mean, geom = "line", colour = "grey20", linewidth = .8) + labs(x = "Boolean complexity", y = "Experiments with the species")
cand <- d %>% filter(evo_id > 0) %>% distinct(eco_id, evo_id) %>% group_by(eco_id, evo_id) %>% group_modify(~ {
  a9 <- S9 %>% filter(eco_id == .y$eco_id); xx <- setdiff(SP, a9$sp); pres <- d %>% filter(eco_id == .y$eco_id, evo_id == .y$evo_id, present)
  data.frame(sp = xx, emerge = as.integer(xx %in% pres$sp), sup0 = sumin(a9$sp, a9$mean_abundance, xx)) }) %>% ungroup() %>% inner_join(spinfo, by = "sp")
f3c <- cand %>% group_by(cx) %>% summarise(p = mean(emerge)) %>% ggplot(aes(factor(cx), p)) + geom_col(fill = C_EVO, width = .65) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .06))) + labs(x = "Boolean complexity", y = "Candidate species that emerged")
f3d <- cand %>% mutate(q = ntile(sup0, 5)) %>% group_by(q) %>% summarise(p = mean(emerge)) %>% ggplot(aes(factor(q), p)) +
  geom_col(fill = C_ANC, width = .65) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .06))) +
  labs(x = "Quintile of potential resource supply", y = "Candidate species that emerged")
F3 <- (f3a | f3b) / (f3c | f3d) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIGURAS, "Figure3.jpg"), F3, width = 10, height = 8.4, dpi = 300)

## ---- Figura 4: propiedades y tasa de replicación
ev <- d %>% filter(evo_id > 0, present) %>% inner_join(chem, by = "eco_id") %>% left_join(anc_sp %>% mutate(anc = TRUE), by = c("eco_id", "sp")) %>%
  mutate(o = factor(case_when(sp == chem ~ "Chemostat inflow", !is.na(anc) ~ "Resident", TRUE ~ "Novel"), names(ORIG)))
f4a <- ggplot(ev, aes(o, mean_abundance, fill = o)) + geom_violin(alpha = .45, colour = NA) + geom_boxplot(width = .14, outlier.size = .3) +
  scale_y_log10(labels = LN) + scale_fill_manual(values = ORIG, guide = "none") + labs(x = NULL, y = "Abundance at steady state")
f4b <- ggplot(ev %>% filter(!is.na(rep)), aes(o, rep, fill = o)) + geom_violin(alpha = .45, colour = NA) + geom_boxplot(width = .14, outlier.size = .3) +
  scale_fill_manual(values = ORIG, guide = "none") + labs(x = NULL, y = "Replication rate")
se <- ev %>% group_by(eco_id, evo_id) %>% mutate(sup = sumin(sp, mean_abundance)) %>% ungroup() %>% filter(o != "Chemostat inflow", sup > 0, !is.na(rep))
f4c <- ggplot(se, aes(sup, mean_abundance, colour = o)) + geom_point(alpha = .3, size = .7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = .7) + scale_x_log10(labels = LN) + scale_y_log10(labels = LN) +
  scale_colour_manual(values = ORIG) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Resource supply received from the species present", y = "Abundance at steady state")
GRP <- c("Chemostat inflow" = C_CHEM, "Direct consumers" = C_ANC, "Other residents" = "#8FA8D0")
pr <- S9 %>% select(eco_id, sp, rep_anc = rep) %>% inner_join(ev %>% select(eco_id, evo_id, sp, rep_evo = rep, o, chem), by = c("eco_id", "sp")) %>%
  filter(!is.na(rep_anc), !is.na(rep_evo)) %>%
  mutate(g = factor(ifelse(o == "Chemostat inflow", "Chemostat inflow", ifelse(U[cbind(chem, sp)] > 0, "Direct consumers", "Other residents")), names(GRP)),
         ratio = rep_evo / rep_anc)
f4d <- ggplot(pr, aes(g, ratio, fill = g)) + geom_hline(yintercept = 1, linetype = 2, colour = GR) + geom_violin(alpha = .45, colour = NA) +
  geom_boxplot(width = .12, outlier.size = .3) + scale_y_log10(labels = LN) + scale_fill_manual(values = GRP, guide = "none") +
  labs(x = NULL, y = "Evolved / ancestral replication rate")
F4 <- (f4a | f4b) / (f4c | f4d) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIGURAS, "Figure4.jpg"), F4, width = 9.5, height = 7.6, dpi = 300)

# =============================================================================
# Figura 5 y Figura S2: topología de las redes (calculada en analisis_evolucion.R)
# =============================================================================
TOPO <- read.csv(file.path(RESULTADOS, "topologia_redes.csv")) %>%
  mutate(fase = factor(ifelse(evo_id == 0, "Ancestral", "Evolved"), c("Ancestral", "Evolved")))
EN <- c(C = "Connectance", LS = "Links per species", vconn = "Vertex connectivity", rec = "Reciprocity",
        cv_in = "In-degree heterogeneity", cv_out = "Out-degree heterogeneity", assort = "Degree assortativity",
        trans = "Transitivity", path = "Mean path length", nodf = "Nestedness", modul = "Modularity",
        scc = "Largest strongly connected component", ffl = "Feed-forward loops", cyc = "Three-species cycles",
        full = "Fully reciprocal triads", capu = "By-product capture", F0 = "Trophic incoherence")
# Cada línea une una comunidad ancestral con la media de sus 10 réplicas; los puntos tenues son las réplicas
panel <- function(vars, ncol, base = 11, strip = 9) {
  L <- TOPO %>% pivot_longer(all_of(vars), names_to = "m") %>% mutate(m = factor(EN[m], EN[vars]))
  Mm <- L %>% group_by(eco_id, fase, m) %>% summarise(value = mean(value), .groups = "drop")
  ggplot(L, aes(fase, value)) +
    geom_jitter(data = L %>% filter(fase == "Evolved"), width = .08, height = 0, colour = C_EVO, alpha = .18, size = .7) +
    geom_line(data = Mm, aes(group = eco_id), colour = "grey55", linewidth = .45) +
    geom_point(data = Mm, aes(colour = fase), size = 1.8) +
    facet_wrap(~ m, scales = "free_y", ncol = ncol, labeller = label_wrap_gen(20)) +
    scale_colour_manual(values = c(Ancestral = C_ANC, Evolved = C_EVO), guide = "none") +
    scale_x_discrete(limits = c("Ancestral", "Evolved")) + labs(x = NULL, y = NULL) + theme_classic(base_size = base) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = strip))
}
## Figura suplementaria S2: las 17 propiedades
ggsave(file.path(FIGURAS, "FigureS2.jpg"), panel(names(EN), ncol = 5), width = 11, height = 9, dpi = 300)

## Figura 5: (A) cinco propiedades y (B) desviaciones frente a los modelos nulos
pA <- panel(c("C", "LS", "capu", "trans", "ffl"), ncol = 5, base = 10, strip = 9)
orden <- c("C", "rec", "cv_in", "cv_out", "assort", "trans", "path", "nodf", "modul", "ffl", "cyc", "full", "capu", "F0")
columnas <- c("Ancestral vs null A" = "Ancestral\nvs null A", "Evolved vs null A" = "Evolved\nvs null A",
              "Evolved vs null B" = "Evolved\nvs null B", "Evolved vs null C" = "Evolved\nvs null C")
zl <- read.csv(file.path(RESULTADOS, "modelos_nulos_resumen.csv")) %>%
  mutate(metric = factor(EN[metrica], rev(EN[orden])), col = factor(columnas[comparacion], columnas),
         sig = p_ajustada < 0.05)
pB <- ggplot(zl, aes(col, metric, fill = z_media)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f%s", z_media, ifelse(sig, "*", ""))), size = 3.2) +
  scale_fill_gradient2(low = C_ANC, mid = "white", high = C_EVO, midpoint = 0, name = "Mean\nstandardized\ndeviation") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), legend.position = "right")
f5 <- (pA / (plot_spacer() + pB + plot_spacer() + plot_layout(widths = c(.08, 1, .08)))) +
  plot_layout(heights = c(.55, 1)) + plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(FIGURAS, "Figure5.jpg"), f5, width = 10, height = 9.2, dpi = 300)

cat("Figuras guardadas en", FIGURAS, "\n")
