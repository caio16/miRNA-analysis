library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(limma)
library(pheatmap)
library(purrr)
library(tibble)
library(knitr)
library(ggplot2)
library(enrichR)


setwd("caminho")

raw <- read_excel("caminho", col_names = FALSE)
mirtarbase <- read_excel("caminho")

sample_cols <- 3:19

groups <- factor(unlist(raw[2, sample_cols]))

expr <- raw[-c(1,2), sample_cols]
expr <- as.data.frame(expr)
expr[] <- lapply(expr, as.numeric)

genes <- unlist(raw[-c(1,2), 1])
rownames(expr) <- genes

get_anova <- function(x, grp) {
  m  <- aov(x ~ grp)
  tb <- summary(m)[[1]]
  c(F.value = tb[1, "F value"],
    p.value = tb[1, "Pr(>F)"])
}

res_mat <- t(apply(expr, 1, get_anova, grp = groups))

res_df <- data.frame(
  Gene = rownames(res_mat),
  F.value = res_mat[, "F.value"],
  p.value = res_mat[, "p.value"],
  stringsAsFactors = FALSE
)
res_df$p.adj <- p.adjust(res_df$p.value, method = "BH")

res_df2 <- res_df
rownames(res_df2) <- res_df2$Gene
res_df2$Gene <- NULL
head(res_df2)

sig <- subset(res_df, p.adj < 0.05)
nrow(sig)
sig

# Conferência - Apenas para conferir os dados gerados, não é necessário para o funcionamento do código
# dim(expr)
# length(groups)
# table(groups)
# 
# sum(res_df$p.value < 0.05)
# sum(res_df$p.adj < 0.05)

# hist(res_df$p.value, breaks = 30, main = "Distribuição p-values", xlab = "p-value")

# Limma
# Cria matriz de indicadores sem intercepto, gerando uma coluna por nível
# do fator (groups)
# Se groups tem níveis "E", "O", "ODM", então design será uma matriz n x 3,
# onde cada linha indica a qual grupo a amostra pertence (1 na coluna do seu grupo
# e 0 nas demais)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Para cada gene (linha de expr), ajusta um modelo linear de intensidade ~ grupos
# usando a matriz design
fit <- lmFit(expr, design)

# Cria matriz de contrastes que diz como formar cada comparação de interesse a partir
# das colunas de design
# Defini aqui 3 contrastes: EvsO, EvsODM, ODMvsO
cont.matrix <- makeContrasts(
  EvsO = E - O,
  EvsODM = E - ODM,
  ODMvsO = ODM - O,
  levels = design
)

# Reestima o modelo (fit), mas em vez de testar cada coeficiente bruto,
# calcula diretamente os constrast estimates e suas variâncias para cada gene
fit2 <- contrasts.fit(fit, cont.matrix)

# Ajusta essas variâncias usando o método empirical Bayes, gera estatísticas
# t moderadas, bem como p-values e estatística b (log-odds)
fit2 <- eBayes(fit2)

res_EvsO <- topTable(fit2, coef="EvsO",   number=Inf, adjust.method="BH")
res_EvsODM <- topTable(fit2, coef="EvsODM", number=Inf, adjust.method="BH")
res_ODMvsO <- topTable(fit2, coef="ODMvsO", number=Inf, adjust.method="BH")

# head(res_EvsO, 10)
# head(res_EvsODM, 10)
# head(res_ODMvsO, 10)

# write.csv(res_EvsO,   "limma_E_vs_O.csv",   row.names=TRUE)
# write.csv(res_EvsODM, "limma_E_vs_ODM.csv", row.names=TRUE)
# write.csv(res_ODMvsO, "limma_ODM_vs_O.csv", row.names=TRUE)

sig_EvsO <- subset(res_EvsO,   adj.P.Val < 0.05)
sig_EvsODM <- subset(res_EvsODM, adj.P.Val < 0.05)
sig_ODMvsO <- subset(res_ODMvsO, adj.P.Val < 0.05)

# Heatmap com os 50 mais significativos. Não é necessário para o funcionamento do código
# top50 <- rownames(res_EvsO)[1:50]
# mat <- expr[top50, ]
# mat_scaled <- t( scale( t(mat) ) )
# annotation_col <- data.frame(
#   Group = groups
# )
# rownames(annotation_col) <- colnames(mat_scaled)

# pheatmap(
#   mat_scaled,
#   annotation_col = annotation_col,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   clustering_distance_rows = "euclidean",
#   clustering_distance_cols = "euclidean",
#   clustering_method = "complete",
#   scale = "none"
# )

# Wilcox e T test - Tanto o Wilcox, quanto o T-test não retornaram adj. p-values significativos.
# pairs <- list(
#   E_vs_O = c("E", "O"),
#   E_vs_ODM = c("E", "ODM"),
#   ODM_vs_O = c("ODM", "O")
# )

# Função que seleciona os valores de expressão, remove NAs (caso tenha),
# executa o teste (wilcox ou t) e retorna só a estatística e o p-value
# run_pair_test <- function(x, grp, levels, test_fn) {
#   g1 <- as.numeric(x[grp == levels[1]])
#   g2 <- as.numeric(x[grp == levels[2]])
#   g1 <- g1[!is.na(g1)]
#   g2 <- g2[!is.na(g2)]
#   tst <- test_fn(g1, g2)
#   c(statistic = unname(tst$statistic),
#     p.value = tst$p.value)
# }
# 
# tt_res_list <- lapply(names(pairs), function(p) {
#   lv <- pairs[[p]]
#   mat <- t(sapply(rownames(expr),
#                    function(g) run_pair_test(expr[g, ], groups, lv, t.test)))
#   df <- data.frame(
#     Gene = rownames(mat),
#     statistic = mat[,"statistic"],
#     p.value = mat[,"p.value"],
#     pair = p,
#     stringsAsFactors = FALSE
#   )
#   df$p.adj <- p.adjust(df$p.value, method="BH")
#   df
# })
# 
# ttests <- do.call(rbind, tt_res_list)
# 
# wilc_res_list <- lapply(names(pairs), function(p) {
#   lv   <- pairs[[p]]
#   mat  <- t(sapply(rownames(expr),
#                    function(g) run_pair_test(expr[g, ], groups, lv, wilcox.test)))
#   df <- data.frame(
#     Gene = rownames(mat),
#     statistic = mat[,"statistic"],
#     p.value = mat[,"p.value"],
#     pair = p,
#     stringsAsFactors = FALSE
#   )
#   df$p.adj <- p.adjust(df$p.value, method="BH")
#   df
# })
# 
# wilcoxon <- do.call(rbind, wilc_res_list)
# 
# sig_t_EO <- subset(ttests,  pair=="E_vs_O"   & p.adj < 0.05)
# sig_w_EO <- subset(wilcoxon, pair=="E_vs_O"   & p.adj < 0.05)
# 
# # E vs ODM
# sig_t_EODM <- subset(ttests,  pair=="E_vs_ODM" & p.adj < 0.05)
# sig_w_EODM <- subset(wilcoxon, pair=="E_vs_ODM" & p.adj < 0.05)
# 
# # ODM vs O
# sig_t_OO <- subset(ttests, pair=="ODM_vs_O" & p.adj < 0.05)
# sig_w_OO <- subset(wilcoxon, pair=="ODM_vs_O" & p.adj < 0.05)

# Gene target de cada miRNA
mirnas <- c("hsa-miR-6514-5p",
            "hsa-miR-216a-3p",
            "hsa-miR-126-3p",
            "hsa-miR-99b-5p",
            "hsa-miR-656-5p")

mt <- mirtarbase %>%
  filter(miRNA %in% mirnas) %>%
  select(miRNA, `Target Gene`) %>%
  distinct() %>%
  arrange(miRNA)

write.csv(mt,   "target_genes.csv",   row.names=TRUE)

target_list <- mt %>%
  group_by(miRNA) %>%
  summarize(targets = list(`Target Gene`)) %>%
  deframe() %>%

pairs <- combn(mirnas, 2, simplify = FALSE)
intersect_df <- lapply(pairs, function(pr) {
  shared <- intersect(target_list[[pr[1]]], target_list[[pr[2]]])
  data.frame(
    miRNA1 = pr[1],
    miRNA2 = pr[2],
    SharedCount = length(shared),
    SharedGenes = paste(shared, collapse = ", "),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Número de targets por miRNA
mt_alvos <- mt %>%
  group_by(miRNA) %>%
  summarise(N_Targets = n()) %>%
  arrange(desc(N_Targets)) %>%
  kable(caption = "Número de targets por miRNA")

# Boxplot
sig_genes <- unique(c(
  rownames(sig_EvsO),
  rownames(sig_EvsODM),
  rownames(sig_ODMvsO)
))

expr_sub <- expr[sig_genes, , drop = FALSE]

df_long <- expr_sub %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  mutate(
    Group = factor(groups[match(Sample, colnames(expr))],
                   levels = c("E","O","ODM"))
  )

boxplot <- ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Expressão dos miRNAs significativos por grupo",
    x = "Grupo",
    y = "Expressão"
  )

# Análise de enriquecimento - Necessário maior aprofundamento na etapa de enriquecimento.
# if (!requireNamespace("enrichR", quietly=TRUE)) install.packages("enrichR")
# library(enrichR)
# 
# setEnrichrSite("https://maayanlab.cloud/Enrichr")
# 
# gene_list <- unique(unlist(target_list))
# 
# dbs <- listEnrichrDbs()
# head(dbs)
# 
# my_dbs <- c("miRTarBase_2017",
#             "KEGG_2021_Human",
#             "GO_Biological_Process_2025"
#             )
# 
# enr <- enrichr(gene_list, my_dbs)
# 
# res_all <- bind_rows(
#   enr$miRTarBase_2017 %>% mutate(Database = "GO_BP"),
#   enr$GO_Biological_Process_2025 %>% mutate(Database = "GO_MF"),
#   enr$KEGG_2021_Human %>% mutate(Database = "KEGG"))
# 
# top20 <- res_all %>%
#   arrange(Adjusted.P.value) %>%
#   head(20)
# print(top20)
# 
# top10_go <- enr$GO_Biological_Process_2025 %>%
#   arrange(Adjusted.P.value) %>%
#   head(10)
# 
# ggplot(top10_go, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(
#     title = "Top 10 - GO Biological Process",
#     x = "",
#     y = "Score combinado"
#   ) +
#   theme_minimal()
