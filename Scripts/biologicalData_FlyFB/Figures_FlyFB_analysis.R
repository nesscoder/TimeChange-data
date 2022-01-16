library(devtools)
library(tidyverse)
library(matrixStats)
library(readr)
library(MASS)
library(circacompare)
library(clusterProfiler)
library(ggnewscale)
# devtools::install_github("YuLab-SMU/enrichplot")
library(enrichplot)
library(GOSemSim)
library(DOSE)
library(org.Dm.eg.db)
library(RcisTarget)
library(DT)
library(ggpubr)
library(viridis)
library(ggupset)
library(RColorBrewer)
library(europepmc)
library(fasano.franceschini.test)
load_all()

load("~/Desktop/TimeChange-data/Results/biologicalData_FlyFB/Results_HC18vsHC25vsJRK_JTK_TimeChange.RData")

# set seed for reproducibility
set.seed(2022)

## Hist of Gene Numbers ----------------------------------------------------
counts <- data.frame(
  cyclingCat = c("No Change", "Differential Cycling", "Cycling at 18°C only", "Cycling at 25°C only"),
  num = c(length(common_nondiffcyclers), length(common_cyclers), length(gainofCycling_25to18), length(lossofCycling_25to18))
)
counts

p1 <- ggplot(data = counts, aes(x = cyclingCat, y = num, fill = cyclingCat, color = cyclingCat, shape = cyclingCat)) +
  geom_bar(stat = "identity") +
  geom_point(aes(y = num + 5), size = 5) +
  labs(title = "Cyling Dynamics from 25°C to 18°C", x = "", y = "Count") +
  scale_x_discrete(limits = c("No Change", "Differential Cycling", "Cycling at 25°C only", "Cycling at 18°C only")) +
  scale_fill_manual("legend", values = c("No Change" = "#868686FF", "Differential Cycling" = "#EFC000FF", "Cycling at 25°C only" = "#CD534CFF", "Cycling at 18°C only" = "#0073C2FF")) +
  scale_colour_manual("legend", values = c("No Change" = "#868686FF", "Differential Cycling" = "#EFC000FF", "Cycling at 25°C only" = "#CD534CFF", "Cycling at 18°C only" = "#0073C2FF")) +
  scale_shape_manual("legend", values = c("No Change" = 24, "Differential Cycling" = 21, "Cycling at 25°C only" = 22, "Cycling at 18°C only" = 23)) +
  coord_flip() +
  theme_linedraw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(size = 3, fill = NA),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )

## Differential Cyclers
common_diffCyclers <- unique(c(gainofCycling_25to18, lossofCycling_25to18, common_cyclers))
group <- c(
  rep("Cyc 18", length(gainofCycling_25to18)),
  rep("Cyc 25", length(lossofCycling_25to18)),
  rep("Diff Cyc", length(common_cyclers)),
  rep("non-Diff Cyc", length(common_nondiffcyclers))
)


rownames(results_HC18vsHC25) <- results_HC18vsHC25$sample1_ID
rownames(HC25vs18_JTK_results) <- HC25vs18_JTK_results$annot
rownames(HC18vs25_JTK_results) <- HC18vs25_JTK_results$annot

diffCyclersPlot <- data.frame(
  x = -log(c(
    HC25vs18_JTK_results[common_diffCyclers, "BH.Q"],
    HC25vs18_JTK_nonDiff_results[common_nondiffcyclers, "BH.Q"]
  )),
  y = -log(c(
    HC18vs25_JTK_results[common_diffCyclers, "BH.Q"],
    HC18vs25_JTK_nonDiff_results[common_nondiffcyclers, "BH.Q"]
  )),
  z = -log(unlist(results_HC18vsHC25[c(common_diffCyclers, common_nondiffcyclers), "FDR"])),
  amp18 = log(c(
    HC18vs25_JTK_results[common_diffCyclers, "AMP"],
    HC18vs25_JTK_nonDiff_results[common_nondiffcyclers, "AMP"]
  )),
  amp25 = log(c(
    HC25vs18_JTK_results[common_diffCyclers, "AMP"],
    HC25vs18_JTK_nonDiff_results[common_nondiffcyclers, "AMP"]
  )),
  phase18 = c(
    HC18vs25_JTK_results[common_diffCyclers, "LAG"],
    HC18vs25_JTK_nonDiff_results[common_nondiffcyclers, "LAG"]
  ),
  phase25 = c(
    HC25vs18_JTK_results[common_diffCyclers, "LAG"],
    HC25vs18_JTK_nonDiff_results[common_nondiffcyclers, "LAG"]
  ),
  group = group
)
rownames(diffCyclersPlot) <- c(common_diffCyclers, common_nondiffcyclers)

p2 <- ggplot(diffCyclersPlot, aes(x = z, fill = group, shape = group)) +
  geom_histogram(breaks=seq(0,17,.23)) +
  ggtitle("TimeChange FDR < 0.1") +
  geom_vline(xintercept = -log(.1), color = "grey", lwd = 1.5, lty = 2) +
  xlab("-log(p.value)") +
  scale_fill_manual("legend", values = c("Diff Cyc" = "#EFC000FF", "Cyc 25" = "#CD534CFF", "Cyc 18" = "#0073C2FF", "non-Diff Cyc" = "#868686FF")) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(size = 3, fill = NA),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )

p3 <- ggplot(diffCyclersPlot, aes(x = x, y = y, color = group, shape = group)) +
  geom_hline(yintercept = -log(.1), color = "grey", lwd = 1.5, lty = 2) +
  geom_vline(xintercept = -log(.1), color = "grey", lwd = 1.5, lty = 2) +
  geom_point(size = 4) +
  ggtitle("JTK_CYCLE FDR < 0.1") +
  xlab("-log(p.value) 25°C") +
  ylab("-log(p.value) 18°C") +
  scale_color_manual("legend", values = c("Diff Cyc" = "#EFC000FF", "Cyc 25" = "#CD534CFF", "Cyc 18" = "#0073C2FF", "non-Diff Cyc" = "#868686FF")) +
  scale_shape_manual("legend", values = c("Diff Cyc" = 16, "Cyc 25" = 15, "Cyc 18" = 18, "non-Diff Cyc" = 17)) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(size = 3, fill = NA),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )

p4 <- ggplot(diffCyclersPlot, aes(x = amp25, y = amp18, color = group, shape = group)) +
  geom_abline(slope = 1, color = "grey", lwd = 1.5, lty = 2) +
  geom_point(size = 4) +
  scale_color_manual("legend", values = c("Diff Cyc" = "#EFC000FF", "non-Diff Cyc" = "#868686FF")) +
  scale_shape_manual("legend", values = c("Diff Cyc" = 16, "non-Diff Cyc" = 17)) +
  ggtitle("Differential vs Non-Differentially Cyclng Genes LogFC AMP") +
  xlab("25°C LogFC AMP") +
  ylab("18°C LogFC AMP") +
  labs(color = "logFC") +
  guides(shape = "none") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(size = 3, fill = "NA"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )

p5 <- ggplot(diffCyclersPlot, aes(x = phase25, y = phase18, color = group, shape = group)) +
  geom_abline(slope = 1, color = "grey", lwd = 1.5, lty = 2) +
  geom_point(size = 4) +
  scale_color_manual("legend", values = c("Diff Cyc" = "#EFC000FF", "non-Diff Cyc" = "#868686FF")) +
  scale_shape_manual("legend", values = c("Diff Cyc" = 16, "non-Diff Cyc" = 17)) +
  ggtitle("Differential vs Non-Differentially Cyclng Genes Phases") +
  xlab("25°C Log Phase (h)") +
  ylab("18°C Log Phase (h)") +
  labs(color = "logFC") +
  guides(shape = "none") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(size = 3, fill = "NA"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )

blank <- ggplot() +
  theme_void()

png(filename = "~/Desktop/TimeChange-data/Results/biologicalData_FlyFB/figs/geneDiffCyc.png",width = 16, height = 8, units = "in", res = 300)
ggarrange(p1, p2, p3, blank, p5, p4, ncol = 3, nrow = 2, labels = c("A", "B", "C", "", "D", "E"))
dev.off()

# Plot Results ----------------------------------------------

# Figure 2A-E----------------------------------------------
fbgn_annotation <- read_table2("~/Desktop/TimeChange-data/Data/biologicalData_FlyFB/fbgn_annotation_ID_fb_2021_06.tsv")
plotExpression <- function(data1, data2, gene, replicates, title, x2 = 1:12) {
  ylim <- range(unlist(unname(data1[gene, ])), unlist(unname(data2[gene, ])))

  sel <- which(fbgn_annotation$primary_FBgn == gene)

  plot(unlist(unname(data1[gene, seq(1, ncol(data1), replicates)])),
       main = fbgn_annotation$gene_symbol[sel],
       cex = 0,
       xlab = "",
       ylab = "",
       ylim = ylim,
       axes = F
  )
  title(xlab = "ZT Time", line = 2)
  title(ylab = "Expression", line = 2)
  mtext(title, side = 2, line = 3, adj = 0.5, font = 2)
  for (i in 1:replicates) {
    points(unlist(unname(data1[gene, seq(i, ncol(data1), replicates)])), type = "o", pch = 15 + i, col = "#0073C2FF")
    points(x = x2, y = unlist(unname(data2[gene, seq(i, ncol(data2), replicates)])), type = "o", pch = 15 + i, col = "#CD534CFF")
  }
  axis(2)
  axis(1, at = seq(1, 12, 1), labels = seq(2, 24, 2))
  box(lwd = 2)
}


# Figure 2F ----------------------------------------------

png(filename = "~/Desktop/TimeChange-data/Results/biologicalData_FlyFB/figs/geneLevelDiffCyc.png",width = 20, height = 12, units = "in", res = 300)

nGenesPlot <- 5
topN_Cyc18 <- diffCyclersPlot %>%
  filter(group == "Cyc 18") %>%
  rownames_to_column(var = "genes") %>%
  group_by(group) %>%
  top_n(n = nGenesPlot, wt = z)

topN_Cyc25 <- diffCyclersPlot %>%
  filter(group == "Cyc 25") %>%
  rownames_to_column(var = "genes") %>%
  group_by(group) %>%
  top_n(n = nGenesPlot, wt = z)

topN_DiffCyc <- diffCyclersPlot %>%
  filter(group == "Diff Cyc") %>%
  rownames_to_column(var = "genes") %>%
  group_by(group) %>%
  top_n(n = nGenesPlot, wt = z)

topN_nonDiffCyc <- diffCyclersPlot %>%
  filter(group == "non-Diff Cyc") %>%
  rownames_to_column(var = "genes") %>%
  group_by(group) %>%
  arrange(desc(z), .by_group = TRUE) %>%
  top_n(n = -nGenesPlot, wt = z)

par(mfrow = c(4, nGenesPlot))
for (idx in 1:nGenesPlot) {
  if (idx == 1) {
    title <- "Cycling 18°C only"
  } else {
    title <- ""
  }

  plotExpression(
    data1 = HC18_data,
    data2 = HC25_data,
    gene = unlist(topN_Cyc18)[idx],
    replicates = 2,
    title = title
  )
}

for (idx in 1:nGenesPlot) {
  if (idx == 1) {
    title <- "Cycling 25°C only"
  } else {
    title <- ""
  }
  plotExpression(
    data1 = HC18_data,
    data2 = HC25_data,
    gene = unlist(topN_Cyc25)[idx],
    replicates = 2,
    title = title
  )
}

for (idx in 1:nGenesPlot) {
  if (idx == 1) {
    title <- "Differential Dynamics"
  } else {
    title <- ""
  }
  plotExpression(
    data1 = HC18_data,
    data2 = HC25_data,
    gene = unlist(topN_DiffCyc)[idx],
    replicates = 2,
    title = title
  )
}

for (idx in 1:nGenesPlot) {
  if (idx == 1) {
    title <- "Non-Differential Dynamics"
  } else {
    title <- ""
  }
  plotExpression(
    data1 = HC18_data,
    data2 = HC25_data,
    gene = unlist(topN_nonDiffCyc)[idx],
    replicates = 2,
    title = title
  )
}
dev.off()

# Figure 3 ----------------------------------------------

rankGeneList <- -log(as.numeric(results_HC18vsHC25$p.value))
names(rankGeneList) <- results_HC18vsHC25$sample1_ID
rankGeneList <- sort(rankGeneList, decreasing = T)

diffDynamicGSEA <- gseGO(
  geneList = rankGeneList,
  OrgDb = org.Dm.eg.db,
  keyType = "FLYBASE",
  scoreType = "pos",
  eps = 0,
  ont = "BP",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

goAnalysisPlotPanel1 <- dotplot(object = diffDynamicGSEA,
                                showCategory = 25,
                                label_format = 100,
                                orderBy = "GeneRatio") +
  ggtitle("GSEA For Go Terms") +
  theme_linedraw() +
  theme(
    panel.border = element_rect(size = 3, fill = "NA"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  )+
  guides(color=guide_colorbar(title = "p.value", reverse = T))

diffDynamicGSEA@result$GeneRatio <- fortify(diffDynamicGSEA, showCategory = nrow(diffDynamicGSEA@result))$GeneRatio
goAnalysisPlotPanel2 <- ridgeplot(diffDynamicGSEA,
                                  showCategory = 25,
                                  label_format = 100,
                                  orderBy = "GeneRatio") +
  xlab("-log(p)") +
  ggtitle("-log(p) Distribution of GSEA Genes") +
  theme_linedraw() +
  theme(
    panel.border = element_rect(size = 3, fill = "NA"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 2),
    axis.text = element_text(face = "bold")
  ) +
  guides(fill=guide_colorbar(title = "p.value", reverse = T))


png(filename = "~/Desktop/TimeChange-data/Results/biologicalData_FlyFB/figs/GoAnalysis.png",width = 16, height = 14, units = "in", res = 300)
ggarrange(goAnalysisPlotPanel1,goAnalysisPlotPanel2, ncol = 2, labels = c("A", "B"))
dev.off()
