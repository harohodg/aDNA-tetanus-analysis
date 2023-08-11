# Benjamin Jean-Marie Tremblay

d <- readr::read_tsv("results/gubbins/core-snippy-genome.per_branch_statistics.csv")
tree <- ape::read.tree("results/gubbins/core-snippy-genome.final_tree.tre")

library(ape)
supp <- readr::read_tsv("data/SuppTable2.txt")
tree <- read.nexus("data/tetani-tree")
tree2 <- read.nexus("data/tetani-fasttree-oct21")
ord <- data.frame(row.names = NULL, Sample = NA_character_, Name = tree2$tip.label)
ss <- data.frame(row.names = NULL, id = supp[[1]], name = supp[[4]])
ss$name[14] <- "'Yámana-Tooth'"
ss$name[15] <- "'Vác-Mummy-Tissue'"
ord$Sample <- structure(ss$id, names = ss$name)[ord$Name]
ord$Final <- ord$Sample
ord$Final[is.na(ord$Sample)] <- paste0("Ref", 1:sum(is.na(ord$Sample)))

d2 <- d[d$Node %in% gsub("'", "", tree$tip.label, fixed = TRUE), ]
d2$acMAG <- grepl("^SAM", unname(structure(ord$Final, names = ord$Name)[d2$Node]))

d3 <- d2[, c("Node", "r/m", "rho/theta", "acMAG")]
d3 <- tidyr::gather(d3, "metric", "value", 2:3)

pval1 <- with(d3[d3$metric == "r/m", ], wilcox.test(value[acMAG], value[!acMAG])$p.value)
pval2 <- with(d3[d3$metric == "rho/theta", ], wilcox.test(value[acMAG], value[!acMAG])$p.value)

with(d3[d3$metric == "r/m", ], length(value[!acMAG]))

p1 <- ggplot(d3[d3$metric == "r/m", ], aes(x = acMAG, y = value, fill = acMAG)) +
  geom_violin(trim = TRUE) +
  scale_x_discrete(name = element_blank(), labels = c("Refs", "acMAGs")) +
  scale_fill_manual(values = palette.colors(4, "R4")[3:4]) +
  ylab("r/m") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid = element_blank()
  )
p2 <- ggplot(d3[d3$metric == "rho/theta", ], aes(x = acMAG, y = value, fill = acMAG)) +
  geom_violin(trim = TRUE) +
  scale_x_discrete(name = element_blank(), labels = c("Refs", "acMAGs")) +
  scale_fill_manual(values = palette.colors(4, "R4")[3:4]) +
  ylab("rho/theta") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid = element_blank()
  )

cowplot::save_plot(
  "figures/recombination_metrics.pdf",
  cowplot::plot_grid(p1, p2, ncol = 2),
  base_height = 2, base_width = 4.2)

