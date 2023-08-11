# Benjamin Jean-Marie Tremblay
# 2021-07-09

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(viridis)
library(cowplot)

sf_use_s2(FALSE)

world <- ne_countries(scale = "medium", returnclass = "sf")

metaNew <- readxl::read_xlsx("data/METADATA-UPDATED.xlsx", skip = 3)

metadata <- metaNew[, c("BioSample ID", "Adjusted-Time", "Date-range-1", "Date-range-2",
  "lat", "lon", "clade", "acMAG damage level")]

dat <- metadata %>%
  group_by_all() %>%
  mutate(count = n()) %>%
  distinct()

colnames(dat) <- c("BioSample", "Midpoint", "Lowpoint", "Highpoint", "Latitude", "Longitude",
  "Clade", "Damage", "count")

dat$Midpoint <- as.numeric(dat$Midpoint)
dat$Lowpoint <- as.numeric(dat$Lowpoint)
dat$Highpoint <- as.numeric(dat$Highpoint)

bioOrder <- dat$Lowpoint
bioOrder[is.na(bioOrder)] <- dat$Midpoint[is.na(bioOrder)]
dat$BioSample <- factor(dat$BioSample, levels = dat$BioSample[order(bioOrder, decreasing = TRUE)])

dat$Midpoint[!is.na(dat$Lowpoint)] <- NA

tosort <- dat$Lowpoint
tosort[is.na(tosort)] <- dat$Midpoint[is.na(tosort)]
dat <- dat[order(tosort), ]

dat2 <- tidyr::gather(dat, "X", "Date", 3:4)

dat2$Date2 <- dat2$Date
dat2$Date2[dat2$X == "Lowpoint"] <- dat2$Date2[dat2$X == "Lowpoint"] + 10
dat2$Date2[dat2$X == "Highpoint"] <- dat2$Date2[dat2$X == "Highpoint"] - 10

dat3 <- dat[order(dat$count), ]

dat4 <- dat3
dat4$Clade[dat4$BioSample=="SAMEA5847426"] <- "1"
dat4$Clade[dat4$BioSample=="SAMEA5847432"] <- "1"
dat4$Clade[dat4$BioSample=="SAMEA5847501"] <- "1"
dat4$Clade[dat4$BioSample=="SAMEA6490841"] <- "1"

ggplot(data = world) +
  geom_sf(fill = "grey99", lwd = 0.2) +
  coord_sf(expand = FALSE) +
  geom_point(aes(Longitude, Latitude, fill = Clade),
    dat4, pch = 21, size = 0.15, alpha = 0.9, colour = "red", fill = "red") +
  ggrepel::geom_label_repel(aes(Longitude, Latitude, label = Clade, fill = Clade),
    size = 1.5, data = dat4, max.overlaps = 1000, segment.size = unit(0.2, "mm"),
    box.padding = 0.1, label.r = 0, label.size = 0, fontface = "bold",
    min.segment.length = 0) +
  scale_fill_manual(values = awtools::ppalette, name = "Clade") +
  theme_bw() +
  guides(fill = guide_legend(keyheight = 0.7)) +
  theme(
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(colour = grey(0.5), linetype = "dashed", size = 0.2),
    panel.background = element_rect(fill = "azure"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )

ggsave("figures/map.pdf", plot = p1v4,  width = 9.99, height = 5.28)

