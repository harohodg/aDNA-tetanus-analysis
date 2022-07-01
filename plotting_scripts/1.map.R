# Benjamin Jean-Marie Tremblay
# 2021-07-09

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(viridis)
library(cowplot)

world <- ne_countries(scale = "medium", returnclass = "sf")

metadata <- readxl::read_xlsx("data/Metadata.xlsx")
metadata <- metadata[, c("FINAL-MP-DATES", "FINAL-DATE-RANGE", "Latitude", "Longitude", "BioSample")]

dat <- metadata %>%
  group_by_all() %>%
  mutate(count = n()) %>%
  distinct()

dat$Midpoint <- gsub("BC", "", dat$`FINAL-MP-DATES`)
dat$Midpoint <- gsub("AD", "", dat$Midpoint)
dat$Midpoint <- as.numeric(dat$Midpoint)
dat$Midpoint[grepl("BC", dat$`FINAL-MP-DATES`)] <- dat$Midpoint[grepl("BC", dat$`FINAL-MP-DATES`)] * -1

d1 <- vapply(strsplit(dat$`FINAL-DATE-RANGE`, "-", TRUE), \(x) x[1], "")
d2 <- vapply(strsplit(dat$`FINAL-DATE-RANGE`, "-", TRUE), \(x) x[2], "")

dat$Lowpoint <- as.numeric(d1)
dat$Highpoint <- gsub("AD", "", d2)
dat$Highpoint <- gsub("BC", "", dat$Highpoint)
dat$Highpoint <- as.numeric(dat$Highpoint)
dat$Lowpoint[grepl("BC", dat$`FINAL-DATE-RANGE`)] <- dat$Lowpoint[grepl("BC", dat$`FINAL-DATE-RANGE`)] * -1
dat$Highpoint[grepl("BC", dat$`FINAL-DATE-RANGE`)] <- dat$Highpoint[grepl("BC", dat$`FINAL-DATE-RANGE`)] * -1

bioOrder <- dat$Lowpoint
bioOrder[is.na(bioOrder)] <- dat$Midpoint[is.na(bioOrder)]
dat$BioSample <- factor(dat$BioSample, levels = dat$BioSample[order(bioOrder, decreasing = TRUE)])

dat$Midpoint[!is.na(dat$Lowpoint)] <- NA

tosort <- dat$Lowpoint
tosort[is.na(tosort)] <- dat$Midpoint[is.na(tosort)]
dat <- dat[order(tosort), ]

calc_plane <- function(x, MidpointBuffer = 50, rerun = 3) {

  plane <- rep(1, nrow(x))

  get_low_high <- function(x) {
    if (is.na(x$Midpoint)) {
      loP <- x$Lowpoint
      hiP <- x$Highpoint
    } else {
      loP <- x$Midpoint - MidpointBuffer
      hiP <- x$Midpoint + MidpointBuffer
    }
    if (anyNA(c(loP, hiP)))
      NA
    else
      list(loP = loP, hiP = hiP)
  }

  for (i in 2:nrow(x)) {
    p <- get_low_high(x[i, ])
    if (is.na(p)) return(plane)
    for (h in 1:rerun) {
      for (j in 1:(i - 1)) {
        pp <- get_low_high(x[j, ])
        if (plane[j] == plane[i] && any(p$loP:p$hiP %in% pp$loP:pp$hiP)) {
          plane[i] <- plane[i] + 1
        }
      }
    }
  }

  plane

}
dat$Plane <- calc_plane(dat)

dat2 <- dat[, -(1:2)] %>%
  tidyr::gather("X", "Date", 6:7)

dat2$Date2 <- dat2$Date
dat2$Date2[dat2$X == "Lowpoint"] <- dat2$Date2[dat2$X == "Lowpoint"] + 10
dat2$Date2[dat2$X == "Highpoint"] <- dat2$Date2[dat2$X == "Highpoint"] - 10

p2 <- ggplot(dat2, aes(Date, Plane, group = BioSample, colour = BioSample)) +
  # geom_path(size=3.2, colour = "black") +
  # geom_path(aes(Date2, Plane, colour = BioSample), size=2.5, alpha = 0.9) +
  geom_path(size=2.3, colour = "black") +
  geom_path(aes(Date2, Plane, colour = BioSample), size=1.6, alpha = 0.9) +
  geom_point(aes(x = Midpoint, y = Plane+.09, group = BioSample, fill = BioSample),
    data = dat,
    # size = 2, pch = 25,
    size = 1.3, pch = 25,
    alpha = 0.9, colour = "black") +
  geom_point(aes(x = 0, y = 5.5), colour = "white") +
  scale_fill_manual(values = c(viridis_pal(option = "C")(34), "grey", "grey", "grey", "grey")) +
  scale_colour_manual(values = c(viridis_pal(option = "C")(34), "grey", "grey", "grey", "grey")) +
  # scale_fill_viridis(option = "C", discrete = TRUE) +
  # scale_colour_viridis(option = "C", discrete = TRUE) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  ylim(c(0.5, NA)) +
  scale_x_continuous(
    breaks = seq(from = -4000, to = 2000, by = 1000),
    limits = c(-4000, 2000),
    expand = c(0, 0),
    labels = c("4000 BCE", "3000 BCE", "2000 BCE", "1000 BCE", "1 CE", "1000 CE", "2000 CE")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "black"),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black", size = 0.4),
    panel.grid = element_blank()
  )

dat3 <- dat[order(dat$count), ]

p1 <- ggplot(data = world) +
  geom_sf(fill = "beige", lwd = 0.2) +
  coord_sf(expand = FALSE) +
  geom_point(aes(Longitude, Latitude, fill = BioSample),
    dat3, pch = 21, size = 2, alpha = 0.9) +
  # scale_fill_viridis(option = "C", discrete = TRUE, name = NULL) +
  scale_fill_manual(name = NULL,
    values = c(viridis_pal(option = "C")(34), "grey20", "grey40", "grey60", "grey80")) +
  theme_bw() +
  guides(fill = guide_legend(keyheight = 0.7)) +
  theme(
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(colour = grey(0.5), linetype = "dashed", size = 0.2),
    panel.background = element_rect(fill = "azure"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )

p <- plot_grid(NULL, p2, NULL, p1, ncol = 1, align = "v", axis = "lr",
  # rel_heights = c(0.01, 0.195, -0.18, 1)
  rel_heights = c(0.034, 0.195, -0.15, 1)
)
# p
ggsave("figures/map.pdf", plot = p,  width = 9.99, height = 5.28)
# ggsave("figures/map.pdf")
