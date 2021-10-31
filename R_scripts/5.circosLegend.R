# Benjamin Jean-Marie Tremblay
# 2021-10-31

library(ggplot2)

PctIDcolours <- rev(viridis::viridis_pal()(12))[-1]

p <- ggplot() +
  geom_segment(aes(x = 89, xend = 90, y = 1, yend = 1), size = 10, colour = PctIDcolours[1]) +
  geom_segment(aes(x = 90, xend = 91, y = 1, yend = 1), size = 10, colour = PctIDcolours[2]) +
  geom_segment(aes(x = 91, xend = 92, y = 1, yend = 1), size = 10, colour = PctIDcolours[3]) +
  geom_segment(aes(x = 92, xend = 93, y = 1, yend = 1), size = 10, colour = PctIDcolours[4]) +
  geom_segment(aes(x = 93, xend = 94, y = 1, yend = 1), size = 10, colour = PctIDcolours[5]) +
  geom_segment(aes(x = 94, xend = 95, y = 1, yend = 1), size = 10, colour = PctIDcolours[6]) +
  geom_segment(aes(x = 95, xend = 96, y = 1, yend = 1), size = 10, colour = PctIDcolours[7]) +
  geom_segment(aes(x = 96, xend = 97, y = 1, yend = 1), size = 10, colour = PctIDcolours[8]) +
  geom_segment(aes(x = 97, xend = 98, y = 1, yend = 1), size = 10, colour = PctIDcolours[9]) +
  geom_segment(aes(x = 98, xend = 99, y = 1, yend = 1), size = 10, colour = PctIDcolours[10]) +
  geom_segment(aes(x = 99, xend = 100, y = 1, yend = 1), size = 10, colour = PctIDcolours[11]) +
  scale_y_continuous(expand = c(0, 0), labels = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(89.01, 90:99, 99.99), labels = paste0(89:100, "%"),
    position = "top") +
  xlab(NULL) + ylab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10)
  )

ggsave("CircosPercentID_legend.pdf", plot = p, width = 0.816, height = 2.62)
