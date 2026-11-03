library(ggplot2)

ure_y <- c(24.455, 26.715, 29.13, 31.63)
ure_x <- c(0, -0.698970004, -1.397940009, -2.096910013)

ure_plot <- ggplot(mapping = aes(y = ure_y, x = ure_x))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#E69F00")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  ylab("Ct")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=24))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-2), y = max(25),
    label = ("R² = 0.999\nEfficiency = 95.9%"),
    hjust = 0, vjust = 1,
    size = 8)
print(ure_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/ure_eff.png", plot = ure_plot, width = 10, height = 8, dpi = 600)

lac_y <- c(24.72, 26.25, 28.98, 31.68)
lac_x <- c(0, -0.698970004, -1.397940009, -2.096910013)

lac_plot <- ggplot(mapping = aes(y = lac_y, x = lac_x))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#009E73")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  ylab("Ct")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=24))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-2), y = max(25),
    label = ("R² = 0.985\nEfficiency = 97.7%"),
    hjust = 0, vjust = 1,
    size = 8)
print(lac_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/lac_eff.png", plot = lac_plot, width = 10, height = 8, dpi = 600)

gap_y <- c(25.54, 28.16, 30.635, 32.715)
gap_x <- c(-0.698970004, -1.397940009, -2.096910013, -2.795880017)

gap_plot <- ggplot(mapping = aes(y = gap_y, x = gap_x))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#56B4E9")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  ylab("Ct")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=24))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-2), y = max(25),
    label = ("R² = 0.992\nEfficiency = 95.5%"),
    hjust = 0, vjust = 1,
    size = 8)
print(gap_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/gap_eff.png", plot = gap_plot, width = 10, height = 8, dpi = 600)

act_y <- c(25.62, 28.32, 30.38, 32.98)
act_x <- c(0, -0.698970004, -1.397940009, -2.096910013)

act_plot <- ggplot(mapping = aes(y = act_y, x = act_x))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#0072B2")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  ylab("Ct")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=24))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-2), y = max(25),
    label = ("R² = 0.978\nEfficiency = 94.7%"),
    hjust = 0, vjust = 1,
    size = 8)
print(act_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/act_eff.png", plot = act_plot, width = 10, height = 8, dpi = 600)

cap_y <- c(25.915, 27.665, 30.4, 33.54)
cap_x <- c(0, -0.698970004, -1.397940009, -2.096910013)

cap_plot <- ggplot(mapping = aes(y = cap_y, x = cap_x))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D55E00")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  ylab("Ct")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=24))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-2), y = max(25),
    label = ("R² = 0.984\nEfficiency = 87.4%"),
    hjust = 0, vjust = 1,
    size = 8)
print(cap_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/cap_eff.png", plot = cap_plot, width = 10, height = 8, dpi = 600)