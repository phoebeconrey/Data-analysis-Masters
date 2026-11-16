library(tidyverse)
library(ggsignif)
library(ggbreak)
library(scales)

urease_raw <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/TM_urease_data.txt",
                          header = TRUE, sep = "\t", dec = ".")

collapsed_urease <- urease_raw %>%
  group_by(Time, Temp, Condition)%>%
  summarise(Abs = mean(Abs, na.rm = TRUE), .groups = "drop")

clean_urease <- collapsed_urease[-25,]

print(clean_urease)

#function ----

orig <- c(0.5, 1, 1.5, 2, 24, 48)
comp <- c(0.5, 1, 1.5, 2, 2.5, 3)   # smaller steps = narrower gaps; adjust to taste

fwd <- function(x) approx(orig, comp, xout = x, rule = 2)$y
inv <- function(y) approx(comp, orig, xout = y, rule = 2)$y
broken <- new_transform("broken", transform = fwd, inverse = inv)

#alternative plot----

cond_cols <- c("YM Broth" =  "#ffb000",
               "MM" = "#dc267f")

plot_urease <- clean_urease%>%
  filter(Time %in% c("0.5", "48"))

urease_plot_simple <- plot_urease %>%
  ggplot(aes(x = factor(Time), y = Abs, fill = Condition))+
  geom_col(position = position_dodge(width = 0.8), width = 0.7)+
  facet_wrap(~ Temp)+
  scale_fill_manual(values = cond_cols)+
  theme_classic()
print(urease_plot_simple)

#Plot Figure

temp_cols <- c("25" =  "#ffb000",
               "37" = "#dc267f")

urease_plot <- clean_urease %>%
  ggplot(aes(x = Time, y = Abs, colour = factor(Temp)))+
  geom_line(linewidth = 1)+
  geom_point(size = 2)+
  facet_grid(~ Condition, scales = "free_x", space = "free_x")+
  scale_x_continuous(
    transform = broken,
    breaks = c(0.5, 1, 2, 24, 48)) +
  annotate("text", x = 13, y = -Inf, label = "//", vjust = 1.2, size = 5) +
  annotate("text", x = 36, y = -Inf, label = "//", vjust = 1.2, size = 5) +
  coord_cartesian(clip = "off")+
  scale_colour_manual(values = temp_cols)+
  theme_classic()+
  theme(strip.text= element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20))+
  labs(y = "Abs (490 nm)", x = "Time (h)", colour = "Temp (°C)")
print(laccase_plot)
ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/TM_Laccase_plot.png", 
       plot = laccase_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)