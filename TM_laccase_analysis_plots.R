library(tidyverse)
library(ggsignif)
library(ggbreak)
library(scales)

laccase_raw <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/TM_laccase_data.txt",
                          header = TRUE, sep = "\t", dec = ".")

collapsed_laccase <- laccase_raw %>%
  group_by(Time, Temp, Condition)%>%
  summarise(Abs = mean(Abs, na.rm = TRUE), .groups = "drop")

clean_laccase <- collapsed_laccase[-25,]

print(clean_laccase, n = 24)



print(shapiro.test(clean_laccase$Abs))

wide <- clean_laccase %>%
  pivot_wider(names_from = Condition, values_from = Abs)

## MM vs YM broth ----

cond_tests <- wide %>%
  group_by(Temp)%>%
  group_modify(~{
    t <- wilcox.test(.x$MM, .x$'YM Broth', paired = TRUE)
    tribble(n_pairs = nrow(.x), V = unname(t$statistic), p = t$p.value)
  })%>%
  ungroup()

print(cond_tests)

#function

condition_levels <- levels(factor(clean_laccase$Condition))
time_levels <- levels(factor(clean_laccase$Time))
n_times <- length(condition_levels)
dodge_w <- 0.8

bar_x <- function(sample, time) {
  match(sample, sample_levels) +
    (match(time, time_levels) - (n_times +1)/2)*dodge_w / n_times
}

sig <- data.frame(
  Gene = "IL-10",
  xmin = bar_x("CN: MØ", "2"),
  xmax = bar_x("CN: MØ", "4"),
  y = 3,
  label = "**"
)

orig <- c(0.5, 1, 1.5, 2, 24, 48)
comp <- c(0.5, 1, 1.5, 2, 2.5, 3)   # smaller steps = narrower gaps; adjust to taste

fwd <- function(x) approx(orig, comp, xout = x, rule = 2)$y
inv <- function(y) approx(comp, orig, xout = y, rule = 2)$y
broken <- new_transform("broken", transform = fwd, inverse = inv)

#alternative plot

cond_cols <- c("YM Broth" =  "#ffb000",
               "MM" = "#dc267f")

plot_laccase <- clean_laccase%>%
  filter(Time %in% c("0.5", "48"))

laccase_plot_simple <- plot_laccase %>%
  ggplot(aes(x = factor(Time), y = Abs, fill = Condition))+
  geom_col(position = position_dodge(width = 0.8), width = 0.7)+
  facet_wrap(~ Temp)+
  scale_fill_manual(values = cond_cols)+
  theme_classic()
print(laccase_plot_simple)

# #Plot Figure
# 
# temp_cols <- c("25" =  "#ffb000",
#                "37" = "#dc267f")
# 
# laccase_plot <- clean_laccase %>%
#   ggplot(aes(x = Time, y = Abs, colour = factor(Temp)))+
#   geom_line(linewidth = 1)+
#   geom_point(size = 2)+
#   facet_grid(~ Condition, scales = "free_x", space = "free_x")+
#   scale_x_continuous(
#     transform = broken,
#     breaks = c(0.5, 1, 2, 24, 48)) +
#   annotate("text", x = 13, y = -Inf, label = "//", vjust = 1.2, size = 5) +
#   annotate("text", x = 36, y = -Inf, label = "//", vjust = 1.2, size = 5) +
#   coord_cartesian(clip = "off")+
#   scale_colour_manual(values = temp_cols)+
#   theme_classic()+
#   theme(strip.text= element_text(size = 12),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         axis.text = element_text(size=14),
#         axis.title = element_text(size=20))+
#   labs(y = "Abs (490 nm)", x = "Time (h)", colour = "Temp (°C)")
# print(laccase_plot)
# ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/TM_Laccase_plot.png", 
#        plot = laccase_plot, 
#        width = 10, 
#        height = 8, 
#        dpi = 600)