library(tidyverse)
library(ggsignif)
library(ggpubr)

uncollapsed_data <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/noncedo_relative_quant-qPCR_MBRT.txt",
                               header = TRUE, sep = "\t")
colnames(uncollapsed_data) <- c("Time", "Rep", "Sample", "Ratio", "Gene", "dCt", "fCh")

excluded <- uncollapsed_data %>%
  filter(!(Time == 2 & Rep == 1))

collapse_replicates_dct <- excluded %>%
  group_by(Gene, Sample, Time, Ratio) %>%
  summarise(dCt = mean (dCt, na.rm = TRUE), .groups = "drop")

collapse_replicates_Fch <- excluded %>%
  group_by(Gene, Sample, Time, Ratio) %>%
  summarise(FCh = mean(fCh, na.rm = TRUE), .groups = "drop")


mean_data <- left_join(collapse_replicates_dct, collapse_replicates_Fch, 
                                     by = join_by("Gene", "Sample", "Time", "Ratio"))

print(mean_data)

#Normalised Expression plot

time_cols <- c("2" = "#dc267f",
               "4" = "#ffb000")


#Plot FCh
log2fCh_expression <- mean_data %>%
  filter(Sample != "PBS: MØ", Gene != "TGF-β", Gene != "CAMP", Ratio != "(1:1)")%>%
  mutate(log2fCh = log2(FCh))

print(log2fCh_expression)

time_levels <- levels(factor(log2fCh_expression$Time))
sample_levels <- levels(factor(log2fCh_expression$Sample))
n_times <- length(time_levels)
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

spss_sd <- tribble(
  ~Gene, ~Time, ~Sample,         ~sd_spss,
  "IL-10", 2, "CN: MØ", 0.29086454,
  "IL-10", 4, "CN: MØ", 0.0452,
  "IL-1β", 2, "CN: MØ", 9.1255,
  "IL-1β", 4, "CN: MØ", 0.2579,
  "TNF-α", 2, "CN: MØ", 0.001315765,
  "TNF-α", 4, "CN: MØ", 0.2486
)

print(spss_sd)

fch_comparison <- mean_data %>%
  filter(Ratio != "(1:1)", !Gene %in% c("CAMP", "TGF-β"))

print(fch_comparison)

plot_data_fch<- fch_comparison %>%
  left_join(spss_sd, by = c("Gene", "Time", "Sample"))

print(plot_data)

plot_data_fch %>% filter(is.na(sd_spss)) %>%
  distinct(Gene, Sample)

sample_cols <- c("CN: MØ" = "#dc267f",
                 "PBS: MØ" = "#ffb000")


plot_data_new <- plot_data_fch %>%
filter(Gene %in% c("IL-10", "TNF-α"))

print(plot_data_new)

cn_bars <- plot_data_new %>%
  filter(Sample == "CN: MØ") %>%
  mutate(Group = paste(Time, "h")) %>%
  select(Gene, Group, FCh, sd_spss)

ctrl <- tibble(Gene = "Control", Group = "Control", FCh = 1, sd_spss = NA_real_)

gene_order <- c("Control", "IL-10", "TNF-α")

bars <- bind_rows(ctrl, cn_bars) %>%
  mutate(Gene  = factor(Gene, levels = gene_order),
         Group = factor(Group, levels = c("Control", "2 h", "4 h")))

fill_cols <- c("Control" = "grey60", "2 h" = "#dc267f", "4 h" = "#ffb000")


fch_plot <- ggplot(bars, aes(x = Gene, y = FCh, fill = Group)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"),
           width = 0.7) +
  geom_errorbar(aes(ymin = FCh - sd_spss, ymax = FCh + sd_spss),
                position = position_dodge2(width = 0.8, preserve = "single",
                                           padding = 0.5),
                width = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = fill_cols, name = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(face = c("plain", "italic", "italic"), size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19),
        legend.text = element_text(size = 14),
        legend.position = c(0.3, 0.8)) +
  labs(y = "Fold Change")
print(fch_plot)
ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/noncedo_qPCR_fch.png", 
       plot = fch_plot, 
       width = 7, 
       height = 8, 
       dpi = 600)


# relative_expression_plot <- log2fCh_expression %>%
#   ggplot(aes(x = Sample, y = log2fCh, fill = factor(Time)))+
#   geom_col(position = position_dodge(width = 0.8), width = 0.7) +
#   geom_hline(yintercept = 0) +
#   geom_signif(data = sig,
#                aes( xmin = xmin, xmax = xmax, annotations = label, y_position = y),
#               manual = TRUE, inherit.aes = FALSE, 
#               tip_length = 0.02, textsize = 8, vjust = 0.4)+
#   facet_grid(~ Gene, scales = "free_y")+
#   scale_fill_manual(name = "Time (h)", values = time_cols)+
#   theme_classic() +
#   theme(strip.text = element_text(face = "italic", size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         axis.text = element_text(size=14),
#         axis.title = element_text(size=19),
#         axis.title.x = element_blank())+
#   labs(y = "Relative Expression")
# 
# print(relative_expression_plot)
# ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/noncedo_qPCR_relative_expression_plot.png", 
#        plot = relative_expression_plot, 
#        width = 10, 
#        height = 8, 
#        dpi = 600)

# normalised_expression_plot <- mean_data %>%
#   ggplot(aes(x = Sample, y = dCt, fill = factor(Time)))+
#   geom_col(position = position_dodge(width = 0.7), width = 0.6) +
#   geom_hline(yintercept = 0) +
#   facet_grid(Gene ~ Ratio)+
#   scale_fill_manual(name = "Time (h)", values = time_cols)+
#   theme_classic() +
#   theme(strip.text.x.top = element_text(size = 12),
#         strip.text.y.right = element_text(face = "italic", size = 12),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         axis.text = element_text(size=14),
#         axis.title = element_text(size=20),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   labs(y = "Normalised Expression")
# 
# print(normalised_expression_plot)
# ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/noncedo_qPCR_normalise_expression_plot.png", 
#        plot = normalised_expression_plot, 
#        width = 10, 
#        height = 8, 
#        dpi = 600)
