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

log2fCh_expression <- mean_data %>%
  filter(Sample != "PBS: MØ", Gene != "TGF-β", Gene != "CAMP", Ratio != "(1:1)")%>%
  mutate(log2fCh = log2(FCh))

print(log2fCh_expression)

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

plot_data_new <- plot_data_fch %>%
  filter(Gene %in% c("IL-10", "TNF-α"), Sample != "PBS: MØ")

print(plot_data_new)

time_levels <- levels(factor(plot_data_new$Time))
sample_levels <- levels(factor(plot_data_new$Sample))
n_times <- length(time_levels)
dodge_w <- 0.8

bar_x <- function(sample, time) {
  match(sample, sample_levels) +
    (match(time, time_levels) - (n_times +1)/2)*dodge_w / n_times
}

sig <- data.frame(
  Gene  = "IL-10",                     # must match your Gene level exactly
  xmin  = bar_x("CN: MØ", "2"),  # replace with the Sample the bars belong to
  xmax  = bar_x("CN: MØ", "4"),
  y     = 2,
  label = "**"
)

fill_cols <- c("2" = "#2E2585", "4" = "#E0DEF6")


fch_plot <- ggplot(plot_data_new, aes(x = Gene, y = FCh, fill = factor(Time))) +
  geom_col(position = position_dodge2(width = 0.8),
           width = 0.7) +
  geom_errorbar(aes(ymin = FCh - sd_spss, ymax = FCh + sd_spss),
                position = position_dodge(width = 0.7),
                width = 0.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    data = sig,
    aes( xmin = xmin, xmax = xmax, annotations = label, y_position = y),
    manual = TRUE, inherit.aes = FALSE,
    tip_length = 0.02, textsize = 6, vjust = 0.4
  )+
  scale_fill_manual(values = fill_cols) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic", size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))+
  labs(y = "Fold Change", fill = "Time (h)")
print(fch_plot)
ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/noncedo_qPCR_fch.png", 
       plot = fch_plot, 
       width = 8, 
       height = 8, 
       dpi = 600)

