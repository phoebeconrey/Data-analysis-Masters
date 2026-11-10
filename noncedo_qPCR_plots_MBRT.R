library(tidyverse)
library(ggsignif)
library(ggh4x)

uncollapsed_data <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/noncedo_relative_quant-qPCR_MBRT.txt",
                               header = TRUE, sep = "\t")
colnames(uncollapsed_data) <- c("Time", "Rep", "Sample", "Ratio", "Gene", "dCt", "fCh")

collapse_replicates_dct <- uncollapsed_data %>%
  group_by(Gene, Sample, Time, Ratio) %>%
  summarise(dCt = mean (dCt, na.rm = TRUE), .groups = "drop")

collapse_replicates_Fch <- uncollapsed_data %>%
  group_by(Gene, Sample, Time, Ratio) %>%
  summarise(FCh = mean(fCh, na.rm = TRUE), .groups = "drop")


mean_data <- left_join(collapse_replicates_dct, collapse_replicates_Fch, 
                                     by = join_by("Gene", "Sample", "Time", "Ratio"))

print(mean_data)

#Normalised Expression plot

time_cols <- c("2" = "#dc267f",
               "4" = "#ffb000")

normalised_expression_plot <- mean_data %>%
  ggplot(aes(x = Sample, y = dCt, fill = factor(Time)))+
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_grid2(Gene ~ Ratio, scales = "free_y")+
  scale_fill_manual(name = "Time (h)", values = time_cols)+
  theme_classic() +
  theme(strip.text.x.top = element_text(size = 12),
        strip.text.y.right = element_text(face = "italic", size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Normalised Expression")

print(normalised_expression_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/noncedo_qPCR_normalise_expression_plot.png", 
       plot = normalised_expression_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)

#Plot FCh

log2fCh_expression <- mean_data %>%
  filter(Sample != "PBS: MØ")%>%
  mutate(log2fCh = log2(FCh))

print(log2fCh_expression)

relative_expression_plot <- log2fCh_expression %>%
  ggplot(aes(x = Ratio, y = log2fCh, fill = factor(Time)))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_grid(~ Gene, scales = "free_y")+
  scale_fill_manual(name = "Time (h)", values = time_cols)+
  theme_classic() +
  theme(strip.text = element_text(face = "italic", size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Expression")

print(relative_expression_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/noncedo_qPCR_relative_expression_plot.png", 
       plot = relative_expression_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)
