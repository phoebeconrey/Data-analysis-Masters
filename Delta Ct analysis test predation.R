library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)

data <- "C:/Users/User/OneDrive/Documents/Research/Data analysis/test_predation_ct_values_r.txt"

data.Ct <- read_Ct_long(path = data, 
                        sep = "\t",
                        dec = ".",
                        skip = 0,
                        add.column.Flag = FALSE,
                        column.Sample = 2,
                        column.Gene = 3,
                        column.Ct = 4,
                        column.Group = 1)

data.CtF.ready <- make_Ct_ready(data = data.Ct,
                                imput.by.mean.within.groups = TRUE)
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         normalise = TRUE,
                         ref = c("GAPDH 3", "Act 1"),
                         transform = TRUE)

results.pn.dCt <- RQ_dCt(data = data.dCt.exp,
                      do.tests = FALSE, 
                      group.study = "Predated",
                      group.ref = "Naive")

results.vn.dCt <- RQ_dCt(data = data.dCt.exp,
                        do.tests = FALSE, 
                        group.study = "Virulence Control",
                        group.ref = "Naive")
results.zero.dCt <- RQ_dCt(data = data.dCt.exp,
                        do.tests = FALSE, 
                        group.study = "Naive",
                        group.ref = "Naive")
print(results.vn.dCt)
print(results.pn.dCt)


naive_comparison <- results.pn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Predated vs Naive")

vnaive_comparison <- results.vn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "37 °C vs Naive")

zero_comparison <- results.zero.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Naive vs Naive")

plot_data <- bind_rows(naive_comparison, vnaive_comparison, zero_comparison)

plot_data <- plot_data %>%
  mutate(Comparison = fct_recode(Comparison,
                                 "Predated" = "Predated vs Naive",
                                 "37 °C" = "37 °C vs Naive",
                                 "Naive" = "Naive vs Naive"))
plot_data <- plot_data %>%
  mutate(Comparison = fct_relevel(Comparison,
                                  "Naive",
                                  "Predated",
                                  "37 °C"))

log10FCh.plot <- ggplot(plot_data, aes(x = Comparison, y = FCh, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  ylab("Relative Gene Expression")+
  xlab("")+
  scale_fill_manual(values = c("Lac 1" = "#5e3c99",
                                "Ure 5" = "#b2abd2"))+
  theme_classic()+
  theme(legend.text = element_text(face = "italic", size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
         axis.title = element_text(size=24))

print(log10FCh.plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/relative_gene_expre.png", 
       plot = log10FCh.plot, 
       width = 10, 
       height = 8, 
       dpi = 600)