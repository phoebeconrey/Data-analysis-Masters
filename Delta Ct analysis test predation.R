library(tidyverse)
library(coin)
library(ggsignif)

ct_data <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/test_predation_ct_values_r_replicates.txt",
                        header = TRUE, sep = "\t") 

colnames(ct_data) <- c("Group", "Sample", "Gene", "Ct")
ct_data <- ct_data %>%
  rename(
    Target = Gene,
    Cq = Ct
  )
ct_data$Target <- gsub(fixed("Act 1"), "Actin", ct_data$Target)

#Gene Expression Analysis function ----
mean_cq <- function(data){
  data %>%
    group_by(Target, Sample, Group) %>%
    summarise(Cq = mean (Cq, na.rm = TRUE), .groups = "drop")
}

average_hk <- function(data, hk_genes) {
  data %>%
    filter(Target %in% hk_genes) %>%
    group_by(Sample) %>%
    summarise(hk_Ct = mean(Cq, na.rm = TRUE), .groups = "drop")
}

deltaCt <- function(data, hk_genes, target_gene){
  hk = average_hk (data, hk_genes)
  
  data %>%
    filter(Target %in% target_gene) %>%
    left_join(hk, by = "Sample") %>%
    mutate(dCt = Cq - hk_Ct)
}

ddct <- function(data, group.ref){
  group.ref.mean <- data %>%
    filter(Group == group.ref) %>%
    group_by(Target) %>%
    summarise(ref_dCt = mean(dCt, na.rm = TRUE), .groups = "drop")
  
  ddct.result <- data %>%
    #filter(Group == group.study) %>%
    left_join(group.ref.mean, by = "Target") %>%
    mutate(ddCt = dCt - ref_dCt)
}

fCh <- function(data){
  data %>%
    mutate(fCh = 2^(-ddCt))
}

log2fCh <- function(data) {
  data %>% 
    mutate(log2fCh = log2(fCh))
}

AP_dCt <- function(data, hk_genes, target_gene, group.ref){
  data %>%
    mean_cq() %>%
    deltaCt(hk_genes, target_gene)%>%
    ddct(group.ref) %>%
    fCh() %>%
    log2fCh()
}

# Delta Ct Manual ----
result_dCt <- AP_dCt(ct_data, 
                     hk_genes = c("Actin", "GAPDH"), 
                     target_gene = c("Lac 1", "Ure 5"),
                     group.ref = "Naive")
print(result_dCt)

write.table(result_dCt, file = "test_predation_gene_analysis_data.txt", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)

#plot normalised count ----

plot_normalised_data <- result_dCt %>%
  dplyr::select(Target, Sample, dCt, Group)

group_cols <- c("Naive" = "#66c2a5",
                "Temp" = "#74a9cf", 
                "Predated" = "#a50f15")

normalised_plot <- plot_normalised_data %>% 
  ggplot(aes(x = Sample, y = dCt, fill = Group))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ Target, scales = "free_y")+
  scale_fill_manual(values = group_cols)+
  scale_x_discrete(labels = c("CNCE(25)" = "Picked", 
                              "CN(25)" = "25 °C",
                              "CN(37)" = "37 °C"))+
  theme_classic() +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=10),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Normalised Expression")

print(normalised_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/normalised_counts_test_qpcr.png", 
       plot = normalised_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)

# plot log fold change ----

plot_log_data <- result_dCt %>%
  filter(Group != "Naive")%>%
  dplyr::select(Target, Sample, Group, log2fCh)

gene_cols <- c("Lac 1" = "#5e3c99",
               "Ure 5" = "#e66101",
               "CAP 10" = "#b2abd2",
               "COX 1" = "#fdb863")

log_plot <- plot_log_data %>%
  ggplot(aes(x = Sample, y = log2fCh, fill = Target))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_grid(Target ~ Group, scales = "free", space = "free_x")+
  scale_fill_manual(values = gene_cols)+
  scale_x_discrete(labels = c("CNCE(25)" = "Picked", 
                              "CN(25)" = "25 °C",
                              "CN(37)" = "37 °C"))+
  theme_classic() +
  theme(strip.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 18, face = "italic"),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=10),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Gene Expression")
print(log_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/relative_gene_expression_test_qpcr.png", 
       plot = log_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)