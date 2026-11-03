library(tidyr)
library(dplyr)
library(drc)
library(ggplot2)


efficiency_data <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/cap_and_plb1_altered_primer_conc_efficiency_r.txt",
                              header = TRUE, sep = "\t")

colnames(efficiency_data) <- c("Conc", "Rep", "Cap", "PLB1")

efficiency_long <- efficiency_data %>% pivot_longer(cols = -c("Conc", "Rep"),
                                                    names_to = "Gene", values_to = "Ct")

cap_data <- efficiency_long %>% 
  filter(Gene == "Cap", Conc != 0) %>%
  group_by(Conc) %>%
  summarise(Ct = mean(Ct))
plb1_data <- efficiency_long %>% 
  filter (Gene == "PLB1", Conc != 0) %>%
  group_by(Conc) %>%
  summarise(Ct = mean(Ct, na.rm = TRUE))


plb1_plot <- ggplot(plb1_data, aes(x = Conc, y = Ct))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#5e3c99")+
  geom_vline(xintercept = 0) +
  xlab(expression(log[10]~"[sample]"))+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.text = element_text(size=12),
              axis.title = element_text(size=20))+
  scale_y_continuous(limits = c(0, NA))+
  annotate(
    "text",
    x = min(-3), y = max(30),
    label = paste0(
      "R² = ", round(get_efficiency(plb1_model)$R2, 3),
      "\nEfficiency = ", round(get_efficiency(plb1_model)$efficiency_pct, 1), "%"
    ),
    hjust = 0, vjust = 1,
    size = 4
  )
print(plb1_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/plb1.png", plb1_plot, width = 6, height = 5, dpi = 600)

plb1_model <- lm(Ct ~ Conc, data = plb1_data)

summary(cap_model)
summary(plb1_model)

get_efficiency <- function(model) {
  slope <- coef(model)[["Conc"]]
  r2 <- summary(model)$r.squared
  efficiency <- (10^(-1 / slope) - 1) * 100
  
  data.frame(slope = slope, R2 = r2, efficiency_pct = efficiency)
}

print(get_efficiency(cap_model))
print(get_efficiency(plb1_model))


