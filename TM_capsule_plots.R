library(tidyverse)

capsule_raw <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/TM_capsule_data.txt",
                         header = TRUE, sep = "\t", dec = ".")

capsule_size <- capsule_raw %>%
  mutate(
    Capsule = Total.Cell - Cell.Body,
  )

collapsed_capsule <- capsule_size %>%
  group_by(Temperature, Medium)%>%
  summarise(
    c_sd = sd(Capsule, na.rm = TRUE),
    n = sum(!is.na(Capsule)),
    Capsule_mean = mean(Capsule, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(c_se = c_sd/sqrt(n))

print(collapsed_capsule)


temp_cols <- c("25" =  "#ffb000",
               "37" = "#dc267f")

capsule_plot <- collapsed_capsule %>%
  ggplot(aes(x = factor(Temperature), y = Capsule_mean, fill = factor(Temperature)))+
  geom_col(position = position_dodge(width = 0.8), width = 0.7)+
  facet_wrap(~ Medium)+
  geom_errorbar(aes(ymin = Capsule_mean - c_se, ymax = Capsule_mean + c_se),
                position = position_dodge(width = 0.9),
                width = 0.2)+
  scale_fill_manual(values = temp_cols)+
  theme_classic()+
  theme(strip.text= element_text(size = 12),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20),
        legend.position = "none")+
  labs(y = "Capsule (μm)", x = "Temperature (°C)", fill = NULL)
print(capsule_plot)
ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/TM_Capsule_plot.png", 
       plot = capsule_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)