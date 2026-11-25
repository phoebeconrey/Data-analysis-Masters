library(tidyverse)

capsule_raw <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/TM_capsule_data.txt",
                         header = TRUE, sep = "\t", dec = ".")

capsule_size <- capsule_raw %>%
  mutate(
    Capsule = Total.Cell - Cell.Body
  )

capsule_flagged <- capsule_size %>%
  group_by(Temperature, Medium) %>%
  mutate(
    lower = quantile(Capsule, 0.25, na.rm = TRUE) - 1.5 * IQR(Capsule, na.rm = TRUE),
    upper = quantile(Capsule, 0.75, na.rm = TRUE) + 1.5 * IQR(Capsule, na.rm = TRUE),
    outlier = Capsule < lower | Capsule > upper
  ) %>%
  ungroup()

capsule_clean <- capsule_flagged %>%
  filter(!outlier)

collapsed_capsule <- capsule_clean %>%
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

#Box plot ----
box_cap_plot <- ggplot(capsule_clean, aes(x = factor(Temperature), y = Capsule, fill = factor(Temperature)))+
  geom_point(aes(colour = factor(Temperature)), position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5)+
  geom_boxplot(outlier.shape = NA, staplewidth = 0.4)+
  scale_fill_manual(values = temp_cols)+
  scale_colour_manual(values = temp_cols)+
  theme_classic()+
  theme(strip.text= element_text(size = 12),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20),
        legend.position = "none")+
  labs(y = "Capsule (μm)", x = "Temperature (°C)", fill = NULL)
 

print(box_cap_plot)
ggsave("C:/Users/pheebsy/OneDrive/Documents/Research/Results/TM_Capsule_boxplot.png", 
       plot = box_cap_plot, 
       width = 5, 
       height = 8, 
       dpi = 600)

#Bar graph ----
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