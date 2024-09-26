library(ggplot2)
library(stringr)  # For str_wrap function

# Your existing ggplot code
ggplot(top20, aes(x = reorder(term_name, -p.adjust), y = Count, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Make plot horizontal
  scale_fill_viridis_c(option = "cividis", direction = -1, name = "Adjusted P-Value") +  # Use cividis color map
  labs(title = header_F, x = "Term Name", y = "Count") +
  theme_minimal(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),  # White plot background
    panel.background = element_rect(fill = "white", colour = NA),  # White panel background
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, margin = margin(0, 0, 10, 0)),  # Adjust margin if needed for wrapped text
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right"
  ) +
  theme(plot.margin = margin(1, 1, 1.5, 3, "cm")) +
  scale_y_continuous(labels = scales::comma) +  # Ensure Y-axis labels are readable
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))  # Wrap X-axis labels
# Save the plot
ggsave(paste0(results_folder, "/",header_F, ".tiff"), width = 12, height = 12, units = 'in', dpi = 300)
