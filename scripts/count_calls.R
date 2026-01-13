
pkgs <- c("qqconf", "qqplotr", "ggplot2", "cowplot")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(qqconf)
library(qqplotr)
library(ggplot2)
library(cowplot)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript count_calls.R <output_directory>")
}
output_dir <- args[1]

# Construct the path to the sample_counts.tsv file
sample_counts_file <- file.path(output_dir, "sample_counts.tsv")

# Read the sample counts file
sample_counts <- read.delim(sample_counts_file)
mean <- mean(sample_counts$Count)
sd <- sd(sample_counts$Count)
sample_counts$outlier <- F
sample_counts$outlier[sample_counts$Count < mean -3*sd | sample_counts$Count > mean+3*sd] <-T

# Plot 1
p1 <- ggplot(sample_counts, aes(x = Count , y = after_stat(density)) )+ geom_histogram(fill = "grey", bins = 100) + 
  theme_bw() + ylab("Density") + 
  geom_vline(xintercept = c(mean, mean+sd, mean+2*sd, mean+3*sd, mean-sd, mean-2*sd, mean-3*sd), 
             lty = 2, lwd = 0.1, color = "black") + xlab("EHdn detected tandem repeats") +
  annotate("text", x = c(mean - 200, mean+sd - 200, mean+2*sd - 200, mean+3*sd - 200,
                         mean-sd - 200, mean-2*sd - 200, mean-3*sd - 200), 
           y =  0.0004, label = 
             c("mean", "+1SD", "+2SD", "+3SD", "-1SD", "-2SD", "-3SD"), angle = 90, cex = 3)

# Plot 2
p2 <- ggplot(sample_counts, aes(sample = Count)) +  stat_qq_band() + stat_qq_point(color = "steelblue") + stat_qq_line() +
  scale_color_manual(values = c("steelblue", "red")) + ylab("Sample Quantiles") + xlab("Theoretical Quantiles") +
  geom_vline(xintercept = c(mean, mean+sd, mean+2*sd, mean+3*sd, mean-sd, mean-2*sd, mean-3*sd), lty = 2, lwd = 0.1, color = "black") + 
  annotate("text", x = c(mean - 200, mean+sd - 200, mean+2*sd - 200, mean+3*sd - 200,
                         mean-sd - 200, mean-2*sd - 200, mean-3*sd - 200), y = 10000, label = 
             c("mean", "+1SD", "+2SD", "+3SD", "-1SD", "-2SD", "-3SD"), angle = 90, cex = 3) 

# Display the plots side by side
gridExtra::grid.arrange(p1, p2, nrow = 1)
ggsave(file.path(output_dir, "ext.fig2a.png"), width = 10, height = 4)

### 2SD
tmp.clean <- sample_counts[sample_counts$Count > mean-2*sd &
                             sample_counts$Count < mean + 2*sd, ]

p1 <- ggplot(tmp.clean, aes(x = Count, y = after_stat(density))) + geom_histogram(fill = "grey", bins = 100) +  xlab("EHdn detected tandem repeats") +
  theme_bw() + ylab("Density")
p2 <- ggplot(tmp.clean, aes(sample = Count)) + stat_qq_band() +
  stat_qq_point(color = "steelblue") + stat_qq_line() +
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles")
plot_grid(p1, p2, nrow = 1, rel_widths = c(1.5, 1))
ggsave(file.path(output_dir, "quality.control2SD.png"), width = 10, height = 4)

### 3SD
tmp.clean <- sample_counts[sample_counts$Count > mean-3*sd &
                             sample_counts$Count < mean + 3*sd, ]

p1 <- ggplot(tmp.clean, aes(x = Count, y = after_stat(density))) + geom_histogram(fill = "grey", bins = 100) +  xlab("EHdn detected tandem repeats") +
  theme_bw() + ylab("Density")
p2 <- ggplot(tmp.clean, aes(sample = Count)) + stat_qq_band() + 
  stat_qq_point(color = "steelblue") + stat_qq_line() +
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles")
plot_grid(p1, p2, nrow = 1, rel_widths = c(1.5, 1))
ggsave(file.path(output_dir, "quality.control3SD.png"), width = 10, height = 4)

write.table(sample_counts, file.path(output_dir, "outliers.tsv"), sep = "\t", row.names = FALSE)

