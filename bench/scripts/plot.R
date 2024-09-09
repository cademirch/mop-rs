library("ggplot2")

# Function to parse the benchmark files and extract the execution time in seconds
parse_time <- function(file_path) {
  lines <- readLines(file_path)
  # Extract the time in seconds from the first line (after the header)
  time_in_seconds <- as.numeric(strsplit(lines[2], '\t')[[1]][1])
  return(time_in_seconds)
}

# Directory containing the benchmark files
benchmark_dir <- 'results/benchmarks'

# Parse mop time
mop_time <- parse_time(file.path(benchmark_dir, 'mop.txt'))

# Parse moprs times
moprs_time <- parse_time(file.path(benchmark_dir, 'moprs.txt'))
merge_d4_time <- parse_time(file.path(benchmark_dir, 'merge_d4.txt'))

# Parse mosdepth times and calculate the mean time
mosdepth_dir <- file.path(benchmark_dir, 'mosdepth')
mosdepth_files <- list.files(mosdepth_dir, pattern = '\\.txt$', full.names = TRUE)
mosdepth_times <- sapply(mosdepth_files, parse_time)
mean_mosdepth_time <- mean(mosdepth_times)

# Prepare data for plotting
data <- data.frame(
  Tool = c('mop', rep('moprs', 3)),
  Component = c('mop', 'mosdepth mean time', 'merge_d4 time', 'moprs time'),
  Time = c(mop_time, mean_mosdepth_time, merge_d4_time, moprs_time)
)

total_time_data <- data.frame(
  Tool = c('mop', 'moprs'),
  Time = c(mop_time, moprs_time + merge_d4_time + mean_mosdepth_time)
)

# Prepare data for moprs components plot
moprs_data <- data.frame(
  Component = c('mosdepth mean time', 'merge_d4 time', 'moprs time'),
  Time = c(mean_mosdepth_time, merge_d4_time, moprs_time)
)

# Create the plot
plot <- ggplot(data, aes(x = Tool, y = Time, fill = Component)) +
  geom_bar(stat = 'identity', position = 'stack') +
  labs(y = 'Time (seconds)', title = 'Execution Time of mop and moprs(+friends)')
total_time_plot <- ggplot(total_time_data, aes(x = Tool, y = Time, fill = Tool)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Time (seconds)', title = 'Total Execution Time Comparison')
  

# Create the moprs components plot
moprs_plot <- ggplot(moprs_data, aes(x = Component, y = Time, fill = Component)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Time (seconds)', title = 'moprs Components Execution Time') +
  scale_fill_manual(values = c('orange', 'green', 'red'))
# Save the plot
ggsave('results/plot.png', plot)
ggsave('results/total_time_plot.png', total_time_plot)
ggsave('results/moprs_components_plot.png', moprs_plot)

