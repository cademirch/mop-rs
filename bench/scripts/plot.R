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

# Parse times for moprs components
moprs_time <- parse_time(file.path(benchmark_dir, 'moprs.txt'))
merge_d4_time <- parse_time(file.path(benchmark_dir, 'merge_d4.txt'))

# Parse mosdepth times and calculate the mean time
mosdepth_dir <- file.path(benchmark_dir, 'mosdepth')
mosdepth_files <- list.files(mosdepth_dir, pattern = '\\.txt$', full.names = TRUE)
mosdepth_times <- sapply(mosdepth_files, parse_time)
max_mosdepth_time <- max(mosdepth_times)

# Calculate the total time for moprs components
total_moprs_time <- max_mosdepth_time + merge_d4_time + moprs_time

# Create a table of moprs components and total time
moprs_components_table <- data.frame(
  Component = c('mosdepth', 'd4tools merge', 'moprs', 'Total'),
  Time = c(max_mosdepth_time, merge_d4_time, moprs_time, total_moprs_time)
)

# Save the table as a CSV file
write.csv(moprs_components_table, file = 'results/moprs_components_table.csv', row.names = FALSE)
