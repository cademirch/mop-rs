import matplotlib.pyplot as plt
import pandas as pd
import os

# Function to parse the benchmark files and extract the execution time in seconds
def parse_time(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Extract the time in seconds from the first line (after the header)
        time_in_seconds = float(lines[1].split('\t')[0])
    return time_in_seconds

# Directory containing the benchmark files
benchmark_dir = 'results/benchmarks'

# Parse mop time
mop_time = parse_time(os.path.join(benchmark_dir, 'mop.txt'))

# Parse moprs times
moprs_time = parse_time(os.path.join(benchmark_dir, 'moprs.txt'))
merge_d4_time = parse_time(os.path.join(benchmark_dir, 'merge_d4.txt'))

# Parse mosdepth times and calculate the mean time
mosdepth_times = []
mosdepth_dir = os.path.join(benchmark_dir, 'mosdepth')
for file in os.listdir(mosdepth_dir):
    if file.endswith('.txt'):
        mosdepth_times.append(parse_time(os.path.join(mosdepth_dir, file)))

mean_mosdepth_time = sum(mosdepth_times) / len(mosdepth_times)

# Prepare data for plotting
moprs_components = [mean_mosdepth_time, merge_d4_time, moprs_time]
labels = ['mop', 'moprs']
moprs_labels = ['mosdepth mean time', 'merge_d4 time', 'moprs time']

# Plotting
fig, ax = plt.subplots()

# Bar for mop
ax.bar(labels[0], mop_time, label='mop time', color='blue')

# Stacked bar for moprs
ax.bar(labels[1], moprs_components[0], label=moprs_labels[0], color='orange')
ax.bar(labels[1], moprs_components[1], bottom=moprs_components[0], label=moprs_labels[1], color='green')
ax.bar(labels[1], moprs_components[2], bottom=sum(moprs_components[:2]), label=moprs_labels[2], color='red')

# Adding labels and title
ax.set_ylabel('Time (seconds)')
ax.set_title('Execution Time Comparison of mop and moprs')
ax.legend()
plt.savefig("results/plot.png")
plt.show()
