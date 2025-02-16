import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker

# Load the CSV file
if len(sys.argv) < 2:
    print("Usage: python script.py <csv_file>")
    sys.exit(1)

# Get the CSV file paths from the command-line arguments
csv_files = sys.argv[1:]

# Ensure the CSV has the correct columns
# Assuming the columns are: "Time", "Case1", "Case2", "Case3"
time_column = "Time"
cases = ["Point_A", "Point_B", "Point_C"]

# Dictionary to hold data for each case
plot_data = {case: [] for case in cases}

benchmark = {}

# Read and process each CSV file
for csv_file in csv_files:
    try:
        # Load the CSV file
        data = pd.read_csv(csv_file)

        # Extract the file name without extension for labeling
        file_label = os.path.splitext(os.path.basename(csv_file))[0]

        # Append data for each case
        for case in cases:
            if "Benchmark" in file_label:
                benchmark[case] = data[case]
            plot_data[case].append((data[time_column], data[case], file_label))

    except FileNotFoundError:
        print(f"Error: File not found: {csv_file}")
    except pd.errors.EmptyDataError:
        print(f"Error: File is empty or invalid: {csv_file}")
    except KeyError as e:
        print(f"Error: Missing expected column in {csv_file}: {e}")


# Loop through the cases and plot each
for case in cases:
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot data from all files for the current case
    for time, values, label in plot_data[case]:
        linestyle = "-"
        if "Benchmark" in label:
            linestyle = "--"
        plt.plot(time, values, label=label, linewidth=2, linestyle=linestyle)


        norm = np.sum(np.square(benchmark[case] - values))
        norm = np.sqrt(norm)

        maxnorm = np.max(np.abs(values - benchmark[case]))
        print(f"Case: {case}, Label: {label}, norm: {norm}, maxnorm: {maxnorm}")

    ax.xaxis.set_major_locator(ticker.MultipleLocator(60))  # Major ticks every 2 units


    # Customize the plot
    plt.xlabel("Time (s)")
    plt.ylabel("Water height perturbation (m)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()

    plt.axhline(0, color='black', linestyle='--', linewidth=1.5)


    # Save the plot as a PNG file
    output_file = f"plots/{case}.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to save memory

    print(f"Plot saved: {output_file}")