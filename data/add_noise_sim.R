# Load necessary library
library(data.table)

# Define file paths
input_path <- "height_sim.phen"
output_path <- "height_sim2.phen"

# Load the data
data <- fread(input_path, header = FALSE)

# Ensure the file has at least three columns
if (ncol(data) < 3) {
  stop("The input file must have at least three columns.")
}

# Extract necessary columns
data_out <- data[, .(V1 = 0, V2 = V2, V3 = V3 + rnorm(.N, mean = 0, sd = 1))]

# Save the modified dataset (tab-separated, no header)
fwrite(data_out, output_path, sep = "\t", col.names = FALSE, quote = FALSE)

cat("Modified file saved to:", output_path, "\n")