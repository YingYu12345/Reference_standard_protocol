#!/usr/bin/env Rscript
# config.R
# The script provides definitions of color and shape 

# Set CRAN mirror for the script
options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

# Function to install a package if it's not already installed
install_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Check if required R packages are installed and prompt the user to install if not
required_packages <- c("RColorBrewer")

# Function to check and install packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Package", pkg, "is not installed.\n")
      cat("Please run the following command to install it:\n")
      cat("install.packages(\"", pkg, "\")\n")
      # Install and load the required packages
      lapply(required_packages, install_if_not_installed)
      stop("Exiting script. Please install the required packages and re-run the script.\n", call.=FALSE)
    }
  }
}

# Run the function to check for required packages
check_and_install_packages(required_packages)

# Now load the required R packages
suppressPackageStartupMessages(library(RColorBrewer))

# Define the default colors (9)
colors <- RColorBrewer::brewer.pal(9, "Set1")
#e41a1c
#377eb8
#4daf4a
#984ea3
#ff7f00
#ffff33
#a65628
#f781bf
#999999

# Define the default shapes (26)
shapes <- c(0:25)
# 0 Square
# 1 Circle
# 2 Triangle pointing up
# 3 Plus
# 4 Cross
# 5 Diamond
# 6 Triangle Pointing Down
# 7 Square Cross
# 8 Star
# 9 Diamond Plus
# 10 Circle Plus
# 11 Triangles pointing up and down
# 12 Square Plus
# 13 Circle cross
# 14 Square and Triangle Down
# 15 Filled Squares
# 16 Filled Circle
# 17 Filled Triangle Pointing Up
# 18 Filled Diamond
# 19 Solid circle
# 20 Bullet (small circle)
# 21 Filled circle blue
# 22 Filled square blue
# 23 Filled Diamond Blue
# 24 Filled Triangle Dot Down Blue
# 25 Filled Triangle Dot Down Blue