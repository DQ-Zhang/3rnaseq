#!/usr/bin/env Rscript
# RNA_analyzer_installer.R
# Install all required packages for 3RNA_analyzer_260104.R

# Function to check and install packages
install_if_missing <- function(packages, type = "CRAN") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      
      if (type == "CRAN") {
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      } else if (type == "BioC") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org/")
        }
        BiocManager::install(pkg, ask = FALSE)
      }
      
      # Verify installation
      if (requireNamespace(pkg, quietly = TRUE)) {
        message(paste("✓", pkg, "installed successfully"))
      } else {
        warning(paste("✗ Failed to install", pkg))
      }
    } else {
      message(paste("✓", pkg, "already installed"))
    }
  }
}

# Function to install species-specific OrgDb packages
install_orgdb_packages <- function() {
  orgdb_packages <- c(
    "org.Hs.eg.db",  # Human
    "org.Mm.eg.db",  # Mouse
    "org.Rn.eg.db"   # Rat
  )
  
  message("\nInstalling species-specific annotation databases...")
  
  # Install BiocManager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  
  # Install each OrgDb package
  for (orgdb in orgdb_packages) {
    if (!requireNamespace(orgdb, quietly = TRUE)) {
      message(paste("Installing", orgdb, "..."))
      BiocManager::install(orgdb, ask = FALSE)
      
      if (requireNamespace(orgdb, quietly = TRUE)) {
        message(paste("✓", orgdb, "installed successfully"))
      } else {
        warning(paste("✗ Failed to install", orgdb))
      }
    } else {
      message(paste("✓", orgdb, "already installed"))
    }
  }
}

# Function to install optional packages (not strictly required but used in script)
install_optional_packages <- function() {
  optional_packages <- c(
    "ggridges",      # For ridge plots in GSVA section
    "ggrepel"        # For label repulsion in PCA plots
  )
  
  message("\nInstalling optional packages...")
  
  for (pkg in optional_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
      message(paste("✓", pkg, "installed"))
    } else {
      message(paste("✓", pkg, "already installed"))
    }
  }
}

# Main installation routine
main <- function() {
  message("=========================================")
  message("RNA Analyzer Package Installer")
  message("=========================================")
  
  # CRAN packages
  cran_packages <- c(
    "dplyr",
    "data.table",
    "ggplot2",
    "pheatmap",
    "corrplot",
    "tidyr",
    "stringr",
    "RColorBrewer",
    "gridExtra",
    "grid",
    "svglite",
    "reshape2",
    "tibble",
    "ggsci"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "DESeq2",
    "limma",
    "clusterProfiler",
    "enrichplot",
    "DOSE",
    "GSVA",
    "GO.db",
    "AnnotationDbi"
  )
  
  message("\n[1/3] Installing CRAN packages...")
  install_if_missing(cran_packages, type = "CRAN")
  
  message("\n[2/3] Installing Bioconductor packages...")
  install_if_missing(bioc_packages, type = "BioC")
  
  # Install OrgDb packages
  install_orgdb_packages()
  
  # Install optional packages
  install_optional_packages()
  
  # Check for EnhancedVolcano (special case - on Bioconductor)
  message("\nChecking for EnhancedVolcano...")
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    message("Installing EnhancedVolcano...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org/")
    }
    BiocManager::install("EnhancedVolcano", ask = FALSE)
  } else {
    message("✓ EnhancedVolcano already installed")
  }
  
  # Final check
  message("\n=========================================")
  message("Verifying installations...")
  message("=========================================")
  
  all_packages <- c(cran_packages, bioc_packages, 
                    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
                    "EnhancedVolcano", "ggridges", "ggrepel")
  
  missing <- c()
  for (pkg in all_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(missing) == 0) {
    message("✓ All packages installed successfully!")
  } else {
    message("✗ Some packages failed to install:")
    for (pkg in missing) {
      message(paste("  -", pkg))
    }
    message("\nYou may need to install these manually:")
    message("For CRAN packages: install.packages('package_name')")
    message("For Bioconductor packages: BiocManager::install('package_name')")
  }
  
  message("\n=========================================")
  message("Installation complete!")
  message("Note: Local databases (MsigDb, collectri) need to be")
  message("placed in 'localdata/' directory separately.")
  message("=========================================")
}

# Run installation
if (interactive()) {
  main()
} else {
  # If run as Rscript
  main()
}