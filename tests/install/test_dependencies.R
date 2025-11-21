# Force empty user libraries and make a temp library for installs
# This simulates a fresh R install with no packages installed
Sys.setenv(R_LIBS_USER = "")
dir.create("temp_lib", showWarnings = FALSE)
.libPaths(c(file.path(getwd(), "temp_lib"), .libPaths()[length(.libPaths())]))

# Use BiocManager to tell R where to look for 'plotgardener'
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# This adds the correct Bioconductor 3.x repos to your list
r <- getOption("repos")
r["CRAN"] <- "https://cloud.r-project.org"
r <- c(r, BiocManager::repositories())
options(repos = r)

cat("Installing dependencies from:\n")
print(getOption("repos"))

# Install
options(Ncpus = 4)
remotes::install_local(".", 
                       lib = "temp_lib", 
                       upgrade = "always", 
                       quiet = FALSE)
