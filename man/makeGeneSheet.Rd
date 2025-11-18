\name{makeGeneSheet}
\alias{makeGeneSheet}
\title{Generate Excel or CSV Output from locusPackRat Project}
\description{
  Creates formatted Excel workbooks or CSV files from a locusPackRat project,
  combining core data with supplementary tables. Supports filtering, highlighting,
  and multi-tab Excel outputs with custom formatting.
}
\usage{
makeGeneSheet(format = c("excel", "csv"), output_file = NULL,
                  filter_expr = NULL, highlight_genes = NULL,
                  include_supplementary = TRUE, split_by = NULL,
                  split_criteria = NULL, project_dir = ".",
                  verbose = TRUE)
}
\arguments{
  \item{format}{Output format, either "excel" or "csv".}
  \item{output_file}{Path for the output file. If NULL, auto-generates based on
    project name and timestamp.}
  \item{filter_expr}{R expression as a string for filtering rows (e.g.,
    "expression > 100").}
  \item{highlight_genes}{Vector of gene symbols to highlight in Excel output
    (yellow background).}
  \item{include_supplementary}{Controls which supplementary tables to include:
    TRUE (all tables), FALSE (none), or character vector of table names.}
  \item{split_by}{For Excel, how to split data across tabs: NULL (single sheet),
    "criteria" (use split_criteria), or a column name.}
  \item{split_criteria}{Named list of filter expressions for creating multiple
    Excel tabs (e.g., list("High" = "value > 10", "Low" = "value <= 10")).}
  \item{project_dir}{Path to the project directory. Defaults to current directory.}
  \item{verbose}{If TRUE, prints progress messages. Default is TRUE.}
}
\value{
  Invisibly returns the path to the created output file. For Excel format,
  creates a formatted workbook with:
  \itemize{
    \item Auto-sized columns
    \item Frozen header row
    \item Bold headers with borders
    \item Yellow highlighting for specified genes
    \item Multiple tabs if split_by is specified
  }
}
\examples{
# Generate simple CSV output
makeGeneSheet(format = "csv", project_dir = "my_project")

# Generate Excel with all supplementary tables
makeGeneSheet(format = "excel", include_supplementary = TRUE,
              project_dir = "my_project")

# Generate Excel with specific tables only
makeGeneSheet(format = "excel",
              include_supplementary = c("expression", "phenotypes"),
                  project_dir = "my_project")

# Generate multi-tab Excel with filtering
makeGeneSheet(
  format = "excel",
  split_by = "criteria",
  split_criteria = list(
    "High_Expression" = "expression > 100",
    "Disease_Associated" = "!is.na(disease)",
    "Significant" = "p_value < 0.05"
  ),
  highlight_genes = c("BRCA1", "TP53"),
  project_dir = "my_project"
)
}
\seealso{
  \code{\link{initPackRat}}, \code{\link{addRatTable}},
  \code{\link{listPackRatTables}}
}
\author{
  Brian Gural, Todd Kimball, Anh Luu, Christoph D. Rau
}