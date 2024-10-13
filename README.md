# Mapping CAFs subtypes convergence using public available CRC spatial transcriptomics

This repository provides the code and detailed instructions required to reproduce the RMarkdown report for the collaborative research between OSR TIGET and IRCC Candiolo.

## Repository Structure

- `/script/`: Contains the R script responsible for data processing, statistical analysis, and result generation. Also contains the RMarkdown file used to produce the final report, incorporating code, visualizations, and explanatory content. 
- `/data/`: Currently with two small files. Users must download the required data files separately (see instructions below).
- `/report/`: Quick access to the PDF version of the Report.

## Instructions

### 1. Download Data

The necessary data files are too large to be hosted directly on GitHub. Please download them from Zenodo using the following link:

[Download Data from Zenodo](https://zenodo.org/records/13922930)

Once downloaded, place the data files in the `/data/` directory of this repository. Do not move any other file in the same directory.

### 2. Run the Analysis

To reproduce the report, follow these steps:

1. **Install Necessary Packages**: Ensure all required R packages are installed. You can do this by executing the following command in R:

   ```r
   install.packages(c("Matrix", "data.table", "Seurat", "harmony", "dplyr", "tidyr", "ggplot2","reshape2", "lisi", "reticulate", "circlize", "Giotto", "ComplexHeatmap", "SeuratObject"))
   ```

2. **Run the Analysis Script**: Execute the analysis script located in the `/scripts/` directory to preprocess the data and generate the results:

   ```r
   source("script/Report_OSR_Candiolo_2024.R")
   ```

   Note that during this phase, the script will generate under the same folder various .Robj files, each one is a ggplot2 figure used in the report. Some of them might be large. 

3. **Render the Report**: Knit the RMarkdown file to produce the final report in your desired format (e.g., HTML, PDF). This can be done within RStudio by opening the `.Rmd` file and clicking the "Knit" button, or by running:

   ```r
   rmarkdown::render("report.Rmd")
   ```

## Citation

If you utilize this repository or the associated data in your research, please cite this collaborative work between OSR TIGET and IRCC Candiolo. Additionally, cite the Zenodo data record linked above.

## Contact

For questions or issues, please contact:

- **Carlo Leonardi (IRCC Candiolo)**
- **Email**: carlo.leonardi\@unito.it (PhD) / carlo.leonardi\@ircc.it (Hospital) / c.leonardi95\@protonmail.com (Personal)

If you encounter any problems or have suggestions for improvement, please open an issue on this GitHub repository.

## License

This repository is made available under the [MIT License](LICENSE). This license allows for free use, modification, and distribution of the code.
