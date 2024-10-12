# Report on the Collaboration between OSR TIGET and IRCC Candiolo

This repository provides the code and detailed instructions required to reproduce the RMarkdown report for the collaborative research between OSR TIGET and IRCC Candiolo.

## Repository Structure

- `/script/`: Contains the R script responsible for data processing, statistical analysis, and result generation.
- `/data/`: Currently empty. Users must download the required data files separately (see instructions below).
- `.Rmd`: The RMarkdown file used to produce the final report, incorporating code, visualizations, and explanatory content.

## Instructions

### 1. Download Data

The necessary data files are too large to be hosted directly on GitHub. Please download them from Zenodo using the following link:

[Download Data from Zenodo](https://zenodo.org/records/13922930)

Once downloaded, place the data files in the `/data/` directory of this repository.

### 2. Run the Analysis

To reproduce the report, follow these steps:

1. **Install Necessary Packages**: Ensure all required R packages are installed. You can do this by executing the following command in R:

   ```r
   install.packages(c("ggplot2", "dplyr", "readr", "rmarkdown"))
   ```

2. **Run the Analysis Script**: Execute the analysis script located in the `/scripts/` directory to preprocess the data and generate the results:

   ```r
   source("script/analysis.R")
   ```

   Note that during this phase, the script will generate under the same folder various .Robj files, each one is a ggplot2 figure used in the report.

3. **Render the Report**: Knit the RMarkdown file to produce the final report in your desired format (e.g., HTML, PDF). This can be done within RStudio by opening the `.Rmd` file and clicking the "Knit" button, or by running:

   ```r
   rmarkdown::render("report.Rmd")
   ```

## Dependencies

This project requires several R packages, including but not limited to:

- **ggplot2**: Used for data visualization.
- **dplyr**: Utilized for data manipulation.
- **readr**: For reading data files.
- **rmarkdown**: For rendering the final report.

Please ensure these packages are installed prior to executing the analysis.

## Citation

If you utilize this repository or the associated data in your research, please cite the collaborative work between OSR TIGET and IRCC Candiolo. Additionally, cite the Zenodo data record linked above.

## Contact

For questions or issues, please contact:

- **Carlo Leonardi (IRCC Candiolo)**
- **Email**: carlo.leonardi\@unito.it (PhD) / carlo.leonardi\@ircc.it (Hospital) / c.leonardi95\@protonmail.com (Personal)

If you encounter any problems or have suggestions for improvement, please open an issue on this GitHub repository.

## License

This repository is made available under the [MIT License](LICENSE). This license allows for free use, modification, and distribution of the code.
