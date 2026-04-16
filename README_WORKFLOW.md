# Gray Whale Abundance Estimation Workflow Documentation

## 1. Three Main Stages
### Data Extraction
In this stage, data related to gray whale sightings and related metrics are collected and preprocessed. It is essential to ensure that the data is clean and formatted correctly for the subsequent analysis. Key functions to use:
- `Extract_Data_All_v2.Rmd`: This R Markdown document serves as the starting point for data extraction, allowing for covariate selection and insight generation based on preliminary data exploration.

### JAGS Analysis
This involves the statistical analysis phase where the data is modeled using JAGS (Just Another Gibbs Sampler). The primary scripts involved in this stage include:
- `data2WinBUGS_input`: This function prepares the data output for JAGS, converting it into the required format for the WinBUGS/JAGS models.
- `Jags_Richards_AllData_n_models.R`: This script runs the JAGS models for estimating the gray whale abundance based on the input data from the previous stage.

### Reporting
Once the analysis is complete, the results will be reported. This stage includes visualizing the outputs and generating reports summarizing the findings. This is crucial for communicating the results to stakeholders effectively.

## 2. Essential Files List
- `Extract_Data_All_v2.Rmd`
- `data2WinBUGS_input`
- `Jags_Richards_AllData_n_models.R`
- `Granite_Canyon_Counts_fcns.R`
- Resulting reports and visualizations

## 3. Step-by-Step Processing Instructions
### Data Extraction
1. Open the `Extract_Data_All_v2.Rmd` file.
2. Source the required libraries and set paths.
3. Run the code chunks sequentially to extract the data.

### JAGS Analysis
1. Ensure the extracted data is in the correct format.
2. Use the `data2WinBUGS_input` function to convert the data.
3. Run the `Jags_Richards_AllData_n_models.R` script and monitor the progress.

### Reporting
1. Compile the results in R Markdown.
2. Generate visualizations using the outputs.
3. Save and share the report with stakeholders.

## 4. Input/Output Specifications
- **Input Data:** Preprocessed data frame from `Extract_Data_All_v2.Rmd`.
- **Output Data:** Summary statistics, JAGS model outputs, and visualizations.

## 5. MCMC Parameter Guidance
- Ensure appropriate burn-in and thinning for MCMC chains.
- Recommended parameters:
  - Number of iterations: 10000
  - Thinning: 10
  - Chains: 3

## 6. Troubleshooting Tips
- If models fail to converge, check input data for errors.
- Ensure that prior distributions are set correctly.
- For any errors in RMarkdown, refer to the console output for clues on where the issue lies.
