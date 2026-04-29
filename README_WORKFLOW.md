# Gray Whale Abundance Estimation Workflow Documentation

## 1. Three Main Stages

### Data Extraction

In this stage, data from the field are reviewed and corrections are made as needed. It is essential to ensure that the data is clean and formatted correctly for the subsequent analysis. Key functions to use: - `Extract_Data_All_v2.Rmd`: This R Markdown document serves as the starting point for data extraction. All data up to the 2025/2026 season have been extracted using this script.

Create a new output directory for the new season. This is done in lines 28 and 29: Line 28: out.dir \<- "RData/V2.1_Feb2026" Line 29: if (!dir.exists(out.dir)) dir.create(out.dir)

It's best to create a new directory, for example V2.1_Feb2026, for the analysis in 2026. This way, old data are preserved for reanalysis if necessary.

V2.1 refers to the version 2.1 for the data extraction method. This is the most recent and should provide necessary data for JAGS and WinBUGS analyses. Data from previous seasons should be copied to this new RData folder.

New data should be extracted by changing Line 35: YEAR \<- 2026 Note that the YEAR is defined as the second year of each season. For example, the 2025/2026 season should be YEAR \<- 2026 and so forth.

Various functions in "Granite_Canyon_Counts_fcns.R" are used for extracting data from the raw data files.

Error messages may return while running this script for the first time. They are mostly from input data file errors. Find which data file is causing the problem, double check where the problem(s) is by looking at the output of get.data on line 85 and get.shift on line 134. The file number (ff) and shift number (k on line 131) can be fixed and run the function line by line within the Granite_Canyon_Counts_fcns.R. The intermediate variables in those functions, e.g., data in get.data or shifts.df in get.shift, usually provide a clue to the error message(s).

Once the data extraction process is completed, make sure all the extracted data files are in the same directory (e.g., RData/V2.1_Feb2026), including those from the previous seasons. These files should have the filename convention of out_YYYY_minMM_Tomo_v2.rds, where YYYY corresponds to the season identification and MM indicates the minimum shift duration, which is defined at line 46:

min_dur \<- shift_dur_min - grace_min

Once this step is completed, the models can be fit to the data using WinBUGS or JAGS.

### Analysis

#### N-mixture model using WinBUGS

This is the analysis described in Durban et al. (2015). It is run in 'WinBUGS Ver2.Rmd'. Edit the lines 33-36 to make sure (1) the data.dir is pointing to the correct path (created in the previous step) and (2) the new season is included in the years object. min.dur can be changed but not necessary. It's been 60 minutes for a few years. For this to be changed, data need to be extracted again with the new minimum shift duration in the previous step.

WinBUGS should be downloaded from https://www.mrc-bsu.cam.ac.uk/software/bugs-project. Please note that WinBUGS is no longer updated. The executable should be placed and the path to the executable noted. The path needs to be specified in line 23:

Line 23: WinBUGS.dir <- paste0(Sys.getenv("HOME"), "/WinBUGS14")

WinBUGS code is 'GW_Nmix_Orig.bugs'. This file should be in the same directory as 'WinBUGS Ver2.Rmd'.

It is a good practice to start with short chains, e.g., n.iter = 200, to test if WinBUGS runs without issues. Furthermore, turn "debug = T" in the bugs function (line 89) so that WinBUGS will stop if it encounters a problem. Some errors that I encountered are listed in the script as well as in 'WinBUGS errors and fixes.txt'. When a short run completes without issues, change the MCMC setup to longer chains (e.g., n.iter = 100,000) and make "debug = F" and run the code again. This may run for many hours.

After WinBUGS completes the analysis, several plots are created to show the results within the .Rmd file. WinBUGS output is saved as a .rds file in the '/RData' directory. This output will be used in writing annual reports. The process for writing an annual report is described below.  

#### Hierarchical State-Space models using JAGS

This is the analysis described in Eguchi (2026). It is run in 'Jags_Richards_AllData_n_models.R'. This script should be edited when new data are added to the analysis. Line 28 defines analysis years, beyond the data in Laake's analysis and Line 29 defines the path to the data files (.rds files). 

max.day on line 30 refers to the assumed number of days in a migration season; 100 days include all years. MCMC setup can be modified in line 38. 

For this analysis, multiple models are fitted to the data. This script uses a function script 'Richards_HSSM_model_definition.R' that creates 8 different models according to the input to the function and runs JAGS to estimate parameters. Differences in the models are the likelihood function (Poisson vs. Negative Binomial) and which parameters of the Richards function are assumed constant over time or season-specific. See Eguchi (2026) for details. I explored other Other parameters to be season-specific but they did not converge well, likely due to the limited information in the data. If more data become available, other models may be developed and fitted.  

Results of the analysis are summarized and visualized using 'Jags_Richards_HSSM_Results.Rmd'. 


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

1.  Open the `Extract_Data_All_v2.Rmd` file.
2.  Source the required libraries and set paths.
3.  Run the code chunks sequentially to extract the data.

### JAGS Analysis

1.  Ensure the extracted data is in the correct format.
2.  Use the `data2WinBUGS_input` function to convert the data.
3.  Run the `Jags_Richards_AllData_n_models.R` script and monitor the progress.

### Reporting

1.  Compile the results in R Markdown.
2.  Generate visualizations using the outputs.
3.  Save and share the report with stakeholders.

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
