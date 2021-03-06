# Computationally-Efficient-Approximations-for-Distributionally-Robust-Optimization

- The available data, results, and MATLAB codes in this project are the exact data, results, and MATLAB codes used in the paper named **"Computationally Efficient Approximations for
Distributionally Robust Optimization under Moment and Wasserstein Ambiguity"** submitted to **"INFORMS Journal on Computing"**.

- The data and results corresponding to each table stored in an Excel file indicated by the table number in the paper. The folder "Problem data and results" includes the excel files of all tables.

- We summarize the lower bounds of the DRO problem with the moment-based ambiguity set on both applications in Tables 1 and 2, while Tables 3 and 4 report the results for the DRO problem with the combined ambiguity set. We report performance of the upper bound (17) in Tables 5 - 6 and report that of (20) and (46) in Tables 7 - 10. The interval results of the DRO problem with the moment-based ambiguity set are summarized in Tables 11 and 12, while Tables 13 and 14 report the results for the DRO problem with the
combined ambiguity set. We show the results related to the benefits of choosing leading components in Table 15. Tables 16  - 18 report the results of sensitivity analysis with respect to γ1, γ2, and R0. We report the performance of lower bound (EC.46) and upper bound (EC.48) of Problem (DRO-C2) in Tables EC.1 and EC.2, respectively. The short description of each table is as follows:
  - Table 1: Lower bound (8) on the production-transportation problem
  - Table 2: Lower bound (8) on the newsvendor problem
  - Table 3: Lower bound (43) on the production-transportation problem
  - Table 4: Lower bound (43) on the newsvendor problem
  - Table 5: Upper bound (17) on the production-transportation problem
  - Table 6: Upper bound (17) on the newsvendor problem
  - Table 7: Upper bound (20) on the production-transportation problem
  - Table 8: Upper bound (20) on the newsvendor problem
  - Table 9: Upper bound (46) on the production-transportation problem
  - Table 10: Upper bound (46) on the newsvendor problem
  - Table 11: Intervals on the production-transportation problem
  - Table 12: Intervals on the newsvendor problem
  - Table 13: Intervals on the production-transportation problem
  - Table 14: Intervals on the newsvendor problem
  - Table 15: Arbitrary vs. leading components on the newsvendor problem
  - Table 16: Sensitivity analysis for lower bound (8) on the newsvendor problem with respect to γ1
  - Table 17: Sensitivity analysis for lower bound (8) on the newsvendor problem with respect to γ2
  - Table 18: Sensitivity analysis for lower bound (43) on the newsvendor problem with respect to R0
  - Table EC.1: Lower bound (EC.46) on the production-transportation problem
  - Table EC.2: Upper bound (EC.48) on the production-transportation problem
  - Table EC.3 Lower bounds (43) and (EC.46) on the production-transportation problem
  - Table EC.4 Upper bounds (46) and (EC.48) on the production-transportation problem
  - Table EC.5 Sensitivity analysis for lower bound (EC.46) with respect to (γ1;R0)
  - Table EC.6 Sensitivity analysis for upper bound (EC.46) with respect to (γ1;R0)

- The folder "MATLAB Codes" includes two folders named "Multi-Product Newsvendor_MATLAB Codes" and "Risk-averse Production-Transportation_MATLAB Codes". The folder "Multi-Product Newsvendor_MATLAB Codes" includes all the MATLAB codes related to the newsvendor problem and the folder  "Risk-averse Production-Transportation_MATLAB Codes" includes all the MATLAB codes related to the production-transportation problem. 
- The folder of each problem has two folders named Experiments and Functions. Each MATLAB file in the folder Experiments is related to at least one experiment (Table) of the paper so that it calls a few MATLAB files from the folder Functions. 
- In the beginning of each MATLAB file of folder Experiments, the datasets of the related experiment (Table) are randomly generated. By running the MATLAB file, those randomly generated datasets appeare in the Workspace of MATLAB as ".mat" files.
- Each MATLAB file has been documented comprehensively. Therefore, one can easily understand different parts of a MATLAB code by reading the comments provided in the MATLAB file. 
- In addition to documentations provided in the MATLAB files, one can use the following guidlines to understand the content of each MATLAB file from its name:
  - UB = Upper Bound
  - LB = Lower Bound
  - moment = Moment-based Ambiguity Set
  - Combined = Combined Ambiguity Set
  - NewsVendor = The newsvendor problem
  - Production = The production-transportation problem
  - Parameter name (e.g., Gamma1, Gamma2) = The experiment related to the sensitivity analysis with respect to that parameter
  - Unsolvable = The experiment related to the interval performance 
  - Arbitrary = The experiment related to the benefits of choosing leading components
- The documentations inside each MATLAB file coupled with the above guidlines about the file names help the interested reader locate the MATLAB file of an experiment or Table  


