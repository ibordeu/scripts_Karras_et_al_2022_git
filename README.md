# scripts_Karras_et_al_Nature_2022_git
Scripts used in the data analysis in the manuscript entitled "A developmental cellular hierarchy in melanoma uncouples growth and metastatic phenotypes", Karras et al. (2022)

To run the scripts:
1) Make sure you have installed Matlab (some toolboxes, such as the Optimization tooldbox might be required) 
2) Download the files in the repository and put them all in the same folder
Note: script_simulation_SP_model.m uses as Input the results table produced by the script_biexponential_fit.m script.
However, this table is already provided, so you can run the scripts in any order.
3) To produce estimate de fitting parameters nbar, nbar_p and p0, and plot the corresponding cumulative distributions, run script_biexponential_fit.m.
Make sure that when you run the code, you Change Folder to the directory where you have stored the files.
The script script_biexponential_fit.m, loads the clone sizes in the clone_sizes.xlsx spreadsheet.
4) Run script_simulation_SP_model.m to perfomr stochastic simulations of the two-compartment model using the parameter in result_table.mat (it also loads the clone_sizes.xlsx spreadsheet).
This script produces plots showing the empirical CDF together with the CDF obtianed from the numerical simulations.
