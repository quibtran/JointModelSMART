# Joint_model_SMART

This repo consists of the R codes to simulate SMART data with survival outcomes, fit the joint model described in "Joint Modeling and Multiple Comparisons with
the Best of data from a SMART with Survival Outcomes" by Qui Tran, Alex Tsodikov, and Kelley M. Kidwell, and plot the predicted survival curves agaisnt other methods by Lunceford et al (2002), Guo & Tsiatis (2005), Tang & Wahed (2015) (available via R package 'DTR').  

In particular:
- simSMART.R: Function simSMART() to simulate 2 stage SMART design II with time-to-response and overall survival outcomes from 4 reigmens.
- calH.R: Function calH() to estimate the 2 baseline hazards (time-to-response and time-to-death) non-parametrically. 
- loglik.R: Function loglik() to calculate the log-likelihood for any set of beta coefficients.
- main.R: Main code files that source simSMART.R, calH.R, and loglik.R
- relative_efficiency.R: Contain procedure to calculate bootstraped standard error (SE) for survival estimates and calculate Relative Efficiency between the proposed joint model and previous methods.
- multiple_comparison_with_the_best.R: Procedure to pick set of best regimens (regimens that has highest survival rates at a particular time of interest) by conducting multiple comparison with the best.
