## Adapt synthetic difference in differences to staggered treatment timing settings
Code to incorporate staggered treatment adoption (based on appendix from Arkhangelsky et al. 2021) into synthdid package
Package can be installed as follows:
```
devtools::install_github(“synth-inference/synthdid”)
```
### This code is developed for estimating the ATT motivated by the idea from Callaway and Sant’Anna (2021).
The main idea of extending to the staggered treatment timing is to iterate through the SynthDiD algorithm on subsets of the data limited to the never-treated units and each treatment cohort individually, and then obtain a weighted average of each treatment cohort’s individual effects to estimate the aggregate group average treatment effect.
### Use property of influence functions for compute variance and covariance estimators.
Standard errors are computed in a more effecient way, based on the influence functions for each observation and treatment cohort. The overall summary parameter’s variance, V_theta, is calculated as follows: V_theta=w’Vw where V is calculated as described in [Kahn (2022)](https://j-kahn.com/files/influencefunctions.pdf).
### Real-world data application to demonstrate the practical use of staggered synthetic DiD estimation
Winnow_waste_data.csv is used for demonstration where 1,200 commercial kitchens adopted AI-powered waste watcher at different times between Jan 2014 and Aug 2021.
Author: Yu Nu (yn292@cornell.edu)
