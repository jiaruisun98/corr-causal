# corr-causal
This project presents code used for the article `Unveiling the Unobservable: Causal Inference on Multiple Derived Outcomes'.

There are 8 documents.

1.Functions.R is the code of the functions for simulation, including the function for data generating, function for proposed method, function for BH procedure.

2.Simulation.R is the code for simulation of the proposed method on correlation. In Simulation.R, `Simulation_block_proposed` is the result of proposed method under block diagonal setting, `Simulation_block_BH` is the result of BH method under block diagonal setting, `Simulation_off_proposed` is the result of proposed method under off diagonal setting, `Simulation_off_BH` is the result of BH method under off diagonal setting. The simulation results are further organized in `simulation-result.xlsx` to produce Figure 1 and Figure 2.

3.Simulation_regression.R is the code for simulation of the proposed method on regression parameters of lasso and debiased-lasso. In Simulation_regression.R,`result` provides the analysis result used to construct Figure S1. The simulation results are further organized in `simulation-result.xlsx` to produce Figure S1.

4.Simulation-result.xlsx is our simulation results which are presented in the article.

5.plot.R uses simulation results in Simulation.R and Simulation_regression.R to produce Figure 1, 2 and S1. In plot.R, `plot_block` produces Figure 1, `plot_off` produces Figure 2, and `plot_lasso` produces Figure S1.

6.RealCaseStudy.R is the code for application. This code produces table 1 and the brain connection in Figure S2. In RealCaseStudy.R, `Connection_diff1` provides the Nerwork in `Table1` and `Figure S2`, `Tau_sig` provides the `Estimated Effect` in `Table1`, `CI_low` and `CI_high` provides the 95% CI in `Table1`, `Ave-trt` and `Ace-cl` provide Ave-trt and Ace-cl in `Table1`.

7.data_dictionary.txt tells the meaning of the data used in this article.

8.Process to request Data.txt shows the process to request the data.
