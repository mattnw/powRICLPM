# powRICLPM_summary() works

    Code
      cat(summary(output, parameter = "wB2~wA1"))
    Output
      
      powRICLPM analysis completed:
       - Monte Carlo replications:  5
       - Sample size(s): 500 550 600
       - Number of time points: 3
       - Proportion(s) random intercept variance: 0.5
       - Target power: 0.5
      
      Number of conditions that reach target power: 2
      
      Suggested next steps:
       - If this is a preliminary powRICLPM analysis, validate these recommendations (or a selection) by rerunning the analysis with an increased number of replications (e.g., `reps = 1000`).
       - Use `plot_powRICLPM()` to visualize results across all experimental conditions.

