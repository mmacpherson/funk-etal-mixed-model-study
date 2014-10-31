funk-etal-mixed-model-study
===========================

Supporting code and data for Funk, Macpherson, and Rakovski study of
phylogenetically-paired designs.

Files:

- **``evaluate_models.R``** -- R script that simulates
  phylogenetically-paired data. Fits and assesses for significance
  Model 1, Model 2, and Model 3 from the paper.

- **``run_simulations.sh``** -- Shell script to run
  ``evaluate_models.R`` over a grid of parameters.

- **``generate_funk_data_comparison_table.sh``** -- Evaluates Model 1,
  Model 2, and Model 3 for each of four traits from Funk & Throop
  (2010).

- **``Funk_2010_herbivory_study.csv``** -- Dataset from Funk & Throop
  (2010).
