# PredSL: Stochastic Cellular Automata for Snow and Glacial Lake Simulation

This repository contains MATLAB code for simulating and forecasting changes in Snow Cover Area (SCA) and Glacial Lake Area (GLA) using a stochastic Cellular Automata-based model named **PredSL**. The model integrates environmental inputs (elevation, surface air temperature, and neighborhood states) with Latin Hypercube Sampling (LHS)-based parameter tuning to optimize transition functions for glacial hazard prediction.

## 🌐 Overview

**PredSL** simulates glacial evolution in mountainous regions by:
- Evaluating environmental influences using a governing function.
- Iteratively minimizing prediction error via MSE.
- Applying thresholding and neighborhood-based convolution to predict snow and lake transitions.
- Comparing results to high-resolution satellite imagery for validation.

## 📂 Structure

- `SCMain_proxy.m` – Takes inputs, normalizes as per requirement.
- `Stochastic_proxy.m` – Starts executing the simulation workflow.
- `generateRaster.m` – Applies thresholding rules to generate discrete class rasters.
- `convolution.m` – Applies spatial refinement through 3×3 box filtering.
- `thresholding2.m` – Computes MSE between simulated and reference rasters during threshold tuning.
- `data/` – Input rasters and ground truth files (not included in this repo).
- `output/` – Final simulated outputs and performance evaluation.
- `matlab_adaptive_threshold_LHS_8.mat` - contains all the necessary computations in the workspace 

## 🧮 Classes and Labels

- **0** – Background (black)
- **1** – Snow (yellow)
- **2** – Glacial Lake (blue)

## 🔁 Workflow Summary

1. **Initialize Parameters** via `xin = [rho, alpha, beta, gamma, pp, qq, rr]`.
2. **Compute Governing Function** using elevation, temperature, and neighborhood state.
3. **Tune Thresholds** over multiple iterations for best MSE using LHS.
4. **Simulate Transitions** using `generateRaster`.
5. **Apply Convolution** for spatial refinement.
6. **Evaluate Accuracy** using MSE across time-steps.
7. **Save Results** to Excel log for analysis.

## 📊 Outputs

- Best snow and lake thresholds.
- MSE before and after convolution.
- Simulated rasters across timestamps.

## 📦 Dependencies

- MATLAB R2021a or newer
- Image Processing Toolbox (for `immse`)
- Your local dataset (not included due to size constraints)

## 📁 Data Requirements

- `data(:,:,t)` – Ground truth 3-class rasters (snow, lake, background)
- `elevation`, `incidence`, `means` – Environmental raster data and stats
- `mask`, `mask2`, `nbors`, `cells` – Preprocessed logical masks and grid structures

> 📝 Example datasets should be placed in the `data/` folder. You can modify the script to point to your local raster files.

## 📈 Sample Visualization

Visualization color scheme:
- Yellow: Snow
- Blue: Lake
- Black: Background
- (Optional: Red overlay for change detection)

## The Dataset can be accessed via the link: https://zenodo.org/records/15337740
** Start the code from SCMain_proxy.m, where the data is loaded and normalized as required. After completing the execution, the Stochastic_proxy.m computes the necessary values, performs all the operations and generates the result as per the requirement.
