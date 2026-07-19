# LiFePO4 Battery SOC Estimation using VFFRLS-MIUKF

Implementation of a hybrid State of Charge (SOC) estimation method for LiFePO4 batteries. Combines Variable Forgetting Factor Recursive Least Squares (VFFRLS) for online parameter identification of a 2nd-order RC Thevenin model with Multi-Innovation Unscented Kalman Filter (MIUKF) for state estimation. Based on a reference paper from the SKPA (Control System and Modeling) course final project.

## Files

- **`main.m`** — Main script. Loads battery test data, computes true SOC via Coulomb counting for reference, then runs the VFFRLS-MIUKF estimation loop. Generates all plots comparing estimated vs. measured terminal voltage, SOC, OCV, and electrical parameters.

- **`VFFRLS.m`** — Variable Forgetting Factor Recursive Least Squares. Identifies discrete-time model parameters (a1, a2, b0, b1, b2) of the battery's transfer function. The forgetting factor λ adapts between 0.90–0.995 based on the MIUKF estimation error — higher error increases forgetting (more weight on new data), lower error stabilizes parameters.

- **`Levenberg_Marquardt.m`** — Converts the discrete-time parameters from VFFRLS into physical electrical parameters (R0, R1, C1, R2, C2) by solving a nonlinear system of equations using MATLAB's `lsqnonlin` with the Levenberg-Marquardt algorithm.

- **`MIUKF.m`** — Multi-Innovation Unscented Kalman Filter. Estimates the battery states (U1, U2, SOC) using sigma-point propagation through the battery model. Uses multi-innovation theory (M=22) — instead of updating with only the latest innovation, it accumulates a window of past innovations weighted by a decay factor, improving convergence over standard UKF.

- **`battery_state_model.m`** — State transition function. Propagates the state (U1, U2, SOC) forward in time using Forward Euler discretization of the 2nd-order RC circuit equations.

- **`battery_output_model.m`** — Measurement function. Computes terminal voltage from the estimated state and parameters: Vt = OCV - U1 - U2 - I × R0.

## Data

- `A1-007-DST-US06-FUDS-40-20120822.xlsx` — Standard battery test cycles (DST, US06, FUDS) from the CALCE Battery Research Group, University of Maryland (<https://calce.umd.edu/battery-data#A123>). A123 LiFePO4 cell data.
- `battery_data.mat` — Preprocessed time, terminal voltage, and current.

## Requirements

MATLAB with Optimization Toolbox (required for `lsqnonlin` in Levenberg-Marquardt)

## Note

Student implementation of a reference paper from the SKPA course final project. Results are preliminary and not validated with quantitative error metrics.
