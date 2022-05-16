# TSdisaggregation
Methods to perform temporal disaggregation and extrapolation.

MAIN FUNCTION: disaggregate.R

Input = 
- Y (n_l x 1 matrix)
- X = (n x p matrix)
- aggMat = 'sum' or 'avg' or 'first' or 'last'
- aggRatio = 4 (for annual-to-quarterly), 3 (for quarterly-to-monthly), etc...
- method = 'Denton'  or 'Denton-Cholette' or 'Chow-Lin' or 'Fernandez' or 'Litterman' or 'spTD' or 'adaptive-spTD'
- Denton = 'first' or 'second' or 'prop' or 'abs' (just for Denton or Denton-Cholette methods)
        
Output = 
- y_Est (n x 1 matrix)
- beta_Est (regression coefficient estimate)
- rho_Est (AR(1) parameter estimate)
- u_l (low-frequency residual)

Extrapolation: if n > aggRatio x n_l then values are extrapolated after aggRatio x n_l using parameter estimates obtained via the fit on the aggregated data (y_l,X_l).

SIMULATION FUNCTION: TempDisaggDGP.R

Simulate a temporal disaggregation type problem. 
