Return codes
============
| Value           | Description |
|:----------------|:------------|
| 1     | Termination, since jacobian matrix became singular |
| 2     | Termination after nItmax iterations (as indicated by the option OPT_NITMAX) |
| 3     | Termination, since damping factor became to small |
| 4     | Warning: Superlinear or quadratic convergence slowed down near the solution. Iteration has been stopped therefore with an approximation of the solution not such accurate as requested by OPT_RTOL because possibly the OPT_RTOL requirement may be too stringent (i.e. the nonlinear problem is ill-conditioned) |
| 5     | Warning: Iteration stopped with termination criterion (using OPT_RTOL as requested precision) satisfied but no superlinear or quadratic convergence has been indicated yet. Therefore, possibly the error estimate for the solution may not match good enough the really achieved accuracy. |
| 21    | Non-positive value for OPT_RTOL supplied |
| 22    | Negative scaling value via vector xScal supplied |
| 30    | One or more fields specified in the options are invalid |
| 80    | Error signalled by linear solver routine nFact. The reason for failure is given by the STATS_IFAIL entry in statistics. The meaning of the value is given in nFact. |
| 81    | Error signalled by linear solver routine nSolv. The reason for failure is given by the STATS_IFAIL entry in statistics. The meaning of the value is given in nSolv. |
| 82    | Error signalled by user function fcn |
| 83    | Error signalled by user function provided in OPT_JACFCN |
