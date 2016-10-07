Statistics
==========
| Name            | Description |
|:----------------|:------------|
| STATS_XSCAL     | After return with retCode >= 0, it contains the latest internal scaling vector used. After return with retCode == -1 in onestep mode it contains a possibly adapted (as described below) user scaling vector: |
|                 | If (xScal[i] <  small) xScal[i] = small <br>If (xScal[i] >  great) xScal[i] = great |
| STATS_RTOL      | Finally achieved (relative) accuracy. The estimated absolute error of component i of x_out is approximately given by abs_err[i] = rTol * xScal_out[i], where (approximately) xScal_out[i] = max(abs(x_out[i]),xScal_in[i]).            |
| STATS_CONV      | The achieved relative accuracy after the latest step |
| STATS_SUMX      | Natural level (not Normx of printouts) of the latest iterate, i.e. sum(dx.^2), where dx = scaled Newton correction. |
| STATS_DLEVF     | Standard level (not Normf of printouts) of the latest iterate, i.e. norm(f(x),2), where f =  nonlinear problem function. |
| STATS_SUBCOND   | Subcondition of the linear system as estimated by the linear solver (only in case of nleq2) |
| STATS_SENS      | Sensitivity of the linear system as estimated by the linear solver (only in case of nleq2) |
| STATS_NITER     | Number of Newton-iterations |
| STATS_NCORR     | Number of corrector steps |
| STATS_NFCN      | Number of fcn-evaluations |
| STATS_NFCNJ     | Number of fcn-evaluations for Jacobian approximation |
| STATS_NJAC      | Number of Jacobian generations or jac-calls |
| STATS_NREJR1    | Number of rejected Newton iteration steps done with a rank-1 approximated Jacobian |
| STATS_NEW       | Count of consecutive rank-1 updates |
| STATS_ICONV     | Current status of of the convergence monitor (only if convergence order monitor is on see OPT_IORMON)<br>=0: No convergence indicated yet<br>=1: Damping factor is 1.0<br>=2: Superlinear convergence in progress<br>=3: Quadratic convergence in progress |
| STATS_IFAIL     | Failure code to be checked in case the return code is 80/81. |
| STATS_PRECISION | The sequence of acheived precisions over the iteration steps. |

Statistics stored only if OPT_STORE = 1
========================================
| Name            | Description |
|:----------------|:------------|
| STATS_XITER     | An array holding all iterates of the Newton iteration run |
| STATS_NATLEVEL  | The sequence of natural levels of the Newton corrections over the iteration steps |
| STATS_SIMLEVEL  | The sequence of natural levels of the simplified Newton corrections over the iteration steps |
| STATS_STDLEVEL  | The sequence of standard levels over the iteration steps |
| STATS_DAMPINGFC | The sequence of accepted damping factors over the iteration steps. |
