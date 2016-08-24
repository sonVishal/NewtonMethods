```@docs
NewtonMethods.nleq1
```
Input parameters of NLEQ1
=========================
| Name       | Type                      | Description|
|:----------:|:-------------------------:|:-----------|
| fcn        | Function                  | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]     | Vector&#123;Float64&#125; | Initial estimate of the solution.                                                        |
| xScal[1:n] | Vector&#123;Float64&#125; | User scaling (lower threshold) of the iteration vector x                                 |
| opt        | OptionsNLEQ               | Options for solving the nonlinear system. Valid options are listed below.                |

Output parameters of NLEQ1
==========================
| Name    | Type                            | Description |
|:-------:|:-------------------------------:|:------------|
| x       | Vector&#123;Float64&#125;       | Solution values (or final values if exit before solution is reached). |
| stats   | Dict&#123;ASCIIString,Any&#125; | A dictionary variable of additional output values. The fields are discussed below. |
| retCode | Int64                           | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |

Statistics
==========
| Name             | Type            | Description |
|:----------------:|:---------------:|:------------|
| STATS_NITER     | Int64 |             |
| STATS_NCORR     | Int64 |             |
| STATS_NFCN      |       |             |
| STATS_NFCNJ     |       |             |
| STATS_NJAC      |       |             |
| STATS_NREJR1    |       |             |
| STATS_XSCAL     |       |             |
| STATS_RTOL      |       |             |
| STATS_XITER     |       |             |
| STATS_NATLEVEL  |       |             |
| STATS_SIMLEVEL  |       |             |
| STATS_STDLEVEL  |       |             |
| STATS_PRECISION |       |             |
| STATS_DAMPINGFC |       |             |
| STATS_NEW       |       |             |
| STATS_ICONV     |       |             |
| STATS_CONV      |       |             |
| STATS_SUMX      |       |             |
| STATS_DLEVF     |       |             |
