```@docs
NewtonMethods.nleq1
```
Input parameters of NLEQ1
=========================
| Name       | Type                      | Description|
|:-----------|:--------------------------|:-----------|
| fcn        | Function                  | Function for which zero is to be found. Should be in the form of fcn(y,x) with y = f(x). |
| x[1:n]     | Vector&#123;Float64&#125; | Initial estimate of the solution.                                                        |
| xScal[1:n] | Vector&#123;Float64&#125; | User scaling (lower threshold) of the iteration vector x                                 |
| opt        | OptionsNLEQ               | Options for solving the nonlinear system. Valid options are listed below.                |

Options
=======
| Name               | Default value   | Meaning |
|:-------------------|:----------------|:------------|
| OPT_RTOL           | 1e-6 |  |
| OPT_QSUCC          |  |  |
| OPT_MODE           |  |  |
| OPT_JACGEN         |  |  |
| OPT_JACFCN         |  |  |
| OPT_MSTOR          |  |  |
| OPT_ML             |  |  |
| OPT_MU             |  |  |
| OPT_ISCAL          |  |  |
| OPT_PRINTWARNING   |  |  |
| OPT_PRINTITERATION |  |  |
| OPT_PRINTIOWARN    |  |  |
| OPT_PRINTIOMON     |  |  |
| OPT_PRINTIOSOL     |  |  |
| OPT_PRINTSOLUTION  |  |  |
| OPT_NONLIN         |  |  |
| OPT_QRANK1         |  |  |
| OPT_QORDI          |  |  |
| OPT_QSIMPL         |  |  |
| OPT_NOROWSCAL      |  |  |
| OPT_BOUNDEDDAMP    |  |  |
| OPT_IORMON         |  |  |
| OPT_NITMAX         |  |  |
| OPT_FCBAND         |  |  |
| OPT_SIGMA          |  |  |
| OPT_SIGMA2         |  |  |
| OPT_AJDEL          |  |  |
| OPT_AJMIN          |  |  |
| OPT_ETADIF         |  |  |
| OPT_ETAINI         |  |  |
| OPT_NBROY          |  |  |
| OPT_FCSTART        |  |  |
| OPT_FCMIN          |  |  |

Output parameters of NLEQ1
==========================
| Name    | Type                            | Description |
|:--------|:--------------------------------|:------------|
| x       | Vector&#123;Float64&#125;       | Solution values (or final values if exit before solution is reached). |
| stats   | Dict&#123;String,Any&#125; | A dictionary variable of additional output values. The fields are discussed below. |
| retCode | Int64                           | An integer value signifying the exit code. The meaning of the exit codes are discussed below. |

Statistics
==========
| Name            | Type            | Description |
|:----------------|:----------------|:------------|
| STATS_XSCAL     | Array{Float64}  |             |
| STATS_XITER     | Vector{Float64} |             |
| STATS_NATLEVEL  | Vector{Float64} |             |
| STATS_SIMLEVEL  | Vector{Float64} |             |
| STATS_STDLEVEL  | Vector{Float64} |             |
| STATS_PRECISION | Vector{Float64} |             |
| STATS_DAMPINGFC | Vector{Float64} |             |
| STATS_RTOL      | Float64         |             |
| STATS_CONV      | Float64         |             |
| STATS_SUMX      | Float64         |             |
| STATS_DLEVF     | Float64         |             |
| STATS_NITER     | Int64           |             |
| STATS_NCORR     | Int64           |             |
| STATS_NFCN      | Int64           |             |
| STATS_NFCNJ     | Int64           |             |
| STATS_NJAC      | Int64           |             |
| STATS_NREJR1    | Int64           |             |
| STATS_NEW       | Int64           |             |
| STATS_ICONV     | Int64           |             |

Return code
===========
