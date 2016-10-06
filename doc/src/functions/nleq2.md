```@docs
NewtonMethods.deccon
```

```@docs
NewtonMethods.solcon
```

```@docs
NewtonMethods.n2prjn
```

```@docs
NewtonMethods.nleq2
```

Options
=======
| Name               | Default value   | Meaning     |
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
| OPT_IRANK          |  |  |
| OPT_COND           |  |  |

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
| STATS_SUBCOND   | Float64         |             |
| STATS_SENS      | Float64         |             |
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

```@docs
NewtonMethods.n2int
```
