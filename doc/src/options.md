Common options
==============
| Name               | Default value | Value | Meaning |
|:-------------------|:--------------|:------|:--------|
| OPT_RTOL           | 1e-6          |   | Required relative precision of solution components.                      |
| OPT_QSUCC          | 0             | 0 | (false) initial call: NLEQ solver is not yet initialized, i.e. this is the first call for this nonlinear system. At successful return with OPT_MODE=1, OPT_QSUCC is set to 1. |
|                    |               | 1 | (true) successive call: NLEQ solver is initialized already and is now called to perform one or more following Newton iteration steps. |
| OPT_MODE           | 0             | 0 | Standard mode initial call: Return when the required accuracy for the iteration vector is reached. User defined parameters are evaluated and checked. |
|                    |               | 0 | Standard mode successive call: If NLEQ solver was called previously with OPT_MODE=1, it performs all remaining iteration steps. |
|                    |               | 1 | Stepwise mode: Return after one Newton iteration step. |
| OPT_JACGEN         | 2             | 1 | User supplied subroutine OPT_JACFCN will be called to generate Jacobian matrix |
|                    |               | 2 | Jacobian approximation by numerical differentation (see subroutines nJac and nJacb) |
|                    |               | 3 | Jacobian approximation by numerical differentation with feedback control (see subroutines nJcf and nJcfb) |
|                    |               | 4 | Jacobian approximation by forward mode automatic differentiation using the ForwardDiff package |
| OPT_JACFCN         | 0             |   | User supplied Jacobian generation function of the form Jac(a,x), where a = Jacobian(x). Only required if OPT_JACGEN = 2 |
| OPT_MSTOR          | 0             | 0 | Jacobian is a dense matrix  |
|                    |               | 1 | Jacobian is a band matrix   |
| OPT_ML             | 0             |   | Lower bandwidth of the Jacobian (excluding the diagonal). Only required if OPT_MSTOR = 1  |
| OPT_MU             | 0             |   | Upper bandwidth of the Jacobian (excluding the diagonal). Only required if OPT_MSTOR = 1  |
| OPT_ISCAL          | 0             | 0 | The user supplied scaling vector XSCAL is used as a (componentwise) lower threshold of the current scaling vector  |
|                    |               | 1 | The vector xScal is always used as the current scaling vector |
| OPT_PRINTWARNING   | 0             | 0 | Do not print any warning or error messages |
|                    |               | 1 | Print warning and error messages |
| OPT_PRINTITERATION | 0             | 0 | Do not print iteration monitor  |
|                    |               | 1 | Standard output |
|                    |               | 2 | Summary iteration monitor additionally |
|                    |               | 3 | Detailed iteration monitor additionally |
|                    |               | 4,5,6 | Outputs with increasing level additional increasing information for code testing purposes. Level 6 produces in general extremely large output! |
| OPT_PRINTSOLUTION  | 0             | 0 | Do not print solutions  |
|                    |               | 1 | Print initial values and solution values |
|                    |               | 2 | Print intermediate values additionally |
| OPT_PRINTIOWARN    | 0             | 0 | Prints the warning and error messages to STDOUT  |
|                    |               | IOStream | Prints the warning and error messages to the file provided. Please look at the examples to understand this option. |
| OPT_PRINTIOMON     | 0             | 0 | Prints the iteration monitor to STDOUT  |
|                    |               | IOStream | Prints the iteration monitor to the file provided. Please look at the examples to understand this option. |
| OPT_PRINTIOSOL     | 0             | 0 | Prints the solutions to STDOUT  |
|                    |               | IOStream | Prints the solutions to the file provided. Please look at the examples to understand this option. |
| OPT_NONLIN         | 3             | 1 | Linear problem.<br>*Warning: If specified, no check will be done, if the problem is really linear, and the NLEQ solver terminates unconditionally after one Newton-iteration step.*  |
|                    |               | 2 | Mildly nonlinear problem |
|                    |               | 3 | Highly nonlinear problem |
|                    |               | 4 | Extremely nonlinear problem |
| OPT_QRANK1         | 0             | 0 | (false) Rank-1 updates by Broyden-approximation are inhibited. |
|                    |               | 1 | (true) Rank-1 updates by Broyden-approximation are allowed. |
| OPT_QORDI          | 0             | 0 | (false) Standard program mode. |
|                    |               | 1 | (true) Special program mode: Ordinary Newton iteration is done, e.g.: No damping strategy and no monotonicity test is applied. |
| OPT_QSIMPL         | 0             | 0 | (false) Standard program mode |
|                    |               | 1 | (true)  Special program mode: Simplified Newton iteration is done, e.g.: The Jacobian computed at the starting point is fixed for all subsequent iteration steps, and no damping strategy and no monotonicity test is applied. |
| OPT_NOROWSCAL      | 0             | 0 | (false) Automatic row scaling of the linear system is active:<br> Rows i=1,...,n will be divided by max j=1,...,n (abs(a[i,j])) |
|                    |               | 1 | (true) No row scaling of the linear system. Recommended only for well scaled nonlinear systems. |
| OPT_BOUNDEDDAMP    | 0             | 0 | The default bounded damping strategy switch takes place, dependent on the setting of OPT_NONLIN :<br>OPT_NONLIN = 0,1,2,3 -> OPT_BOUNDEDDAMP = off,<br>OPT_NONLIN = 4 -> OPT_BOUNDEDDAMP = on |
|                    |               | 1 | means always OPT_BOUNDEDDAMP = on |
|                    |               | 2 | means always OPT_BOUNDEDDAMP = off |
| OPT_IORMON         | 2             | 1 | Convergence order is not checked, the iteration will be always proceeded until the solution has the required precision OPT_RTOL (or some error condition occured) |
|                    |               | 2 | Use additional 'weak stop' criterion: Convergence order is monitored and termination due to slowdown of the convergence may occur. |
|                    |               | 3 | Use additional 'hard stop' criterion: Convergence order is monitored and termination due to superlinear convergence slowdown may occur. |
|                    |               |   | In case of termination due to convergence slowdown, the warning code retCode=4 will be set. In cases, where the Newton iteration converges but superlinear convergence order has never been detected, the warning code retCode=5 is returned. |
| OPT_NITMAX         | 50            |   | Maximum number of permitted iteration steps |
| OPT_FCBAND         | 10.0          |   | Bounded damping strategy restriction factor |
| OPT_SIGMA          |  |  | Broyden-approximation decision parameter. Required choice: OPT_SIGMA >= 1. Increasing this parameter make it less probable that the algorithm performs rank-1 updates. Rank-1 updates are inhibited, if OPT_SIGMA > 1/OPT_FCMIN is set. (see note 4) |
| OPT_SIGMA2         |  |  | Decision parameter about increasing damping factor to corrector if predictor is small. Required choice: OPT_SIGMA2 > 1. Increasing this parameter make it less probable that the algorithm performs rank-1 updates.  |
| OPT_AJDEL          | sqrt(epMach*10) |  | Jacobian approximation without feedback: Relative pertubation for components. epMach: relative machine precision) |
| OPT_AJMIN          | 0.0 |  | Jacobian approximation without feedback: Threshold value. The absolute pertubation for component k is computed by delx := OPT_AJDEL*max(abs(x[k]),OPT_AJMIN) |
| OPT_ETADIF         | 1e-6 |  | Jacobian approximation with feedback: Target value for relative pertubation eta of x |
| OPT_ETAINI         | 1e-6 |  | Jacobian approximation with feedback: Initial value for denominator differences |
| OPT_NBROY          | |  | Maximum number of possible consecutive iterative Broyden steps. Default is n if OPT_MSTOR=0, and OPT_ML+OPT_MU+1 if OPT_MSTOR=1 (but minimum is always 10) provided that Broyden is allowed. If Broyden is inhibited, NBROY is always set to zero.|
| OPT_FCSTART        |  |  | Damping factor for first Newton iteration: overrides option OPT_NONLIN, if set (see note 4) |
| OPT_FCMIN          |  |  | Minimal allowed damping factor (see note 4) |

Options specific to NLEQ2
=========================
| Name               | Meaning       |
|:-------------------|:--------------|
| OPT_IRANK          |Initially proposed (in) and final (out) rank of Jacobian |
| OPT_COND           |Maximum permitted subcondition for rank-decision by linear solver. |
