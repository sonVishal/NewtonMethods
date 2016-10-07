### Note 1 :
The machine dependent values small, great and epMach are
set as global variables over here. As delivered, this function is adapted
to use constants suitable for all machines with IEEE arithmetic.
If you use another type of machine, you have to change these variables suitable
for your machine.

### Note 2 :
If NLEQ solver terminates with retCode=2 (maximum iterations)
or retCode=3 (small damping factor), you may try to continue
the iteration by increasing OPT_NITMAX or decreasing OPT_FCMIN
and setting OPT_QSUCC to 1.

### Note 3 :
Storage of user supplied banded Jacobiann the band matrix case, the following lines may build
up the analytic Jacobian A; Here AFL denotes the quadratic matrix A in dense form,
and ABD the rectangular matrix A in banded form :
```
ML = getOption(opt,OPT_ML,0)
MU = getOption(opt,OPT_MU,0)
MH = MU+1
for J = 1:N
    I1 = MAX(1,J-MU)
    I2 = MIN(N,J+ML)
    for I = I1:I2
        K = I-J+MH
        ABD[K,J] = AFL[I,J]
    end
end
```
The total number of rows needed in ABD is ML+MU+1 .
The MU by MU upper left triangle and the ML by ML lower right
triangle are not referenced.

### Note 4 :
The default values of the internal parameters may be obtained
from the monitor output with at least the field OPT_PRINTIOMON set to 2
and by initializing the corresponding options to zero.

### Note 5 :
in case of failure:
- use non-standard options
- or turn to Newton-algorithm with rank strategy
- use another initial guess
- or reformulate model
- or apply continuation techniques
