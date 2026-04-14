# EES25 tableau coefficients in 2N low-storage form.
#
# The 3-stage recurrence is:
#   tmp = A[i] * tmp + dt * f(u, p, t + C[i]*dt)
#   u   = u + B[i] * tmp
#
# where A[1] = 0 implicitly (first stage has no history term).

const EES25_A2end = (-0.5, -2.0)
const EES25_B1    = 0.5
const EES25_B2end = (1.0, 0.25)
const EES25_C2end = (0.5, 1.0)
