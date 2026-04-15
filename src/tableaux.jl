# Standard explicit RK tableaux supplied for the effectively-symmetric methods.

const _SQRT2 = sqrt(2.0)

const EES25_A21 = 1 / 3
const EES25_A31 = -5 / 48
const EES25_A32 = 15 / 16
const EES25_B = (1 / 10, 1 / 2, 2 / 5)
const EES25_C = (0.0, 1 / 3, 5 / 6)

const EES27_A21 = (2 - _SQRT2) / 3
const EES27_A31 = (-4 + _SQRT2) / 24
const EES27_A32 = (4 + _SQRT2) / 8
const EES27_A41 = (-176 + 145 * _SQRT2) / 168
const EES27_A42 = 3 * (8 - 5 * _SQRT2) / 56
const EES27_A43 = 3 * (3 - _SQRT2) / 7
const EES27_B = (
    (5 - 3 * _SQRT2) / 14,
    (3 + _SQRT2) / 14,
    3 * (-1 + 2 * _SQRT2) / 14,
    (9 - 4 * _SQRT2) / 14,
)
const EES27_C = (
    0.0,
    (2 - _SQRT2) / 3,
    (2 + _SQRT2) / 6,
    (4 + _SQRT2) / 6,
)

# Existing Williamson 2N / commutator-free EES25 coefficients retained for the
# dedicated low-storage and Lie-group backends.
const EES25_A2end = (-0.5, -2.0)
const EES25_B1    = 0.5
const EES25_B2end = (1.0, 0.25)
const EES25_C2end = (0.5, 1.0)
