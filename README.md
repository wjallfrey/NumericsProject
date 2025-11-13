# Numerics project - Will Allfrey
Requires numpy and matplotlib

This is a report for simulating solutions to the 1-D Viscous Burgers' equation on a periodic domain using either FTCS (with NT time steps, N grid spacings, up to time T, with the visoscity parameter nu) or a Fast Fourier Transform Spectral method (with NT timesteps, N modes, up to time T, visoscity parameter nu). The initial condition in both cases is sin^2(2x) on [0,pi].

It can test for stability of the FTCS scheme on certain space-time grids and parameter values. For the FTCS scheme it can plot at solutions of the scheme at the final time T (using plotting==True), or can plot at 10 evenly spaced time intervals between 0 and T on the same grid (using multiplot==True).

In the FTCS scheme, if testing for instability then ensure the "stabilityTest" argument is True, otherwise the function will break and return NaNs. If stabilityTest==True then the output will be a True if unstable, and False if the scheme is stable up to T.

