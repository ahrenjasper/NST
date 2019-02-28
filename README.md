
NST: A simple minimum-on-the-seam-of-crossings (MSX) optimizer and 
nonadiabatic statistical theory (NST) flux calculator

A. Jasper, M. Pfeifle and S. J. Klippenstein
February 28, 2019

This code can be used to
(1) refine the geometry of a good MSX (often called MECP) guess,
(2) compute Landau-Zener NST (often called NA TST) rates, and 
(3) generate the requisite data files for use in MESS and VariFlex.

I. References

The MSX optimization strategy is from
M. J. Bearpark, M. A. Robb, and H. B. Schlegel, Chem. Phys. Lett. 223, 
269 (1994).

The step-size control in the MSX optimizer is adopted from
J. Barzilai and J. M. Borwein, IMA J. Numer. Anal. 8, 141 (1988).

The effective Hessian used in the NST calculation is from
J. N. Harvey and M. Aschi, Phys. Chem. Chem. Phys. 1, 5555 (1999).

NST (often called NA TST) has been described in
J. N. Harvey, Phys. Chem. Chem. Phys. 9, 331 (2007) and
A. W. Jasper, J. Phys. Chem. A 119, 7339 (2015).

The refined expression for the curve-crossing probability (Airy formula) is from
J. B. Delos, J. Chem. Phys. 59, 2365 (1973).

