
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

II. Notes

The optimizer is not robust, and a good guess must be provided.
By default, the Barzilai-Borwein step-size control is used in the 
gradient descent search. This usually leads to faster convergence, 
but may fail in some cases. Appending the command-line argument 
noBB to the nst.x call turns off this algorithm.

The distributed code is set up to call Molpro.

III. Description of the Distribution

Compilation and execution 
  src/makefile          ! type make to create ../exe/nst.x
  exe/nst.x             ! executable: ./nst.x < [std input] > [std output]

Source files
  src/nst.f             ! driver, performs the MSX search
  src/msxfreq.f         ! computes the effective hessian, creates nej_lz.dat
                          (nej_delos.dat) and ne_lz.dat (ne_delos.dat) files
  src/dd.f              ! direct dynamics interface for gaussian; change 
                        ! to molpro calls via command-line argument (see above)
  src/proj.f            ! projects out translations, rotations, and the
                        ! gradient of the gap

Input files - edit these
Note: These should all be in the run time directory.
(see, e.g., examples/C2H4O/*)
  input (std input)     ! main input file, see next section for formatting
  qc.1                  ! gaussian template for state 1 (e.g., a triplet)
  qc.3                  ! gaussian template for state 2 (e.g., a singlet)
  g.x                   ! system call to Gaussian
  m.x                   ! system call to Molpro

Output files
(see, e.g., examples/C2H4O/*)
  output                ! standard output, prints code progress, the 
                        ! optimized MSX geometry, rotational constants 
                        ! and frequencies
  fort.80               ! a Molden-style optimization movie
  nej_lz.dat            ! NST state counts in 2D VariFlex format, using Landau-Zener transition probabilities
  nej_delos.dat         ! NST state counts in 2D VariFlex format, using Airy formula by Delos
  ne_lz.dat             ! NST state counts in 1D MESS format, using Landau-Zener transition probabilities
  ne_delos.dat          ! NST state counts in 1D MESS format, using Airy formula by Delos

