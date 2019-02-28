
NST: A simple minimum-on-the-seam-of-crossings (MSX) optimizer and 
nonadiabatic statistical theory (NST) flux calculator

A. Jasper, M. Pfeifle and S. J. Klippenstein
June 27, 2016

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

The distributed code is set up to call Gaussian, but it can be made to 
call Molpro by appending the argument M to the command-line call of 
nst.x.

The code assumes a nonlinear molecule when counting states.

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
(see, e.g., runs/c2h4o/*)
  input (std input)     ! main input file, see next section for formatting
  qc.1                  ! gaussian template for state 1 (e.g., a triplet)
  qc.3                  ! gaussian template for state 2 (e.g., a singlet)
  g.x                   ! system call to gaussian
  m.x                   ! system call to Molpro, in case "nst.x M" was invoked

Output files
(see, e.g., runs/c2h4o/*)
  output                ! standard output, prints code progress, the 
                        ! optimized MSX geometry, rotational constants 
                        ! and frequencies
  fort.80               ! a Molden-style optimization movie
  nej_lz.dat            ! NST state counts in 2D VariFlex format, using Landau-Zener transition probabilities
  nej_delos.dat         ! NST state counts in 2D VariFlex format, using Airy formula by Delos
  ne_lz.dat             ! NST state counts in 1D MESS format, using Landau-Zener transition probabilities
  ne_delos.dat          ! NST state counts in 1D MESS format, using Airy formula by Delos

IV. Standard Input Formatting
The following example for O+C2H4 is provided with the distribution in runs/c2h4o/* .

C2H4O m062x/dz S0/T0    ! header, not used
7 -153.6288736          ! number of atoms, zero of energy in Hartrees
O 15.9949      1.140494420833375    -0.280674562604505    -0.034835946592096  ! atom label, mas
s in amu, and Cartesian coordinates in A
C 12.0000      0.020624135447152     0.503088095375588     0.035030638055883
C 12.0000     -1.265276715371149    -0.242197552521405    -0.018073232163525
H 1.007825     0.057467001267133     1.332071729395067    -0.702477480127337
H 1.007825     0.142204307870437     1.029502775942943     1.020002581489454
H 1.007825    -1.293043764132468    -1.265522467200460     0.355766721259008
H 1.007825    -2.187220160488727     0.252073923366285    -0.322329491425112
T                       ! T = calculate Hessians and write nej.dat and ne.dat; F = don't, just optimize the MSX
100. 50000.             ! Energy grid spacing , maximum grid energy in cm-1
6 300                   ! Angular momentum grid spacing, maximum grid angular momentum in hbar
29. 3.                  ! Spin-orbit coupling in cm-1, and scaling factor for the state counts in ne_*.dat and nej_*.dat

