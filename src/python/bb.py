
# Set global keys
ETOL = 5.e-5      # energy gap convergence tolerance in kcal/mol
GTOL = 1.e-5      # gradient norm convergence tolerance in Hartree/bohr
MAXSTEP = 1000    # maximum number of steps
STEPSCALE = 1.    # scales stepsize
N1=1              # PES call index for surface 1
N2=2              # PES call index for surface 2
EXCTOLE = 10.0e0  # energy tolerance for excessive step identification (kcal/mol)
EXCTOLG = 0.3e0   # gradient norm tolerance for excessive step id. (Hartree/bohr)
EXCSCAL = 0.5e0   # scale down step by this amount if excessive step was detected


def calc_grad_and_ene(geom, method, basis, prog)
    """ calls elstruct to get the energy and gradinet
    """

    return grad, energy

# diabatic energies
pemd = np.array([n1n1, n1n2],
                [n2n1, n2n2])


#####




