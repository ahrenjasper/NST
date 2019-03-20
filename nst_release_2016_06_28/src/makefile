
FC=gfortran
LIB= -llapack

nstmake: dd.f nst.f msxfreq.f proj.f airy.f
	$(FC) -O3 dd.f nst.f msxfreq.f proj.f airy.f -o ../exe/nst.x $(LIB)


