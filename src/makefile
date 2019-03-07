
FC=gfortran
LIB= -llapack

nstmake: dd.f nst.f msxfreq.f proj.f airy.f opt.f
	$(FC) -g -O3 dd.f nst.f msxfreq.f proj.f airy.f opt.f -o ./nst.x $(LIB)


