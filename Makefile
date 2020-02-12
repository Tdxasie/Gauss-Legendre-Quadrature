quad: precision_mod.o functions_mod.o quad_mod.o quad.o
	gfortran -o quad precision_mod.o functions_mod.o quad_mod.o quad.o

precision_mod.o: precision_mod.f90
	gfortran -c precision_mod.f90

functions_mod.o: precision_mod.o functions_mod.f90
	gfortran -c functions_mod.f90

quad_mod.o: precision_mod.o functions_mod.o quad_mod.f90
	gfortran -c quad_mod.f90

quad.o: quad_mod.o precision_mod.o functions_mod.o quad.f90
	gfortran -c quad.f90
