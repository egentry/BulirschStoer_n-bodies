compiler=gfortran

non_axisymmetric_potential: non_axisymmetric_potential.f08 ../bulirsch_stoer.f08 physics_non_axisymmetric.f08
	$(compiler) -c physics_non_axisymmetric.f08
	$(compiler) -c ../bulirsch_stoer.f08 
	$(compiler) -o non_axisymmetric_potential  physics_non_axisymmetric.f08 ../bulirsch_stoer.f08 non_axisymmetric_potential.f08
	rm *.mod
	rm *.o
	./non_axisymmetric_potential