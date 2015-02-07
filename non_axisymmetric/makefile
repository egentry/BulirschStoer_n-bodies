compiler=gfortran

non_axisymmetric_potential: non_axisymmetric_potential.f08 ../bulirsch_stoer.f08 physics_non_axisymmetric.f08
	$(compiler) -O2 -c physics_non_axisymmetric.f08
	$(compiler) -O2 -c ../bulirsch_stoer.f08 
	$(compiler) -O2 -o non_axisymmetric_potential  physics_non_axisymmetric.f08 ../bulirsch_stoer.f08 non_axisymmetric_potential.f08

clean:
	rm physics.mod physics_non_axisymmetric.o bulirsch_stoer.mod bulirsch_stoer.o
	rm non_axisymmetric_potential
	rm *.dat
	rm surface_of_section.eps surface_of_section.pdf
