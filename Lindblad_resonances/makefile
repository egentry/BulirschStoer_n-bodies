compiler=gfortran

Lindblad_resonances: Lindblad_resonances.f08 ../bulirsch_stoer.f08 physics_Lindblad_resonances.f08
	$(compiler) -O2 -c physics_Lindblad_resonances.f08
	$(compiler) -O2 -c ../bulirsch_stoer.f08 
	$(compiler) -O2 -o Lindblad_resonances  physics_Lindblad_resonances.f08 ../bulirsch_stoer.f08 Lindblad_resonances.f08

clean:
	rm physics.mod physics_Lindblad_resonances.o bulirsch_stoer.mod bulirsch_stoer.o
	rm Lindblad_resonances
	rm *.dat

debug:
	$(compiler) -g -c physics_Lindblad_resonances.f08
	$(compiler) -g -c ../bulirsch_stoer.f08 
	$(compiler) -g -o Lindblad_resonances  physics_Lindblad_resonances.f08 ../bulirsch_stoer.f08 Lindblad_resonances.f08