compiler=gfortran

LidovKozai: LidovKozai.f08 ../bulirsch_stoer.f08 physics_LidovKozai.f08
	$(compiler) -O2 -c physics_LidovKozai.f08
	$(compiler) -O2 -c ../bulirsch_stoer.f08 
	$(compiler) -O2 -o LidovKozai  physics_LidovKozai.f08 ../bulirsch_stoer.f08 LidovKozai.f08

clean:
	rm physics.mod physics_LidovKozai.o bulirsch_stoer.mod bulirsch_stoer.o
	rm LidovKozai
	rm *.dat