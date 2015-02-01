compiler=gfortran

pythagorean_3body: pythagorean_3body.f08 ../bulirsch_stoer.f08 physics_pythagorean.f08
	$(compiler) -c physics_pythagorean.f08
	$(compiler) -c ../bulirsch_stoer.f08 
	$(compiler) -o pythagorean_3body  physics_pythagorean.f08 ../bulirsch_stoer.f08 pythagorean_3body.f08
	rm *.mod
	rm *.o
	./pythagorean_3body