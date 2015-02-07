compiler=gfortran

pythagorean_3body: pythagorean_3body.f08 ../bulirsch_stoer.f08 physics_pythagorean.f08
	$(compiler) -O2 -c physics_pythagorean.f08
	$(compiler) -O2 -c ../bulirsch_stoer.f08 
	$(compiler) -O2 -o pythagorean_3body  physics_pythagorean.f08 ../bulirsch_stoer.f08 pythagorean_3body.f08

clean:
	rm physics.mod physics_pythagorean.o bulirsch_stoer.mod bulirsch_stoer.o
	rm pythagorean_3body
	rm *.dat
