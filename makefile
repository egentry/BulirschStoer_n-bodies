cc=gfortran

pythagorean_3body: physics_pythagorean.f08 bulirsch_stoer.f08 pythagorean_3body.f08 
	$(cc) -o pythagorean_3body physics_pythagorean.f08 bulirsch_stoer.f08 pythagorean_3body.f08
	./pythagorean_3body