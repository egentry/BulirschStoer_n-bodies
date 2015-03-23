FC=gfortran
FLAGS        = -O2
FLAGS_DEBUG  = -ggdb -g -Wall
EXE_FILE = pythagorean_3body
OBJS =  pythagorean_3body.o physics_pythagorean.o bulirsch_stoer.o 
 
############# SUFFIX RULES ########################################
$(EXE_FILE) : $(OBJS) ../bulirsch_stoer.f08
	$(FC) $(FLAGS) $(OBJS) -o $(EXE_FILE)
	@echo "linking..."

.SUFFIXES :             # clear all defaults
.SUFFIXES : .f08 .o     # and replace

.f08.o:
	$(FC) $(FLAGS) -c  $<

############# MAKE RULES ########################################
clean:
	@echo "cleaning build files and integration data"
	@rm -f *.o *.mod $(EXE_FILE)
	@rm -f bulirsch_stoer.f08
	@rm -f *.dat

debug: FLAGS = $(FLAGS_DEBUG) 
debug: clean 			 # is there a better way to ensure all .o's are debug friendly?
debug: $(EXE_FILE)

############# DEPENDENCIES ########################################
bulirsch_stoer.f08: ../bulirsch_stoer.f08
	cp ../bulirsch_stoer.f08 bulirsch_stoer.f08  # there's probably a better way to do this...

pythagorean_3body.o: physics_pythagorean.o bulirsch_stoer.o

bulirsch_stoer.o: physics_pythagorean.o