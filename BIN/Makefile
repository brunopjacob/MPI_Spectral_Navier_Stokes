FOR = mpif90
FLUID_RES_DIR = ../RESULTS/FLUID_FIELD
STATS_RES_DIR = ../RESULTS/STATISTICS
OBJ_DIR	= .

EXEC1 = SPECTRAL3DP


SRCS =	../SRC/declare_variables.f90\
	../SRC/read_input.f90\
	../SRC/constants.f90\
	../SRC/incompressible_solver.f90\
	../SRC/vtk.f90\
	../FFTE-6/pzfft3dv.f\
	../FFTE-6/fft235.f\
	../FFTE-6/kernel.f\
	../FFTE-6/factor.f\
	../SRC/main.f90\



CLEAR = clear

#OPT =  -O2  -stand f03  -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone

OPT =  -O2 -ffree-line-length-512

#OPT =  -fpe0 -traceback -heap-arrays
 
LIBS = -lm

$(EXEC1) : $(SRCS) 
	 $(FOR) $(OPT) -o $(EXEC1) $(SRCS) $(LIBS) $(OPT) 
	rm -rf $(OBJ_DIR)/*.mod
clean : 
	rm -f *.o core
	rm -f $(EXEC1)	
	rm -rf $(OBJ_DIR)/pbs.out
	rm -rf $(OBJ_DIR)/pbs.err
	rm -rf $(OBJ_DIR)/*.txt
	rm -rf $(OBJ_DIR)/out
	rm -rf $(OBJ_DIR)/log_error
	rm -rf $(OBJ_DIR)/log_out
	rm -rf $(OBJ_DIR)/*~
	rm -rf $(OBJ_DIR)/*.mod

clean_results : 
	rm -rf $(FLUID_RES_DIR)/*.*
	rm -rf $(STATS_RES_DIR)/*.*

run :	$(EXEC1) 
	clear
	mpirun -np 4 $(EXEC1)
