# intel and gnu compilers are supported
compiler = intel
flag = -O3

# user does not have to take care of following variables
gnumkl = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

libFoptim.a: show_time.o linalg.o \
Wolfe_1st.o strong_Wolfe_1st.o steepest_descent.o steepest_descent_verbose.o CGDY.o CGDY_verbose.o CGPR.o CGPR_verbose.o \
Wolfe_2nd.o strong_Wolfe_2nd.o NewtonRaphson.o BFGS.o \
trust_region.o trust_region_verbose.o Gauss_BFGS.o \
ALagrangian_NewtonRaphson.o ALagrangian_BFGS.o
ifeq ($(compiler),intel)
	xiar rcs $@ $^
else
	ar rcs $@ $^
endif

# utilities
%.o: source/utility/%.f90
ifeq ($(compiler),intel)
	ifort -fpp -parallel -mkl -static-intel -ipo $(flag) -c $<
else
	gfortran -cpp -ffree-line-length-0 -fno-range-check -I${MKLROOT}/include $(flag) -c $<
endif
# 1st-order line search
%.o: source/line-search_1st/%.f90
ifeq ($(compiler),intel)
	ifort -fpp -parallel -mkl -static-intel -ipo $(flag) -c $<
else
	gfortran -cpp -ffree-line-length-0 -fno-range-check -I${MKLROOT}/include $(flag) -c $<
endif
# 2nd-order line search
%.o: source/line-search_2nd/%.f90
ifeq ($(compiler),intel)
	ifort -fpp -parallel -mkl -static-intel -ipo $(flag) -c $<
else
	gfortran -cpp -ffree-line-length-0 -fno-range-check -I${MKLROOT}/include $(flag) -c $<
endif
# least square
%.o: source/least-square/%.f90
ifeq ($(compiler),intel)
	ifort -fpp -parallel -mkl -static-intel -ipo $(flag) -c $<
else
	gfortran -cpp -ffree-line-length-0 -fno-range-check -I${MKLROOT}/include $(flag) -c $<
endif
# equality constraint
%.o: source/constraint/%.f90
ifeq ($(compiler),intel)
	ifort -fpp -parallel -mkl -static-intel -ipo $(flag) -c $<
else
	gfortran -cpp -ffree-line-length-0 -fno-range-check -I${MKLROOT}/include $(flag) -c $<
endif

.PHONY: install
install: | lib
	mv *.a lib

lib:
	mkdir lib

.PHONY: test
test: test/f90.exe test/cpp.exe
	./test/f90.exe > test/f90.log
	./test/cpp.exe > test/cpp.log

test/f90.exe: test/main.f90 lib/libFoptim.a
	ifort -parallel -mkl -static-intel -ipo $(flag) $^ -o $@

test/cpp.exe: test/main.cpp lib/libFoptim.a
	icpc  -parallel -mkl -static-intel -ipo $(flag) -Iinclude/ $^ -lifcore -o $@

.PHONY: clean
clean:
	rm -f *.mod *.o *.a *.exe
