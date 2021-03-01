# intel and gnu compilers are supported
compiler = intel
flag = -O3

# User does not have to take care of following variables
gnumkl = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

libFoptim.a: trust_region.o strong_Wolfe.o linalg.o BFGS.o
ifeq ($(compiler),intel)
	xiar rcs $@ $^
else
	ar rcs $@ $^
endif

%.o: source/%.f90
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
test: f90.exe cpp.exe
	./f90.exe
	./cpp.exe

f90.exe: test/main.f90
	ifort -parallel -mkl -static-intel -ipo $(flag) $^ lib/libFoptim.a -o $@

cpp.exe: test/main.cpp
	icpc  -parallel -mkl -static-intel -ipo $(flag) -Iinclude/ $^ lib/libFoptim.a -lifcore -o $@

.PHONY: clean
clean:
	rm -f *.mod *.o *.a *.exe
