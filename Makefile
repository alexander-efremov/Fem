all: icc openmp cuda
    
icc:
	make -f makefile.icc

openmp:
	make -f makefile.openmp
    
cuda:
	make -f makefile.cuda

clean:
	make -f makefile.icc clean
	make -f makefile.openmp clean
	make -f makefile.cuda clean
	
clobber: clean