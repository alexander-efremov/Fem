BUILD := build

all: icc openmp cuda
    
icc:
	make -f $(BUILD)/makefile.icpc

openmp:
	make -f $(BUILD)/makefile.openmp
    
cuda:
	make -f $(BUILD)/makefile.cuda

clean:
	make -f $(BUILD)/makefile.icpc clean
	make -f $(BUILD)/makefile.openmp clean
	make -f $(BUILD)/makefile.cuda clean
	
clobber: clean