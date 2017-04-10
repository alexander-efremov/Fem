BUILD := build

all: icpc openmp cuda main

main:
	make -f $(BUILD)/makefile.main

icpc:
	make -f $(BUILD)/makefile.icpc

openmp:
	make -f $(BUILD)/makefile.openmp
    
cuda:
	make -f $(BUILD)/makefile.cuda

clean:
	make -f $(BUILD)/makefile.icpc clean
	make -f $(BUILD)/makefile.openmp clean
	make -f $(BUILD)/makefile.cuda clean
	make -f $(BUILD)/makefile.main clean
	
clobber: clean