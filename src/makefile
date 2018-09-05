all: BAckPay.so Generate.so 

BAckPay.so: simulateddata1.c mainFunction.c utils.c myfunction.c myfunction.h utils.h 
	R CMD SHLIB simulateddata1.c mainFunction.c utils.c myfunction.c -lgsl -lgslcblas -o BAckPay.so 

Generate.so: simulateddata1.c  utils.c  myfunction.h utils.h 
	R CMD SHLIB simulateddata1.c  utils.c  -lgsl -lgslcblas -o Generate.so 



clean:
	rm -rf *.o *.so



