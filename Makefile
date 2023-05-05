# This Makefile builds the cofact program. It assumes that the GMP library has been installed at the system level,
# and that links to the following gwnum files are created in the local directory:
#	gwnum.ld gwnum.a giants.h gwcommon.h gwnum.h gwthread.h
 
cofact: cofact.o gwnum.a
	gcc cofact.o gwnum.a gwnum.ld -lm -lgmp -ldl -lpthread -lstdc++ -o cofact

cofact.o: cofact.c
	gcc -c -O2 -m64 -Wall -funroll-loops -fno-inline cofact.c

clean:
	rm -f *.o cofact

