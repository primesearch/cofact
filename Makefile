# This Makefile builds the cofact program. It assumes that the GMP library has been 
# installed at the system level, and requires symbolic links to the following gwnum and 
# Prime95/mprime files from the source:
#
# p95v30 directory:
#	common.h exponentiate.c exponentiate.h proof_hash.c sha3.c sha3.h
# gwnum  directory:
#	giants.h gwcommon.h gwini.h gwnum.a gwnum.h gwthread.h gwutil.h
#
# To make this easy, either unzip cofact inside the gwnum directory, or cd p95v30.../gwnum 
# and run:
#	git clone --branch cxc https://github.com/primesearch/cofact.git
#
# From there, start with:
#	make -f [makemac|make64|makemw64|...]
# to compile gwnum.a for your 64bit [Mac|Linux|Win|...] system.
#
# Once that is done, you may proceed to:
#	cd cofact; bash cmake.sh
#
# The shell script creates symbolic links to the gwnum files above, then runs make on this Makefile.
#
# Later than version v30.8 of the Prime95 source is desirable, to avoid also requiring gwnum.ld
#
# If you have any suggestions or issues please raise these at:
#	https://github.com/primesearch/cofact/issues/

cofact: cofact.o gwnum.a
	gcc cofact.o gwnum.a -lm -lgmp -lpthread -lstdc++ -o cofact
cofact.o: cofact.c
	gcc -c -O2 -m64 -Wall -funroll-loops -fno-inline cofact.c
clean:
	rm -f *.o cofact
