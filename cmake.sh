#!/bin/bash

# Shell script for creating symbolic links for cofact to enable building the executable.
# This script requires a path to the Prime95 source directory, linking to the required
# files in the Prime95 source (p95v3019b...) and gwnum directories.
#
# You should have built the gwnum.a library already, using the appropriate Makefile in the
# gwnum directory (please let us know if cofact works on various flavours of devices!):
# 
# make -f make64            (Linux 64bit)
# make -f makebsd64         (FreeBSD 12.3 64bit)
# make -f makefile          (Linux, FreeBSD 32bit)
# make -f haiku             (Haiku 32bit)
# make -f makemac           (macOS)
# make -f makemsys          (Mingw/Msys)
# make -f makemw64          (Mingw/MSys 64bit)
# make -f makeos2           (OS/2 32bit)
#
# After making cofact, install to your binary directory e.g., sudo cp cofact /usr/local/bin
#
# If you are using an older version of the Prime95 source you may also need:
# ln -s $1/gwnum/gwnum.ld gwnum.ld

ln -s $1/common.h common.h
ln -s $1/exponentiate.c exponentiate.c
ln -s $1/exponentiate.h exponentiate.h
ln -s $1/gwnum/giants.h giants.h
ln -s $1/gwnum/gwcommon.h gwcommon.h
ln -s $1/gwnum/gwini.h gwini.h
ln -s $1/gwnum/gwnum.a gwnum.a
ln -s $1/gwnum/gwnum.h gwnum.h
ln -s $1/gwnum/gwthread.h gwthread.h
ln -s $1/gwnum/gwutil.h gwutil.h
ln -s $1/proof_hash.c proof_hash.c
ln -s $1/sha3.c sha3.c
ln -s $1/sha3.h sha3.h
make
rm *.h exponentiate.c gwnum.a proof_hash.c sha3.c
./cofact -h
