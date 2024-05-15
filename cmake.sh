#!/bin/bash

# Shell script for creating symbolic links for cofact to enable building the executable.
# This script assumes that the cofact directory is contained within the Prime95 source's
# gwnum directory, and creates links to the required files in the directories above.
#
# You should have built the gwnum.a library already, using the appropriate Makefile in the
# gwnum directory:
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
# ln -s ../gwnum.ld gwnum.ld

ln -s ../../common.h common.h
ln -s ../../exponentiate.c exponentiate.c
ln -s ../../exponentiate.h exponentiate.h
ln -s ../giants.h giants.h
ln -s ../gwcommon.h gwcommon.h
ln -s ../gwini.h gwini.h
ln -s ../gwnum.a gwnum.a
ln -s ../gwnum.h gwnum.h
ln -s ../gwthread.h gwthread.h
ln -s ../gwutil.h gwutil.h
ln -s ../../proof_hash.c proof_hash.c
ln -s ../../sha3.c sha3.c
ln -s ../../sha3.h sha3.h
make
rm *.h exponentiate.c gwnum.a proof_hash.c sha3.c
./cofact -h
