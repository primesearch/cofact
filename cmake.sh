#!/bin/bash

# Shell script for creating symbolic links for cofact to enable building the executable
# This script assumes that the cofact directory is contained within the Prime95 source's
# gwnum directory, and creates links to the required files in the directories above.
#
# The gwnum.a file should be built using the appropriate Makefile in the gwnum directory,
# e.g.
# make -f makemac           (macOS)

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
./cofact -h
