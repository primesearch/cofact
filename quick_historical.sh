# Run cofact on each Fermat number using the best method for that number, "best" meaning reasonably 
# fast.
#
# This run should take less than an hour, as we exclude several mode 2 -cpr tests (i.e. running 
# the full Pepin test) for three reasonably large Fermat numbers (F20, F22, and F24).
#
# To reduce overall runtime, the F22 run has been changed to a mode 3 -upr test (so that cofact 
# only evaluates the Suyama test, reading from a proof file). Using a Pepin test here with the 
# interim flag -i would allow verification of the results of Trevisan and Carvalho (1995).
# 
# To run the Pepin test, simply change the -upr flag to -cpr. (This holds if you want to 
# directly run Pepin tests on any Fermat number greater than F17.)
#
# The Suyama test is not possible for F20 and F24 as they have no known factors, so cofact
# can display the Suyama A residue and no more.
#
# This script assumes that mprime has already generated proof files for F17 to F29 (the F30 test is 
# commented out).
#
# (Proof files for F14 to F29 can be downloaded from https://64ordle.au/fermat/ and should be 
# located in the same directory as this script.)
#
# We suggest copying the compiled cofact executable to a directory for binaries in your path, and 
# redirecting the output of this script to a file in the background, with something like:
#
#       nice ./quick_historical.sh > quick_history.txt &
#
# If you have a problem executing this script try:      sudo chmod +x quick_historical.sh
#
# First, we print the verbose help as explanation of the remainder of the output.

cofact -h -sv

# For F0 - F16, mprime does not generate proof files. So just have cofact perform the Pepin and 
# Suyama tests.

echo "1640: Pierre de Fermat found F0 to F4 to be prime:
"
cofact 0
cofact 1
cofact 2
cofact -i 3
cofact -i -sep 4

echo "Example of Suyama test on F5 factor, 641 (Euler, 1732):
"
cofact -amv -sep 5 641
echo "Example of Suyama test on F5 factor, 6700417 (Euler, 1732):
"
cofact -m -sep 5 6700417
cofact -imv -sep 6 274177
cofact -m -sep 6 67280421310721

echo "Pepin test on F7 (1905); note Morehead's residue = P7 + 1, while Western's is 3^2^119 mod F7:
"
cofact -iov -sep 7 59649589127497217
cofact -m -sep 7 5704689200685129054721

echo "Pepin test on F8 (1909); note Morehead & Western's result is the interim residue at iteration 246:
"
cofact -iov -sep 8 1238926361552897
cofact -m -sep 8 93461639715357977769163558199606896584051237541638188580280321

echo "Suyama test on 148-digit cofactor of F9 (Brillhart, 1967; 1 known factor):
"
cofact -iv -sep 9 2424833
echo "F9 is completely factored with 3 factors (largest factor with 99 digits):
"
cofact -m -sep 9 2424833 7455602825647884208337395736200454918783366342657

echo "Pepin test on F10 (Robinson, 1952, Lehmer, 1953):"
echo "Suyama test on 291-digit cofactor of F10 (Brillhart, 1967; 2 known factors):
"
cofact -iov -sep 10 45592577 6487031809
echo "F10 is completely factored with 4 factors (largest factor with 252 digits):
"
cofact -m -sep 10 45592577 6487031809 4659775785220018543264560743076778192897

echo "Suyama test on 606-digit cofactor of F11 (Wagstaff, 2 known factors):
"
cofact -iv -sep 11 319489 974849
echo "F11 is completely factored with 5 factors (largest factor with 564 digits):
"
cofact -m -sep 11 319489 974849 167988556341760475137 3560841906445833920513

echo "Suyama test on 1,202-digit cofactor of F12 (Wagstaff, 4 known factors):
"
cofact -iv -sep 12 114689 26017793 63766529 190274191361
echo "Suyama test on 1,187-digit cofactor of F12 (Baillie, 1986, 5 known factors):
"
cofact -sep 12 114689 26017793 63766529 190274191361 1256132134125569
echo "Suyama test on 1,133-digit cofactor of F12 (Vang, Batalov et al., 2010, 6 known factors):
"
cofact -sep 12 114689 26017793 63766529 190274191361 1256132134125569 568630647535356955169033410940867804839360742060818433

echo "Pepin test on F13 (Paxson, 1960):"
echo "Suyama test on 2,454-digit cofactor of F13 (Wagstaff, 1 known factor):
"
cofact -iov -sep 13 2710954639361 
echo "Suyama test on 2,436-digit cofactor of F13 (Lenstra, Keller, 1991, 2 known factors):
"
cofact -sep 13 2710954639361 2663848877152141313
echo "Suyama test on 2,417-digit cofactor of F13 (Crandall, 1991, 3 known factors):
"
cofact -sep 13 2710954639361 2663848877152141313 3603109844542291969
echo "Suyama test on 2,391-digit cofactor of F13 (Brent, 1995, 4 known factors):
"
cofact -sep 13 2710954639361 2663848877152141313 3603109844542291969 319546020820551643220672513

echo "Pepin test on F14 (Hurwitz & Selfridge, 1961):"
echo "Suyama test on 4,880-digit cofactor of F14 (Lipp, Silverman, Rajala, Moore, 2010, 1 known factor):
"
cofact -iov -sep 14 116928085873074369829035993834596371340386703423373313

echo "Suyama test on 9,856-digit cofactor of F15 (Hiromi Suyama, 1984, 1 known factor):
"
cofact -iv -sep 15 1214251009
echo "Suyama test on 9,840-digit cofactor of F15 (Suyama, Baillie, 1987, 2 known factors):
"
cofact -sep 15 1214251009 2327042503868417
echo "Suyama test on 9,808-digit cofactor of F15 (Brent & Crandall, 1997, 3 known factors):
"
cofact -sep 15 1214251009 2327042503868417 168768817029516972383024127016961

echo "Suyama test on 19,720-digit cofactor of F16 (Baillie, 1987, 1 known factor):
"
cofact -iv -sep 16 825753601
echo "Suyama test on 19,694-digit cofactor of F16 (Brent & Crandall, 1996, 2 known factors):
"
cofact -sep 16 825753601 188981757975021318420037633 

# For F17 - F24, cofact is fairly fast. So normally we would have it perform the Pepin and Suyama 
# tests and then check the mprime proof file A residue.
# However, we don't really need Pepin results here as most of the interest is in the cofactors, so 
# only F17 runs the Pepin test, once (mode 2, -cpr option).

echo "Selfridge-Hurwitz (1964) published interim residue of Pepin test at iteration 20:"
echo "Suyama test on 39,444-digit cofactor of F17 (Baillie, 1987, 1 known factor):
"
cofact -ciov F17.proof -sep 17 31065037602817
echo "Suyama test on 39,395-digit cofactor of F17 (Chia, Hoeglund, Sorbera, 2011, 2 known factors):
"
cofact -u F17.proof -sep 17 31065037602817 7751061099802522589358967058392886922693580423169

echo "Suyama test on 78,907-digit cofactor of F18 (Chudnovsky twins, 1990, 1 known factor):
"
cofact -u F18.proof -sep 18 13631489
echo "Suyama test on 78,884-digit cofactor of F18 (Crandall, 1999, 2 known factors):
"
cofact -u F18.proof -sep 18 13631489 81274690703860512587777

echo "Suyama test on 157,804-digit cofactor of F19 (Crandall, Doenias et al., 1993, 2 known factors):
"
cofact -mu F19.proof -sep 19 70525124609 646730219521
echo "Suyama test on 157,770-digit cofactor of F19 (JRK, Kruppa et al., 2009, 3 known factors):
"
cofact -u F19.proof -sep 19 70525124609 646730219521 37590055514133754286524446080499713

echo "Omitting Pepin test on F20 (Young & Buell, 1987), we show Suyama A residue from proof:
"
#./cofact -ov -cpr F20.proof -sep 20
cofact -o -upr F20.proof -sep 20
# N.B. You should comment out either the -upr (quick) or -cpr (slow) line here (and perhaps remove "Omitting").

echo "Suyama test on 631,294-digit cofactor of F21 (Crandall, Doenias et al., 1993, 1 known factor):
"
cofact -mu F21.proof -sep 21 4485296422913

#echo "Pepin tests on F22 (Crandall, Doenias et al., Trevisan & Carvalho, both 1993):"
#echo "Trevisan & Carvalho also published (1995) 13 interim residues at various values:"
echo "Suyama test on 1,262,577-digit cofactor of F22 (Domanov, Yamada, 2010, 1 known factor):
"
#./cofact -iv -cpr F22.proof -sep -t 4 22 64658705994591851009055774868504577
cofact -upr F22.proof -sep 22 64658705994591851009055774868504577
# N.B. You should comment out either the -upr (quick) or -cpr (slow) line and add the comments.

echo "Suyama test on 2,525,215-digit cofactor of F23 (Crandall, Mayer et al., 2000, 1 known factor):
"
cofact -u F23.proof -sep 23 167772161

echo "Omitting Pepin test on F24 (Crandall, Mayer, Papadopoulos, 1999), we show Suyama A residue from proof:
"
#./cofact -v -cpr F24.proof -sep -t 8 24
cofact -upr F24.proof -sep 24
# N.B. You should comment out either the -upr (quick) or -cpr (slow) line here (and perhaps remove "Omitting").

# Note the F22 and F24 Pepin tests above will try to use 4 or 8 threads if possible.
#
# For F25 - F30, mprime is much faster than cofact. So here we have cofact use the A residue from the 
# mprime proof files and perform the Suyama test.
# In this mode, the gwnum library is not used. So specifying multiple threads would have no effect on speed.

echo "Suyama test on 10,100,842-digit cofactor of F25 (Yamada, Hoeglund, 2009, 3 known factors):
"
cofact -u F25.proof -sep 25 25991531462657 204393464266227713 2170072644496392193

echo "Suyama test on 20,201,768-digit cofactor of F26 (Hoeglund, 2009, 1 known factor):
"
cofact -u F26.proof -sep 26 76861124116481

echo "Suyama test on 40,403,531-digit cofactor of F27 (Hoeglund, 2010, 2 known factors):
"
cofact -u F27.proof -sep 27 151413703311361 231292694251438081

echo "Suyama test on 80,807,103-digit cofactor of F28 (Mayer, 2022, 1 known factor):
"
cofact -u F28.proof -sep 28 1766730974551267606529

echo "Suyama test on 161,614,233-digit cofactor of F29 (Mayer, 2022, 1 known factor):
"
cofact -u F29.proof -sep 29 2405286912458753

# Feel free to try running the following if you have a proof of F30. It isn't efficient to use 
# cofact to run Pepin tests at this size.

# echo "Suyama test on 323,228,467-digit cofactor of F30 (Mayer, 2022, 2 known factors):"
# ./cofact -u F30.proof -sep 30 640126220763137 1095981164658689
# # N.B. F30 proof file not generated yet; mprime requires AVX-512 to run on exponents this large

# Fermat numbers F31 and beyond are not supported by gwnum; Mlucas does not yet generate proofs

