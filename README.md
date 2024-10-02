# cofact
## A utility for testing the character of a Fermat or Mersenne cofactor
`cofact` is a companion program for George Woltman’s `mprime` that tests a small Fermat number for 
primality, or optionally in version 0.9, a small Mersenne number for probable primality (PRP). 
(In this context, ‘small’ means, up to about 320 million digits in size!)

If known factors are given, then the cofactor is likewise tested for probable primality.

`cofact` reports typical hexadecimal Res64 and Selfridge–Hurwitz residues for the Pépin or 
Fermat-PRP tests on Fermat and Mersenne numbers respectively, and the Suyama cofactor test, 
for comparison with the residues reported by other programs.

## Requirements
* Intel or compatible x86-family processor (AVX-512 required for the largest exponents)
* GNU Multiple Precision Arithmetic Library
* Copy of the Prime95/mprime source code, including the `gwnum` large number library. 

## Installation
You should either unzip the `cofact` download as a directory  within the Prime95 source library 
directory `gwnum`, or run `git clone` in that directory, as `cofact` uses a number of files from 
the Prime95 source. First, run the `gwnum` makefile for your platform, e.g. for \[macOS|Linux|Win\]:
```bash
cd p95v3019b17/gwnum; make -f [makemac|make64|makemw64]
```
Then for `cofact`:
```bash
git clone --branch cxc https://github.com/primesearch/cofact.git
cd cofact; bash cmake.sh
```
Finish the installation by copying the `cofact` executable to your usual directory for binaries.
### Basic operation
To test a Fermat number, type `cofact` followed by an exponent (a non-negative integer up to 30), 
optionally followed by any factors of that Fermat number. For instance, to test the fifth Fermat 
number $F_5$ and the factor $641$:
```bash
cofact 5 641
```
The `cofact` distribution includes a verifiable delay function (VDF) proof for a small Fermat 
number $F_{14}$, so to check the file has a correct value stored, type:
```bash
cofact -cpr F14.proof 14
```
This computes a primality test for $F_{14}$ and compares the result with the proof file. If you are 
satisfied you have a valid proof file the check may be skipped; to proceed directly to test the 
cofactor, by supplying known factors, type:
```bash
cofact -upr F14.proof 14 116928085873074369829035993834596371340386703423373313
```
You may check or use any proof file with or without factors.
## Basic command line options
The initial command line options have burgeoned since version 0.6. The basic options are:

Command line option | Function
--------------------|------------------------------
-cpr _filename_     | Run a primality test and then check the result against a proof file; if factors are supplied, run the cofactor test
-upr _filename_     | Skip primality testing and use the proof file to immediately run the cofactor test
-p _n_              | Print a progress report every _n_ iterations
-t _n_              | Run `cofact` with _n_ threads
-h                  | Print basic help
-v                  | Print verbose test information

The full set of menu options available in later versions are discussed [below](#full-list-of-features--command-line-options).

## Authors
Gary B. Gostin, versions 0.2 to 0.8.1 (the original, `main` branch)

Catherine X. Cowie, versions 0.6 to 0.9 (the `cxc` branch)

Using `git` you may switch between the two versions, prior to the build instructions above:
```bash
cd cofact; git switch [main|cxc]
```
## Copyright
This program is © 2023–2024 Gostin and Cowie under the Creative Commons Zero (CC-0) licence, we 
will not be responsible for whether you find this work useful, or whether it opens a portal to 
another hellish dystopian dimension where your computer used to be.

The `gwnum` library, and the proof validation module is © 2002–24 Mersenne Research, Inc, used 
with permission. All rights reserved.

The GMP library is © 1991, 1993–2016, 2018–2024 Free Software Foundation, Inc.

## Rationale
When running very long calculations such as testing a number for primality, one approach 
to ensure accuracy is to perform the calculation twice using two different programs running 
on two different computers, and then compare various calculation residues. This helps protect 
against both program design errors and hardware errors. Historically, part of the process of 
investigating untested Fermat numbers was to double-check previous Pépin test results obtained 
by previous researchers, which helps verify that the mathematics involved in carrying out the 
test is being correctly calculated.

In 1999, Ernst Mayer (with Richard Crandall and Jason Papadopoulos) had run the Pépin test on 
then-smallest known Fermat number lacking known prime factors, $F_{24}$, finding it composite. 
Beginning in 2013, he added the capability to run the Pépin test to his 
[`Mlucas` program](https://github.com/primesearch/Mlucas) as part of 
a decade-long programme to run the Pépin and Suyama cofactor tests on each of the next six 
Fermat numbers up to $F_{30}$, proving each of the cofactors to be composite.

Complete `Mlucas` results of the Pépin tests starting with $F_{13}$ up to $F_{30}$, and the Suyama 
test for each of those except for $F_{13}$, $F_{20}$, and $F_{24}$, have been 
[repeated](https://www.mersenneforum.org/node/17112?p=889067#post889067), 
and to meet the ideal of hardware and software-independent verification Wilfrid Keller and 
Gary Gostin began a follow-up project to verify Mayer’s `Mlucas` results using Woltman’s `mprime`.

`mprime` does not produce Selfridge–Hurwitz residues, however `cofact` uses mprime’s `gwnum` 
maths library for running the Pépin test on Fermat numbers, and results including 
Selfridge–Hurwitz residues from `cofact` are practical up to about $F_{26}$. `mprime` can 
generate Verifiable Delay Function proofs for Fermat-PRP tests (which constitutes the first 
and largest calculation in the Suyama test, described as the $A$ value), so for each 
composite Fermat number from $F_{12}$ up to $F_{29}$, VDF proofs have been 
[created](https://64ordle.au/fermat/) using `mprime` to permit the Suyama test to be run 
with `cofact`.

With the exception of $F_{30}$, and the Selfridge–Hurwitz residues of Pépin tests of 
$F_{28}$ and $F_{29}$, all of Mayer’s results obtained in 2013–22 using `Mlucas`, have been 
replicated independently using `mprime` and `cofact`.

## Mathematics
The Fermat numbers have the form $2^{2^m}+1$, $m \ge 0$, and the Mersenne numbers have the form 
$2^p - 1$, where $p$ is a prime number. The five smallest Fermat numbers are prime, and of the 
first ten million prime numbers, 51 prime exponents $p$ are currently known to give rise to a 
Mersenne prime. This utility can run a primality or probable primality test on any of these 
numbers in such a way that the test allows further testing of the cofactor, if the Fermat or 
Mersenne number is composite and has one or more known factors.

The practical computational limit for this sort of testing is reached around $F_{30}$ for Fermat 
numbers, and $M_{1,073,741,789}$ for Mersenne numbers.

If $F$ is a Fermat number, we are interested in obtaining the result $P$ of a Pépin test, where 
$P \equiv b^{\frac 1 2 (F - 1)} \equiv -1$ (mod $F$) if and only if $F$ is prime. Since 
$F = 2^{2^m} + 1$ is the result of $m$ squarings of the number $2$ and adding $1$, $P$ may be 
obtained by $2^m-1$ squarings of the base $b$, modulo $F$. Pépin’s theorem is a definitive 
primality test for bases such as 3, 5, 6, and so on. One further squaring then follows to prepare 
for the cofactor test below.

A similar test for Mersenne numbers is _not_ definitive (a different test, the Lucas–Lehmer test is 
used for that) but serves to determine whether the Mersenne number is _probably_ prime, which is 
sufficient for our purposes. Here, $R \equiv b^{M - 1} \equiv 1$ (mod $M$) is used to 
establish probable primality, using Fermat’s little theorem. Since $M = 2^p - 1$, evaluating $R$ 
requires $p$ squarings of the base $b$ followed by a modular division of $b^2$ to reach the 
Fermat-PRP result, also described as value $A$ of the cofactor test, which was devised by Hiromi 
Suyama in 1984 and subsequently improved upon by Hendrik W. Lenstra, Jr.

### Suyama cofactor test
Suyama’s method also utilises Fermat’s little theorem to test the cofactor. If a Fermat number $F$ 
or Mersenne number $M$ is composite and equal to $Q \times C$, where $Q$ is the product of known 
factors and $C$ is the cofactor, then having calculated either $A \equiv b^{F-1}$ (mod $F$) or 
$A \equiv b^{2^{p}-2}$ (mod $M$), we also calculate $B \equiv b^{Q-1}$ modulo $F$ or $M$ 
respectively; this then allows a comparison by simple subtraction, modulo the cofactor. If 
$$A - B \equiv 0\ (\text{mod}\ C)$$
then the cofactor is probably prime to the base $b^Q$; otherwise it is composite. Taking the 
greatest common divisor $\text{gcd}(A - B, C)$ also tests whether the cofactor is a 
prime power, divisible by the gcd.

## Computation
Most of the smaller calculations in `cofact` use the GNU Multiple Precision (GMP) arithmetic 
library, however for the heavy lifting of modular squarings required by the Pépin test or 
Fermat-PRP test, the `gwnum` library provides better, multi-threaded performance. `cofact` 
does not save interim results however, so for any Mersenne exponent greater than a million, 
a dedicated piece of software such as Prime95 or GPUowl will provide save and restart 
functionality, and may generate a proof which cofact can then utilise whenever new factors
are discovered.
### Proof files
Prime95 is capable of generating verifiable delay function proof files for any Mersenne 
exponent larger than $105,000$ (though the source code readily permits files to be generated 
for smaller exponents if desired). For Fermat exponents, proof files have been furnished for 
$F_{14}$ up to $F_{29}$, available from the co-author’s [website](https://64ordle.au/fermat/). 
If you wish to use Prime95 to generate a proof for $F_{30}$, we would be 
[most interested](https://www.mersenneforum.org/showthread.php?t=22668&page=6)
in knowing about it.
### Residues
The default output prints a hexadecimal residue modulo $2^{64}$, along with the triplet of 
smaller residues devised by Alexander Hurwitz and John Selfridge in 1964, which by default 
are given in decimal. Historical references or more modern software may reproduce these as 
octal or hexadecimal, both of which are supported (by `-o` and `-x`). A variety of options 
are available for printing interim residues while the program is computing the primality or 
PRP test, or the Suyama test (the `-i`, `-p`, and `-v` options may be of use).

The JSON reporting option for Mersenne numbers also includes the lower two kilobits of the 
residue as hexadecimal (modulo $2^{2,048}$).
## Full list of features / command line options

Long command line options use the same format for the parameter, if one is required (usually 
a filename, number, or string). The short, single-letter options may be combined in one flag, 
e.g. `-imo`, to enable multiple options `-i`, `-m`, and `-o`, but only one such option with 
a trailing parameter is permitted per flag.

Short option   | Long option          | Function
---------------|----------------------|---------
-a             |--all-residues        | Print residues for every modular squaring.
-b             |--binary              | Output final residues in binary.
-c _filename_  |--check-proof         | Check a VDF proof by computing the Fermat-PRP/Suyama $A$ residue for direct comparison. (This can be lengthy for large exponents; most often you will want to use `-u` or `--use-proof` below.)
-d             |--debug               | Print debug information.
-h             |--help                | Print basic help (`-hv` and `-h -sv` are increasingly verbose).
-i             |--interim-residues    | Print interim residues at various points.
-j             |--report-json         | Print a `JSON` report string for a Mersenne cofactor test, to allow submission of cofactor results. PrimeNet user and computer names may be entered into the string (see `-w` and `-q`). If a proof file is used, it will be verified to ensure the final residue can be correctly generated from the file.
-k             |--known-factors       | Use known factors of Fermat numbers (as of 2012) in place of supplying them after the exponent. When testing a Mersenne number in combination with checking or using a VDF proof, `cofact` can read the proof file’s description to import any known factors saved in the file header.
-m             |--mod-c               | Reduce Suyama $A$ and $B$ values, modulo $C$ and print residues.
-o             |--octal               | Print Selfridge–Hurwitz residues in octal as well as decimal.
-p _iterations_|--iterations          | Display progress every _iterations_ modular squarings (the default is 10% of a total run in excess of 100,000). If `-i` is also specified, residues will be printed at each progress point.
-q _string_    |--computer            | Supplies a PrimeNet computer name for reporting results (see `-j`).
-sep           |--separator           | Draw a horizontal line after a test.
-t _threads_   |--threads             | Use multi-threaded `gwnum` by specifying the number of _threads_ (the default is single-threaded).
-u _filename_  |--use-proof           | Use a VDF proof without double-checking the residue.
-v             |--verbose             | Print verbose test information.
-w _string_    |--user                | Supplies a PrimeNet username for reporting results (see `-j`).
-x             |--hex or --hexadecimal| Print Selfridge–Hurwitz residues in hexadecimal as well as decimal.
-y             |--mersenne            | Specify that the exponent is for a Mersenne number. This command line flag must immediately precede the exponent.
-z _base_      |--base                | Use a different base for primality testing (the default is 3).

## Future feature list
Some future nice things to consider adding: full Gerbicz error checking, proof generation, and giving everyone a unicorn.
