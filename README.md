# cofact
## A utility for testing the character of a Fermat or Mersenne cofactor
cofact tests a small Fermat number for primality, or optionally in version 0.9, a small Mersenne 
number for probable primality. (In this context, ‘small’ means, up to about 320 million digits in 
size!)

If known factors are given, then the cofactor is likewise tested for probable primality.

## Requirements
* Intel or compatible x86-family processor (AVX-512 required for the largest exponents)
* GNU Multiple Precision Arithmetic Library
* Copy of the Prime95/mprime source, including the `gwnum` large number library. 

## Installation
Unzip the `cofact` download as a directory  within the Prime95 source library directory `gwnum`, 
as `cofact` uses several files at the higher directory hierachies of the Prime95 source. First, 
run the `gwnum` make file for your platform, e.g. for macOS:
```bash
cd p95v3019b17/gwnum; make -f makemac
```
Then for `cofact`:
```bash
cd cofact; make
```
Install by simply copying cofact to your usual directory for binaries.
### Basic operation
To test a Fermat number, type `cofact` followed by an exponent (a non-negative integer up to 30), 
optionally followed by any factors of that Fermat number. For instance, to test the fifth Fermat 
number $F_5$ and the factor $641$:
```bash
cofact 5 641
```
The `cofact` distribution includes a verifiable delay function (VDF) proof for a middling Fermat 
number $F_{17}$, so to check the file has a correct value stored, type:
```bash
cofact -cpr F17.proof 17
```
This computes a primality test for $F_{17}$ and compares the result with the proof file. If you are 
satisfied you have a valid proof file the check may be skipped; to proceed directly to test the 
cofactor, by supplying known factors, type:
```bash
cofact -upr F17.proof 17 31065037602817 7751061099802522589358967058392886922693580423169
```
You may check or use any proof file with or without factors.
## Command line options
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
Gary B. Gostin, versions 0.2 to 0.8.1

Catherine X. Cowie, versions 0.7 to 0.9

## Copyright
This program is © 2023–2024 Gostin and Cowie under the Creative Commons Zero (CC-0) licence, we 
will not be responsible for whether you find this work useful, or whether it opens a portal to 
another hellish dystopian dimension where your computer used to be.

The gwnum library, and the proof validation module is © 2002–24 Mersenne Research, Inc, used 
with permission. All rights reserved.

The GMP library is © 1991, 1993–2016, 2018–2024 Free Software Foundation, Inc.

## Mathematics
The Fermat numbers have the form $2^{2^m} +1$, $m \ge 0$, and the Mersenne numbers have the form 
$2^p - 1$, where $p$ is a prime number. The five smallest Fermat numbers are prime, and of the 
first ten million prime numbers, 51 prime exponents _p_ are currently known to give rise to a 
Mersenne prime. This utility can run a primality or probable primality test on any of these 
numbers in such a way that the test allows further testing of the cofactor, if the Fermat or 
Mersenne number is composite and has one or more known factors.

The limit of this sort of testing is reached around $F_{30}$ for the Fermat numbers and 
$M_{1,073,741,789}$ for Mersenne numbers.

If $F$ is a Fermat number, we are interested in obtaining the result $P$ of a Pépin test, where 
$P \equiv b^{\frac 1 2 (F - 1)}$ (mod $F$) $\equiv -1$ (mod $F$) if and only if $F$ is prime. 
Since $F = 2^{2^m} + 1$ is the result of $m$ squarings of the number $2$ and adding $1$, $P$ may 
be obtained by squaring the base $2^m-1$ times, modulo $F$. Pépin’s theorem is a definitive 
primality test for bases such as 3, 5, 6, and so on. One further squaring follows to prepare for 
the cofactor test below.

A similar test for Mersenne numbers is not definitive (a different test, the Lucas–Lehmer test is 
used for that) but serves to determine whether the Mersenne number is _probably_ prime, which is 
sufficient for our purposes. Here, $R \equiv b^{M - 1}$ (mod $M$) $\equiv 1$ (mod $M$) is used to 
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
greatest common divisor $\textrm{gcd}(A - B, C)$ also tests whether the cofactor is a 
prime power, divisible by the gcd.

## Computation
Most of the smaller calculations in `cofact` use the GNU Multiple Precision (GMP) arithmetic 
library, however for the heavy lifting of modular squarings required by the Pépin test or 
Fermat-PRP test, the `gwnum` library provides better, multi-threaded performance. `cofact` 
does not save interim results however, so for any Mersenne exponent greater than a million, 
a dedicated piece of software such as Prime95 or GPUowl will provide save and restart 
functionality, and may generate a proof which cofact can then utilise.

## Full list of features / command line options

Long command line options use the same format for the parameter, if one is required (usually a filename, number, or string).

Short option   | Long option          | Function
---------------|----------------------|---------
-a             |--all-residues        | Print residues for every modular squaring.
-b             |--binary              | Output final residues in binary.
-c _filename_  |--check-proof         | Check a VDF proof by computing the Fermat-PRP/Suyama A residue.
-d             |--debug               | Print debug information.
-h             |--help                | Print basic help (-hv and -h -sv are increasingly verbose).
-i             |--interim-residues    | Print interim residues at various points.
-j             |--report-json         | Print a `JSON` report string for a Mersenne cofactor test. User and computer names can be entered into the string (see `-w` and `-q`).
-k             |--known-factors       | Use known factors of Fermat numbers (as of 2012) in place of supplying them after the exponent.
-m             |--mod-c               | Reduce Suyama $A$ and $B$ values, modulo $C$ and print residues.
-o             |--octal               | Print Selfridge–Hurwitz residues in octal as well as decimal.
-p _iterations_|--iterations          | Display progress every _iterations_ modular squarings (the default is 10% of the total run).
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