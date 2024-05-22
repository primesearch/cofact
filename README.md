# cofact

cofact is a companion program for the mprime / Prime95 program created by George Woltman. mprime and cofact can be used together in several ways to test a small Fermat number for primality using Pepin's test and, if known factors are provided, to test the Fermat cofactor for probable primality using the Suyama test. cofact reports the RES64 and Selfridge-Hurwitz residues for the final residue of Pepin's test and for each step of the Suyama test to enable comparison with the residues reported by other programs.

Starting in cofact version 0.9, cofact can perform the same tests on a small Mersenne number and the Mersenne cofactor.

## Motivation and History

When running very long calculations such as testing a number for primality, one approach to ensure accuracy is to perform the calculation twice using two different programs running on two different computers and then compare various calculation residues. This helps protect against both program design errors and hardware errors. 

Starting in 2013, Ernst Mayer added the ability to his Mlucas program to test a Fermat number for primality using Pepin's test, saving the final Pepin residue to a file. Then, if the Fermat number has known factors, Mlucas can use the final Pepin residue to test the Fermat cofactor for probable primality using the Suyama test. From 2013 to 2022 he tested each of the Fermat numbers from F25 through F30 for primality, a huge effort particularly for F30. Once complete, he ran the Suyama test on each of the Fermat cofactors, proving each of them to be composite. For each Fermat number, Mlucas reported the RES64 residue and the Selfridge-Hurwitz residues (see Residues below) for the final Pepin residue and for each step of the Suyama test: A, B and (A-B) mod C (see Mathematics below).

Similar to Mlucas, mprime / Prime95 has the ability to test a Fermat cofactor for primality using the Pepin and Suyama tests. Unlike Mlucas, mprime reports only the RES64 residue for the Suyama A step. So it was not possible with mprime alone to verify all the residues reported by Mlucas. Fortunately, mprime stores the full Suyama A step residue in the generated proof file. So, to verify all of the Mlucas residues, cofact was created to read the full Suyama A residue from the mprime proof file, perform the Suyama test, and then report all the same residues as Mlucas. As of November 2023, mprime + cofact have been used to verify all Mlucas Suyama residues for $F_{21}$ through $F_{23}$ and $F_{25}$ through $F_{29}$. 

## Repository Branches
There are two permanent branches in this cofact repository:

* The "main" branch is the original cofact program used with mprime to verify the Mlucas residues. Versions through 0.6 were written by Gary Gostin, while versions 0.7 and beyond were in colaboration with Catherine Cowie. This branch is currently at version 0.8.2 and is considered a "finished product" (no future changes are currently planned).

* The "cxc" branch is an evolution of cofact starting from version 0.6 created by Catherine Cowie. It includes many extensions to the original cofact including the ability to test Mersenne numbers. Development of this branch is still ongoing.

Using `git`, you may easily switch between the two branches using: `git switch [main|cxc]`.

The remainder of this README.md file will document the "main" branch. For the "cxc" branch, see the README.md file in that branch for documentation. 


## Requirements
* Intel or compatible x86-family processor. AVX-512 is required to perform Pepin's test on $F_{30}$.
* GNU Multiple Precision Arithmetic Library.
* The mprime/Prime95 source code, including the `gwnum` large number library.

## Installation

Note: The "main" branch of cofact assumes the program is being built on a Linux system. To build for other OSes, please see the "cxc" branch. 

First, download and install the GMP library following the documentation at gmplib.org. The cofact Makefile assumes that "make install" was run for GMP.

Next, download the mprime/Prime95 source code from mersenne.org/download. Version 30.8 or beyond is required for cofact. Then, either build the entire mprime program or build only the gwnum library for your platform \[macOS|Linux|Win\] using:
```
cd <mprime root directory>/gwnum
make -f [makemac|make64|makemw64]
```

Finally, either download and unzip cofact as a directory or, if you have git installed, create a "cofact" directory then run the following in that directory:
```
git clone --branch main https://github.com/primesearch/cofact.git
```

The cofact Makefile uses a number of files from the mprime gwnum directory. To create links to those files and then build cofact, run:
```
link_gwnum <mprime root directory>
make
```

## Basic operation

cofact can be run in one of three modes:

1. Test a Fermat number for primality using Pepin's test. Then, if known factors are provided, use the Pepin residue to perform the Suyama probable primality (PRP) test on the cofactor. This mode is selected if neither -cpr or -upr are specified on the command line.
2. Perform all steps in mode 1. Also compare the Suyama A residue calculated by cofact to the A residue read from the proof file generated by mprime when testing the same cofactor. This mode is selected via -cpr on the command line.
3. Read the Suyama A residue for a Fermat number from the mprime proof file (assumes it is correct), then perform the Suyama PRP test on the cofactor. This mode is selected via -upr on the command line.

To test a Fermat number in mode 1, type `cofact` followed by a Fermat exponent (a non-negative integer up through 30), optionally followed by any factors of that Fermat number. For instance, to test the fifth Fermat number $F_5$ using the known factor 641:
```bash
cofact 5 641
```
cofact will report that $F_5$ is composite and the $F_5$ cofactor is probably prime.

Mode 2 can be used to verify the Suyama A residue stored in a verifiable delay function (VDF) proof file generated by mprime. For example, the cofact distribution includes a proof file for the small Fermat number $F_{14}$. To check that the file has a correct value stored, type:
```bash
cofact -cpr F14.proof 14
```
This performs a Pepin primality test for $F_{14}$, then calculates the Suyama A residue and compares the with the residue stored in the proof file. If no residue mismatch is reported, then the proof file has a valid residue. You may check any proof file with or without factors.

Mode 3 can be used to test a Fermat cofactor using the Suyama A residue stored in the corresponding proof file. For example, to test the $F_{14}$ cofactor using the $F_{14}$ proof file, type:
```bash
cofact -upr F14.proof 14 116928085873074369829035993834596371340386703423373313
```
Mode 3 avoids the lengthy Pepin calculation, allowing a cofactor to be tested in a matter of minutes even for Fermat numbers as large as $F_{29}$. To enable yourself to test the resulting cofactor after the next factor of $F_{14}$ through $F_{29}$ is discovered, download the $F_{14}$ through $F_{29}$ proof files from Catherine's excellent [website](https://64ordle.au/fermat/).

The cofact distribution includes a script called `run_all` that will run cofact on each Fermat number from $F_0$ through $F_{29}$ using the best method for that number, with "best" meaning reasonably fast. Before running the script, download proof files for $F_{17}$ through $F_{29}$ into the same directory as cofact.

## Command line options
The following command line options are supported by cofact (main branch):
Command line option | Function
--------------------|------------------------------
-cpr _file_         | Read the Suyama A residue from the mprime proof file and compare it to the A residue calculated by cofact (mode 2)
-d                  | Print debug information
-h                  | Print this help and exit
-p _iter_           | Print progress every iter iterations, instead of the default of every 10% of total iterations for longer runs
-sep                | Print a separator at the end of the run to better see multiple run's output in a single output file
-t _threads_        | Specifies the number of threads to use in the gwnum library. Defaults to 1.
-upr _file_         | Read the Suyama A residue from the mprime proof file and use it to complete the Suyama test (mode 3)
-v                  | Print more verbose information

## Authors
Gary B. Gostin (gary641), versions 0.2 to 0.8.2 (the original, `main` branch)

Catherine X. Cowie (xanthe-cat), versions 0.6 to 0.9 (the `cxc` branch)


## Copyright
This program is © 2023–2024 Gostin and Cowie under the Creative Commons Zero (CC-0) licence. We will not be responsible for whether you find this work useful, or whether it opens a portal to another hellish dystopian dimension where your computer used to be.

The gwnum library is © 2002–24 Mersenne Research, Inc, used with permission. All rights reserved.

The GMP library is © 1991, 1993–2016, 2018–2024 Free Software Foundation, Inc.

## Mathematics
A Fermat number is an integer of the form $2^{2^n}+1$, where $n \ge 0$. The five smallest Fermat numbers, $F_0$ through $F_4$, are prime. $F_5$ through $F_{11}$ have been fully factored. $F_{20}$ and $F_{24}$ have been proven composite but not factored, while many other Fermat numbers have known factors.

A Fermat cofactor is the integer $F_n / P$, where $P$ is the product of all known factors of $F_n$. 

### Pépin's test
Pépin's test is a definitive primality test that can be used to determine if a Fermat number is prime. If $F_n$ is the $n$th Fermat number, then $F_n$ is prime if and only if $3^{(F_n - 1)/2} \equiv -1$ (mod $F$). The term $3^{(F_n - 1)/2}$ (mod $F$) may be calculated by $2^n-1$ squarings of the base $3$, modulo $F_n$. Pépin’s test is also a definitive primality test when performed with bases other than 3, such as 5, 6, and so on.

### Suyama cofactor test
The Suyama test utilises Fermat’s little theorem to test a Fermat number's cofactor for probably primality. 
It was devised by Hiromi Suyama in 1984 and subsequently improved upon by Hendrik W. Lenstra, Jr.
If a Fermat number $F_n$ is composite and equal to $P \times C$, where $P$ is the product of known factors and $C$ is the cofactor, then the Suyama test is performed as follows:

* Calculate $A = b^{F_n-1}$ (mod $F_n$) for some base $b$. If base = 3 then this is the square mod $F_n$ of the final Pepin residue.
* Calculate $B = b^{P-1}$ (mod $F_n$).
* Calculate $R = (A - B)$ (mod $C$). If $R = 0$ then the cofactor $C$ is probably prime; otherwise it is composite. 

Proving the cofactor composite requires finding only one base for which $R \neq 0$. cofact (main branch) only supports base $b = 3$ for the Suyama test. So far, this has been sufficient since each Fermat cofactor from $F_{12}$ through $F_{29}$ is currently composite. 

Finally, the cofactor $C$ can be tested to determine if it is a prime power by calculating the greatest common divisor $G = \text{gcd}(A - B, C)$. If $G = 1$ then the cofactor is not a prime power. If $G \neq 1$ then the cofactor is a prime power and is divisible by $G$. 

### Computation
The less performance-critical calculations in cofact use the GNU Multiple Precision (GMP) arithmetic library. However for the heavy lifting of modular squarings required by Pépin's test, the `gwnum` library is used since it is multi-threaded and therefore provides much higher performance. Even so, cofact (main branch) does not have Gerbitz error checking of calculations or save interim results to a file. So, for a Fermat cofactor test that is expected to run more than a few hours, it is preferable to first either generate the proof file using mprime / Prime95 (which has Gerbitz error checking and checkpoint / restart capability) or download the proof file from Catherine's [website](https://64ordle.au/fermat/). Once a Fermat number's proof file is in hand, `cofact -upr` can be used to test the new cofactor whenever a new factor of the Fermat number is discovered.

### Proof files
mprime / Prime95 is capable of generating a verifiable delay function proof files for any Fermat number through $F_{30}$. Proof files have been collected for $F_{14}$ through $F_{29}$ and are available from Catherine’s [website](https://64ordle.au/fermat/). If you wish to use mprime / Prime95 to generate a proof file for $F_{30}$ (requires a system with AVX512 support), we would be [most interested](https://www.mersenneforum.org/showthread.php?t=22668&page=6) in knowing about it!

### Residues
For each full residue from Pepin's test and each step of the Suyama test, cofact prints the residue mod $2^{64}$ in hexadecimal along with the triplet of 
smaller residues (mod $2^{35}-1$, mod $2^{36}-1$ and mod $2^{36}$) devised by Alexander Hurwitz and John Selfridge in 1964. These Selfridge-Hurwitz residues are printed in decimal and octal, to enable easier comparison with residues reported in historical references.
