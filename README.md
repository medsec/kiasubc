# Code for validating boomerang probabilities on Kiasu-BC

## Description
Contains three implementations for validating the 6- and 7-round boomerang
distinguishers on Kiasu-BC and a downscaled version Mini-Kiasu-BC and 
round-reduced versions thereof.

Original publication:

https://eprint.iacr.org/2016/1170

## Working Manner
Since the 7-round distinguisher is infeasible to test with Kiasu-BC (and slow
for the 6-round distinguisher), we defined Mini-Kiasu-BC is as a 
nibblewise-operating variant of Kiasu-BC that employs the same high-level
structure as Kiasu-BC in downscaled manner, i.e., the same number and order of
operations, equal number of rounds and key schedule, the same ShiftRows, AddKey,
and AddTweak operations, the same MDS matrix, though, with multiplica- tions in
GF(2^4) under the irreducible polynomial x^4 + x + 1, operating on nibbles
instead of bytes, and with the Piccolo S-box instead of that of the AES.

For the full number of rounds, the implementations simply follow the trails and
build structures. For round-reduced versions of the trails, they start from the
state after the first round.

The tests 
(1) chose a random base plaintext, 
(2) create plaintexts by iterating over the values of the first column and all
possible values of the first tweak byte,
(3) encrypt the resulting plaintext-tweak pairs over the remaining rounds, 
(4) apply a delta-shift, and 
(5) decrypt 3, 4, 5, or 6 rounds backwards. 

The base plaintexts of each structure and the keys is chosen from /dev/urandom.

## Building
Simply run the make file in the individual subdirectories. There are several 
targets, among the important ones are:

- `make/make boomerang-differentials`
  Compile boomerang test program

- `make test-kiasu`
  Compiles test code for the Kiasu-BC/Mini-Kiasu-BC implementations.

- `make clean`
  Cleans executables

## Dependencies
- clang
- make
