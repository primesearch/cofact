// Copyright 2020 Mihai Preda and George Woltman.

// This routine, written by Mihai Preda, calculates the small roots of unity of a Mersenne number.
// This is used to thwart a roots-of-unity attack on a PRP proof where the final result is multiplied
// by a root to disguise the discovery of a new Mersenne prime.

typedef unsigned char u8;
typedef unsigned int u32;

// "Small primes up to 2^15": gaps between successive primes.
u8 gaps[] = {
2,1,2,2,4,2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,2,4,14,4,6,2,10,2,6,6,4,6,6,2,10,2,4,2,12,12,4,2,4,6,2,10,6,6,6,2,6,4,2,10,14,4,2,4,14,6,10,2,4,6,8,6,6,4,6,8,4,8,10,2,10,2,6,4,6,8,4,2,4,12,8,4,8,4,6,12,2,18,6,10,6,6,2,6,10,6,6,2,6,6,4,2,12,10,2,4,6,6,
2,12,4,6,8,10,8,10,8,6,6,4,8,6,4,8,4,14,10,12,2,10,2,4,2,10,14,4,2,4,14,4,2,4,20,4,8,10,8,4,6,6,14,4,6,6,8,6,12,4,6,2,10,2,6,10,2,10,2,6,18,4,2,4,6,6,8,6,6,22,2,10,8,10,6,6,8,12,4,6,6,2,6,12,10,18,2,4,6,2,6,4,2,4,12,2,6,34,6,6,8,18,10,14,4,2,4,6,8,4,2,6,12,
10,2,4,2,4,6,12,12,8,12,6,4,6,8,4,8,4,14,4,6,2,4,6,2,6,10,20,6,4,2,24,4,2,10,12,2,10,8,6,6,6,18,6,4,2,12,10,12,8,16,14,6,4,2,4,2,10,12,6,6,18,2,16,2,22,6,8,6,4,2,4,8,6,10,2,10,14,10,6,12,2,4,2,10,12,2,16,2,6,4,2,10,8,18,24,4,6,8,16,2,4,8,16,2,4,8,6,6,4,12,
2,22,6,2,6,4,6,14,6,4,2,6,4,6,12,6,6,14,4,6,12,8,6,4,26,18,10,8,4,6,2,6,22,12,2,16,8,4,12,14,10,2,4,8,6,6,4,2,4,6,8,4,2,6,10,2,10,8,4,14,10,12,2,6,4,2,16,14,4,6,8,6,4,18,8,10,6,6,8,10,12,14,4,6,6,2,28,2,10,8,4,14,4,8,12,6,12,4,6,20,10,2,16,26,4,2,12,6,4,12,
6,8,4,8,22,2,4,2,12,28,2,6,6,6,4,6,2,12,4,12,2,10,2,16,2,16,6,20,16,8,4,2,4,2,22,8,12,6,10,2,4,6,2,6,10,2,12,10,2,10,14,6,4,6,8,6,6,16,12,2,4,14,6,4,8,10,8,6,6,22,6,2,10,14,4,6,18,2,10,14,4,2,10,14,4,8,18,4,6,2,4,6,2,12,4,20,22,12,2,4,6,6,2,6,22,2,6,16,6,
12,2,6,12,16,2,4,6,14,4,2,18,24,10,6,2,10,2,10,2,10,6,2,10,2,10,6,8,30,10,2,10,8,6,10,18,6,12,12,2,18,6,4,6,6,18,2,10,14,6,4,2,4,24,2,12,6,16,8,6,6,18,16,2,4,6,2,6,6,10,6,12,12,18,2,6,4,18,8,24,4,2,4,6,2,12,4,14,30,10,6,12,14,6,10,12,2,4,6,8,6,10,2,4,14,6,6,
4,6,2,10,2,16,12,8,18,4,6,12,2,6,6,6,28,6,14,4,8,10,8,12,18,4,2,4,24,12,6,2,16,6,6,14,10,14,4,30,6,6,6,8,6,4,2,12,6,4,2,6,22,6,2,4,18,2,4,12,2,6,4,26,6,6,4,8,10,32,16,2,6,4,2,4,2,10,14,6,4,8,10,6,20,4,2,6,30,4,8,10,6,6,8,6,12,4,6,2,6,4,6,2,10,2,16,6,20,4,12,
14,28,6,20,4,18,8,6,4,6,14,6,6,10,2,10,12,8,10,2,10,8,12,10,24,2,4,8,6,4,8,18,10,6,6,2,6,10,12,2,10,6,6,6,8,6,10,6,2,6,6,6,10,8,24,6,22,2,18,4,8,10,30,8,18,4,2,10,6,2,6,4,18,8,12,18,16,6,2,12,6,10,2,10,2,6,10,14,4,24,2,16,2,10,2,10,20,4,2,4,8,16,6,6,2,12,16,
8,4,6,30,2,10,2,6,4,6,6,8,6,4,12,6,8,12,4,14,12,10,24,6,12,6,2,22,8,18,10,6,14,4,2,6,10,8,6,4,6,30,14,10,2,12,10,2,16,2,18,24,18,6,16,18,6,2,18,4,6,2,10,8,10,6,6,8,4,6,2,10,2,12,4,6,6,2,12,4,14,18,4,6,20,4,8,6,4,8,4,14,6,4,14,12,4,2,30,4,24,6,6,12,12,14,6,4,
2,4,18,6,12,8,6,4,12,2,12,30,16,2,6,22,14,6,10,12,6,2,4,8,10,6,6,24,14,6,4,8,12,18,10,2,10,2,4,6,20,6,4,14,4,2,4,14,6,12,24,10,6,8,10,2,30,4,6,2,12,4,14,6,34,12,8,6,10,2,4,20,10,8,16,2,10,14,4,2,12,6,16,6,8,4,8,4,6,8,6,6,12,6,4,6,6,8,18,4,20,4,12,2,10,6,2,
10,12,2,4,20,6,30,6,4,8,10,12,6,2,28,2,6,4,2,16,12,2,6,10,8,24,12,6,18,6,4,14,6,4,12,8,6,12,4,6,12,6,12,2,16,20,4,2,10,18,8,4,14,4,2,6,22,6,14,6,6,10,6,2,10,2,4,2,22,2,4,6,6,12,6,14,10,12,6,8,4,36,14,12,6,4,6,2,12,6,12,16,2,10,8,22,2,12,6,4,6,18,2,12,6,4,12,
8,6,12,4,6,12,6,2,12,12,4,14,6,16,6,2,10,8,18,6,34,2,28,2,22,6,2,10,12,2,6,4,8,22,6,2,10,8,4,6,8,4,12,18,12,20,4,6,6,8,4,2,16,12,2,10,8,10,2,4,6,14,12,22,8,28,2,4,20,4,2,4,14,10,12,2,12,16,2,28,8,22,8,4,6,6,14,4,8,12,6,6,4,20,4,18,2,12,6,4,6,14,18,10,8,10,
32,6,10,6,6,2,6,16,6,2,12,6,28,2,10,8,16,6,8,6,10,24,20,10,2,10,2,12,4,6,20,4,2,12,18,10,2,10,2,4,20,16,26,4,8,6,4,12,6,8,12,12,6,4,8,22,2,16,14,10,6,12,12,14,6,4,20,4,12,6,2,6,6,16,8,22,2,28,8,6,4,20,4,12,24,20,4,8,10,2,16,2,12,12,34,2,4,6,12,6,6,8,6,4,2,6,
24,4,20,10,6,6,14,4,6,6,2,12,6,10,2,10,6,20,4,26,4,2,6,22,2,24,4,6,2,4,6,24,6,8,4,2,34,6,8,16,12,2,10,2,10,6,8,4,8,12,22,6,14,4,26,4,2,12,10,8,4,8,12,4,14,6,16,6,8,4,6,6,8,6,10,12,2,6,6,16,8,6,6,12,10,2,6,18,4,6,6,6,12,18,8,6,10,8,18,4,14,6,18,10,8,10,12,2,
6,12,12,36,4,6,8,4,6,2,4,18,12,6,8,6,6,4,18,2,4,2,24,4,6,6,14,30,6,4,6,12,6,20,4,8,4,8,6,6,4,30,2,10,12,8,10,8,24,6,12,4,14,4,6,2,28,14,16,2,12,6,4,20,10,6,6,6,8,10,12,14,10,14,16,14,10,14,6,16,6,8,6,16,20,10,2,6,4,2,4,12,2,10,2,6,22,6,2,4,18,8,10,8,22,2,10,
18,14,4,2,4,18,2,4,6,8,10,2,30,4,30,2,10,2,18,4,18,6,14,10,2,4,20,36,6,4,6,14,4,20,10,14,22,6,2,30,12,10,18,2,4,14,6,22,18,2,12,6,4,8,4,8,6,10,2,12,18,10,14,16,14,4,6,6,2,6,4,2,28,2,28,6,2,4,6,14,4,12,14,16,14,4,6,8,6,4,6,6,6,8,4,8,4,14,16,8,6,4,12,8,16,2,
10,8,4,6,26,6,10,8,4,6,12,14,30,4,14,22,8,12,4,6,8,10,6,14,10,6,2,10,12,12,14,6,6,18,10,6,8,18,4,6,2,6,10,2,10,8,6,6,10,2,18,10,2,12,4,6,8,10,12,14,12,4,8,10,6,6,20,4,14,16,14,10,8,10,12,2,18,6,12,10,12,2,4,2,12,6,4,8,4,44,4,2,4,2,10,12,6,6,14,4,6,6,6,8,6,
36,18,4,6,2,12,6,6,6,4,14,22,12,2,18,10,6,26,24,4,2,4,2,4,14,4,6,6,8,16,12,2,42,4,2,4,24,6,6,2,18,4,14,6,28,18,14,6,10,12,2,6,12,30,6,4,6,6,14,4,2,24,4,6,6,26,10,18,6,8,6,6,30,4,12,12,2,16,2,6,4,12,18,2,6,4,26,12,6,12,4,24,24,12,6,2,12,28,8,4,6,12,2,18,6,4,
6,6,20,16,2,6,6,18,10,6,2,4,8,6,6,24,16,6,8,10,6,14,22,8,16,6,2,12,4,2,22,8,18,34,2,6,18,4,6,6,8,10,8,18,6,4,2,4,8,16,2,12,12,6,18,4,6,6,6,2,6,12,10,20,12,18,4,6,2,16,2,10,14,4,30,2,10,12,2,24,6,16,8,10,2,12,22,6,2,16,20,10,2,12,12,18,10,12,6,2,10,2,6,10,18,
2,12,6,4,6,2,24,28,2,4,2,10,2,16,12,8,22,2,6,4,2,10,6,20,12,10,8,12,6,6,6,4,18,2,4,12,18,2,12,6,4,2,16,12,12,14,4,8,18,4,12,14,6,6,4,8,6,4,20,12,10,14,4,2,16,2,12,30,4,6,24,20,24,10,8,12,10,12,6,12,12,6,8,16,14,6,4,6,36,20,10,30,12,2,4,2,28,12,14,6,22,8,4,
18,6,14,18,4,6,2,6,34,18,2,16,6,18,2,24,4,2,6,12,6,12,10,8,6,16,12,8,10,14,40,6,2,6,4,12,14,4,2,4,2,4,8,6,10,6,6,2,6,6,6,12,6,24,10,2,10,6,12,6,6,14,6,6,52,20,6,10,2,10,8,10,12,12,2,6,4,14,16,8,12,6,22,2,10,8,6,22,2,22,6,8,10,12,12,2,10,6,12,2,4,14,10,2,6,18,
4,12,8,18,12,6,6,4,6,6,14,4,2,12,12,4,6,18,18,12,2,16,12,8,18,10,26,4,6,8,6,6,4,2,10,20,4,6,8,4,20,10,2,34,2,4,24,2,12,12,10,6,2,12,30,6,12,16,12,2,22,18,12,14,10,2,12,12,4,2,4,6,12,2,16,18,2,40,8,16,6,8,10,2,4,18,8,10,8,12,4,18,2,18,10,2,4,2,4,8,28,2,6,22,12,
6,14,18,4,6,8,6,6,10,8,4,2,18,10,6,20,22,8,6,30,4,2,4,18,6,30,2,4,8,6,4,6,12,14,34,14,6,4,2,6,4,14,4,2,6,28,2,4,6,8,10,2,10,2,10,2,4,30,2,12,12,10,18,12,14,10,2,12,6,10,6,14,12,4,14,4,18,2,10,8,4,8,10,12,18,18,8,6,18,16,14,6,6,10,14,4,6,2,12,12,4,6,6,12,2,16,
2,12,6,4,14,6,4,2,12,18,4,36,18,12,12,2,4,2,4,8,12,4,36,6,18,2,12,10,6,12,24,8,6,6,16,12,2,18,10,20,10,2,6,18,4,2,40,6,2,16,2,4,8,18,10,12,6,2,10,8,4,6,12,2,10,18,8,6,4,20,4,6,36,6,2,10,6,24,6,14,16,6,18,2,10,20,10,8,6,4,6,2,10,2,12,4,2,4,8,10,6,12,18,14,12,
16,8,6,16,8,4,2,6,18,24,18,10,12,2,4,14,10,6,6,6,18,12,2,28,18,14,16,12,14,24,12,22,6,2,10,8,4,2,4,14,12,6,4,6,14,4,2,4,30,6,2,6,10,2,30,22,2,4,6,8,6,6,16,12,12,6,8,4,2,24,12,4,6,8,6,6,10,2,6,12,28,14,6,4,12,8,6,12,4,6,14,6,12,10,6,6,8,6,6,4,2,4,8,12,4,14,18,
10,2,16,6,20,6,10,8,4,30,36,12,8,22,12,2,6,12,16,6,6,2,18,4,26,4,8,18,10,8,10,6,14,4,20,22,18,12,8,28,12,6,6,8,6,12,24,16,14,4,14,12,6,10,12,20,6,4,8,18,12,18,10,2,4,20,10,14,4,6,2,10,24,18,2,4,20,16,14,10,14,6,4,6,20,6,10,6,2,12,6,30,10,8,6,4,6,8,40,2,4,2,12,
18,4,6,8,10,6,18,18,2,12,16,8,6,4,6,6,2,52,14,4,20,16,2,4,6,12,2,6,12,12,6,4,14,10,6,6,14,10,14,16,8,6,12,4,8,22,6,2,18,22,6,2,18,6,16,14,10,6,12,2,6,4,8,18,12,16,2,4,14,4,8,12,12,30,16,8,4,2,6,22,12,8,10,6,6,6,14,6,18,10,12,2,10,2,4,26,4,12,8,4,18,8,10,14,16,
6,6,8,10,6,8,6,12,10,20,10,8,4,12,26,18,4,12,18,6,30,6,8,6,22,12,2,4,6,6,2,10,2,4,6,6,2,6,22,18,6,18,12,8,12,6,10,12,2,16,2,10,2,10,18,6,20,4,2,6,22,6,6,18,6,14,12,16,2,6,6,4,14,12,4,2,18,16,36,12,6,14,28,2,12,6,12,6,4,2,16,30,8,24,6,30,10,2,18,4,6,12,8,22,2,6,
22,18,2,10,2,10,30,2,28,6,14,16,6,20,16,2,6,4,32,4,2,4,6,2,12,4,6,6,12,2,6,4,6,8,6,4,20,4,32,10,8,16,2,22,2,4,6,8,6,16,14,4,18,8,4,20,6,12,12,6,10,2,10,2,12,28,12,18,2,18,10,8,10,48,2,4,6,8,10,2,10,30,2,36,6,10,6,2,18,4,6,8,16,14,16,6,14,4,20,4,6,2,10,12,2,
6,12,6,6,4,12,2,6,4,12,6,8,4,2,6,18,10,6,8,12,6,22,2,6,12,18,4,14,6,4,20,6,16,8,4,8,22,8,12,6,6,16,12,18,30,8,4,2,4,6,26,4,14,24,22,6,2,6,10,6,14,6,6,12,10,6,2,12,10,12,8,18,18,10,6,8,16,6,6,8,16,20,4,2,10,2,10,12,6,8,6,10,20,10,18,26,4,6,30,2,4,8,6,12,12,
18,4,8,22,6,2,12,34,6,18,12,6,2,28,14,16,14,4,14,12,4,6,6,2,36,4,6,20,12,24,6,22,2,16,18,12,12,18,2,6,6,6,4,6,14,4,2,22,8,12,6,10,6,8,12,18,12,6,10,2,22,14,6,6,4,18,6,20,22,2,12,24,4,18,18,2,22,2,4,12,8,12,10,14,4,2,18,16,38,6,6,6,12,10,6,12,8,6,4,6,14,30,6,
10,8,22,6,8,12,10,2,10,2,6,10,2,10,12,18,20,6,4,8,22,6,6,30,6,14,6,12,12,6,10,2,10,30,2,16,8,4,2,6,18,4,2,6,4,26,4,8,6,10,2,4,6,8,4,6,30,12,2,6,6,4,20,22,8,4,2,4,72,8,4,8,22,2,4,14,10,2,4,20,6,10,18,6,20,16,6,8,6,4,20,12,22,2,4,2,12,10,18,2,22,6,18,30,2,10,
14,10,8,16,50,6,10,8,10,12,6,18,2,22,6,2,4,6,8,6,6,10,18,2,22,2,16,14,10,6,2,12,10,20,4,14,6,4,36,2,4,6,12,2,4,14,12,6,4,6,2,6,4,20,10,2,10,6,12,2,24,12,12,6,6,4,24,2,4,24,2,6,4,6,8,16,6,2,10,12,14,6,34,6,14,6,4,2,30};

const u32 GAPS_SIZE = sizeof(gaps) / sizeof(gaps[0]);

// This is a compacted table computed in pari-gp of z(p) = znorder(Mod(2,p)), associating for each "z" entry
// the set of primes "p" that have that z(p) value.
// The keys "z" are stored as negative values to distinguish them from the primes that follow each z.
int merged[] = {
#include "rootsdata.txt"
};

const u32 MERGED_SIZE = sizeof(merged) / sizeof(merged[0]);

// find the first negative entry under pos.
int locateZ(u32 pos) {
  assert(pos < MERGED_SIZE);
  while (merged[pos] >= 0) { --pos; }
  assert(merged[pos] < 0);
  return pos;
}

// Search entry for "z" in the "merged" table. Binary search, the table is ordered by "z".
int searchZ(u32 z) {
  u32 a = 0;
  u32 b = locateZ(MERGED_SIZE - 1);

  assert(z >= -merged[a]);
  if (z == -merged[a]) { return a; }
  if (z > -merged[b]) { return -1; }
  if (z == -merged[b]) { return b; }

  assert(a < b && z > -merged[a] && z < -merged[b]);

  // binary search
  while (1) {
    u32 mid = locateZ((a + b) / 2);
    if (mid == a) {
      mid = locateZ(b - 1);
      if (mid == a) { return -1; }
    }
    assert(a < mid && mid < b);
    if (z < -merged[mid]) {
      b = mid;
    } else if (z > -merged[mid]) {
      a = mid;
    } else {
      return mid;
    }
  }
}

// 2400 divisors should be enough for everybody :)
#define DIV_SIZE 2401
u32 divisors[DIV_SIZE] = {1};

u32 findDivisors(u32 n) {
  u32 nDivs = 1;

  u32 prime = 0;
  for (u8 *it = gaps, *end = gaps + GAPS_SIZE; n > 1 && it < end; ++it) {
    prime += *it;

    // as soon as we pass sqrt(n), treat the rest as one final factor.
    if (prime * prime > n) { prime = n; }

    u32 multiplicity = 0;
    while (n % prime == 0) {
      ++multiplicity;      
      n /= prime;
    }

    // printf("%u %u\n", prime, multiplicity);

    u32 oldNDivs = nDivs;
    for (u32 m=0, p=prime; m < multiplicity; ++m, p *= prime) {
      for (int i = 0; i < oldNDivs; ++i) {
        assert(nDivs < DIV_SIZE);
        if (nDivs < DIV_SIZE) { divisors[nDivs++] = divisors[i] * p; }
      }
    }
  }
  return nDivs;
}

// small prime factors of 2^n - 1 with multiplicity.
void mersenne_roots(u32 n, mpz_t accumulator) {
  u32 nDivs = findDivisors(n);
  int i, j, k;
  double bits = 0, bitsSmall = 0;

  // skip the first divisor that's always "1".
  for (i = 1; i < nDivs; ++i) {
    u32 d = divisors[i];
    int zPos = searchZ(d);
    if (zPos >= 0) {
      for (j = zPos + 1; j < MERGED_SIZE && merged[j] >= 0; ++j) {
        u32 prime = merged[j];
        u32 m = 1;
        u32 left = n / d;
        while (left % prime == 0) {
          ++m;
          left /= prime;
        }

	// HERE! accumulate the factors with multiplicity into a huge exponent.
	for (k = 0; k < m; k++) mpz_mul_ui (accumulator, accumulator, prime);
      }
    }
  }
}

