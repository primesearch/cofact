/*----------------------------------------------------------------------
| Copyright 2020-2021 Mersenne Research, Inc.  All rights reserved
|
| This program verifies a PRP proof.
+---------------------------------------------------------------------*/

/* Additional modifications 2024 by Catherine X. Cowie to use proof verification in cofact. 
|  Basically, whenever you see lines of code doubled with one commented out, that will 
|  usually be one of my modifications. 
*/

//static const char JUNK[]="Copyright 2020 Mersenne Research, Inc. All rights reserved";

//#define NO_HWLOC
//#include "common.h"
//#include "time.h"
//#include "exponentiate.c"
//#include "proof_hash.c"
//#include "roots.c"

// Guess a partial proof's file size.  Proof files start as all zeros, so lop off trailing zeros in the proof file to get the file size.
uint64_t guess_filesize (
	FILE	*fd)
{
	uint64_t last_non_zero, file_offset, filesize;
	char	buf[4096];

	last_non_zero = 0;
	file_offset = 0;
	for ( ; ; ) {
		size_t i, bytes_read;
		bytes_read = fread (buf, 1, sizeof(buf), fd);
		if (bytes_read == 0) break;
		//scan for last zero
		for (i = 1; i <= bytes_read; i++) {
			if (buf[bytes_read-i]) {
				last_non_zero = file_offset + bytes_read - i + 1;	// last_non_zero is one-based
				break;
			}
		}
		file_offset += bytes_read;
	}

	// If more than the last 10 bytes are zero, return last_non_zero
	fseek (fd, 0, SEEK_END);
	filesize = ftell (fd);
	fseek (fd, 0, SEEK_SET);
	if (last_non_zero < filesize - 10) return (last_non_zero);
	else return (filesize);
}

// Read next residue from the proof file
void readResidue (
	gwhandle *gwdata,
	FILE	*fd,
	gwnum	x,			/* Returned residue as a gwnum */
	mpz_t	*mpz_res,		/* Returned residue as a GMP mpz_t */
	char	*res64)			/* Returned 64-bit residue */
{
	int	i, residue_size;
	uint32_t *array;		/* Array to contain the binary value */
	uint32_t arraylen;		/* Size of the array */

	// Allocate an array for binary value
	residue_size = divide_rounding_up ((int) ceil(gwdata->bit_length), 8);
	arraylen = divide_rounding_up (residue_size, 4);
	array = (uint32_t *) malloc (arraylen * sizeof(uint32_t));
	if (array == NULL) {
    	printf ("Malloc error\n");
		return;
	}
	array[arraylen-1] = 0;

	// Read the residue
	if (fread (array, 1, residue_size, fd) != residue_size) {
		printf ("Error reading residue from proof file\n");
		free (array);
	}

	// Convert from an array of bytes (LSB to MSB) to an array of uint32_t
	for (i = 0; i < (int) arraylen; i++) {
		uint32_t val;
		val = ((unsigned char *)&array[i])[3];
		val = (val << 8) + ((unsigned char *)&array[i])[2];
		val = (val << 8) + ((unsigned char *)&array[i])[1];
		val = (val << 8) + ((unsigned char *)&array[i])[0];
		array[i] = val;
	}

	// Convert binary to mpz
	if (mpz_res != NULL) {
		mpz_init (*mpz_res);
		mpz_import (*mpz_res, arraylen, -1, sizeof (array[0]), 0, 0, array);
	}

	// Convert binary to 64-bit residue
	if (res64 != NULL) {
		sprintf (res64, "%02X%02X%02X%02X", array[1] >> 24, (array[1] >> 16) & 0xFF, (array[1] >> 8) & 0xFF, array[1] & 0xFF);
		sprintf (res64+8, "%02X%02X%02X%02X", array[0] >> 24, (array[0] >> 16) & 0xFF, (array[0] >> 8) & 0xFF, array[0] & 0xFF);
	}

	// Convert binary to gwnum
	binarytogw (gwdata, array, arraylen, x);
	free (array);
}

// Convert a giant to a 64-bit residue
void giantToRes64 (
	giant	g,
	char	*res64)
{
	uint64_t tmp;

	tmp = (g->sign > 1) ? g->n[1] : 0;
	tmp <<= 32;
	tmp += (g->sign > 0) ? g->n[0] : 0;
	sprintf (res64, "%016" PRIX64, tmp);
}

// Generate next prime > 32 bits, < 49 bits
uint64_t nextprime48 (uint64_t r) {
	mpz_t	a;
	r = r & 0xFFFFFFFFFFFFULL;	// 48 bits
	if (r < 0x100000000ULL) r = r | 0x100000000ULL;
	r = r | 1;
	mpz_init (a);
	mpz_set_d (a, (double) r);
	mpz_nextprime (a, a);
	r = (uint64_t) mpz_get_d (a);
	mpz_clear (a);
	return (r);
}

// Mihai's explanation for how a proof verification works:
//	The verifier obtains the proof consisting of:
//		{kbnc, N, topK, B, M[0 .. N-1]}
//		(i.e.: the number k*b^n+c, the proof power N, the upper bound of the proof topK,
//		"B" the final residue at iteration topK, and a vector of middles containing N middles)
//	The verifier sets:
//		A0 = 3 (the PRP base), B0 = B
//		rootHash = hash(B)
//		And iteratively computes:
//		h0 = hash(rootH, M0)
//		A1 = A0^h0 * M0, B1 = M0^h0 * B0
//		h1 = hash(h0, M1)
//		A2 = A1^h1 * M1, B2 = M1^h1 * B1
//	And continues thus for the full set of middles, obtaining A[N] and B[N].
//	At this point, topK/(2^N) PRP iterations on A[N] should match B[N].  That is:
//		A[N] ^ (2 ^ (topK/(2^N))) == B[N]
//	If this equality holds, the proof is valid, otherwise it isn't.
//
// However, Pavel Atnashev enhanced the process by:
//		1) set r=random number
//		2) compute A[N]^r, B[N]^r
//		3) save hash(B[N]^r)
//		4) send to any random user the value x=A[N]^r and task them with computing x^(2^(topK/2^N))
//		5) user returns the hash of the result which must match our saved hash
//
// One final step, we must insure the original tester did not misreport a PRP as a composite.
// Assume the iteration that checks for a Fermat-PRP 3^((k*b^n+c)-1) is 1.  Run the PRP test
// forward to iteration topK.  If the result matches B, then the user PRP test result was misreported.

//int main (int argc, char** argv)
int verify (char *filename, int verbose)
{
	FILE	*fd;
	int	version, i, power, power_mult, prp_base, hashlen, excess_squarings, isPRP, retcode, topK, T, pm;
	uint64_t h, r;
	hash256_t rooth, *prevh, thish;
	gwhandle gwdata;
	gwnum	A, B, M, next_starting_A, saved_output_A, saved_output_B;
	hash256_t hashA, hashB;
	mpz_t	mpz_known_factors;
	char	number[2048], newline[2], buf[1024], type3_res64[17], type5_res64[17];
	double k; unsigned long b; unsigned long n; signed long c;	// k*b^n+c inputs to gwsetup
	int	square_y = TRUE;
	int	run_partial = FALSE;
	uint64_t filesize;			// Partial proof file size
	int	num_residues;			// Num residues in proof file.  Is power+1 for a complete proof file.
	int	partial_power;			// Proof power we will use because of a partial proof file
	int fft_length;
	fft_length = 0;				// Initialise as zero

// Make sure we have an input argument

//	if (argc <= 1) goto noarg;

// Get optional input argument:  -partial
// This is used to process partial proof files

//	if (strcmp (argv[1], "-partial") == 0) {
//		run_partial = TRUE;
//		argv[1] = argv[2];
// 	}

// Open the PRP proof file and parse the header

	fd = fopen (filename, "rb");
//	fd = fopen (argv[1], "rb");
	if (fd == NULL) {
		if (verbose) printf ("Error opening proof file %s\n", filename);
// 		printf ("Error opening proof file %s\n", argv[1]);
		goto end;
	}

// Get the partial PRP proof file size

	if (run_partial) filesize = guess_filesize (fd);

// Parse the PRP proof file header

	fscanf (fd, "PRP PROOF\n");
	if (fscanf (fd, "VERSION=%d\n", &version) != 1 || (version != 1 && version != 2)) {
		if (verbose) printf ("Error getting version number from proof header\n");
//		printf ("Error getting version number from proof header\n");
		goto end;
	}
	if (fscanf (fd, "HASHSIZE=%d\n", &hashlen) != 1 || hashlen < 32 || hashlen > 64) {
		if (verbose) printf ("Error getting hash size from proof header\n");
//		printf ("Error getting hash size from proof header\n");
		goto end;
	}
	if (fscanf (fd, "POWER=%d\n", &power) != 1 || power <= 0 || power >= 16) {
		if (verbose) printf ("Error getting power from proof header\n");
//		printf ("Error getting power from proof header\n");
		goto end;
	}
	if (fscanf (fd, "x%d\n", &power_mult) != 1) power_mult = 1;		// Power multiplier is an optional prime95-only feature
	if (fscanf (fd, "BASE=%d\n", &prp_base) != 1) prp_base = 3;		// BASE is an optional prime95-only construct
	if (fscanf (fd, "NUMBER=%2047[^\n]%1[\n]", number, newline) != 2) {
		if (verbose) printf ("Error getting number from proof header\n");
//		printf ("Error getting number from proof header\n");
		goto end;
	}

// Parse the number.  Handle various number formats:
//	Mn
//	[k*]2^n+c
//	([k*]2^n+c)
// Followed by 0 or more "/factor"

	if (sscanf (number, "M%ld", &n) == 1) { k = 1.0; b = 2; c = -1; }
	else if (sscanf (number, "(M%ld)", &n) == 1) { k = 1.0; b = 2; c = -1; }
	else if (sscanf (number, "F%ld", &n) == 1) { k = 1.0; b = 2; n = 1 << n; c = 1; } // plus Fermat numbers - CXC
	else if (sscanf (number, "(F%ld)", &n) == 1) { k = 1.0; b = 2; n = 1 << n; c = 1; }
	else if (sscanf (number, "2^%ld%ld", &n, &c) == 2) { k = 1.0; b = 2; }
	else if (sscanf (number, "(2^%ld%ld)", &n, &c) == 2) { k = 1.0; b = 2; }
	else if (sscanf (number, "%lf*2^%ld%ld", &k, &n, &c) == 3) { b = 2; }
	else if (sscanf (number, "(%lf*2^%ld%ld)", &k, &n, &c) == 3) { b = 2; }
	else {
		if (verbose) printf ("Error parsing number: %s\n", number);
//		printf ("Error parsing number: %s\n", number);
		goto end;
	}

	mpz_init_set_ui (mpz_known_factors, 1);
	for (char *p = strchr (number, '/'); p != NULL; p = strchr (p+1, '/')) {
		mpz_t	mpz_factor;
		char	factor[2048];
		sscanf (p, "/%2047[0-9]", &factor);
		mpz_init_set_str (mpz_factor, p+1, 10);
		mpz_mul (mpz_known_factors, mpz_known_factors, mpz_factor);
		mpz_clear (mpz_factor);
	}

	if (version == 1)
		topK = round_up_to_multiple_of (n, power_mult << power);
	else
		topK = n;

// Init the gwnum library

	gwinit (&gwdata);
	gwset_safety_margin (&gwdata, 0.1);	// Just to be safe, allow 0.1 fewer bits per word
	retcode = gwsetup(&gwdata, k, b, n, c);
	if (retcode) goto init_failed;
	gwsetnormroutine (&gwdata, 0, 1, 0);	// Error checking on

// Check the actual file size against the expected file size
// Get the PRP proof file size

	if (run_partial) {
		uint64_t residue_size;
		residue_size = divide_rounding_up ((int) ceil(gwdata.bit_length), 8);
		num_residues = (filesize - ftell (fd)) / residue_size;
		// Turn off run_partial flag if we have a complete proof file
		if (num_residues == power_mult * (power + 1)) run_partial = FALSE;
	}

// Output some info on the proof

	sprintf (buf, "Proof power = %d, ", power);
	if (power_mult > 1) sprintf (buf+strlen(buf), " power multiplier = %d, ", power_mult);
	if (version == 1) sprintf (buf+strlen(buf), " TopK = %d, ", (int) topK);
	sprintf (buf+strlen(buf), " hash length = %d", hashlen);
	if (verbose) printf ("Verifying proof for %s.  %s\n", number, buf);
//	printf ("Verifying proof for %s.  %s\n", number, buf);
	if (run_partial) {
		partial_power = num_residues - (power_mult - 1) * (power + 1) - 1;
		if (partial_power <= 0) {			// Need at least 2 residues in last power_mult proof
			if (verbose) printf ("Partial proof file contains %d residues.  Too few, exiting.\n", num_residues);
//			printf ("Partial proof file contains %d residues.  Too few, exiting.\n", num_residues);
			goto end;
		}
		if (verbose) printf ("Partial proof file contains %d residues.\n", num_residues);
		if (verbose) printf ("Attempting proof power = %d\n", partial_power);
//		printf ("Partial proof file contains %d residues.\n", num_residues);
//		printf ("Attempting proof power = %d\n", partial_power);
	}

	gwfft_description (&gwdata, buf);
	if (verbose) printf("Using %s\n", buf);
//	printf("Using %s\n", buf);

	// Calculate the excess squarings performed during the original PRP test
	excess_squarings = topK - n;
	if (excess_squarings && verbose) printf ("Prover did %d excess squarings\n", excess_squarings);
//	if (excess_squarings) printf ("Prover did %d excess squarings\n", excess_squarings);

// Start with A = 3^k (plus enough to make n a multiple of proof power multiplier)

	A = gwalloc (&gwdata);
	if (A == NULL) goto oom;
	dbltogw (&gwdata, (double) prp_base, A);
	exponentiate (&gwdata, A, (uint64_t) k);
	for (i = 0; i < topK % power_mult; i++) gwsquare (&gwdata, A);

// Allocate more memory

	B = gwalloc (&gwdata);
	if (B == NULL) goto oom;
	M = gwalloc (&gwdata);
	if (M == NULL) goto oom;
	if (power_mult > 1) {
		next_starting_A = gwalloc (&gwdata);
		if (next_starting_A == NULL) goto oom;
		saved_output_A = gwalloc (&gwdata);
		if (saved_output_A == NULL) goto oom;
		saved_output_B = gwalloc (&gwdata);
		if (saved_output_B == NULL) goto oom;
	}

// Loop over the number of proofs in proof file

	for (pm = 1; pm <= power_mult; pm++) {

		// Read the final proof residue and if excess_squarings is zero and this is the last proof in the proof file,
		// then the returned 64-bit residue is the type-3 residue.
		readResidue (&gwdata, fd, B, NULL, type3_res64);

// In version 1 of proofs, we chose topK to be round_up_to_multiple_of(n,2^proof_power).
// In version 2 of proofs, we Pavel showed us how to handle odd n values, so no rounding or excess squarings are needed.

// See if this is the residue one would get for a positive PRP test
// The comments below show how we derive the PRP residue we are looking for.
// We start with a Fermat PRP test where we residue is one and change to
// a Gerbicz / proof generating test where the exponent end in many zero bits.
//
// Test of a Mersenne number: 2^n-1
// Fermat PRP test is 3 ^ ((2^n - 1) - 1)		PRP residue: 1
// Instead calc 3 ^ (2^n - 1)				PRP residue: 3
// Instead calc 3 ^ (2^n)				PRP residue: 9
// Instead calc 3 ^ (2^(n + x))				PRP residue: 9 ^ (2^x)	where x is round_up_to_multiple_of(n,2^proof_power)-n
//
// Test of k * 2^n + c
// Fermat PRP test is 3 ^ ((k * 2^n + c) - 1)		PRP residue: 1
// Instead calc 3 ^ (k * 2^n + c)			PRP residue: 3
// Instead calc 3 ^ (k * 2^n)				PRP residue: 3 * 3^-c = 3^(1-c)
// Instead calc (3^k) ^ (2^n)				PRP residue: 3^(1-c)
// Instead calc (3^k) ^ (2^(n + x))			PRP residue: (3^(1-c)) ^ (2^x)	where x is round_up_to_multiple_of(n,2^proof_power)-n
//
// Test of (k * 2^n + c) / kf
// Fermat PRP test is 3 ^ ((k * 2^n + c) / kf - 1)	PRP residue: 1
// Instead calc 3 ^ ((k * 2^n + c) / kf)		PRP residue: 3
// Instead calc 3 ^ (k * 2^n + c)			PRP residue: 3^kf
// Instead calc 3 ^ (k * 2^n)				PRP residue: 3^kf * 3^-c = 3^(kf-c)
// Instead calc (3^k) ^ (2^n)				PRP residue: 3^(kf-c)
// Instead calc (3^k) ^ (2^(n + x))			PRP residue: (3^(kf-c)) ^ (2^x)	where x is round_up_to_multiple_of(n,2^proof_power)-n

		// The final residue of the last proof is also the final PRP residue for us to check for a PRP
		if (pm == power_mult) {
			mpz_t	mpz_modulus, mpz_inverse, mpz_tmp, mpz_roots;
			giant	g;
			gwnum	hardened_M;
			int	gmp_squarings;

			// Allocate a giant to use later
			g = popg (&gwdata.gdata, (n >> 5) + 5);
			if (g == NULL) goto oom;

			// Modulus = 2^n-1.  Later we will divide by the known factors.  Type-5 residue recovery requires the undivided modulus.
			mpz_init_set_d (mpz_modulus, k);
			mpz_mul_2exp (mpz_modulus, mpz_modulus, n);
			mpz_add_si (mpz_modulus, mpz_modulus, c);

			// Handle version 1 proof files which did excess squarings
			if (excess_squarings) {

				// Calculate the result of a successful PRP test: 3^(kf-c)
				mpz_init (mpz_tmp);
				mpz_sub_si (mpz_tmp, mpz_known_factors, c);
				if (mpz_eq_ui (mpz_tmp, 2) && prp_base < 65536) {	// mpz_powm is slow, bypass it for the common Mersenne case
					mpz_init_set_ui (mpz_inverse, prp_base * prp_base);
				} else {
					mpz_init_set_ui (mpz_inverse, prp_base);
					mpz_powm (mpz_inverse, mpz_inverse, mpz_tmp, mpz_modulus);
				}
				mpz_clear (mpz_tmp);

				// Our goal is to apply the inverse of 3^(kf-c) ^ (2^excess_squarings) to B.  The result will be one for a PRP.
				// If we invert first and then do the squarings, we run into an issue where the FFT bit pattern is not random and
				// gwsquare_carefully would be required.  Thus, and for better performance, we do some of the excess_squarings
				// in GMP before calculating the inverse.
				for (gmp_squarings = 0; gmp_squarings < excess_squarings && mpz_sizeinbase (mpz_inverse, 2) <= 256; gmp_squarings++) {
					mpz_mul (mpz_inverse, mpz_inverse, mpz_inverse);
				}
				mpz_invert (mpz_inverse, mpz_inverse, mpz_modulus);

				// Convert the inverse to a gwnum for the excess squarings
				mpztog (mpz_inverse, g);
				gianttogw (&gwdata, g, M);
				mpz_clear (mpz_inverse);

				// Do the remaining excess squarings
				gwstartnextfft (&gwdata, 1);
				for (i = gmp_squarings; i < excess_squarings; i++) {
					gwsquare (&gwdata, M);
				}
				gwstartnextfft (&gwdata, 0);

				// Mul B by the inverse^(2^excess_squarings).  This will be one if we have a PRP.
				// If excess_squarings is small, M may be a non-random FFT bit pattern.  Multiply carefully in that case.
				if (excess_squarings <= 4) {
					gwmul_carefully (&gwdata, B, M);
				} else {
					gwsafemul (&gwdata, B, M);
				}
			}

			// Handle the version 2 proof files, no excess squarings and recovery of res64 is possible
			else {
				if (verbose) printf ("Type-3 res64: %s\n", type3_res64);
//				printf ("Type-3 res64: %s\n", type3_res64);

				// Recover the type-1 (kf == 1) or type-5 (kf != 1) residue by multiplying B by the inverse of 3^(1-c)
				// NOTE: Good optimizations are available for c == 1 and c > 1 cases.
				if (c == -1 && prp_base < 65536) {	// mpz_powm is slow, bypass it for the common Mersenne case
					mpz_init_set_ui (mpz_inverse, prp_base * prp_base);
				} else {
					mpz_init_set_si (mpz_tmp, 1 - c);
					mpz_init_set_ui (mpz_inverse, prp_base);
					mpz_powm (mpz_inverse, mpz_inverse, mpz_tmp, mpz_modulus);
					mpz_clear (mpz_tmp);
				}
				mpz_invert (mpz_inverse, mpz_inverse, mpz_modulus);

				// Convert the inverse to a gwnum for the residue recovery
				mpztog (mpz_inverse, g);
				gianttogw (&gwdata, g, M);
				mpz_clear (mpz_inverse);

				// Mul B by the inverse.  Multiply carefully is required.
				gwmul_carefully (&gwdata, B, M);

				// Generate the 64-bit residue
				retcode = gwtogiant (&gwdata, M, g);
				if (retcode < 0) goto oom;				/* Should not happen */
				giantToRes64 (g, type5_res64);
				if (verbose) printf ("Type-5 res64: %s\n", type5_res64);
//				printf ("Type-5 res64: %s\n", type5_res64);

				// Our end goal is to apply the inverse of 3^(kf-c) to B.  The result will be one for a PRP.
				// If there are no known factors we've already accomplished this.  Otherwise, more work is required.
				// Mul B by inverse of 3^(kf-c).
				if (!mpz_eq_ui (mpz_known_factors, 1)) {
					mpz_init (mpz_tmp);
					mpz_sub_si (mpz_tmp, mpz_known_factors, c);
					mpz_init_set_ui (mpz_inverse, prp_base);
					mpz_powm (mpz_inverse, mpz_inverse, mpz_tmp, mpz_modulus);
					mpz_clear (mpz_tmp);
					mpz_invert (mpz_inverse, mpz_inverse, mpz_modulus);

					// Convert the inverse to a gwnum
					mpztog (mpz_inverse, g);
					gianttogw (&gwdata, g, M);
					mpz_clear (mpz_inverse);

					// Mul B by the inverse.  This will be one if we have a PRP.
					gwsafemul (&gwdata, B, M);
				}
			}

			// Pavel Atnashev and Ravi Fernando discovered a root-of-unity attack where a user can modify the proof
			// of a PRP by multiplying B (and the middle residue) by a root of (k*b^n+c)-1.  We can thwart this bizarre
			// madman by exponentiating with all the small roots (in essense converting his root-of-unity back to one).
			mpz_init_set_ui (mpz_roots, 2);				// Start with a root of 2 since (k*b^n+c)-1 is even

			// Special code path for finding roots of Mersenne numbers
			if (k == 1.0 && b == 2 && c == -1) {
				mersenne_roots (n - 1, mpz_roots);		// Now pick up the roots of 2^(n-1)-1
			}

			// Do the general case here
			else {
				// BUG/NOTE:  We have not coded hardening against the root-of-unity attack in the general case
			}

			// Raise PRP result to product-of-roots power
			double	start_fft_count = gwdata.fft_count;
			hardened_M = gwalloc (&gwdata);
			if (hardened_M == NULL) goto oom;
			gwcopy (&gwdata, M, hardened_M);
			exponentiate_mpz (&gwdata, hardened_M, mpz_roots);
			if (verbose) printf ("Proof hardening product-of-roots %d bits\n", (int) mpz_sizeinbase (mpz_roots, 2));
			if (verbose) printf ("Proof hardening cost %d squarings\n", (int) ceil ((gwdata.fft_count - start_fft_count) / 2.0));
//			printf ("Proof hardening product-of-roots %d bits\n", (int) mpz_sizeinbase (mpz_roots, 2));
//			printf ("Proof hardening cost %d squarings\n", (int) ceil ((gwdata.fft_count - start_fft_count) / 2.0));
			mpz_clear (mpz_roots);

			// We now need the divided modulus
			mpz_divexact (mpz_modulus, mpz_modulus, mpz_known_factors);

			// Convert the result back to GMP mpz_t.  If the result is one we have a PRP!
			retcode = gwtogiant (&gwdata, hardened_M, g);
			if (retcode < 0) goto oom;				/* Should not happen */
			mpz_init (mpz_tmp);
			gtompz (g, mpz_tmp);
			if (mpz_cmpabs (mpz_tmp, mpz_modulus) >= 0) mpz_mod (mpz_tmp, mpz_tmp, mpz_modulus);

			// See if we have a PRP.  If we do, see if it was one hidden by a root-of-unity attack.
			if (mpz_eq_ui (mpz_tmp, 1)) {
				retcode = gwtogiant (&gwdata, M, g);
				if (retcode < 0) goto oom;			/* Should not happen */
				gtompz (g, mpz_tmp);
				mpz_mod (mpz_tmp, mpz_tmp, mpz_modulus);
				isPRP = (mpz_eq_ui (mpz_tmp, 1) ? 1 : 2);
			} else
				isPRP = 0;

			// Cleanup
			fft_length = gwfftlen (&gwdata);
			gwfree (&gwdata, hardened_M);
			pushg (&gwdata.gdata, 1);
			mpz_clear (mpz_modulus);
			mpz_clear (mpz_tmp);

// Print out our PRP determination
		  if (verbose) { // CXC was here, too
			if (isPRP == 2) printf ("Verifying a root-of-unity obscured positive PRP test\n");
			else if (isPRP) printf ("Verifying a positive PRP test\n");
			else printf ("Verifying a negative PRP test\n");
		  }
		}

// If there are multiple proofs in the proof file, each final proof residue is the starting A value for next proof

		if (pm < power_mult) gwcopy (&gwdata, B, next_starting_A);

// Start with B = final residue, roothash = hash(B)

		roothash (&gwdata, B, &rooth);
		if (verbose) printf ("Root hash = %s\n", hash_to_string (rooth));
//		printf ("Root hash = %s\n", hash_to_string (rooth));

// Iteratively computes:
//	h0 = hash(rootH, M0)
//	A1 = A0^h0 * M0, B1 = M0^h0 * B0
//	h1 = hash(h0, M1)
//	A2 = A1^h1 * M1, B2 = M1^h1 * B1
//
// The original Pietrzak scheme handles odd n by squaring B values and setting next n to (n+1)/2.
// Pavel Atnashev comments that one can also handle odd n by squaring A values and setting next n to (n-1)/2.
// Note:  The original paper uses different variables.  Mihai's A is x and Mihai's B is y.  Here are two
// examples showing odd n values in a power 3 proof.

// A power=3, distance8=25, squaring y proof							(this col. is log of value on left with hash omitted)
// x_0 = 3													2^0
// y_0 = x^(distance8)												2^25
// h0 = hash(y_0)
// distance8 is odd, so distance4 = (distance8 + 1) / 2, y_0 = y_0^2			distance4=13		2^26
//
// u1 = x_0^(distance4)												2^13
// x_1 = x_0^h0 * u1 = x^h0 * x^distance4									2^0+2^13
// y_1 = u1^h0 * y_0 = x^(distance4*h0) * x^(distance8)								2^13+2^26
// h1 = hash(y_1)
// distance4 is odd, so distance2 = (distance4 + 1) / 2, y_1 = y_1^2			distance2=7		2^14+2^27
//
// u2 = x_1^(distance2) = x^(distance2*h0) * x^(distance6)							2^7+2^20
// x_2 = x_1^h1 * u2 = x^(h1*h0) * x^(distance4*h1) * x^(distance2*h0) * x^(distance6)				2^0+2^13+2^7+2^20
// y_2 = u2^h1 * y_1 = x^(distance2*h1*h0) * x^(distance6*h1) * x^(distance4*h0) * x^(distance8)		2^7+2^20+2^14+2^27
// if distance2 is odd, so distance = (distance2 + 1) / 2, y_2 = y_2^2			distance=4		2^8+2^21+2^15+2^28
// else distance = distance2 / 2
//
// u3 = x_2^distance = x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7)			2^4+2^17+2^11+2^24
// x_3 = x_2 * u3 = x^(h1*h0) * x^(distance4*h1) * x^(distance2*h0) * x^(distance6) *				2^0+2^13+2^7+2^20+
//			x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7)			2^4+2^17+2^11+2^24
// y_3 = u3 * y_2 = x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7) *			2^4+2^17+2^11+2^24+
//			x^(distance2*h1*h0) * x^(distance6*h1) * x^(distance4*h0) * x^(distance8)			2^8+2^21+2^15+2^28
//
// Certifier shows x_3^distance == y_3!
// 
// u1,u2,u3,y_0 are in the proof file
// verifier does the x and y calculations
// PRPer does the u calculations via recursive descent using residues at iteration 4 7 11 13 17 20 24

// A power=3, distance8=27, squaring x proof							(this col. is log of value on left with hash omitted)
// x_0 = 3													2^0
// y_0 = x^(distance8)												2^27
// h0 = hash(y_0)
// distance8 is odd, so distance4 = (distance8 - 1) / 2, x_0 = x_0^2,		distance4=13			2^1
//
// u1 = x_0^(distance4)												2^14
// x_1 = x_0^h0 * u1 = x^h0 * x^distance4									2^1+2^14
// y_1 = u1^h0 * y_0 = x^(distance4*h0) * x^(distance8)								2^14+2^27
// h1 = hash(y_1)
// distance4 is odd, so distance2 = (distance4 - 1) / 2, x_1 = x_1^2		distance2=6			2^2+2^15
//
// u2 = x_1^(distance2) = x^(distance2*h0)	* x^(distance6)							2^8+2^21
// x_2 = x_1^h1 * u2 = x^(h1*h0) * x^(distance4*h1) * x^(distance2*h0) * x^(distance6)				2^2+2^15+2^8+2^21
// y_2 = u2^h1 * y_1 = x^(distance2*h1*h0) * x^(distance6*h1) * x^(distance4*h0) * x^(distance8)		2^8+2^21+2^14+2^27
// if distance2 is even, so distance = distance2 / 2				distance=3
//
// u3 = x_2^distance = x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7)			2^5+2^18+2^11+2^24
// x_3 = x_2 * u3 = x^(h1*h0) * x^(distance4*h1) * x^(distance2*h0) * x^(distance6) *				2^2+2^15+2^8+2^21+
//			x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7)			2^5+2^18+2^11+2^24
// y_3 = u3 * y_2 = x^(distance*h1*h0) * x^(distance5*h1) * x^(distance3*h0) * x^(distance7) *			2^5+2^18+2^11+2^24+
//			x^(distance2*h1*h0) * x^(distance6*h1) * x^(distance4*h0) * x^(distance8)			2^8+2^21+2^14+2^27
//
// x_3^distance == y_3   CHECK!
//
// PRPer does the u calculations via recursive descent using residues at iteration 5 8 11 14 18 21 24

		for (i = 0, prevh = &rooth, T = topK / power_mult; i < (run_partial ? partial_power : power); i++, prevh = &thish) {
			if ((T & 1) == 0)
				T = T / 2;
			else if (square_y) {
				gwsquare (&gwdata, B);
				T = (T + 1) / 2;
			} else {
				gwsquare (&gwdata, A);
				T = (T - 1) / 2;
			}
			readResidue (&gwdata, fd, M, NULL, NULL);
			if (i != power - 1) {	// Proof generator only used power-1 hashes
				hash (&gwdata, prevh, M, &thish);
				h = truncate_hash (thish, hashlen);
				printf ("h%d = %016" PRIX64 "\n", i, h);
				exponentiate (&gwdata, A, h);
				gwsafemul (&gwdata, M, A);
				exponentiate (&gwdata, M, h);
				gwmul (&gwdata, M, B);
			} else {		// Apply last middle value without a hash
				gwmul (&gwdata, M, A);
				gwfftmul (&gwdata, M, B);
			}
		}

// Skip the residues we won't use because this is a partial power_mult proof

		for ( ; i < power; i++) {
			readResidue (&gwdata, fd, M, NULL, NULL);
		}

// If there are multiple proofs in the proof file, multiply all the output A and B values together
// Also set the starting A value for the next proof.

		if (power_mult != 1) {
			if (pm > 1) {
				gwmul (&gwdata, saved_output_A, A);
				gwmul (&gwdata, saved_output_B, B);
			}
			if (pm < power_mult) {
				gwswap (saved_output_A, A);
				gwswap (saved_output_B, B);
				gwswap (A, next_starting_A);
			}
		}
	}
	fclose (fd);

// Check if an error occurred before publishing our A[n]^random and hashB
// In production code, switch to bigger FFT size if round off error was too high

	if (gw_test_for_error (&gwdata) && verbose) printf ("Some kind of gwerror happened.\n");
//	if (gw_test_for_error (&gwdata)) printf ("Some kind of gwerror happened.\n");
#ifdef NOT_CODED_YET
	if (gw_get_maxerr (&gwdata) > 0.4) goto somewhere;
#endif

	if (verbose) printf ("Total server processing cost would be %d squarings\n", (int) ceil (gwdata.fft_count / 2.0));
	if (verbose) printf ("Certification cost is %d squarings\n", T);
//	printf ("Total server processing cost would be %d squarings\n", (int) ceil (gwdata.fft_count / 2.0));
//	printf ("Certification cost is %d squarings\n", T);

// Hash final B value

	sha3_gwnum (&gwdata, B, &hashB);
	if (verbose) printf ("Hash B = %s\n", hash_to_string (hashB));
//	printf ("Hash B = %s\n", hash_to_string (hashB));

// Do the step a random user would be tasked with

	for (i = 0; i < T; i++) {
		gwsquare (&gwdata, A);
	}
	sha3_gwnum (&gwdata, A, &hashA);
	if (verbose) printf ("Hash A = %s\n", hash_to_string (hashA));
//	printf ("Hash A = %s\n", hash_to_string (hashA));

// Check if an error occurred during the simulated user task

	if (gw_test_for_error (&gwdata)) printf ("Some kind of gwerror happened.\n");
	if (verbose) printf ("Max roe during calculations: %g\n", gw_get_maxerr (&gwdata));
//	printf ("Max roe during calculations: %g\n", gw_get_maxerr (&gwdata));

// Compare hashes
  if (verbose) {
	if (memcmp (hashA, hashB, sizeof (hash256_t)) == 0) printf ("PROOF VERIFICATION SUCCESS!\n");
	else printf("FAILED PROOF VERIFICATION\n");
  } else {
	if (memcmp (hashA, hashB, sizeof (hash256_t)) == 0) printf ("Proof successfully verified, FFT length %d\n", fft_length);
	else printf ("Proof verification failed; try re-running with --verbose or --debug for more information.\n");
  }

// End program

	return (fft_length);

//noarg:	printf ("Proof file name not given\n");
//	return (1);
init_failed:
	if (verbose) printf ("gwsetup failed: %d\n", retcode);
//	printf ("gwsetup failed: %d\n", retcode);
	return (1);
oom:	if (verbose) printf ("Out of memory\n");
//oom:	printf ("Out of memory\n");
end:	return (1);
}

