/*
 * Check the Fermat number cofactor to determine if it is prime or composite.
 *
 * Values: F is the Fermat number (2^2^n + 1), P is the product of known factors, C is the remaining cofactor.
 * Prints Res64 on each step to compare with other programs. Specifically Res64 for:
 *	R = 3^((F-1)/2) mod F			Pepin test: if R == -1 mod F then F is prime
 *	A = R^2 mod F = 3^(F-1) mod F		Prime95/mprime type 5 residue
 *	B = 2^(P-1) mod F
 *	(A - B) mod C
 *
 * Can also read an mprime proof file to get the residue of 3^(F-1) mod F, where F is a Fermat number.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <errno.h>
#include <gmp.h>

#include "gwnum.h"

#define CMDLEN 1024		// Length of the command line string
#define N_FACT 10		// Number of Fermat factors supported
#define TIME_STRING_SIZE 64

#define tv_secs(tv) (tv.tv_sec + tv.tv_usec / 1000000.0)

const char *prog_name  = "cofact";
const char *prog_vers  = "0.1";
const char *build_date = __DATE__;
const char *build_time = __TIME__;

mpz_t mask64;			// 64 bit mask for print_res64
mpz_t r64;			// 64 LSBs of a residue

// Print the 64 LSBs of the number in hex with a label
void print_res64 (mpz_t n, char *name) {
    unsigned long val;

    mpz_and (r64, n, mask64);
    val = mpz_get_ui (r64);
    printf ("%s Res64: 0x%016lX\n", name, val);
}

// Print an mpz_t number with a label. Used for debug only.
void print_mpz (mpz_t n, int base, char *name) {

    printf ("%s = ", name);
    mpz_out_str (stdout, base, n);
    printf ("\n");
}

int main (int argc, char **argv) {
    int n;			// N of the Fermat number number
    int n_fact;			// The number of factors entered
    unsigned long k;		// Always 1 for a Fermat number
    unsigned long exp;		// The fermat exponent: 2^n
    unsigned long x;		// Number of Pepin test squaring operations
    unsigned long a;		// The base to test: always 3
    unsigned long m;
    unsigned long m_progress;
    int digits;			// Number of digits int the cofactor
    char *cof_s;		// The cofactor as a string
    int cores;
    int i;

    float wall_time;
    static struct timeval tv_start, tv_stop, tv_progress;
    int wall_hours, wall_mins, wall_secs;

    time_t	current_time;
    struct tm	*time_block;
    char	time_string[TIME_STRING_SIZE];

    gwhandle gwdata;
    gwnum r_gw;
    int gwerr;
    double maxerr;

    size_t r_bin_buf_len;	// Number of longs in r_bin buffer
    unsigned long *r_bin;	// Binary array for tranfer of residue from GWNUM to GMP
    int len;

    mpz_t fact[N_FACT];		// The known factors of the Fermat number
    mpz_t Fm1;			// The Fermat number - 1
    mpz_t F;			// The Fermat number
    mpz_t R;			// The Pepin residue, later the Suyama residue
    mpz_t A;			// The A residue
    mpz_t B;			// The B residue
    mpz_t P;			// The product of the known factors
    mpz_t Pm1;			// The product of the known factors minus one
    mpz_t C;			// The remaining cofactor
    mpz_t three;		// The base 3


    // Start the wall time timer
    (void) gettimeofday(&tv_start, (struct timezone *) NULL);

    // Print compile time program information
    printf ("Program to check the primality of Fermat cofactors: %s, Version %s (%s %s)\n", prog_name, prog_vers, build_date, build_time);

    // Print date and time the program was started
    current_time = time(NULL);
    time_block = localtime(&current_time);
    strftime(time_string, TIME_STRING_SIZE, "%A %d %B %Y %X", time_block);
    printf ("Run started:     %s\n", time_string);
    printf ("\n");

/*
    // Print the command line
    char cmdline[CMDLEN];
    int argi;
    bzero(cmdline, sizeof(cmdline));
    for (argi=0; argi<argc; argi++) {
	if ((strlen(cmdline) + strlen(argv[argi])) > (CMDLEN-1)) break;
	strcat (cmdline, argv[argi]);
	strcat (cmdline," ");
    }
    printf ("Command line: %s\n\n", cmdline);
*/
    // Initialize GMP variables
    for (i=0; i<=N_FACT; i++) mpz_init (fact[i]);
    mpz_init (Fm1);
    mpz_init (F);
    mpz_init (R);
    mpz_init (A);
    mpz_init (B);
    mpz_init (P);
    mpz_init (Pm1);
    mpz_init (C);
    mpz_init (mask64);
    mpz_init (r64);
    mpz_init (three);

    mpz_set_ui (mask64, 0xffffffffffffffffLL);
    mpz_set_ui (three, 3LL);

    // Parse command line parameters
    if (argc < 2) {
    	printf ("Error: must specify the Fermat number\n");
    	printf ("Usage: cofact n factor1 factor2 ...\n");
	exit(1);
    }

    n = atoi (argv[1]);
    n_fact = argc - 2;

    for (i = 0; i < n_fact; i++) {
    	if (mpz_set_str(fact[i], argv[2+i], 10) != 0) {
	    printf ("Error: cannot parse factor: %s\n", argv[2+i]);
	    exit(1);
	}
    }

    cores = 16;			// Number of cores to use; may make this a parameter later

    printf ("Testing F%d for primality using the Pepin test\n", n);
    printf ("Using %d cores in gwnum library\n", cores);
    fflush (stdout);

    k = 1;			// K for modulo value
    exp = 1 << n;		// Exponent of 2 for modulo value
    x = exp - 1;		// Number of Pepin test square mod steps: x = 2^n - 1
    a = 3;			// Base for Pepin test

    mpz_set_ui (Fm1, 1LL);					// Fm1 = 2^2^n
    mpz_mul_2exp (Fm1, Fm1, exp);				// "
    mpz_add_ui (F, Fm1, 1LL);					// F = Fm1 + 1

    gwinit (&gwdata);						// Initialize the gwnum handle
    gwset_num_threads (&gwdata, cores);				// Set number of cores to use
    gwset_safety_margin (&gwdata, 3);				// Had to set this to 3 in pmfs to prevent calc errors with very large N

    gwerr = gwsetup (&gwdata, (double) k, 2LL, exp, 1LL);	// Setup to use modulo F = 2^2^n + 1
								// Note that K is double, so only values <= 53 bits can be represented. GWNUM checks for this.
    if (gwerr) {
	printf ("gwsetup error = %d\n", gwerr);
	exit(1);
    }
    gwsetnormroutine (&gwdata, 0, 1, 0);			// Set flag to enable round-off error checking
								// This call must be AFTER gwsetup
/*
    char line[1024];
    gwfft_description (&gwdata, line);
    printf ("fft_description: %s\n", line);
    printf ("fftlen = %ld\n", gwfftlen (&gwdata));
    printf ("near_fft_limit = %d\n", gwnear_fft_limit (&gwdata, (double)3.0));
    printf ("\n");
*/
    r_gw = gwalloc (&gwdata);					// Allocate a GW number for the residue
    if (r_gw == NULL) {
	printf ("gwalloc for r_gw failed\n");
	exit(1);
    }

    binary64togw (&gwdata, &a, 1, r_gw);			// Initialize r_gw = a
    gw_clear_maxerr (&gwdata);

    // Create buffer for transfer of residues from GWNUM to GMP
    r_bin_buf_len = exp / 64 + 1;
    r_bin = (unsigned long *) calloc (r_bin_buf_len, sizeof (unsigned long));
/*
    len = gwtobinary64 (&gwdata, r_gw, r_bin, r_bin_buf_len);
    printf ("a = %2ld, %016lx %016lx\n", a, r_bin[1], r_bin[0]);
*/

    m_progress = x / 10;

    // Almost all the runtime is in the following loop
    for (m = 1; m <= x; m++) {
	if (m < 24) {						// FIXME Good for n <= 2^24? Could this be set more intelligently?
	    gwsquare2_carefully (&gwdata, r_gw, r_gw);		// r_gw = (r_gw ^ 2) mod F
	} else {
	    gwsquare2 (&gwdata, r_gw, r_gw);			// r_gw = (r_gw ^ 2) mod F
	}
/*
	if (m > (x-8)) {
	    len = gwtobinary64 (&gwdata, r_gw, r_bin, r_bin_buf_len);
	    printf ("m = %ld, %016lx %016lx\n", m, r_bin[1], r_bin[0]);
	}
*/
	maxerr = gw_get_maxerr (&gwdata);
	if (maxerr >= 0.45) {
	    printf ("Roundoff warning: k = %ld, n = %d, m = %ld, maxerr = %22.20lf\n", k, n, m, maxerr);
	    gw_clear_maxerr (&gwdata);
	}

	if (m > m_progress) {
	    (void) gettimeofday(&tv_progress, (struct timezone *) 0);
	    wall_time = tv_secs(tv_progress) - tv_secs(tv_start);
	    wall_hours = wall_time / 3600;
	    wall_mins = (wall_time - (wall_hours * 3600)) / 60;
	    wall_secs = (wall_time - (wall_hours * 3600) - (wall_mins * 60));
	    printf ("%9ld steps (%5.1f%%), Wall time = %d:%02d:%02d (HH:MM:SS)\n", m, 100.0 * m / x, wall_hours, wall_mins, wall_secs);
	    fflush (stdout);
	    m_progress += x / 10;
	}
    }

    // Check for errors
    gwerr =  gw_test_for_error (&gwdata);
    if (gwerr) {
	printf ("Error: gw_test_for_error = %d\n", gwerr);
	exit(1);
    }

    // Convert Pepin residue r_gw to r_bin to R
    len = gwtobinary64 (&gwdata, r_gw, r_bin, r_bin_buf_len);
    mpz_import (R, len, -1, 8, 0, 0, r_bin);

    printf ("Pepin residue:  len = %d, %016lx ... %016lx %016lx\n", len, r_bin[len-1], r_bin[1], r_bin[0]);

    print_res64 (R, "Pepin");

    if (mpz_cmp (R, Fm1) == 0) {
	printf ("F%d is prime!\n", n);
    } else {
	printf ("F%d is composite\n", n);
    }
    printf ("\n");

    // If known factors were provided, perform the Suyama steps to determine whether the remaining cofactor C is prime or composite
    if (n_fact > 0) {
	printf ("Testing the F%d cofactor for primality using the following known factors: ", n);
	for (i = 0; i < n_fact; i++) {
	    mpz_out_str (stdout, 10, fact[i]);
	    printf (" ");
	}
	printf ("\n");

	// Calculate P = product of the known factors
	mpz_set_ui (P, 1LL);
	for (i = 0; i < n_fact; i++) {
	    mpz_mul (P, P, fact[i]);
	}
//	print_mpz (P, 10, "P");

	// Calculate the cofactor C = F / P
	mpz_div (C, F, P);
	digits = mpz_sizeinbase (C, 10);
	cof_s = malloc (digits + 2);
	mpz_get_str (cof_s, 10, C);
	digits = strlen (cof_s);

	if (digits < 600) {
	    printf ("Cofactor (%d digits): ", digits);
	    mpz_out_str (stdout, 10, C);
	    printf ("\n");
	} else {
	    printf ("Cofactor is %d digits long\n", digits);
	}

	// Square/mod one more time to get A. This is the mprime type 5 residue.
	gwsquare2 (&gwdata, r_gw, r_gw);
	len = gwtobinary64 (&gwdata, r_gw, r_bin, r_bin_buf_len);
	mpz_import (A, len, -1, 8, 0, 0, r_bin);

	printf ("Type 5 residue: len = %d, %016lx ... %016lx %016lx\n", len, r_bin[len-1], r_bin[1], r_bin[0]);
	print_res64 (A, "A");

	// Calculate B = 3^(P-1) mod F
	mpz_sub_ui (Pm1, P, 1LL);
	mpz_powm (B, three, Pm1, F);

	print_res64 (B, "B");

        // Calculate R = (A - B) mod C
	mpz_sub (R, A, B);			// R = A - B
	mpz_mod (R, R, C);			// R = (A - B) mod C

	print_res64 (R, "Suyama");

	if (mpz_cmp_ui (R, 0LL) == 0) {
	    printf ("F%d cofactor is prime!\n", n);
	} else {
	    printf ("F%d cofactor is composite\n", n);
	}
	printf ("\n");
    }

    gwfree (&gwdata, r_gw);				// Free the GW number: GW docs do not make it clear when this is needed
    gwdone (&gwdata);					// Free all data

    // Print the date and time the program ended and the total wall time
    current_time = time(NULL);
    time_block = localtime(&current_time);
    strftime(time_string, TIME_STRING_SIZE, "%A %d %B %Y %X", time_block);

    (void) gettimeofday(&tv_stop, (struct timezone *) 0);
    wall_time = tv_secs(tv_stop) - tv_secs(tv_start);
    wall_hours = wall_time / 3600;
    wall_mins = (wall_time - (wall_hours * 3600)) / 60;
    wall_secs = (wall_time - (wall_hours * 3600) - (wall_mins * 60));

    printf ("Run ended: %s, Wall time = %d:%02d:%02d (HH:MM:SS)\n\n", time_string, wall_hours, wall_mins, wall_secs);
}

