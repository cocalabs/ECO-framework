#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//GERADOR MERSENNE TWISTER
enum { N = 624 };        // length of state vector
enum { M = 397 };        // period parameter
unsigned long state[N];  // internal state
unsigned long *pNext;    // next value to get from state
int left;          	 // number of values left before reload needed
unsigned long MT_randInt(unsigned long n);
unsigned long randInt();
void reload();
unsigned long twist(unsigned long m, unsigned long s0, unsigned long s1);
unsigned long hiBit(unsigned long u);
unsigned long loBit(unsigned long u);
unsigned long loBits(unsigned long u);
unsigned long mixBits(unsigned long u, unsigned long v );
void MT_seed();
unsigned long MT_hash(time_t t, clock_t c);
void MT_seedfinal(unsigned long oneSeed);
void MT_initialize(unsigned long seed);
float MT_randfloat();
double MT_randExc(const double  *n );

