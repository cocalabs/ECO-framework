
/*===============================================================================================================
			         G E R A D O R      M E R S E N N E      T W I S T E R
===============================================================================================================*/
/*-------------------------------------------------------------------------------------------------------------
Nome:		MT_randInt
Descricao:	gerador Mersenne Twister
entrada:	inteiro que limita o intervalo de numeros
saida:		numero inteiro
Data de criacao: 08/08/08
-------------------------------------------------------------------------------------------------------------*/
#include "mersenne.h"

unsigned long MT_randInt(unsigned long n)      // inteiro entre [0,n] para n < 2^32
{
	unsigned long used = n;
	used |= used >> 1;
	used |= used >> 2;
	used |= used >> 4;
	used |= used >> 8;
	used |= used >> 16;
	unsigned long i;
	do{
		i = randInt() & used;
	}while( i > n );
	return i;
}
unsigned long randInt()
{
	register unsigned long s1;

	if( left == 0 ) reload();
	--left;

	s1 = *pNext++;
	s1 ^= (s1 >> 11);
	s1 ^= (s1 <<  7) & 0x9d2c5680UL;
	s1 ^= (s1 << 15) & 0xefc60000UL;

	return ( s1 ^ (s1 >> 18) );
}

void reload()
{
	register unsigned long *p = state;
	register int i;
	for( i = N - M; i--; ++p )
		*p = twist( p[M], p[0], p[1] );
	for( i = M; --i; ++p )
		*p = twist( p[M-N], p[0], p[1] );
	*p = twist( p[M-N], p[0], state[0] );

	left = N, pNext = state;
}
unsigned long twist(unsigned long m, unsigned long s0, unsigned long s1 )
{
	return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL);
}


void MT_seed()
{
	FILE* urandom = fopen( "/dev/urandom", "rb" );
	/*if( urandom )
	{
		unsigned long bigSeed[N];
		register unsigned long *s = bigSeed;
		register int i=N;
		register bool success = true;
		while( success && i-- )
			success = fread( s++, sizeof(unsigned long), 1, urandom );
		fclose(urandom);
		if( success ) { seed( bigSeed, N );  return; }
	}*/
	MT_seedfinal( MT_hash( time(NULL), clock() ) );
}

unsigned long MT_hash(time_t t, clock_t c)
{
	size_t i, j;
	static unsigned long  differ = 0;

	unsigned long  h1 = 0;
	unsigned char *p = (unsigned char *) &t;
	for(i = 0; i < sizeof(t); ++i)

	{
		h1 *= UCHAR_MAX + 2U;
		h1 += p[i];
	}
	unsigned long  h2 = 0;
	p = (unsigned char *) &c;
	for(j = 0; j < sizeof(c); ++j)
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += p[j];
	}
	return ( h1 + differ++ ) ^ h2;
}

void MT_seedfinal(unsigned long oneSeed)
{
	MT_initialize(oneSeed);
	reload();
}

void MT_initialize(unsigned long seed)
{
	register unsigned long *s = state;
	register unsigned long *r = state;
	register int i = 1;
	*s++ = seed & 0xffffffffUL;
	for( ; i < N; ++i )
	{
		*s++ = ( 1812433253UL * ( *r ^ (*r >> 30) ) + i ) & 0xffffffffUL;
		r++;
	}
}

float MT_randfloat()
{
	return (float)(randInt()) * (1.0/4294967295.0);
}

double MT_rand()
    { return (double) (randInt()) * (1.0/4294967296.0);
    }

double MT_randExc(const double  *n )
    { return MT_rand() * *n;
    }

unsigned long hiBit(unsigned long u) { return u & 0x80000000UL; }
unsigned long loBit(unsigned long u) { return u & 0x00000001UL; }
unsigned long loBits(unsigned long u){ return u & 0x7fffffffUL; }
unsigned long mixBits(unsigned long u, unsigned long v ) { return hiBit(u) | loBits(v); }

