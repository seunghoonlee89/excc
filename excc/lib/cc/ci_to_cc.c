#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

int Sc(int i, int a, int nocc) { return nocc * ( a + 1 ) - ( i + 1 ); }
int Dc(int i, int j, int a, int b, int nocc2)
{ return (int) (nocc2 * ( b*(b-1)/2 + a + 1 ) - ( j*(j-1)/2 + i + 1 )); }
int Tc(int i, int j, int k, int a, int b, int c, int nocc3)
{ return (int) (nocc3 * ( c*(c-1)*(c-2)/6 + b*(b-1)/2 + a + 1 ) - ( k*(k-1)*(k-2)/6 + j*(j-1)/2 + i + 1 )); }
int DSc(int i, int j, int k, int a, int b, int c, int nocc, int nvir, int nocc2)
{ return Dc(i, j, a, b, nocc2) * nocc * nvir + Sc(k, c, nocc); }
int S(int i, int a, int nvir) { return i*nvir+a; } 
int D(int i, int j, int a, int b, int nocc, int nvir)
{ return ((i*nocc+j)*nvir+a)*nvir+b; } 
int Dab(int i, int j, int a, int b, int nvira, int noccb, int nvirb)
{ return ((i*noccb+j)*nvira+a)*nvirb+b; } 
size_t T(int i, int j, int k, int a, int b, int c, int nocc, int nvir)
{ return ((((i*nocc+(size_t)(j))*nocc+k)*nvir+a)*nvir+b)*nvir+c; } 
size_t Taab(int i, int j, int k, int a, int b, int c, int nocca, int nvira, int noccb, int nvirb)
{ return ((((i*nocca+(size_t)(j))*noccb+k)*nvira+a)*nvira+b)*nvirb+c; } 
size_t Tabb(int i, int j, int k, int a, int b, int c, int nvira, int noccb, int nvirb)
{ return ((((i*noccb+(size_t)(j))*noccb+k)*nvira+a)*nvirb+b)*nvirb+c; } 
int64_t Q(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir)
{ return (((((((int64_t)i*nocc+j)*nocc+k)*nocc+l)*nvir+a)*nvir+b)*nvir+c)*nvir+d; } 

double t1xt1aa(int i, int j, int a, int b, int nocc, int nvir, double *t1)
{
    double t1xt1 = 0.0;
    t1xt1 += t1[S(i, a, nvir)] * t1[S(j, b, nvir)];
    t1xt1 -= t1[S(i, b, nvir)] * t1[S(j, a, nvir)];
    return t1xt1;
}

double t1xt1ab(int i, int j, int a, int b, int nocc, int nvir, double *t1)
{
    return t1[S(i, a, nvir)] * t1[S(j, b, nvir)];
}

double t1xt1xt1aaa(int i, int j, int k, int a, int b, int c, int nocc, int nvir, double *t1)
{
    double t1xt1xt1 = 0.0;
    t1xt1xt1 += t1[S(i, a, nvir)] * t1[S(j, b, nvir)] * t1[S(k, c, nvir)];
    t1xt1xt1 -= t1[S(i, a, nvir)] * t1[S(j, c, nvir)] * t1[S(k, b, nvir)];
    t1xt1xt1 -= t1[S(i, b, nvir)] * t1[S(j, a, nvir)] * t1[S(k, c, nvir)];
    t1xt1xt1 += t1[S(i, b, nvir)] * t1[S(j, c, nvir)] * t1[S(k, a, nvir)];
    t1xt1xt1 += t1[S(i, c, nvir)] * t1[S(j, a, nvir)] * t1[S(k, b, nvir)];
    t1xt1xt1 -= t1[S(i, c, nvir)] * t1[S(j, b, nvir)] * t1[S(k, a, nvir)];
    return t1xt1xt1;
}

double t1xt2aaa(int i, int j, int k, int a, int b, int c, int nocc, int nvir, double *t1, double *t2aa)
{
    double t1xt2 = 0.0;
    t1xt2 += t1[S(i, a, nvir)] * t2aa[D(j, k, b, c, nocc, nvir)];
    t1xt2 -= t1[S(i, b, nvir)] * t2aa[D(j, k, a, c, nocc, nvir)];
    t1xt2 += t1[S(i, c, nvir)] * t2aa[D(j, k, a, b, nocc, nvir)];
    t1xt2 -= t1[S(j, a, nvir)] * t2aa[D(i, k, b, c, nocc, nvir)];
    t1xt2 += t1[S(j, b, nvir)] * t2aa[D(i, k, a, c, nocc, nvir)];
    t1xt2 -= t1[S(j, c, nvir)] * t2aa[D(i, k, a, b, nocc, nvir)];
    t1xt2 += t1[S(k, a, nvir)] * t2aa[D(i, j, b, c, nocc, nvir)];
    t1xt2 -= t1[S(k, b, nvir)] * t2aa[D(i, j, a, c, nocc, nvir)];
    t1xt2 += t1[S(k, c, nvir)] * t2aa[D(i, j, a, b, nocc, nvir)];
    return t1xt2;
}

double t1xt2aab(int i, int j, int k, int a, int b, int c, int nocc, int nvir, double *t1, double *t2aa, double *t2ab)
{
    double t1xt2 = 0.0;
    t1xt2 += t1[S(i, a, nvir)] * t2ab[D(j, k, b, c, nocc, nvir)];
    t1xt2 -= t1[S(i, b, nvir)] * t2ab[D(j, k, a, c, nocc, nvir)];
    t1xt2 -= t1[S(j, a, nvir)] * t2ab[D(i, k, b, c, nocc, nvir)];
    t1xt2 += t1[S(j, b, nvir)] * t2ab[D(i, k, a, c, nocc, nvir)];
    t1xt2 += t1[S(k, c, nvir)] * t2aa[D(i, j, a, b, nocc, nvir)];
    return t1xt2;
}

double t1xt1xt1aab(int i, int j, int k, int a, int b, int c, int nocc, int nvir, double *t1)
{
    double t1xt1xt1 = 0.0;
    t1xt1xt1 += t1[S(i, a, nvir)] * t1[S(j, b, nvir)] * t1[S(k, c, nvir)];
    t1xt1xt1 -= t1[S(i, b, nvir)] * t1[S(j, a, nvir)] * t1[S(k, c, nvir)];
    return t1xt1xt1;
}

double t1xt3aaaa(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1a, double *t3aaa)
{
    double t1xt3 = 0.0;
    t1xt3 += t1a[S(i, a, nvir)] * t3aaa[T(j, k, l, b, c, d, nocc, nvir)];
    t1xt3 -= t1a[S(j, a, nvir)] * t3aaa[T(i, k, l, b, c, d, nocc, nvir)];
    t1xt3 += t1a[S(k, a, nvir)] * t3aaa[T(i, j, l, b, c, d, nocc, nvir)];
    t1xt3 -= t1a[S(l, a, nvir)] * t3aaa[T(i, j, k, b, c, d, nocc, nvir)];
    t1xt3 -= t1a[S(i, b, nvir)] * t3aaa[T(j, k, l, a, c, d, nocc, nvir)];
    t1xt3 += t1a[S(j, b, nvir)] * t3aaa[T(i, k, l, a, c, d, nocc, nvir)];
    t1xt3 -= t1a[S(k, b, nvir)] * t3aaa[T(i, j, l, a, c, d, nocc, nvir)];
    t1xt3 += t1a[S(l, b, nvir)] * t3aaa[T(i, j, k, a, c, d, nocc, nvir)];
    t1xt3 += t1a[S(i, c, nvir)] * t3aaa[T(j, k, l, a, b, d, nocc, nvir)];
    t1xt3 -= t1a[S(j, c, nvir)] * t3aaa[T(i, k, l, a, b, d, nocc, nvir)];
    t1xt3 += t1a[S(k, c, nvir)] * t3aaa[T(i, j, l, a, b, d, nocc, nvir)];
    t1xt3 -= t1a[S(l, c, nvir)] * t3aaa[T(i, j, k, a, b, d, nocc, nvir)];
    t1xt3 -= t1a[S(i, d, nvir)] * t3aaa[T(j, k, l, a, b, c, nocc, nvir)];
    t1xt3 += t1a[S(j, d, nvir)] * t3aaa[T(i, k, l, a, b, c, nocc, nvir)];
    t1xt3 -= t1a[S(k, d, nvir)] * t3aaa[T(i, j, l, a, b, c, nocc, nvir)];
    t1xt3 += t1a[S(l, d, nvir)] * t3aaa[T(i, j, k, a, b, c, nocc, nvir)];
    return t1xt3;
}

double t1xt3aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1, double *t3aaa, double *t3aab)
{
    double t1xt3 = 0.0;
    t1xt3 += t1[S(i, a, nvir)] * t3aab[T(j, k, l, b, c, d, nocc, nvir)];
    t1xt3 -= t1[S(i, b, nvir)] * t3aab[T(j, k, l, a, c, d, nocc, nvir)];
    t1xt3 += t1[S(i, c, nvir)] * t3aab[T(j, k, l, a, b, d, nocc, nvir)];
    t1xt3 -= t1[S(j, a, nvir)] * t3aab[T(i, k, l, b, c, d, nocc, nvir)];
    t1xt3 += t1[S(j, b, nvir)] * t3aab[T(i, k, l, a, c, d, nocc, nvir)];
    t1xt3 -= t1[S(j, c, nvir)] * t3aab[T(i, k, l, a, b, d, nocc, nvir)];
    t1xt3 += t1[S(k, a, nvir)] * t3aab[T(i, j, l, b, c, d, nocc, nvir)];
    t1xt3 -= t1[S(k, b, nvir)] * t3aab[T(i, j, l, a, c, d, nocc, nvir)];
    t1xt3 += t1[S(k, c, nvir)] * t3aab[T(i, j, l, a, b, d, nocc, nvir)];
    t1xt3 += t1[S(l, d, nvir)] * t3aaa[T(i, j, k, a, b, c, nocc, nvir)];
    return t1xt3;
}

double t1xt3aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1, double *t3aab)
{
    double t1xt3 = 0.0;
    t1xt3 += t1[S(i, a, nvir)] * t3aab[T(k, l, j, c, d, b, nocc, nvir)];
    t1xt3 -= t1[S(i, b, nvir)] * t3aab[T(k, l, j, c, d, a, nocc, nvir)];
    t1xt3 -= t1[S(j, a, nvir)] * t3aab[T(k, l, i, c, d, b, nocc, nvir)];
    t1xt3 += t1[S(j, b, nvir)] * t3aab[T(k, l, i, c, d, a, nocc, nvir)];
    t1xt3 += t1[S(k, c, nvir)] * t3aab[T(i, j, l, a, b, d, nocc, nvir)];
    t1xt3 -= t1[S(k, d, nvir)] * t3aab[T(i, j, l, a, b, c, nocc, nvir)];
    t1xt3 -= t1[S(l, c, nvir)] * t3aab[T(i, j, k, a, b, d, nocc, nvir)];
    t1xt3 += t1[S(l, d, nvir)] * t3aab[T(i, j, k, a, b, c, nocc, nvir)];
    return t1xt3;
}

double ut1xt1ab(int i, int j, int a, int b, int nvira, int nvirb, double *t1a, double *t1b)
{
    return t1a[S(i, a, nvira)] * t1b[S(j, b, nvirb)];
}

double ut1xt2aab(int i, int j, int k, int a, int b, int c, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t2aa, double *t2ab)
{
    double t1xt2 = 0.0;
    t1xt2 += t1a[S(i, a, nvira)] * t2ab[Dab(j, k, b, c, nvira, noccb, nvirb)];
    t1xt2 -= t1a[S(i, b, nvira)] * t2ab[Dab(j, k, a, c, nvira, noccb, nvirb)];
    t1xt2 -= t1a[S(j, a, nvira)] * t2ab[Dab(i, k, b, c, nvira, noccb, nvirb)];
    t1xt2 += t1a[S(j, b, nvira)] * t2ab[Dab(i, k, a, c, nvira, noccb, nvirb)];
    t1xt2 += t1b[S(k, c, nvirb)] * t2aa[D(i, j, a, b, nocca, nvira)];
    return t1xt2;
}

double ut1xt2abb(int i, int j, int k, int a, int b, int c, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t2ab, double *t2bb)
{
    double t1xt2 = 0.0;
    t1xt2 += t1a[S(i, a, nvira)] * t2bb[D(j, k, b, c, noccb, nvirb)];
    t1xt2 += t2ab[Dab(i, j, a, b, nvira, noccb, nvirb)] * t1b[S(k, c, nvirb)];
    t1xt2 -= t2ab[Dab(i, j, a, c, nvira, noccb, nvirb)] * t1b[S(k, b, nvirb)];
    t1xt2 -= t2ab[Dab(i, k, a, b, nvira, noccb, nvirb)] * t1b[S(j, c, nvirb)];
    t1xt2 += t2ab[Dab(i, k, a, c, nvira, noccb, nvirb)] * t1b[S(j, b, nvirb)];
    return t1xt2;
}

double ut1xt1xt1aab(int i, int j, int k, int a, int b, int c, int nvira, int nvirb, double *t1a, double *t1b)
{
    double t1xt1xt1 = 0.0;
    t1xt1xt1 += t1a[S(i, a, nvira)] * t1a[S(j, b, nvira)] * t1b[S(k, c, nvirb)];
    t1xt1xt1 -= t1a[S(i, b, nvira)] * t1a[S(j, a, nvira)] * t1b[S(k, c, nvirb)];
    return t1xt1xt1;
}

double ut1xt3aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t3aaa, double *t3aab)
{
    double t1xt3 = 0.0;
    t1xt3 += t1a[S(i, a, nvira)] * t3aab[Taab(j, k, l, b, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(i, b, nvira)] * t3aab[Taab(j, k, l, a, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1a[S(i, c, nvira)] * t3aab[Taab(j, k, l, a, b, d, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(j, a, nvira)] * t3aab[Taab(i, k, l, b, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1a[S(j, b, nvira)] * t3aab[Taab(i, k, l, a, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(j, c, nvira)] * t3aab[Taab(i, k, l, a, b, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1a[S(k, a, nvira)] * t3aab[Taab(i, j, l, b, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(k, b, nvira)] * t3aab[Taab(i, j, l, a, c, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1a[S(k, c, nvira)] * t3aab[Taab(i, j, l, a, b, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1b[S(l, d, nvirb)] * t3aaa[T(i, j, k, a, b, c, nocca, nvira)];
    return t1xt3;
}

double ut1xt3aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t3aab, double *t3abb)
{
    double t1xt3 = 0.0;
    t1xt3 += t1a[S(i, a, nvira)] * t3abb[Tabb(j, k, l, b, c, d, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(i, b, nvira)] * t3abb[Tabb(j, k, l, a, c, d, nvira, noccb, nvirb)];
    t1xt3 -= t1a[S(j, a, nvira)] * t3abb[Tabb(i, k, l, b, c, d, nvira, noccb, nvirb)];
    t1xt3 += t1a[S(j, b, nvira)] * t3abb[Tabb(i, k, l, a, c, d, nvira, noccb, nvirb)];
    t1xt3 += t1b[S(k, c, nvirb)] * t3aab[Taab(i, j, l, a, b, d, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1b[S(k, d, nvirb)] * t3aab[Taab(i, j, l, a, b, c, nocca, nvira, noccb, nvirb)];
    t1xt3 -= t1b[S(l, c, nvirb)] * t3aab[Taab(i, j, k, a, b, d, nocca, nvira, noccb, nvirb)];
    t1xt3 += t1b[S(l, d, nvirb)] * t3aab[Taab(i, j, k, a, b, c, nocca, nvira, noccb, nvirb)];
    return t1xt3;
}

double ut1xt3abbb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t3abb, double *t3bbb)
{
    double t1xt3 = 0.0;
    t1xt3 += t1a[S(i, a, nvira)] * t3bbb[T(j, k, l, b, c, d, noccb, nvirb)];
    t1xt3 += t3abb[Tabb(i, j, k, a, b, c, nvira, noccb, nvirb)] * t1b[S(l, d, nvirb)];
    t1xt3 -= t3abb[Tabb(i, j, k, a, b, d, nvira, noccb, nvirb)] * t1b[S(l, c, nvirb)];
    t1xt3 += t3abb[Tabb(i, j, k, a, c, d, nvira, noccb, nvirb)] * t1b[S(l, b, nvirb)];
    t1xt3 -= t3abb[Tabb(i, j, l, a, b, c, nvira, noccb, nvirb)] * t1b[S(k, d, nvirb)];
    t1xt3 += t3abb[Tabb(i, j, l, a, b, d, nvira, noccb, nvirb)] * t1b[S(k, c, nvirb)];
    t1xt3 -= t3abb[Tabb(i, j, l, a, c, d, nvira, noccb, nvirb)] * t1b[S(k, b, nvirb)];
    t1xt3 += t3abb[Tabb(i, k, l, a, b, c, nvira, noccb, nvirb)] * t1b[S(j, d, nvirb)];
    t1xt3 -= t3abb[Tabb(i, k, l, a, b, d, nvira, noccb, nvirb)] * t1b[S(j, c, nvirb)];
    t1xt3 += t3abb[Tabb(i, k, l, a, c, d, nvira, noccb, nvirb)] * t1b[S(j, b, nvirb)];
    return t1xt3;
}

double c3tot3aab(int i, int j, int k, int a, int b, int c, int nocc, int nocc2, int nvir, double *t1, double *t2aa, double *t2ab, double *c3aab, double c0)
{
    double t3 = 0.0;
    t3 = c3aab[DSc(i, j, k, a, b, c, nocc, nvir, nocc2)] / c0;
    t3-= t1xt2aab(i, j, k, a, b, c, nocc, nvir, t1, t2aa, t2ab); 
    t3-= t1xt1xt1aab(i, j, k, a, b, c, nocc, nvir, t1); 
    return t3;
}

double c3tot3aaa(int i, int j, int k, int a, int b, int c, int nocc, int nocc3, int nvir, double *t1, double *t2aa, double *c3aaa, double c0)
{
    double t3 = 0.0;
    t3 = c3aaa[Tc(i, j, k, a, b, c, nocc3)] / c0;
    t3-= t1xt2aaa (i, j, k, a, b, c, nocc, nvir, t1, t2aa); 
    t3-= t1xt1xt1aaa (i, j, k, a, b, c, nocc, nvir, t1); 
    return t3;
}

double t1xc3aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nocc2, int nocc3, int nvir, double *t1, double *t2aa, double *t2ab, double *c3aaa, double *c3aab, double c0)
{
    double t1xt3 = 0.0;
    t1xt3 += t1[S(i, a, nvir)] * c3tot3aab(j, k, l, b, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(i, b, nvir)] * c3tot3aab(j, k, l, a, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(i, c, nvir)] * c3tot3aab(j, k, l, a, b, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(j, a, nvir)] * c3tot3aab(i, k, l, b, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(j, b, nvir)] * c3tot3aab(i, k, l, a, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(j, c, nvir)] * c3tot3aab(i, k, l, a, b, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(k, a, nvir)] * c3tot3aab(i, j, l, b, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(k, b, nvir)] * c3tot3aab(i, j, l, a, c, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(k, c, nvir)] * c3tot3aab(i, j, l, a, b, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(l, d, nvir)] * c3tot3aaa(i, j, k, a, b, c, nocc, nocc3, nvir, t1, t2aa, c3aaa, c0);
    return t1xt3;
}

double t1xc3aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nocc2, int nvir, double *t1, double *t2aa, double *t2ab, double *c3aab, double c0)
{
    double t1xt3 = 0.0;
    t1xt3 += t1[S(i, a, nvir)] * c3tot3aab(k, l, j, c, d, b, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(i, b, nvir)] * c3tot3aab(k, l, j, c, d, a, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(j, a, nvir)] * c3tot3aab(k, l, i, c, d, b, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(j, b, nvir)] * c3tot3aab(k, l, i, c, d, a, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(k, c, nvir)] * c3tot3aab(i, j, l, a, b, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(k, d, nvir)] * c3tot3aab(i, j, l, a, b, c, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 -= t1[S(l, c, nvir)] * c3tot3aab(i, j, k, a, b, d, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    t1xt3 += t1[S(l, d, nvir)] * c3tot3aab(i, j, k, a, b, c, nocc, nocc2, nvir, t1, t2aa, t2ab, c3aab, c0);
    return t1xt3;
}

double ut2xt2aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t2aa, double *t2ab, double *t2bb)
{
    double t2xt2 = 0.0;
    t2xt2 += t2aa[D(i, j, a, b, nocca, nvira)] * t2bb[D(k, l, c, d, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, k, a, c, nvira, noccb, nvirb)] * t2ab[Dab(j, l, b, d, nvira, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, k, a, d, nvira, noccb, nvirb)] * t2ab[Dab(j, l, b, c, nvira, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, k, b, c, nvira, noccb, nvirb)] * t2ab[Dab(j, l, a, d, nvira, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, k, b, d, nvira, noccb, nvirb)] * t2ab[Dab(j, l, a, c, nvira, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, l, a, c, nvira, noccb, nvirb)] * t2ab[Dab(j, k, b, d, nvira, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, l, a, d, nvira, noccb, nvirb)] * t2ab[Dab(j, k, b, c, nvira, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, l, b, c, nvira, noccb, nvirb)] * t2ab[Dab(j, k, a, d, nvira, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, l, b, d, nvira, noccb, nvirb)] * t2ab[Dab(j, k, a, c, nvira, noccb, nvirb)];
    return t2xt2;
}

double ut2xt2abbb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t2ab, double *t2bb)
{
    double t2xt2 = 0.0;
    t2xt2 += t2ab[Dab(i, j, a, b, nvira, noccb, nvirb)] * t2bb[D(k, l, c, d, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, j, a, c, nvira, noccb, nvirb)] * t2bb[D(k, l, b, d, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, j, a, d, nvira, noccb, nvirb)] * t2bb[D(k, l, b, c, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, k, a, b, nvira, noccb, nvirb)] * t2bb[D(j, l, c, d, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, k, a, c, nvira, noccb, nvirb)] * t2bb[D(j, l, b, d, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, k, a, d, nvira, noccb, nvirb)] * t2bb[D(j, l, b, c, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, l, a, b, nvira, noccb, nvirb)] * t2bb[D(j, k, c, d, noccb, nvirb)];
    t2xt2 -= t2ab[Dab(i, l, a, c, nvira, noccb, nvirb)] * t2bb[D(j, k, b, d, noccb, nvirb)];
    t2xt2 += t2ab[Dab(i, l, a, d, nvira, noccb, nvirb)] * t2bb[D(j, k, b, c, noccb, nvirb)];
    return t2xt2;
}

double t2xt2aaaa(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t2aa)
{
    double t2xt2 = 0.0;
    t2xt2 += t2aa[D(i, j, a, b, nocc, nvir)] * t2aa[D(k, l, c, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, j, a, c, nocc, nvir)] * t2aa[D(k, l, b, d, nocc, nvir)];
    t2xt2 += t2aa[D(i, j, a, d, nocc, nvir)] * t2aa[D(k, l, b, c, nocc, nvir)];
    t2xt2 -= t2aa[D(i, k, a, b, nocc, nvir)] * t2aa[D(j, l, c, d, nocc, nvir)];
    t2xt2 += t2aa[D(i, k, a, c, nocc, nvir)] * t2aa[D(j, l, b, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, k, a, d, nocc, nvir)] * t2aa[D(j, l, b, c, nocc, nvir)];
    t2xt2 += t2aa[D(i, l, a, b, nocc, nvir)] * t2aa[D(j, k, c, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, l, a, c, nocc, nvir)] * t2aa[D(j, k, b, d, nocc, nvir)];
    t2xt2 += t2aa[D(i, l, a, d, nocc, nvir)] * t2aa[D(j, k, b, c, nocc, nvir)];
    return t2xt2;
}

double t2xt2aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t2aa, double *t2ab)
{
    double t2xt2 = 0.0;
    t2xt2 += t2aa[D(i, j, a, b, nocc, nvir)] * t2ab[D(k, l, c, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, j, a, c, nocc, nvir)] * t2ab[D(k, l, b, d, nocc, nvir)];
    t2xt2 += t2aa[D(i, j, b, c, nocc, nvir)] * t2ab[D(k, l, a, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, k, a, b, nocc, nvir)] * t2ab[D(j, l, c, d, nocc, nvir)];
    t2xt2 += t2aa[D(i, k, a, c, nocc, nvir)] * t2ab[D(j, l, b, d, nocc, nvir)];
    t2xt2 -= t2aa[D(i, k, b, c, nocc, nvir)] * t2ab[D(j, l, a, d, nocc, nvir)];
    t2xt2 += t2aa[D(j, k, b, c, nocc, nvir)] * t2ab[D(i, l, a, d, nocc, nvir)];
    t2xt2 -= t2aa[D(j, k, a, c, nocc, nvir)] * t2ab[D(i, l, b, d, nocc, nvir)];
    t2xt2 += t2aa[D(j, k, a, b, nocc, nvir)] * t2ab[D(i, l, c, d, nocc, nvir)];
    return t2xt2;
}

double ut2xt2aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t2aa, double *t2ab)
{
    double t2xt2 = 0.0;
    t2xt2 += t2aa[D(i, j, a, b, nocca, nvira)] * t2ab[Dab(k, l, c, d, nvira, noccb, nvirb)];
    t2xt2 -= t2aa[D(i, j, a, c, nocca, nvira)] * t2ab[Dab(k, l, b, d, nvira, noccb, nvirb)];
    t2xt2 += t2aa[D(i, j, b, c, nocca, nvira)] * t2ab[Dab(k, l, a, d, nvira, noccb, nvirb)];
    t2xt2 -= t2aa[D(i, k, a, b, nocca, nvira)] * t2ab[Dab(j, l, c, d, nvira, noccb, nvirb)];
    t2xt2 += t2aa[D(i, k, a, c, nocca, nvira)] * t2ab[Dab(j, l, b, d, nvira, noccb, nvirb)];
    t2xt2 -= t2aa[D(i, k, b, c, nocca, nvira)] * t2ab[Dab(j, l, a, d, nvira, noccb, nvirb)];
    t2xt2 += t2aa[D(j, k, b, c, nocca, nvira)] * t2ab[Dab(i, l, a, d, nvira, noccb, nvirb)];
    t2xt2 -= t2aa[D(j, k, a, c, nocca, nvira)] * t2ab[Dab(i, l, b, d, nvira, noccb, nvirb)];
    t2xt2 += t2aa[D(j, k, a, b, nocca, nvira)] * t2ab[Dab(i, l, c, d, nvira, noccb, nvirb)];
    return t2xt2;
}

double t2xt2aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t2aa, double *t2ab)
{
    double t2xt2 = 0.0;
    t2xt2 += t2aa[D(i, j, a, b, nocc, nvir)] * t2aa[D(k, l, c, d, nocc, nvir)];
    t2xt2 += t2ab[D(i, k, a, c, nocc, nvir)] * t2ab[D(j, l, b, d, nocc, nvir)];
    t2xt2 -= t2ab[D(i, k, a, d, nocc, nvir)] * t2ab[D(j, l, b, c, nocc, nvir)];
    t2xt2 -= t2ab[D(i, k, b, c, nocc, nvir)] * t2ab[D(j, l, a, d, nocc, nvir)];
    t2xt2 += t2ab[D(i, k, b, d, nocc, nvir)] * t2ab[D(j, l, a, c, nocc, nvir)];
    t2xt2 -= t2ab[D(i, l, a, c, nocc, nvir)] * t2ab[D(j, k, b, d, nocc, nvir)];
    t2xt2 += t2ab[D(i, l, a, d, nocc, nvir)] * t2ab[D(j, k, b, c, nocc, nvir)];
    t2xt2 += t2ab[D(i, l, b, c, nocc, nvir)] * t2ab[D(j, k, a, d, nocc, nvir)];
    t2xt2 -= t2ab[D(i, l, b, d, nocc, nvir)] * t2ab[D(j, k, a, c, nocc, nvir)];
    return t2xt2;
}


double t1xt1xt2aaaa(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1a, double *t2aa)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1a[S(i,a,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,b,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,a,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,c,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,a,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(k,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,d,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(k,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,b,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,c,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,b,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(k,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,d,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(k,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,c,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(k,l,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,d,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(k,l,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(j,a,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,b,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,a,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,c,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,a,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(k,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,d,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(k,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,b,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,c,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,b,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(k,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,d,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(k,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,c,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(k,l,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,d,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(k,l,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(i,a,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,b,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,a,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,c,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,a,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(j,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,d,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(j,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,b,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,c,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,b,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(j,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,d,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(j,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,c,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(j,l,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,d,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(j,l,a,b,nocc,nvir)];

    t1xt1xt2 += t1a[S(k,a,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,b,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,a,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,c,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,a,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,d,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,b,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,c,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,b,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,d,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,c,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,l,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,d,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,l,a,b,nocc,nvir)];

    t1xt1xt2 += t1a[S(i,a,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(j,k,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,b,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(j,k,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,a,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(j,k,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,c,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(j,k,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,a,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,d,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,b,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(j,k,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,c,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(j,k,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,b,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,d,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(i,c,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(j,k,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(i,d,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(j,k,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(l,a,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,k,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,b,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,k,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,a,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,k,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,c,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,k,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,a,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,d,nvir)] * t1a[S(i,a,nvir)] * t2aa[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,b,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,k,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,c,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,k,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,b,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,d,nvir)] * t1a[S(i,b,nvir)] * t2aa[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,c,nvir)] * t1a[S(i,d,nvir)] * t2aa[D(j,k,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,d,nvir)] * t1a[S(i,c,nvir)] * t2aa[D(j,k,a,b,nocc,nvir)];

    t1xt1xt2 += t1a[S(j,a,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,b,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,a,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,c,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,a,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,d,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,b,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,c,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,b,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,d,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,c,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,l,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,d,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,l,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(k,a,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,b,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,a,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,c,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,a,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,d,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,b,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,c,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,b,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,d,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,c,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,l,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,d,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,l,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(j,a,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,k,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,b,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,k,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,a,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,k,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,c,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,k,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,a,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,d,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,b,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,k,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,c,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,k,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,b,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,d,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(j,c,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,k,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(j,d,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,k,a,b,nocc,nvir)];

    t1xt1xt2 += t1a[S(l,a,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,k,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,b,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,k,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,a,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,k,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,c,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,k,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,a,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,d,nvir)] * t1a[S(j,a,nvir)] * t2aa[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,b,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,k,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,c,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,k,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,b,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,k,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,d,nvir)] * t1a[S(j,b,nvir)] * t2aa[D(i,k,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,c,nvir)] * t1a[S(j,d,nvir)] * t2aa[D(i,k,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,d,nvir)] * t1a[S(j,c,nvir)] * t2aa[D(i,k,a,b,nocc,nvir)];

    t1xt1xt2 += t1a[S(k,a,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,j,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,b,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,j,c,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,a,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,j,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,c,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,j,b,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,a,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,j,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,d,nvir)] * t1a[S(l,a,nvir)] * t2aa[D(i,j,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,b,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,j,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,c,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,j,a,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,b,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,j,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,d,nvir)] * t1a[S(l,b,nvir)] * t2aa[D(i,j,a,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(k,c,nvir)] * t1a[S(l,d,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    t1xt1xt2 -= t1a[S(k,d,nvir)] * t1a[S(l,c,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];

    t1xt1xt2 -= t1a[S(l,a,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,j,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,b,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,j,c,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,a,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,j,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,c,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,j,b,d,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,a,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,j,b,c,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,d,nvir)] * t1a[S(k,a,nvir)] * t2aa[D(i,j,b,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,b,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,j,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,c,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,j,a,d,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,b,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,j,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,d,nvir)] * t1a[S(k,b,nvir)] * t2aa[D(i,j,a,c,nocc,nvir)];
    t1xt1xt2 -= t1a[S(l,c,nvir)] * t1a[S(k,d,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    t1xt1xt2 += t1a[S(l,d,nvir)] * t1a[S(k,c,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    return t1xt1xt2;
}

double t1xt1xt2aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1, double *t2aa, double *t2ab)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(j,b,nvir)] * t2ab[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,a,nvir)] * t1[S(j,c,nvir)] * t2ab[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1[S(i,b,nvir)] * t1[S(j,c,nvir)] * t2ab[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,a,nvir)] * t1[S(k,b,nvir)] * t2ab[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(k,c,nvir)] * t2ab[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,b,nvir)] * t1[S(k,c,nvir)] * t2ab[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(l,d,nvir)] * t2aa[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,b,nvir)] * t1[S(l,d,nvir)] * t2aa[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 += t1[S(i,c,nvir)] * t1[S(l,d,nvir)] * t2aa[D(j,k,a,b,nocc,nvir)];
    t1xt1xt2 += t1[S(j,a,nvir)] * t1[S(k,b,nvir)] * t2ab[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(k,c,nvir)] * t2ab[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1[S(j,b,nvir)] * t1[S(k,c,nvir)] * t2ab[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1[S(j,b,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,c,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,k,a,b,nocc,nvir)];
    t1xt1xt2 += t1[S(k,a,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,j,b,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(k,b,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,j,a,c,nocc,nvir)];
    t1xt1xt2 += t1[S(k,c,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(i,b,nvir)] * t2ab[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1[S(j,a,nvir)] * t1[S(i,c,nvir)] * t2ab[D(k,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,b,nvir)] * t1[S(i,c,nvir)] * t2ab[D(k,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1[S(k,a,nvir)] * t1[S(i,b,nvir)] * t2ab[D(j,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(k,a,nvir)] * t1[S(i,c,nvir)] * t2ab[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1[S(k,b,nvir)] * t1[S(i,c,nvir)] * t2ab[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(k,a,nvir)] * t1[S(j,b,nvir)] * t2ab[D(i,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1[S(k,a,nvir)] * t1[S(j,c,nvir)] * t2ab[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(k,b,nvir)] * t1[S(j,c,nvir)] * t2ab[D(i,l,a,d,nocc,nvir)];
    return t1xt1xt2;
}

double ut1xt1xt2aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t2aa, double *t2ab)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1a[S(j,b,nvira)] * t2ab[Dab(k,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1a[S(j,c,nvira)] * t2ab[Dab(k,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,b,nvira)] * t1a[S(j,c,nvira)] * t2ab[Dab(k,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1a[S(k,b,nvira)] * t2ab[Dab(j,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1a[S(k,c,nvira)] * t2ab[Dab(j,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,b,nvira)] * t1a[S(k,c,nvira)] * t2ab[Dab(j,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,a,nvira)] * t1a[S(k,b,nvira)] * t2ab[Dab(i,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1a[S(k,c,nvira)] * t2ab[Dab(i,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,b,nvira)] * t1a[S(k,c,nvira)] * t2ab[Dab(i,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1a[S(i,b,nvira)] * t2ab[Dab(k,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,a,nvira)] * t1a[S(i,c,nvira)] * t2ab[Dab(k,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,b,nvira)] * t1a[S(i,c,nvira)] * t2ab[Dab(k,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(k,a,nvira)] * t1a[S(i,b,nvira)] * t2ab[Dab(j,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(k,a,nvira)] * t1a[S(i,c,nvira)] * t2ab[Dab(j,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(k,b,nvira)] * t1a[S(i,c,nvira)] * t2ab[Dab(j,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(k,a,nvira)] * t1a[S(j,b,nvira)] * t2ab[Dab(i,l,c,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(k,a,nvira)] * t1a[S(j,c,nvira)] * t2ab[Dab(i,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(k,b,nvira)] * t1a[S(j,c,nvira)] * t2ab[Dab(i,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(j,k,b,c,nocca,nvira)];
    t1xt1xt2 -= t1a[S(i,b,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(j,k,a,c,nocca,nvira)];
    t1xt1xt2 += t1a[S(i,c,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(j,k,a,b,nocca,nvira)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,k,b,c,nocca,nvira)];
    t1xt1xt2 += t1a[S(j,b,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,k,a,c,nocca,nvira)];
    t1xt1xt2 -= t1a[S(j,c,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,k,a,b,nocca,nvira)];
    t1xt1xt2 += t1a[S(k,a,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,j,b,c,nocca,nvira)];
    t1xt1xt2 -= t1a[S(k,b,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,j,a,c,nocca,nvira)];
    t1xt1xt2 += t1a[S(k,c,nvira)] * t1b[S(l,d,nvirb)] * t2aa[D(i,j,a,b,nocca,nvira)];
    return t1xt1xt2;
}

double t1xt1xt2aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1, double *t2aa, double *t2ab)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(j,b,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 += t1[S(k,c,nvir)] * t1[S(l,d,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(i,b,nvir)] * t2aa[D(k,l,c,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(l,c,nvir)] * t1[S(k,d,nvir)] * t2aa[D(i,j,a,b,nocc,nvir)];
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(k,c,nvir)] * t2ab[D(j,l,b,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,a,nvir)] * t1[S(k,d,nvir)] * t2ab[D(j,l,b,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,b,nvir)] * t1[S(k,c,nvir)] * t2ab[D(j,l,a,d,nocc,nvir)];
    t1xt1xt2 += t1[S(i,b,nvir)] * t1[S(k,d,nvir)] * t2ab[D(j,l,a,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,a,nvir)] * t1[S(l,c,nvir)] * t2ab[D(j,k,b,d,nocc,nvir)];
    t1xt1xt2 += t1[S(i,a,nvir)] * t1[S(l,d,nvir)] * t2ab[D(j,k,b,c,nocc,nvir)];
    t1xt1xt2 += t1[S(i,b,nvir)] * t1[S(l,c,nvir)] * t2ab[D(j,k,a,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(i,b,nvir)] * t1[S(l,d,nvir)] * t2ab[D(j,k,a,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(k,c,nvir)] * t2ab[D(i,l,b,d,nocc,nvir)];
    t1xt1xt2 += t1[S(j,a,nvir)] * t1[S(k,d,nvir)] * t2ab[D(i,l,b,c,nocc,nvir)];
    t1xt1xt2 += t1[S(j,b,nvir)] * t1[S(k,c,nvir)] * t2ab[D(i,l,a,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,b,nvir)] * t1[S(k,d,nvir)] * t2ab[D(i,l,a,c,nocc,nvir)];
    t1xt1xt2 += t1[S(j,a,nvir)] * t1[S(l,c,nvir)] * t2ab[D(i,k,b,d,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,a,nvir)] * t1[S(l,d,nvir)] * t2ab[D(i,k,b,c,nocc,nvir)];
    t1xt1xt2 -= t1[S(j,b,nvir)] * t1[S(l,c,nvir)] * t2ab[D(i,k,a,d,nocc,nvir)];
    t1xt1xt2 += t1[S(j,b,nvir)] * t1[S(l,d,nvir)] * t2ab[D(i,k,a,c,nocc,nvir)];
    return t1xt1xt2;
}

double ut1xt1xt2aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t2aa, double *t2ab, double *t2bb)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1b[S(k,c,nvirb)] * t1b[S(l,d,nvirb)] * t2aa[D(i,j,a,b,nocca,nvira)];
    t1xt1xt2 -= t1b[S(l,c,nvirb)] * t1b[S(k,d,nvirb)] * t2aa[D(i,j,a,b,nocca,nvira)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1a[S(j,b,nvira)] * t2bb[D(k,l,c,d,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1a[S(i,b,nvira)] * t2bb[D(k,l,c,d,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(k,c,nvirb)] * t2ab[Dab(j,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(k,d,nvirb)] * t2ab[Dab(j,l,b,c,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,b,nvira)] * t1b[S(k,c,nvirb)] * t2ab[Dab(j,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,b,nvira)] * t1b[S(k,d,nvirb)] * t2ab[Dab(j,l,a,c,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(l,c,nvirb)] * t2ab[Dab(j,k,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(l,d,nvirb)] * t2ab[Dab(j,k,b,c,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,b,nvira)] * t1b[S(l,c,nvirb)] * t2ab[Dab(j,k,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,b,nvira)] * t1b[S(l,d,nvirb)] * t2ab[Dab(j,k,a,c,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1b[S(k,c,nvirb)] * t2ab[Dab(i,l,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,a,nvira)] * t1b[S(k,d,nvirb)] * t2ab[Dab(i,l,b,c,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,b,nvira)] * t1b[S(k,c,nvirb)] * t2ab[Dab(i,l,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,b,nvira)] * t1b[S(k,d,nvirb)] * t2ab[Dab(i,l,a,c,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,a,nvira)] * t1b[S(l,c,nvirb)] * t2ab[Dab(i,k,b,d,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,a,nvira)] * t1b[S(l,d,nvirb)] * t2ab[Dab(i,k,b,c,nvira,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(j,b,nvira)] * t1b[S(l,c,nvirb)] * t2ab[Dab(i,k,a,d,nvira,noccb,nvirb)];
    t1xt1xt2 += t1a[S(j,b,nvira)] * t1b[S(l,d,nvirb)] * t2ab[Dab(i,k,a,c,nvira,noccb,nvirb)];
    return t1xt1xt2;
}

double ut1xt1xt2abbb(int i, int j, int k, int l, int a, int b, int c, int d, int nocca, int nvira, int noccb, int nvirb, double *t1a, double *t1b, double *t2ab, double *t2bb)
{
    double t1xt1xt2 = 0.0;
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(j,b,nvirb)] * t2bb[D(k,l,c,d,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(j,c,nvirb)] * t2bb[D(k,l,b,d,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(j,d,nvirb)] * t2bb[D(k,l,b,c,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(k,b,nvirb)] * t2bb[D(j,l,c,d,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(k,c,nvirb)] * t2bb[D(j,l,b,d,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(k,d,nvirb)] * t2bb[D(j,l,b,c,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(l,b,nvirb)] * t2bb[D(j,k,c,d,noccb,nvirb)];
    t1xt1xt2 -= t1a[S(i,a,nvira)] * t1b[S(l,c,nvirb)] * t2bb[D(j,k,b,d,noccb,nvirb)];
    t1xt1xt2 += t1a[S(i,a,nvira)] * t1b[S(l,d,nvirb)] * t2bb[D(j,k,b,c,noccb,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,j,a,b,nvira,noccb,nvirb)] * t1b[S(k,c,nvirb)] * t1b[S(l,d,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,j,a,b,nvira,noccb,nvirb)] * t1b[S(k,d,nvirb)] * t1b[S(l,c,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,k,a,b,nvira,noccb,nvirb)] * t1b[S(j,c,nvirb)] * t1b[S(l,d,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,k,a,b,nvira,noccb,nvirb)] * t1b[S(j,d,nvirb)] * t1b[S(l,c,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,l,a,b,nvira,noccb,nvirb)] * t1b[S(j,c,nvirb)] * t1b[S(k,d,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,l,a,b,nvira,noccb,nvirb)] * t1b[S(j,d,nvirb)] * t1b[S(k,c,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,j,a,c,nvira,noccb,nvirb)] * t1b[S(k,b,nvirb)] * t1b[S(l,d,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,j,a,c,nvira,noccb,nvirb)] * t1b[S(k,d,nvirb)] * t1b[S(l,b,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,k,a,c,nvira,noccb,nvirb)] * t1b[S(j,b,nvirb)] * t1b[S(l,d,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,k,a,c,nvira,noccb,nvirb)] * t1b[S(j,d,nvirb)] * t1b[S(l,b,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,l,a,c,nvira,noccb,nvirb)] * t1b[S(j,b,nvirb)] * t1b[S(k,d,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,l,a,c,nvira,noccb,nvirb)] * t1b[S(j,d,nvirb)] * t1b[S(k,b,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,j,a,d,nvira,noccb,nvirb)] * t1b[S(k,b,nvirb)] * t1b[S(l,c,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,j,a,d,nvira,noccb,nvirb)] * t1b[S(k,c,nvirb)] * t1b[S(l,b,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,k,a,d,nvira,noccb,nvirb)] * t1b[S(j,b,nvirb)] * t1b[S(l,c,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,k,a,d,nvira,noccb,nvirb)] * t1b[S(j,c,nvirb)] * t1b[S(l,b,nvirb)];
    t1xt1xt2 += t2ab[Dab(i,l,a,d,nvira,noccb,nvirb)] * t1b[S(j,b,nvirb)] * t1b[S(k,c,nvirb)];
    t1xt1xt2 -= t2ab[Dab(i,l,a,d,nvira,noccb,nvirb)] * t1b[S(j,c,nvirb)] * t1b[S(k,b,nvirb)];
    return t1xt1xt2;
}

double t1xt1xt1xt1aaaa(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1a)
{
    double t1xt1xt1xt1 = 0.0;
    t1xt1xt1xt1 += t1a[S(i, a, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, a, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, a, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 += t1a[S(i, a, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, b, nvir)];
    t1xt1xt1xt1 += t1a[S(i, a, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, a, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, b, nvir)];

    t1xt1xt1xt1 -= t1a[S(i, b, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 += t1a[S(i, b, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 += t1a[S(i, b, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, b, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, a, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, b, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 += t1a[S(i, b, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, a, nvir)];

    t1xt1xt1xt1 += t1a[S(i, c, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, c, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, a, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, c, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, d, nvir)];
    t1xt1xt1xt1 += t1a[S(i, c, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, d, nvir)] * t1a[S(l, b, nvir)];
    t1xt1xt1xt1 += t1a[S(i, c, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, a, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, c, nvir)] * t1a[S(j, d, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, b, nvir)];

    t1xt1xt1xt1 -= t1a[S(i, d, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, a, nvir)];
    t1xt1xt1xt1 += t1a[S(i, d, nvir)] * t1a[S(j, b, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 += t1a[S(i, d, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, a, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, d, nvir)] * t1a[S(j, c, nvir)] * t1a[S(k, a, nvir)] * t1a[S(l, b, nvir)];
    t1xt1xt1xt1 -= t1a[S(i, d, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, b, nvir)] * t1a[S(l, c, nvir)];
    t1xt1xt1xt1 += t1a[S(i, d, nvir)] * t1a[S(j, a, nvir)] * t1a[S(k, c, nvir)] * t1a[S(l, b, nvir)];
    return t1xt1xt1xt1;
}

double t1xt1xt1xt1aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1)
{
    double t1xt1xt1xt1 = 0.0;
    t1xt1xt1xt1 += t1[S(i, a, nvir)] * t1[S(j, b, nvir)] * t1[S(k, c, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1[S(i, a, nvir)] * t1[S(j, c, nvir)] * t1[S(k, b, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1[S(i, b, nvir)] * t1[S(j, a, nvir)] * t1[S(k, c, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 += t1[S(i, b, nvir)] * t1[S(j, c, nvir)] * t1[S(k, a, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 += t1[S(i, c, nvir)] * t1[S(j, a, nvir)] * t1[S(k, b, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1[S(i, c, nvir)] * t1[S(j, b, nvir)] * t1[S(k, a, nvir)] * t1[S(l, d, nvir)];
    return t1xt1xt1xt1;
}

double ut1xt1xt1xt1aaab(int i, int j, int k, int l, int a, int b, int c, int d, int nvira, int nvirb, double *t1a, double *t1b)
{
    double t1xt1xt1xt1 = 0.0;
    t1xt1xt1xt1 += t1a[S(i, a, nvira)] * t1a[S(j, b, nvira)] * t1a[S(k, c, nvira)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 -= t1a[S(i, a, nvira)] * t1a[S(j, c, nvira)] * t1a[S(k, b, nvira)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 -= t1a[S(i, b, nvira)] * t1a[S(j, a, nvira)] * t1a[S(k, c, nvira)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 += t1a[S(i, b, nvira)] * t1a[S(j, c, nvira)] * t1a[S(k, a, nvira)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 += t1a[S(i, c, nvira)] * t1a[S(j, a, nvira)] * t1a[S(k, b, nvira)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 -= t1a[S(i, c, nvira)] * t1a[S(j, b, nvira)] * t1a[S(k, a, nvira)] * t1b[S(l, d, nvirb)];
    return t1xt1xt1xt1;
}

// need dbg?
double t1xt1xt1xt1aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nocc, int nvir, double *t1)
{
    double t1xt1xt1xt1 = 0.0;
    t1xt1xt1xt1 += t1[S(i, a, nvir)] * t1[S(j, b, nvir)] * t1[S(k, c, nvir)] * t1[S(l, d, nvir)];
    t1xt1xt1xt1 -= t1[S(i, a, nvir)] * t1[S(j, b, nvir)] * t1[S(k, d, nvir)] * t1[S(l, c, nvir)];
    return t1xt1xt1xt1;
}

double ut1xt1xt1xt1aabb(int i, int j, int k, int l, int a, int b, int c, int d, int nvira, int nvirb, double *t1a, double *t1b)
{
    double t1xt1xt1xt1 = 0.0;
    t1xt1xt1xt1 += t1a[S(i, a, nvira)] * t1a[S(j, b, nvira)] * t1b[S(k, c, nvirb)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 -= t1a[S(i, a, nvira)] * t1a[S(j, b, nvira)] * t1b[S(k, d, nvirb)] * t1b[S(l, c, nvirb)];
    t1xt1xt1xt1 -= t1a[S(i, b, nvira)] * t1a[S(j, a, nvira)] * t1b[S(k, c, nvirb)] * t1b[S(l, d, nvirb)];
    t1xt1xt1xt1 += t1a[S(i, b, nvira)] * t1a[S(j, a, nvira)] * t1b[S(k, d, nvirb)] * t1b[S(l, c, nvirb)];
    return t1xt1xt1xt1;
}

void c1_to_t1(double *t1, double *c1, int nocc, int nvir) 
{
    int i, a, ia_c, ia_t;
    ia_c = -1;
    for (a = 0; a < nvir; a++) {
    for (i = nocc-1; i > -1; i--) {
        ia_c += 1;
        ia_t  = i * nvir + a;
        t1[ia_t] = c1[ia_c];
    }
    }
}

void c2_to_t2(double *t2aa, double *t2ab, double *c2aa, double *c2ab, double *t1, int nocc, int nvir, double numzero) 
{
    int i, j, a, b, ijab_c, ijab_t1, ijab_t2, ijab_t3, ijab_t4;
    int ia, jb, iajb_c, ijab_t;
    double tmp;

    ijab_c = -1;
    for (b = 1; b < nvir; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocc-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab_c += 1;
        ijab_t1 = ((i*nocc+j)*nvir+a)*nvir+b;
        ijab_t2 = ((i*nocc+j)*nvir+b)*nvir+a;
        ijab_t3 = ((j*nocc+i)*nvir+a)*nvir+b;
        ijab_t4 = ((j*nocc+i)*nvir+b)*nvir+a;

        tmp = c2aa[ijab_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp -= t1xt1aa (i, j, a, b, nocc, nvir, t1); 
            t2aa[ijab_t1] =  tmp;
            t2aa[ijab_t2] = -tmp;
            t2aa[ijab_t3] = -tmp;
            t2aa[ijab_t4] =  tmp;
        }
    }
    }
    }
    }

    ia = -1;
    for (a = 0; a < nvir; a++) {
    for (i = nocc-1; i > -1; i--) {
        ia += 1;
        jb  =-1;
        for (b = 0; b < nvir; b++) {
        for (j = nocc-1; j > -1; j--) {
            jb += 1;
            iajb_c = ia * nocc*nvir + jb;
            ijab_t = ((i*nocc+j)*nvir+a)*nvir+b;

            tmp = c2ab[iajb_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp -= t1xt1ab (i, j, a, b, nocc, nvir, t1); 
                t2ab[ijab_t] = tmp;
            } 
        }
        }
    }
    }
}

void c3_to_t3(double *t3aaa, double *t3aab, double *c3aaa, double *c3aab, double *t1, double *t2aa, double *t2ab, int nocc, int nvir, double numzero) 
{
    int i, j, k, a, b, c;
    size_t ijkabc_t11, ijkabc_t21, ijkabc_t31, ijkabc_t41, ijkabc_t51, ijkabc_t61;
    size_t ijkabc_t12, ijkabc_t22, ijkabc_t32, ijkabc_t42, ijkabc_t52, ijkabc_t62;
    size_t ijkabc_t13, ijkabc_t23, ijkabc_t33, ijkabc_t43, ijkabc_t53, ijkabc_t63;
    size_t ijkabc_t14, ijkabc_t24, ijkabc_t34, ijkabc_t44, ijkabc_t54, ijkabc_t64;
    size_t ijkabc_t15, ijkabc_t25, ijkabc_t35, ijkabc_t45, ijkabc_t55, ijkabc_t65;
    size_t ijkabc_t16, ijkabc_t26, ijkabc_t36, ijkabc_t46, ijkabc_t56, ijkabc_t66;
    size_t ijab, kc, ijabkc_c, ijkabc_c;

    double tmp, tmp2;
    ijkabc_c = -1;
    for (c = 2; c < nvir; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (k = nocc-1; k > 1; k--) {
    for (j = k-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijkabc_c += 1;
        tmp = c3aaa[ijkabc_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp2 = t1xt2aaa (i, j, k, a, b, c, nocc, nvir, t1, t2aa); 
            tmp2+= t1xt1xt1aaa (i, j, k, a, b, c, nocc, nvir, t1); 
            tmp -= tmp2; 
            ijkabc_t11 = T(i, j, k, a, b, c, nocc, nvir);
            ijkabc_t12 = T(i, j, k, b, c, a, nocc, nvir);
            ijkabc_t13 = T(i, j, k, c, a, b, nocc, nvir);
            ijkabc_t14 = T(i, j, k, a, c, b, nocc, nvir);
            ijkabc_t15 = T(i, j, k, b, a, c, nocc, nvir);
            ijkabc_t16 = T(i, j, k, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t11] =  tmp;
            t3aaa[ijkabc_t12] =  tmp;
            t3aaa[ijkabc_t13] =  tmp;
            t3aaa[ijkabc_t14] = -tmp;
            t3aaa[ijkabc_t15] = -tmp;
            t3aaa[ijkabc_t16] = -tmp;
    
            ijkabc_t21 = T(j, k, i, a, b, c, nocc, nvir);
            ijkabc_t22 = T(j, k, i, b, c, a, nocc, nvir);
            ijkabc_t23 = T(j, k, i, c, a, b, nocc, nvir);
            ijkabc_t24 = T(j, k, i, a, c, b, nocc, nvir);
            ijkabc_t25 = T(j, k, i, b, a, c, nocc, nvir);
            ijkabc_t26 = T(j, k, i, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t21] =  tmp;
            t3aaa[ijkabc_t22] =  tmp;
            t3aaa[ijkabc_t23] =  tmp;
            t3aaa[ijkabc_t24] = -tmp;
            t3aaa[ijkabc_t25] = -tmp;
            t3aaa[ijkabc_t26] = -tmp;
    
            ijkabc_t31 = T(k, i, j, a, b, c, nocc, nvir);
            ijkabc_t32 = T(k, i, j, b, c, a, nocc, nvir);
            ijkabc_t33 = T(k, i, j, c, a, b, nocc, nvir);
            ijkabc_t34 = T(k, i, j, a, c, b, nocc, nvir);
            ijkabc_t35 = T(k, i, j, b, a, c, nocc, nvir);
            ijkabc_t36 = T(k, i, j, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t31] =  tmp;
            t3aaa[ijkabc_t32] =  tmp;
            t3aaa[ijkabc_t33] =  tmp;
            t3aaa[ijkabc_t34] = -tmp;
            t3aaa[ijkabc_t35] = -tmp;
            t3aaa[ijkabc_t36] = -tmp;
    
            ijkabc_t41 = T(i, k, j, a, b, c, nocc, nvir);
            ijkabc_t42 = T(i, k, j, b, c, a, nocc, nvir);
            ijkabc_t43 = T(i, k, j, c, a, b, nocc, nvir);
            ijkabc_t44 = T(i, k, j, a, c, b, nocc, nvir);
            ijkabc_t45 = T(i, k, j, b, a, c, nocc, nvir);
            ijkabc_t46 = T(i, k, j, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t41] = -tmp;
            t3aaa[ijkabc_t42] = -tmp;
            t3aaa[ijkabc_t43] = -tmp;
            t3aaa[ijkabc_t44] =  tmp;
            t3aaa[ijkabc_t45] =  tmp;
            t3aaa[ijkabc_t46] =  tmp;
    
            ijkabc_t51 = T(j, i, k, a, b, c, nocc, nvir);
            ijkabc_t52 = T(j, i, k, b, c, a, nocc, nvir);
            ijkabc_t53 = T(j, i, k, c, a, b, nocc, nvir);
            ijkabc_t54 = T(j, i, k, a, c, b, nocc, nvir);
            ijkabc_t55 = T(j, i, k, b, a, c, nocc, nvir);
            ijkabc_t56 = T(j, i, k, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t51] = -tmp;
            t3aaa[ijkabc_t52] = -tmp;
            t3aaa[ijkabc_t53] = -tmp;
            t3aaa[ijkabc_t54] =  tmp;
            t3aaa[ijkabc_t55] =  tmp;
            t3aaa[ijkabc_t56] =  tmp;
    
            ijkabc_t61 = T(k, j, i, a, b, c, nocc, nvir);
            ijkabc_t62 = T(k, j, i, b, c, a, nocc, nvir);
            ijkabc_t63 = T(k, j, i, c, a, b, nocc, nvir);
            ijkabc_t64 = T(k, j, i, a, c, b, nocc, nvir);
            ijkabc_t65 = T(k, j, i, b, a, c, nocc, nvir);
            ijkabc_t66 = T(k, j, i, c, b, a, nocc, nvir);
    
            t3aaa[ijkabc_t61] = -tmp;
            t3aaa[ijkabc_t62] = -tmp;
            t3aaa[ijkabc_t63] = -tmp;
            t3aaa[ijkabc_t64] =  tmp;
            t3aaa[ijkabc_t65] =  tmp;
            t3aaa[ijkabc_t66] =  tmp;
        }
    }
    }
    }
    }
    }
    }

    ijab = -1;
    for (b = 1; b < nvir; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocc-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab += 1;
        kc  =-1;
        for (c = 0; c < nvir; c++) {
        for (k = nocc-1; k > -1; k--) {
            kc += 1;
            ijabkc_c = ijab * nocc*nvir + kc;

            tmp = c3aab[ijabkc_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = t1xt2aab(i, j, k, a, b, c, nocc, nvir, t1, t2aa, t2ab); 
                tmp2+= t1xt1xt1aab(i, j, k, a, b, c, nocc, nvir, t1); 

                tmp -= tmp2;    
                ijkabc_t11 = T(i, j, k, a, b, c, nocc, nvir);
                ijkabc_t12 = T(i, j, k, b, a, c, nocc, nvir);
        
                t3aab[ijkabc_t11] =  tmp;
                t3aab[ijkabc_t12] = -tmp;
       
                ijkabc_t21 = T(j, i, k, a, b, c, nocc, nvir);
                ijkabc_t22 = T(j, i, k, b, a, c, nocc, nvir);
       
                t3aab[ijkabc_t21] = -tmp;
                t3aab[ijkabc_t22] =  tmp;
            }
        }
        }
    }
    }
    }
    }
}

void c4_to_t4(double *t4aaab, double *t4aabb, double *c4aaab, double *c4aabb, double *t1, double *t2aa, double *t2ab, double *t3aaa, double *t3aab, int nocc, int nvir, double numzero) 
{
    int i, j, k, l, a, b, c, d, m_ijab;
    int ijkabc, ld, ijkabcld_c;
    int64_t ijklabcd_t11, ijklabcd_t21, ijklabcd_t31, ijklabcd_t41, ijklabcd_t51, ijklabcd_t61;
    int64_t ijklabcd_t12, ijklabcd_t22, ijklabcd_t32, ijklabcd_t42, ijklabcd_t52, ijklabcd_t62;
    int64_t ijklabcd_t13, ijklabcd_t23, ijklabcd_t33, ijklabcd_t43, ijklabcd_t53, ijklabcd_t63;
    int64_t ijklabcd_t14, ijklabcd_t24, ijklabcd_t34, ijklabcd_t44, ijklabcd_t54, ijklabcd_t64;
    int64_t ijklabcd_t15, ijklabcd_t25, ijklabcd_t35, ijklabcd_t45, ijklabcd_t55, ijklabcd_t65;
    int64_t ijklabcd_t16, ijklabcd_t26, ijklabcd_t36, ijklabcd_t46, ijklabcd_t56, ijklabcd_t66;
    int ijab, klcd, ijabklcd_c;

    double tmp, tmp2;
    ijkabc = -1;
    for (c = 2; c < nvir; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (k = nocc-1; k > 1; k--) {
    for (j = k-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijkabc += 1;
        ld = -1;
        for (d = 0; d < nvir; d++) {
        for (l = nocc-1; l > -1; l--) {
            ld += 1;
            ijkabcld_c = ijkabc * nocc*nvir + ld;
            tmp = c4aaab[ijkabcld_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = t1xt3aaab (i, j, k, l, a, b, c, d, nocc, nvir, t1, t3aaa, t3aab);
                tmp2+= t2xt2aaab (i, j, k, l, a, b, c, d, nocc, nvir, t2aa, t2ab);        
                tmp2+= t1xt1xt2aaab (i, j, k, l, a, b, c, d, nocc, nvir, t1, t2aa, t2ab); 
                tmp2+= t1xt1xt1xt1aaab (i, j, k, l, a, b, c, d, nocc, nvir, t1);         

                tmp -= tmp2; 
                ijklabcd_t11 = Q(i, j, k, l, a, b, c, d, nocc, nvir);
                ijklabcd_t12 = Q(i, j, k, l, b, c, a, d, nocc, nvir);
                ijklabcd_t13 = Q(i, j, k, l, c, a, b, d, nocc, nvir);
                ijklabcd_t14 = Q(i, j, k, l, a, c, b, d, nocc, nvir);
                ijklabcd_t15 = Q(i, j, k, l, b, a, c, d, nocc, nvir);
                ijklabcd_t16 = Q(i, j, k, l, c, b, a, d, nocc, nvir);

                t4aaab[ijklabcd_t11] =  tmp;
                t4aaab[ijklabcd_t12] =  tmp;
                t4aaab[ijklabcd_t13] =  tmp;
                t4aaab[ijklabcd_t14] = -tmp;
                t4aaab[ijklabcd_t15] = -tmp;
                t4aaab[ijklabcd_t16] = -tmp;
        
                ijklabcd_t21 = Q(j, k, i, l, a, b, c, d, nocc, nvir);
                ijklabcd_t22 = Q(j, k, i, l, b, c, a, d, nocc, nvir);
                ijklabcd_t23 = Q(j, k, i, l, c, a, b, d, nocc, nvir);
                ijklabcd_t24 = Q(j, k, i, l, a, c, b, d, nocc, nvir);
                ijklabcd_t25 = Q(j, k, i, l, b, a, c, d, nocc, nvir);
                ijklabcd_t26 = Q(j, k, i, l, c, b, a, d, nocc, nvir);
        
                t4aaab[ijklabcd_t21] =  tmp;
                t4aaab[ijklabcd_t22] =  tmp;
                t4aaab[ijklabcd_t23] =  tmp;
                t4aaab[ijklabcd_t24] = -tmp;
                t4aaab[ijklabcd_t25] = -tmp;
                t4aaab[ijklabcd_t26] = -tmp;
        
                ijklabcd_t31 = Q(k, i, j, l, a, b, c, d, nocc, nvir);
                ijklabcd_t32 = Q(k, i, j, l, b, c, a, d, nocc, nvir);
                ijklabcd_t33 = Q(k, i, j, l, c, a, b, d, nocc, nvir);
                ijklabcd_t34 = Q(k, i, j, l, a, c, b, d, nocc, nvir);
                ijklabcd_t35 = Q(k, i, j, l, b, a, c, d, nocc, nvir);
                ijklabcd_t36 = Q(k, i, j, l, c, b, a, d, nocc, nvir);
        
                t4aaab[ijklabcd_t31] =  tmp;
                t4aaab[ijklabcd_t32] =  tmp;
                t4aaab[ijklabcd_t33] =  tmp;
                t4aaab[ijklabcd_t34] = -tmp;
                t4aaab[ijklabcd_t35] = -tmp;
                t4aaab[ijklabcd_t36] = -tmp;
        
                ijklabcd_t41 = Q(i, k, j, l, a, b, c, d, nocc, nvir);
                ijklabcd_t42 = Q(i, k, j, l, b, c, a, d, nocc, nvir);
                ijklabcd_t43 = Q(i, k, j, l, c, a, b, d, nocc, nvir);
                ijklabcd_t44 = Q(i, k, j, l, a, c, b, d, nocc, nvir);
                ijklabcd_t45 = Q(i, k, j, l, b, a, c, d, nocc, nvir);
                ijklabcd_t46 = Q(i, k, j, l, c, b, a, d, nocc, nvir);
        
                t4aaab[ijklabcd_t41] = -tmp;
                t4aaab[ijklabcd_t42] = -tmp;
                t4aaab[ijklabcd_t43] = -tmp;
                t4aaab[ijklabcd_t44] =  tmp;
                t4aaab[ijklabcd_t45] =  tmp;
                t4aaab[ijklabcd_t46] =  tmp;
        
                ijklabcd_t51 = Q(j, i, k, l, a, b, c, d, nocc, nvir);
                ijklabcd_t52 = Q(j, i, k, l, b, c, a, d, nocc, nvir);
                ijklabcd_t53 = Q(j, i, k, l, c, a, b, d, nocc, nvir);
                ijklabcd_t54 = Q(j, i, k, l, a, c, b, d, nocc, nvir);
                ijklabcd_t55 = Q(j, i, k, l, b, a, c, d, nocc, nvir);
                ijklabcd_t56 = Q(j, i, k, l, c, b, a, d, nocc, nvir);
        
                t4aaab[ijklabcd_t51] = -tmp;
                t4aaab[ijklabcd_t52] = -tmp;
                t4aaab[ijklabcd_t53] = -tmp;
                t4aaab[ijklabcd_t54] =  tmp;
                t4aaab[ijklabcd_t55] =  tmp;
                t4aaab[ijklabcd_t56] =  tmp;
        
                ijklabcd_t61 = Q(k, j, i, l, a, b, c, d, nocc, nvir);
                ijklabcd_t62 = Q(k, j, i, l, b, c, a, d, nocc, nvir);
                ijklabcd_t63 = Q(k, j, i, l, c, a, b, d, nocc, nvir);
                ijklabcd_t64 = Q(k, j, i, l, a, c, b, d, nocc, nvir);
                ijklabcd_t65 = Q(k, j, i, l, b, a, c, d, nocc, nvir);
                ijklabcd_t66 = Q(k, j, i, l, c, b, a, d, nocc, nvir);
        
                t4aaab[ijklabcd_t61] = -tmp;
                t4aaab[ijklabcd_t62] = -tmp;
                t4aaab[ijklabcd_t63] = -tmp;
                t4aaab[ijklabcd_t64] =  tmp;
                t4aaab[ijklabcd_t65] =  tmp;
                t4aaab[ijklabcd_t66] =  tmp;
            }
        }
        }
    }
    }
    }
    }
    }
    }

    m_ijab = nocc*(nocc-1)/2 * nvir*(nvir-1)/2;
    ijab = -1;
    for (b = 1; b < nvir; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocc-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab += 1;
        klcd  =-1;
        for (d = 1; d < nvir; d++) {
        for (c = 0; c < d; c++) {
        for (l = nocc-1; l > 0; l--) {
        for (k = l-1; k > -1; k--) {
            klcd += 1;
            ijabklcd_c = ijab * m_ijab + klcd;
            tmp = c4aabb[ijabklcd_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = 0.0;
                tmp2 = t1xt3aabb(i, j, k, l, a, b, c, d, nocc, nvir, t1, t3aab); 
                tmp2+= t2xt2aabb(i, j, k, l, a, b, c, d, nocc, nvir, t2aa, t2ab); 
                tmp2+= t1xt1xt2aabb(i, j, k, l, a, b, c, d, nocc, nvir, t1, t2aa, t2ab); 
                tmp2+= t1xt1xt1xt1aabb(i, j, k, l, a, b, c, d, nocc, nvir, t1); 

                tmp -= tmp2; 
                ijklabcd_t11 = Q(i, j, k, l, a, b, c, d, nocc, nvir);
                ijklabcd_t12 = Q(j, i, k, l, b, a, c, d, nocc, nvir);
                ijklabcd_t13 = Q(i, j, k, l, b, a, c, d, nocc, nvir);
                ijklabcd_t14 = Q(j, i, k, l, a, b, c, d, nocc, nvir);

                //printf("ijklabcd11,%" PRId64 "\n", ijklabcd_t11);
                //printf("ijklabcd12,%" PRId64 "\n", ijklabcd_t12);
                //printf("ijklabcd13,%" PRId64 "\n", ijklabcd_t13);
                //printf("ijklabcd14,%" PRId64 "\n", ijklabcd_t14);
 
                t4aabb[ijklabcd_t11] =  tmp;
                t4aabb[ijklabcd_t12] =  tmp;
                t4aabb[ijklabcd_t13] = -tmp;
                t4aabb[ijklabcd_t14] = -tmp;  
    
                ijklabcd_t21 = Q(i, j, l, k, a, b, d, c, nocc, nvir);
                ijklabcd_t22 = Q(j, i, l, k, b, a, d, c, nocc, nvir);
                ijklabcd_t23 = Q(i, j, l, k, b, a, d, c, nocc, nvir);
                ijklabcd_t24 = Q(j, i, l, k, a, b, d, c, nocc, nvir);
        
                t4aabb[ijklabcd_t21] =  tmp;
                t4aabb[ijklabcd_t22] =  tmp;
                t4aabb[ijklabcd_t23] = -tmp;
                t4aabb[ijklabcd_t24] = -tmp;  
    
                ijklabcd_t31 = Q(i, j, k, l, a, b, d, c, nocc, nvir);
                ijklabcd_t32 = Q(j, i, k, l, b, a, d, c, nocc, nvir);
                ijklabcd_t33 = Q(i, j, k, l, b, a, d, c, nocc, nvir);
                ijklabcd_t34 = Q(j, i, k, l, a, b, d, c, nocc, nvir);
        
                t4aabb[ijklabcd_t31] = -tmp;
                t4aabb[ijklabcd_t32] = -tmp;
                t4aabb[ijklabcd_t33] =  tmp;
                t4aabb[ijklabcd_t34] =  tmp;  
    
                ijklabcd_t41 = Q(i, j, l, k, a, b, c, d, nocc, nvir);
                ijklabcd_t42 = Q(j, i, l, k, b, a, c, d, nocc, nvir);
                ijklabcd_t43 = Q(i, j, l, k, b, a, c, d, nocc, nvir);
                ijklabcd_t44 = Q(j, i, l, k, a, b, c, d, nocc, nvir);
        
                t4aabb[ijklabcd_t41] = -tmp;
                t4aabb[ijklabcd_t42] = -tmp;
                t4aabb[ijklabcd_t43] =  tmp;
                t4aabb[ijklabcd_t44] =  tmp;  
            }
        }
        }
        }
        }
    }
    }
    }
    }
}

void c2_to_t2u(double *t2aa, double *t2ab, double *t2bb, double *c2aa, double *c2ab, double *c2bb, double *t1a, double *t1b, int nocca, int nvira, int noccb, int nvirb, double numzero) 
{
    int i, j, a, b, ijab_c;
    int ia, jb, iajb_c;
    int m_jb;
    double tmp;

    ijab_c = -1;
    for (b = 1; b < nvira; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocca-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab_c += 1;
        tmp = c2aa[ijab_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp -= t1xt1aa (i, j, a, b, nocca, nvira, t1a); 
            t2aa[(((i*nocca+j)*nvira+a)*nvira+b)] =  tmp;
            t2aa[(((i*nocca+j)*nvira+b)*nvira+a)] = -tmp;
            t2aa[(((j*nocca+i)*nvira+a)*nvira+b)] = -tmp;
            t2aa[(((j*nocca+i)*nvira+b)*nvira+a)] =  tmp;
        }
    }
    }
    }
    }

    ijab_c = -1;
    for (b = 1; b < nvirb; b++) {
    for (a = 0; a < b; a++) {
    for (j = noccb-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab_c += 1;
        tmp = c2bb[ijab_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp -= t1xt1aa (i, j, a, b, noccb, nvirb, t1b); 
            t2bb[(((i*noccb+j)*nvirb+a)*nvirb+b)] =  tmp;
            t2bb[(((i*noccb+j)*nvirb+b)*nvirb+a)] = -tmp;
            t2bb[(((j*noccb+i)*nvirb+a)*nvirb+b)] = -tmp;
            t2bb[(((j*noccb+i)*nvirb+b)*nvirb+a)] =  tmp;
        }
    }
    }
    }
    }

    m_jb = noccb*nvirb;
    ia = -1;
    for (a = 0; a < nvira; a++) {
    for (i = nocca-1; i > -1; i--) {
        ia += 1;
        jb  =-1;
        for (b = 0; b < nvirb; b++) {
        for (j = noccb-1; j > -1; j--) {
            jb += 1;
            iajb_c = ia * m_jb + jb;
            tmp = c2ab[iajb_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp -= ut1xt1ab (i, j, a, b, nvira, nvirb, t1a, t1b);
                t2ab[(((i*noccb+j)*nvira+a)*nvirb+b)] =  tmp;
            } 
        }
        }
    }
    }
}

void c3_to_t3u(double *t3aaa, double *t3aab, double *t3abb, double *t3bbb, double *c3aaa, double *c3aab, double *c3abb, double *c3bbb, double *t1a, double *t1b, double *t2aa, double *t2ab, double *t2bb, int nocca, int nvira, int noccb, int nvirb, double numzero) 
{
    int i, j, k, a, b, c, m_jkbc, m_kc;
    size_t ijab, kc, ijabkc_c, ijkabc_c, ia, jkbc, iajkbc_c;

    double tmp, tmp2;
    ijkabc_c = -1;
    for (c = 2; c < nvira; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (k = nocca-1; k > 1; k--) {
    for (j = k-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijkabc_c += 1;
        tmp = c3aaa[ijkabc_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp2 = t1xt2aaa (i, j, k, a, b, c, nocca, nvira, t1a, t2aa); 
            tmp2+= t1xt1xt1aaa (i, j, k, a, b, c, nocca, nvira, t1a); 
            tmp -= tmp2; 

            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t3aaa[(((((i*nocca+j)*nocca+k)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t3aaa[(((((j*nocca+k)*nocca+i)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t3aaa[(((((k*nocca+i)*nocca+j)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t3aaa[(((((i*nocca+k)*nocca+j)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t3aaa[(((((k*nocca+j)*nocca+i)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t3aaa[(((((j*nocca+i)*nocca+k)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
        }
    }
    }
    }
    }
    }
    }

    ijkabc_c = -1;
    for (c = 2; c < nvirb; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (k = noccb-1; k > 1; k--) {
    for (j = k-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijkabc_c += 1;
        tmp = c3bbb[ijkabc_c]; 
        if(fabs(tmp) > numzero) 
        {
            tmp2 = t1xt2aaa (i, j, k, a, b, c, noccb, nvirb, t1b, t2bb); 
            tmp2+= t1xt1xt1aaa (i, j, k, a, b, c, noccb, nvirb, t1b); 
            tmp -= tmp2; 
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t3bbb[(((((i*noccb+j)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t3bbb[(((((j*noccb+k)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t3bbb[(((((k*noccb+i)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t3bbb[(((((i*noccb+k)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t3bbb[(((((k*noccb+j)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t3bbb[(((((j*noccb+i)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
        }
    }
    }
    }
    }
    }
    }

    m_kc = noccb * nvirb;
    ijab = -1;
    for (b = 1; b < nvira; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocca-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab += 1;
        kc = -1;
        for (c = 0; c < nvirb; c++) {
        for (k = noccb-1; k > -1; k--) {
            kc += 1;
            ijabkc_c = ijab * m_kc + kc;

            tmp = c3aab[ijabkc_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = ut1xt2aab(i, j, k, a, b, c, nocca, nvira, noccb, nvirb, t1a, t1b, t2aa, t2ab); 
                tmp2+= ut1xt1xt1aab(i, j, k, a, b, c, nvira, nvirb, t1a ,t1b); 
                tmp -= tmp2;    

                t3aab[(((((i*nocca+j)*noccb+k)*nvira+a)*nvira+b)*nvirb+c)] =  tmp;
                t3aab[(((((i*nocca+j)*noccb+k)*nvira+b)*nvira+a)*nvirb+c)] = -tmp;
                t3aab[(((((j*nocca+i)*noccb+k)*nvira+a)*nvira+b)*nvirb+c)] = -tmp;
                t3aab[(((((j*nocca+i)*noccb+k)*nvira+b)*nvira+a)*nvirb+c)] =  tmp;
            }
        }
        }
    }
    }
    }
    }

    m_jkbc = noccb * (noccb - 1) / 2 * nvirb * (nvirb - 1) / 2;
    ia = -1;
    for (a = 0; a < nvira; a++) {
    for (i = nocca-1; i > -1; i--) {
        ia += 1;
        jkbc = -1;
        for (c = 1; c < nvirb; c++) {
        for (b = 0; b < c; b++) {
        for (k = noccb-1; k > 0; k--) {
        for (j = k-1; j > -1; j--) {
            jkbc += 1;
            iajkbc_c = ia * m_jkbc + jkbc;
            tmp = c3abb[iajkbc_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = ut1xt2abb(i, j, k, a, b, c, nvira, noccb, nvirb, t1a, t1b, t2ab, t2bb); 
                tmp2+= ut1xt1xt1aab(j, k, i, b, c, a, nvirb, nvira, t1b ,t1a); 
                tmp -= tmp2;    

                t3abb[(((((i*noccb+j)*noccb+k)*nvira+a)*nvirb+b)*nvirb+c)] =  tmp;
                t3abb[(((((i*noccb+j)*noccb+k)*nvira+a)*nvirb+c)*nvirb+b)] = -tmp;
                t3abb[(((((i*noccb+k)*noccb+j)*nvira+a)*nvirb+b)*nvirb+c)] = -tmp;
                t3abb[(((((i*noccb+k)*noccb+j)*nvira+a)*nvirb+c)*nvirb+b)] =  tmp;
            }
        }
        }
        }
        }
    }
    }
}

void c4_to_t4u(double *t4aaaa, double *t4aaab, double *t4aabb, double *t4abbb, double *t4bbbb, double *c4aaaa, double *c4aaab, double *c4aabb, double *c4abbb, double *c4bbbb, double *t1a, double *t1b, double *t2aa, double *t2ab, double *t2bb, double *t3aaa, double *t3aab, double *t3abb, double *t3bbb,int nocca, int nvira, int noccb, int nvirb, double numzero) 
{
    int i, j, k, l, a, b, c, d, m_klcd, m_jklbcd, m_ld;
    int ijkabc, ld, ijkabcld_c, ijklabcd;
    int ijab, klcd, ijabklcd_c, ia, jklbcd, iajklbcd_c;
    double tmp, tmp2;

    ijklabcd = -1;
    for (d = 3; d < nvira; d++) {
    for (c = 2; c < d; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (l = nocca-1; l > 1; l--) {
    for (k = l-1; k > 0; k--) {
    for (j = k-1; j > -1; j--) {
    for (i = j-1; i > -2; i--) {
        ijklabcd += 1;
        tmp = c4aaaa[ijklabcd]; 
        if(fabs(tmp) > numzero) 
        {
            tmp2 = t1xt3aaaa (i, j, k, l, a, b, c, d, nocca, nvira, t1a, t3aaa);
            tmp2+= t2xt2aaaa (i, j, k, l, a, b, c, d, nocca, nvira, t2aa);
            tmp2+= t1xt1xt2aaaa (i, j, k, l, a, b, c, d, nocca, nvira, t1a, t2aa); 
            tmp2+= t1xt1xt1xt1aaaa (i, j, k, l, a, b, c, d, nocca, nvira, t1a);

            tmp -= tmp2; 

            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+k)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+l)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+i)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+j)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+l)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+l)*nocca+j)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+i)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+k)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+j)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+j)*nocca+k)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+k)*nocca+i)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+l)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+j)*nocca+l)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+k)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+i)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+i)*nocca+j)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+l)*nocca+k)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+k)*nocca+j)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+i)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+i)*nocca+l)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)i*nocca+k)*nocca+j)*nocca+l)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)k*nocca+j)*nocca+l)*nocca+i)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)j*nocca+l)*nocca+i)*nocca+k)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+b)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+c)*nvira+d)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+d)*nvira+a)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+a)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+c)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+d)*nvira+b)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+b)*nvira+a)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+a)*nvira+c)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+d)*nvira+b)*nvira+c)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+b)*nvira+c)*nvira+a)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+c)*nvira+a)*nvira+d)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+a)*nvira+d)*nvira+b)] = -tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+b)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+d)*nvira+c)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+c)*nvira+a)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+a)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+d)*nvira+c)*nvira+b)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+c)*nvira+b)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+b)*nvira+a)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+a)*nvira+d)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+a)*nvira+c)*nvira+b)*nvira+d)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+c)*nvira+b)*nvira+d)*nvira+a)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+b)*nvira+d)*nvira+a)*nvira+c)] =  tmp;
            t4aaaa[((((((((int64_t)l*nocca+i)*nocca+k)*nocca+j)*nvira+d)*nvira+a)*nvira+c)*nvira+b)] =  tmp;
        }
    }
    }
    }
    }
    }
    }
    }
    }

    ijklabcd = -1;
    for (d = 3; d < nvirb; d++) {
    for (c = 2; c < d; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (l = noccb-1; l > 1; l--) {
    for (k = l-1; k > 0; k--) {
    for (j = k-1; j > -1; j--) {
    for (i = j-1; i > -2; i--) {
        ijklabcd += 1;
        tmp = c4bbbb[ijklabcd]; 
        if(fabs(tmp) > numzero) 
        {
            tmp2 = t1xt3aaaa (i, j, k, l, a, b, c, d, noccb, nvirb, t1b, t3bbb);
            tmp2+= t2xt2aaaa (i, j, k, l, a, b, c, d, noccb, nvirb, t2bb);
            tmp2+= t1xt1xt2aaaa (i, j, k, l, a, b, c, d, noccb, nvirb, t1b, t2bb); 
            tmp2+= t1xt1xt1xt1aaaa (i, j, k, l, a, b, c, d, noccb, nvirb, t1b);

            tmp -= tmp2; 

            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+l)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+i)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+j)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+l)*noccb+j)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+i)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+k)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+j)*noccb+k)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+k)*noccb+i)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+l)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+k)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+i)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+i)*noccb+j)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+k)*noccb+j)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+i)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+i)*noccb+l)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)k*noccb+j)*noccb+l)*noccb+i)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)j*noccb+l)*noccb+i)*noccb+k)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+d)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+a)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+d)*nvirb+b)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+a)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+c)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+b)*nvirb+c)*nvirb+a)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+c)*nvirb+a)*nvirb+d)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+d)*nvirb+b)] = -tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+c)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+a)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+a)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+c)*nvirb+b)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+a)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+a)*nvirb+d)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+c)*nvirb+b)*nvirb+d)*nvirb+a)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+b)*nvirb+d)*nvirb+a)*nvirb+c)] =  tmp;
            t4bbbb[((((((((int64_t)l*noccb+i)*noccb+k)*noccb+j)*nvirb+d)*nvirb+a)*nvirb+c)*nvirb+b)] =  tmp;
        }
    }
    }
    }
    }
    }
    }
    }
    }

    m_ld = noccb * nvirb;
    ijkabc = -1;
    for (c = 2; c < nvira; c++) {
    for (b = 1; b < c; b++) {
    for (a = 0; a < b; a++) {
    for (k = nocca-1; k > 1; k--) {
    for (j = k-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijkabc += 1;
        ld = -1;
        for (d = 0; d < nvirb; d++) {
        for (l = noccb-1; l > -1; l--) {
            ld += 1;
            ijkabcld_c = ijkabc * m_ld + ld;
            tmp = c4aaab[ijkabcld_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = ut1xt3aaab (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t3aaa, t3aab);
                tmp2+= ut2xt2aaab (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t2aa, t2ab);
                tmp2+= ut1xt1xt2aaab (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t2aa, t2ab); 
                tmp2+= ut1xt1xt1xt1aaab (i, j, k, l, a, b, c, d, nvira, nvirb, t1a, t1b);
                tmp -= tmp2; 

                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+j)*nocca+k)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+k)*nocca+i)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+i)*nocca+j)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)i*nocca+k)*nocca+j)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)k*nocca+j)*nocca+i)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+a)*nvira+b)*nvira+c)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+b)*nvira+c)*nvira+a)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+c)*nvira+a)*nvira+b)*nvirb+d)] = -tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+a)*nvira+c)*nvira+b)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+c)*nvira+b)*nvira+a)*nvirb+d)] =  tmp;
                t4aaab[((((((((int64_t)j*nocca+i)*nocca+k)*noccb+l)*nvira+b)*nvira+a)*nvira+c)*nvirb+d)] =  tmp;
            }
        }
        }
    }
    }
    }
    }
    }
    }

    m_jklbcd = noccb*(noccb-1)*(noccb-2)/6 * nvirb*(nvirb-1)*(nvirb-2)/6;
    ia = -1;
    for (a = 0; a < nvira; a++) {
    for (i = nocca-1; i > -1; i--) {
        ia += 1;
        jklbcd = -1;
        for (d = 2; d < nvirb; d++) {
        for (c = 1; c < d; c++) {
        for (b = 0; b < c; b++) {
        for (l = noccb-1; l > 1; l--) {
        for (k = l-1; k > 0; k--) {
        for (j = k-1; j > -1; j--) {
            jklbcd += 1;
            iajklbcd_c = ia * m_jklbcd + jklbcd;
            tmp = c4abbb[iajklbcd_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = ut1xt3abbb (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t3abb, t3bbb);
                tmp2+= ut2xt2abbb (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t2ab, t2bb);
                tmp2+= ut1xt1xt2abbb (i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t2ab, t2bb); 
                tmp2+= ut1xt1xt1xt1aaab (j, k, l, i, b, c, d, a, nvirb, nvira, t1b, t1a);
                tmp -= tmp2; 

                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+k)*noccb+l)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+l)*noccb+j)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+j)*noccb+k)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+j)*noccb+l)*noccb+k)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+l)*noccb+k)*noccb+j)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+b)*nvirb+c)*nvirb+d)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+c)*nvirb+d)*nvirb+b)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+d)*nvirb+b)*nvirb+c)] = -tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+b)*nvirb+d)*nvirb+c)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+d)*nvirb+c)*nvirb+b)] =  tmp;
                t4abbb[((((((((int64_t)i*noccb+k)*noccb+j)*noccb+l)*nvira+a)*nvirb+c)*nvirb+b)*nvirb+d)] =  tmp;
            }
        }
        }
        }
        }
        }
        }
    }
    }


    m_klcd = noccb*(noccb-1)/2 * nvirb*(nvirb-1)/2;
    ijab = -1;
    for (b = 1; b < nvira; b++) {
    for (a = 0; a < b; a++) {
    for (j = nocca-1; j > 0; j--) {
    for (i = j-1; i > -1; i--) {
        ijab += 1;
        klcd  =-1;
        for (d = 1; d < nvirb; d++) {
        for (c = 0; c < d; c++) {
        for (l = noccb-1; l > 0; l--) {
        for (k = l-1; k > -1; k--) {
            klcd += 1;
            ijabklcd_c = ijab * m_klcd + klcd;
            tmp = c4aabb[ijabklcd_c]; 
            if(fabs(tmp) > numzero) 
            {
                tmp2 = ut1xt3aabb(i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t3aab, t3abb); 
                tmp2+= ut2xt2aabb(i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t2aa, t2ab, t2bb); 
                tmp2+= ut1xt1xt2aabb(i, j, k, l, a, b, c, d, nocca, nvira, noccb, nvirb, t1a, t1b, t2aa, t2ab, t2bb); 
                tmp2+= ut1xt1xt1xt1aabb(i, j, k, l, a, b, c, d, nvira, nvirb, t1a, t1b); 
                tmp -= tmp2; 

                t4aabb[((((((((int64_t)i*nocca+j)*noccb+k)*noccb+l)*nvira+a)*nvira+b)*nvirb+c)*nvirb+d)] =  tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+k)*noccb+l)*nvira+a)*nvira+b)*nvirb+d)*nvirb+c)] = -tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+l)*noccb+k)*nvira+a)*nvira+b)*nvirb+c)*nvirb+d)] = -tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+l)*noccb+k)*nvira+a)*nvira+b)*nvirb+d)*nvirb+c)] =  tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+k)*noccb+l)*nvira+b)*nvira+a)*nvirb+c)*nvirb+d)] = -tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+k)*noccb+l)*nvira+b)*nvira+a)*nvirb+d)*nvirb+c)] =  tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+l)*noccb+k)*nvira+b)*nvira+a)*nvirb+c)*nvirb+d)] =  tmp;
                t4aabb[((((((((int64_t)i*nocca+j)*noccb+l)*noccb+k)*nvira+b)*nvira+a)*nvirb+d)*nvirb+c)] = -tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+k)*noccb+l)*nvira+a)*nvira+b)*nvirb+c)*nvirb+d)] = -tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+k)*noccb+l)*nvira+a)*nvira+b)*nvirb+d)*nvirb+c)] =  tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+l)*noccb+k)*nvira+a)*nvira+b)*nvirb+c)*nvirb+d)] =  tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+l)*noccb+k)*nvira+a)*nvira+b)*nvirb+d)*nvirb+c)] = -tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+k)*noccb+l)*nvira+b)*nvira+a)*nvirb+c)*nvirb+d)] =  tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+k)*noccb+l)*nvira+b)*nvira+a)*nvirb+d)*nvirb+c)] = -tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+l)*noccb+k)*nvira+b)*nvira+a)*nvirb+c)*nvirb+d)] = -tmp;
                t4aabb[((((((((int64_t)j*nocca+i)*noccb+l)*noccb+k)*nvira+b)*nvira+a)*nvirb+d)*nvirb+c)] =  tmp; 
            }
        }
        }
        }
        }
    }
    }
    }
    }
}
