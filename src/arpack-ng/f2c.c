#include "f2c.h"

double d_sign(double* a, double* b)
{
	double x = fabs(*a);
	return *b >= 0 ? x : -x;
}

double r_sign(float* a, float* b)
{
	double x = (*a >= 0 ? *a : -*a);
	return(*b >= 0 ? x : -x);
}

void r_cnjg(complex* r, complex* z)
{
	float zi = z->i;
	r->r = z->r;
	r->i = -zi;
}

void d_cnjg(doublecomplex* r, doublecomplex* z)
{
	double zi = z->i;
	r->r = z->r;
	r->i = -zi;
}

void c_div(complex* c, complex* a, complex* b)
{
	double ratio, den;
	double abr, abi, cr;

	if ((abr = b->r) < 0.)
		abr = -abr;
	if ((abi = b->i) < 0.)
		abi = -abi;
	if (abr <= abi)
	{
		if (abi == 0)
		{
			float af, bf;
			af = bf = abr;
			if (a->i != 0 || a->r != 0)
				af = 1.;
			c->i = c->r = af / bf;
			return;
		}
		ratio = (double)b->r / b->i;
		den = b->i * (1 + ratio * ratio);
		cr = (a->r * ratio + a->i) / den;
		c->i = (a->i * ratio - a->r) / den;
	}

	else
	{
		ratio = (double)b->i / b->r;
		den = b->r * (1 + ratio * ratio);
		cr = (a->r + a->i * ratio) / den;
		c->i = (a->i - a->r * ratio) / den;
	}
	c->r = cr;
}

void z_div(doublecomplex* c, doublecomplex* a, doublecomplex* b)
{
	double ratio, den;
	double abr, abi, cr;

	if ((abr = b->r) < 0.)
		abr = -abr;
	if ((abi = b->i) < 0.)
		abi = -abi;
	if (abr <= abi)
	{
		if (abi == 0) {
			if (a->i != 0 || a->r != 0)
				abi = 1.;
			c->i = c->r = abi / abr;
			return;
		}
		ratio = b->r / b->i;
		den = b->i * (1 + ratio * ratio);
		cr = (a->r * ratio + a->i) / den;
		c->i = (a->i * ratio - a->r) / den;
	}

	else
	{
		ratio = b->i / b->r;
		den = b->r * (1 + ratio * ratio);
		cr = (a->r + a->i * ratio) / den;
		c->i = (a->i - a->r * ratio) / den;
	}
	c->r = cr;
}
