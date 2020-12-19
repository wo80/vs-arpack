#pragma once

/*
 * NOTE:
 *
 * Putting the constants in a header means that each translation unit (c file)
 * including the header will get its own set of constants.
 * 
 */

static const logical c_true = TRUE_;
static const logical c_false = FALSE_;

static const integer c__0 = 0;
static const integer c__1 = 1;
static const integer c__2 = 2;
static const integer c__3 = 3;
static const integer c__4 = 4;

static const complex c_one = { 1.0f, 0.0f };
static const complex c_zero = { 0.0f, 0.0f };

static const real s_one = 1.0f;
static const real s_zero = 0.0f;
static const real s_m1 = -1.0f;

static const doublereal d_one = 1.0;
static const doublereal d_zero = 0.0;
static const doublereal d_m1 = -1.0;
static const doublereal d_23 = 0.66666666666666663;

static const doublecomplex z_one = { 1.0, 0.0 };
static const doublecomplex z_zero = { 0.0, 0.0 };

/* used in examples */

static const integer c__5 = 5;
static const integer c__6 = 6;
static const integer c__9 = 9;
static const integer c__25 = 25;
static const integer c__30 = 30;
static const integer c__250 = 250;
static const integer c__256 = 256;
static const integer c_n6 = -6;

static const complex c_two = { 2.0f, 0.0f };
static const complex c_four = { 4.0f, 0.0f };
static const complex c_six = { 6.0f, 0.0f };
static const complex c_ten = { 10.0f, 0.0f };

static const doublecomplex z_two = { 2.0, 0.0 };
static const doublecomplex z_four = { 4.0, 0.0 };
static const doublecomplex z_six = { 6.0, 0.0 };
static const doublecomplex z_ten = { 10.0, 0.0 };
