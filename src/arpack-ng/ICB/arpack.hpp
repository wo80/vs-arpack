#ifndef __ARPACK_HPP__
#define __ARPACK_HPP__

#include "arpack_internal.h"

#include <complex.h>
#include <complex>

namespace arpack {
enum class which : int {
  /// 'LA' - compute the NEV largest (algebraic) eigenvalues
  largest_algebraic,
  /// 'SA' - compute the NEV smallest (algebraic) eigenvalues.
  smallest_algebraic,
  /// 'LM' - compute the NEV largest (in magnitude) eigenvalues.
  largest_magnitude,
  /// 'SM' - compute the NEV smallest (in magnitude) eigenvalues.
  smallest_magnitude,
  /// 'LR' - compute the NEV largest (real part) eigenvalues.
  largest_real,
  /// 'SR' - compute the NEV smallest (real part) eigenvalues.
  smallest_real,
  /// 'LI' - compute the NEV largest (imaginary part) eigenvalues.
  largest_imaginary,
  /// 'SI' - compute the NEV smallest (imaginary part) eigenvalues.
  smallest_imaginary,
  /// 'BE' - compute NEV eigenvalues, half from each end of the
  /// spectrum.  When NEV is odd, compute one more from the
  /// high end than from the low end.
  both_ends
};

enum class bmat : int {
  /// B = 'I' -> standard eigenvalue problem A*x = lambda*x
  identity,
  /// B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
  generalized
};

enum class howmny : int {
  /// 'A' Compute NEV Ritz vectors
  ritz_vectors,
  /// 'P' Compute NEV Schur vectors;
  schur_vectors,
  /// 'S' compute some of the Ritz vectors, specified by the logical array
  /// SELECT.
  ritz_specified
};

namespace internal {
#include "arpack.h"

inline char const* convert_to_char(which const option) {
  switch (option) {
    case which::largest_algebraic: {
      return "LA";
      break;
    }
    case which::smallest_algebraic: {
      return "SA";
      break;
    }
    case which::largest_magnitude: {
      return "LM";
      break;
    }
    case which::smallest_magnitude: {
      return "SM";
      break;
    }
    case which::largest_real: {
      return "LR";
      break;
    }
    case which::smallest_real: {
      return "SR";
      break;
    }
    case which::largest_imaginary: {
      return "LI";
      break;
    }
    case which::smallest_imaginary: {
      return "SI";
      break;
    }
    case which::both_ends: {
      return "BE";
      break;
    }
  }
  return "LM";
}

inline char const* convert_to_char(bmat const option) {
  return option == bmat::identity ? "I" : "G";
}

inline char const* convert_to_char(howmny const option) {
  switch (option) {
    case howmny::ritz_vectors: {
      return "A";
      break;
    }
    case howmny::schur_vectors: {
      return "P";
      break;
    }
    case howmny::ritz_specified: {
      return "S";
      break;
    }
  }
  return "A";
}
}  // namespace internal

inline void saupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::ssaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(int rvec, howmny const howmny_option, int* select, float* d,
                  float* z, int ldz, float sigma, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::sseupd_c(rvec, internal::convert_to_char(howmny_option), select, d,
                     z, ldz, sigma, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void saupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dsaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(int rvec, howmny const howmny_option, int* select, double* d,
                  double* z, int ldz, double sigma, bmat const bmat_option,
                  int n, which const ritz_option, int nev, double tol,
                  double* resid, int ncv, double* v, int ldv, int* iparam,
                  int* ipntr, double* workd, double* workl, int lworkl,
                  int& info) {
  internal::dseupd_c(rvec, internal::convert_to_char(howmny_option), select, d,
                     z, ldz, sigma, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::snaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(int rvec, howmny const howmny_option, int* select, float* dr,
                  float* di, float* z, int ldz,
                  float sigmar, float sigmai, float * workev,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, float tol, float* resid, int ncv, float* v, int ldv,
                  int* iparam, int* ipntr, float* workd, float* workl,
                  int lworkl, int& info) {
  internal::sneupd_c(rvec, internal::convert_to_char(howmny_option), select, dr,
                     di, z, ldz, sigmar, sigmai, workev,
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dnaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(int rvec, howmny const howmny_option, int* select,
                  double* dr, double* di, double* z, int ldz,
                  double sigmar, double sigmai, double * workev,
                  bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dneupd_c(rvec, internal::convert_to_char(howmny_option), select, dr,
                     di, z, ldz, sigmar, sigmai, workev,
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol,
                  std::complex<float>* resid, int ncv, std::complex<float>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, int lworkl,
                  float* rwork, int& info) {
  internal::cnaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<a_fcomplex*>(resid), ncv,
                     reinterpret_cast<a_fcomplex*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<a_fcomplex*>(workd),
                     reinterpret_cast<a_fcomplex*>(workl), lworkl,
                     rwork, &info);
}

inline void neupd(int rvec, howmny const howmny_option, int* select,
                  std::complex<float>* d, std::complex<float>* z, int ldz,
                  std::complex<float> sigma, std::complex<float>* workev,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, float tol, std::complex<float>* resid, int ncv,
                  std::complex<float>* v, int ldv, int* iparam, int* ipntr,
                  std::complex<float>* workd, std::complex<float>* workl,
                  int lworkl, float* rwork, int& info) {
  std::complex<float> sigma2 = sigma;
  internal::cneupd_c(rvec, internal::convert_to_char(howmny_option), select,
                     reinterpret_cast<a_fcomplex*>(d),
                     reinterpret_cast<a_fcomplex*>(z), ldz,
                     *reinterpret_cast<a_fcomplex*>(&sigma2),
                     reinterpret_cast<a_fcomplex*>(workev),
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<a_fcomplex*>(resid), ncv,
                     reinterpret_cast<a_fcomplex*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<a_fcomplex*>(workd),
                     reinterpret_cast<a_fcomplex*>(workl), lworkl,
                     rwork, &info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol,
                  std::complex<double>* resid, int ncv, std::complex<double>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, int lworkl,
                  double* rwork, int& info) {
  internal::znaupd_c(&ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<a_dcomplex*>(resid), ncv,
                     reinterpret_cast<a_dcomplex*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<a_dcomplex*>(workd),
                     reinterpret_cast<a_dcomplex*>(workl), lworkl,
                     rwork, &info);
}

inline void neupd(int rvec, howmny const howmny_option, int* select,
                  std::complex<double>* d, std::complex<double>* z, int ldz,
                  std::complex<double> sigma, std::complex<double>* workev,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, double tol, std::complex<double>* resid, int ncv,
                  std::complex<double>* v, int ldv, int* iparam, int* ipntr,
                  std::complex<double>* workd, std::complex<double>* workl,
                  int lworkl, double* rwork, int& info) {
  std::complex<double> sigma2 = sigma;
  internal::zneupd_c(rvec, internal::convert_to_char(howmny_option), select,
                     reinterpret_cast<a_dcomplex*>(d),
                     reinterpret_cast<a_dcomplex*>(z), ldz,
                     *reinterpret_cast<a_dcomplex*>(&sigma2),
                     reinterpret_cast<a_dcomplex*>(workev),
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<a_dcomplex*>(resid), ncv,
                     reinterpret_cast<a_dcomplex*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<a_dcomplex*>(workd),
                     reinterpret_cast<a_dcomplex*>(workl), lworkl,
                     rwork, &info);
}
}  // namespace arpack

#endif
