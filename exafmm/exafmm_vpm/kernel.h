#ifndef kernel_h
#define kernel_h
#include "exafmm.h"
#include "surface.hpp"
// #include "thirdparty/eigen"  // ???????
#include <Eigen/Core>
#include "rapidobj/rapidobj.hpp" // ???????
#include <iostream>
namespace exafmm {
const complex_t        I(0., 1.); //!< Imaginary unit
constexpr double                 EPS = 1e-16;
//!< L2 norm of vector X
template <class T>
inline T norm(T *X)
{
  return X[0] * X[0] + X[1] * X[1] + X[2] * X[2]; // L2 norm
}
//! Odd or even
inline int oddOrEven(int n)
{
  return (((n) & 1) == 1) ? -1 : 1;
}

//! i^2n
inline int ipow2n(int n)
{
  return (n >= 0) ? 1 : oddOrEven(n); // i^2n
}

//! Get r,theta,phi from x,y,z
template <class T>
void cart2sph(T *dX, T &r, T &theta, T &phi)
{
  r     = sqrt(norm(dX));                            // r = sqrt(x^2 + y^2 + z^2)
  theta = real_t(r) == 0.0 ? T(0) : acos(dX[2] / r); // theta = acos(z / r)
  phi   = atan2(dX[1], dX[0]);                       // phi = atan(y / x)
}

//! Spherical to cartesian coordinates
template <class T1, class T2>
void sph2cart(T1 r, T1 theta, T1 phi, T2 *spherical, T2 *cartesian)
{
  cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
                 + std::cos(theta) * std::cos(phi) / r * spherical[1] - std::sin(phi) / r / std::sin(theta) * spherical[2];
  cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
                 + std::cos(theta) * std::sin(phi) / r * spherical[1] + std::cos(phi) / r / std::sin(theta) * spherical[2];
  cartesian[2] = std::cos(theta) * spherical[0] // z component (not z itself)
                 - std::sin(theta) / r * spherical[1];
}

//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t *Ynm, complex_t *YnmTheta)
{
  real_t    x    = std::cos(alpha);                                             // x = cos(alpha)
  real_t    y    = std::sin(alpha);                                             // y = sin(alpha)
  real_t    invY = y == 0 ? 0 : 1 / y;                                          // 1 / y
  real_t    fact = 1;                                                           // Initialize 2 * m + 1
  real_t    pn   = 1;                                                           // Initialize Legendre polynomial Pn
  real_t    rhom = 1;                                                           // Initialize rho^m
  complex_t ei   = std::exp(I * beta);                                          // exp(i * beta)
  complex_t eim  = 1.0;                                                         // Initialize exp(i * m * beta)
  for (int m = 0; m < P; m++) {                                                 // Loop over m in Ynm
    real_t p      = pn;                                                         //  Associated Legendre polynomial Pnm
    int    npn    = m * m + 2 * m;                                              //  Index of Ynm for m > 0
    int    nmn    = m * m;                                                      //  Index of Ynm for m < 0
    Ynm[npn]      = rhom * p * eim;                                             //  rho^m * Ynm for m > 0
    Ynm[nmn]      = std::conj(Ynm[npn]);                                        //  Use conjugate relation for m < 0
    real_t p1     = p;                                                          //  Pnm-1
    p             = x * (2 * m + 1) * p1;                                       //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim;                 //  theta derivative of r^n * Ynm
    rhom *= rho;                                                                //  rho^m
    real_t rhon = rhom;                                                         //  rho^n
    for (int n = m + 1; n < P; n++) {                                           //  Loop over n in Ynm
      int npm = n * n + n + m;                                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                                  //   Index of Ynm for m < 0
      rhon /= -(n + m);                                                         //   Update factorial
      Ynm[npm]      = rhon * p * eim;                                           //   rho^n * Ynm
      Ynm[nmm]      = std::conj(Ynm[npm]);                                      //   Use conjugate relation for m < 0
      real_t p2     = p1;                                                       //   Pnm-2
      p1            = p;                                                        //   Pnm-1
      p             = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);      //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim; // theta derivative
      rhon *= rho;                                                              //   Update rho^n
    }                                                                           //  End loop over n in Ynm
    rhom /= -(2 * m + 2) * (2 * m + 1);                                         //  Update factorial
    pn = -pn * fact * y;                                                        //  Pn
    fact += 2;                                                                  //  2 * m + 1
    eim *= ei;                                                                  //  Update exp(i * m * beta)
  }                                                                             // End loop over m in Ynm
}

void evalMultipole(cplx rho, cplx alpha, cplx beta, multicomplex *Ynm, multicomplex *YnmTheta)
{
  cplx x    = cos(alpha);
  cplx y    = sin(alpha);
  cplx invY = (real_t(y) <= EPS) ? cplx(0, 0) : 1 / y;
  cplx fact = 1;
  cplx pn   = 1;
  cplx rhom = 1;
  // complex_t ei = std::exp(I * beta);
  // complex_t eim = 1.0;
  multicomplex ei  = multi_exp(beta);
  multicomplex eim = init_multicomplex(1.0, 0.0);
  for (int m = 0; m < P; m++) {
    cplx p   = pn;
    int  npn = m * m + 2 * m;
    int  nmn = m * m;
    // Ynm[npn] = rhom * p * eim;
    // Ynm[nmn] = std::conj(Ynm[npn]);
    Ynm[npn] = product(eim, rhom * p);
    Ynm[nmn] = conjugate(Ynm[npn]);
    cplx p1  = p;
    p        = x * (2 * m + 1) * p1;
    // YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim;
    YnmTheta[npn] = product(eim, rhom * (p - (m + 1) * x * p1) * invY);
    rhom *= rho;
    cplx rhon = rhom;
    for (int n = m + 1; n < P; n++) {
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      rhon /= -(n + m);
      // Ynm[npm] = rhon * p * eim;
      // Ynm[nmm] = std::conj(Ynm[npm]);
      Ynm[npm] = product(eim, rhon * p);
      Ynm[nmm] = conjugate(Ynm[npm]);
      cplx p2  = p1;
      p1       = p;
      p        = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      // YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;
      YnmTheta[npm] = product(eim, rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY);
      rhon *= rho;
    }
    rhom /= -(2 * m + 2) * (2 * m + 1);
    pn = -pn * fact * y;
    fact += 2;
    // eim *= ei;
    eim = product(eim, ei);
  }
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t *Ynm)
{
  real_t    x    = std::cos(alpha);                                    // x = cos(alpha)
  real_t    y    = std::sin(alpha);                                    // y = sin(alpha)
  real_t    fact = 1;                                                  // Initialize 2 * m + 1
  real_t    pn   = 1;                                                  // Initialize Legendre polynomial Pn
  real_t    invR = -1.0 / rho;                                         // - 1 / rho
  real_t    rhom = -invR;                                              // Initialize rho^(-m-1)
  complex_t ei   = std::exp(I * beta);                                 // exp(i * beta)
  complex_t eim  = 1.0;                                                // Initialize exp(i * m * beta)
  for (int m = 0; m < P; m++) {                                        // Loop over m in Ynm
    real_t p   = pn;                                                   //  Associated Legendre polynomial Pnm
    int    npn = m * m + 2 * m;                                        //  Index of Ynm for m > 0
    int    nmn = m * m;                                                //  Index of Ynm for m < 0
    Ynm[npn]   = rhom * p * eim;                                       //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn]   = std::conj(Ynm[npn]);                                  //  Use conjugate relation for m < 0
    real_t p1  = p;                                                    //  Pnm-1
    p          = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
    rhom *= invR;                                                      //  rho^(-m-1)
    real_t rhon = rhom;                                                //  rho^(-n-1)
    for (int n = m + 1; n < P; n++) {                                  //  Loop over n in Ynm
      int npm   = n * n + n + m;                                       //   Index of Ynm for m > 0
      int nmm   = n * n + n - m;                                       //   Index of Ynm for m < 0
      Ynm[npm]  = rhon * p * eim;                                      //   rho^n * Ynm for m > 0
      Ynm[nmm]  = std::conj(Ynm[npm]);                                 //   Use conjugate relation for m < 0
      real_t p2 = p1;                                                  //   Pnm-2
      p1        = p;                                                   //   Pnm-1
      p         = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
      rhon *= invR * (n - m + 1);                                      //   rho^(-n-1)
    }                                                                  //  End loop over n in Ynm
    pn = -pn * fact * y;                                               //  Pn
    fact += 2;                                                         //  2 * m + 1
    eim *= ei;                                                         //  Update exp(i * m * beta)
  }                                                                    // End loop over m in Ynm
}

void initKernel()
{
  NTERM = 3 * (P * (P + 1) / 2);
}

void P2PS(Cell *Ci, Cell *Cj)
{
  Body *Bi = Ci->BODY;
  Body *Bj = Cj->BODY;
  int   ni = Ci->NBODY;
  int   nj = Cj->NBODY;
  for (int i = 0; i < ni; i++) {
    // real_t pot = 0;
    real_t ax = 0;
    real_t ay = 0;
    real_t az = 0;
    for (int j = 0; j < nj; j++) {
      for (int d = 0; d < 3; d++)
        dX[d] = Bi[i].X[d] - Bj[j].X[d];
      real_t R2 = norm(dX);
      if (R2 != 0) {
        real_t S2     = 2 * Bj[j].radius * Bj[j].radius;                       //   2 * simga^2
        real_t RS     = R2 / S2;                                               //   R^2 / (2 * sigma^2)
        real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf(std::sqrt(RS)) // cutoff function for first term
                                                            - std::sqrt(4 / M_PI * RS) * exp(-RS));

        ax += (Bi[i].alpha[1] * Bj[j].alpha[2] - Bi[i].alpha[2] * Bj[j].alpha[1]) * cutoff; // x component of first term

        ay += (Bi[i].alpha[2] * Bj[j].alpha[0] - Bi[i].alpha[0] * Bj[j].alpha[2]) * cutoff; // y component of first term

        az += (Bi[i].alpha[0] * Bj[j].alpha[1] - Bi[i].alpha[1] * Bj[j].alpha[0]) * cutoff; // z component of first term
        cutoff = 0.25 / M_PI / R2 / R2 / std::sqrt(R2) * (3 * erf(std::sqrt(RS))            // cutoff function for second term
                                                          - (2 * RS + 3) * std::sqrt(4 / M_PI * RS) * exp(-RS)) *
                 (Bi[i].alpha[0] * dX[0] + Bi[i].alpha[1] * dX[1] + Bi[i].alpha[2] * dX[2]);

        ax += (Bj[j].alpha[1] * dX[2] - Bj[j].alpha[2] * dX[1]) * cutoff; // x component of second term

        ay += (Bj[j].alpha[2] * dX[0] - Bj[j].alpha[0] * dX[2]) * cutoff; // y component of second term

        az += (Bj[j].alpha[0] * dX[1] - Bj[j].alpha[1] * dX[0]) * cutoff; // z component of second term
      }
    }
    // Bi[i].p += pot;
    Bi[i].dadt[0] -= ax;
    Bi[i].dadt[1] -= ay;
    Bi[i].dadt[2] -= az;
  }
}


void P2P(Cell *Ci, Cell *Cj)
{
  Body *Bi = Ci->BODY;
  Body *Bj = Cj->BODY;
  int   ni = Ci->NBODY;
  int   nj = Cj->NBODY;
  for (int i = 0; i < ni; i++) {
    // real_t pot = 0;
    real_t ax = 0;
    real_t ay = 0;
    real_t az = 0;
    for (int j = 0; j < nj; j++) {
      for (int d = 0; d < 3; d++)
        dX[d] = Bi[i].X[d] - Bj[j].X[d];
      real_t R2 = norm(dX);
      if (R2 != 0) {
        real_t S2     = 2 * Bj[j].radius * Bj[j].radius; //    2 * sigma^2
        real_t RS     = R2 / S2;                         //    R^2 / (2 * simga^2)
        real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf(std::sqrt(RS)) - std::sqrt(4 / M_PI * RS) * exp(-RS));
        // real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2);
        ax += (dX[1] * Bj[j].alpha[2] - dX[2] * Bj[j].alpha[1]) * cutoff; // x component of curl G * cutoff
        ay += (dX[2] * Bj[j].alpha[0] - dX[0] * Bj[j].alpha[2]) * cutoff; // y component of curl G * cutoff
        az += (dX[0] * Bj[j].alpha[1] - dX[1] * Bj[j].alpha[0]) * cutoff; // z component of curl G * cutof
      }
    }
    // Bi[i].p += pot;
    Bi[i].velocity[0] -= ax;
    Bi[i].velocity[1] -= ay;
    Bi[i].velocity[2] -= az;
  }
}

Eigen::Vector3d P2cell( const double x,const double y,const double z , const Bodies& bodies) //
{
  Eigen::Vector3d Velocity;
  Velocity.setZero();

    // real_t pot = 0;
    real_t ax = 0;
    real_t ay = 0;
    real_t az = 0;
    for (int j = 0; j < bodies.size(); j++) {

      dX[0] = x- bodies[j].X[0];
      dX[1] = y- bodies[j].X[1];
      dX[2] = z- bodies[j].X[2];
      real_t R2 = norm(dX);
      //std::cout<<"for body "<<j<<" the norm is "<< R2 <<std::endl;
      if (R2 != 0) {
        //real_t S2     = 2 * Bj[j].radius * Bj[j].radius; //    2 * sigma^2
        real_t S2     = 2 * bodies[j].radius * bodies[j].radius; 
        real_t RS     = R2 / S2;                         //    R^2 / (2 * simga^2)
        real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf(std::sqrt(RS)) - std::sqrt(4 / M_PI * RS) * exp(-RS));
        // real_t cutoff = 0.25 / M_PI / R2 / std::sqrt(R2);
        ax += (dX[1] * bodies[j].alpha[2] - dX[2] * bodies[j].alpha[1]) * cutoff; // x component of curl G * cutoff
        ay += (dX[2] * bodies[j].alpha[0] - dX[0] * bodies[j].alpha[2]) * cutoff; // y component of curl G * cutoff
        az += (dX[0] * bodies[j].alpha[1] - dX[1] * bodies[j].alpha[0]) * cutoff; // z component of curl G * cutof
      }
    }



    Velocity[0] =- ax;
    Velocity[1] =- ay;
    Velocity[2] =- az;

    return Velocity;
  
}

void P2M(Cell *C)
{
  complex_t Ynm[P * P], YnmTheta[P * P];
  for (Body *B = C->BODY; B != C->BODY + C->NBODY; B++) {
    {
      for (int d = 0; d < 3; d++)
        dX[d] = B->X[d] - C->X[d];
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
      for (int n = 0; n < P; n++) {
        for (int m = 0; m <= n; m++) {
          int nm  = n * n + n + m;
          int nms = n * (n + 1) / 2 + m;
          for (int d = 0; d < 3; d++) {
            C->M[3 * nms + d] += B->alpha[d] * Ynm[nm];
          }
        }
      }
    }
  }
}

void M2M(Cell *Ci)
{
  complex_t Ynm[P * P], YnmTheta[P * P];
  for (Cell *Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++) {
    for (int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int j = 0; j < P; j++) {
      for (int k = 0; k <= j; k++) {
        int       jks = j * (j + 1) / 2 + k;
        complex_t M[3]{0, 0, 0};
        for (int n = 0; n <= j; n++) {
          for (int m = std::max(-n, -j + k + n); m <= std::min(k - 1, n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            int nm    = n * n + n - m;
            for (int d = 0; d < 3; d++)
              M[d] += Cj->M[3 * jnkms + d] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
          }
          for (int m = k; m <= std::min(n, j + k - n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            int nm    = n * n + n - m;
            for (int d = 0; d < 3; d++)
              M[d] += std::conj(Cj->M[3 * jnkms + d]) * Ynm[nm] * real_t(oddOrEven(k + n + m));
          }
        }
        for (int d = 0; d < 3; d++)
          Ci->M[3 * jks + d] += M[d];
      }
    }
  }
}

void M2L(Cell *Ci, Cell *Cj)
{
  complex_t Ynm2[4 * P * P];
  for (int d = 0; d < 3; d++)
    dX[d] = Ci->X[d] - Cj->X[d];
  real_t rho, alpha, beta;
  cart2sph(dX, rho, alpha, beta);
  evalLocal(rho, alpha, beta, Ynm2);
  for (int j = 0; j < P; j++) {
    real_t Cnm = oddOrEven(j);
    for (int k = 0; k <= j; k++) {
      int       jks = j * (j + 1) / 2 + k;
      complex_t L[3]{0, 0, 0};
      for (int n = 0; n < P; n++) {
        for (int m = -n; m < 0; m++) {
          int nms  = n * (n + 1) / 2 - m;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          for (int d = 0; d < 3; d++)
            L[d] += std::conj(Cj->M[3 * nms + d]) * Cnm * Ynm2[jnkm];
        }
        for (int m = 0; m <= n; m++) {
          int    nms  = n * (n + 1) / 2 + m;
          int    jnkm = (j + n) * (j + n) + j + n + m - k;
          real_t Cnm2 = Cnm * oddOrEven((k - m) * (k < m) + m);
          for (int d = 0; d < 3; d++)
            L[d] += Cj->M[3 * nms + d] * Cnm2 * Ynm2[jnkm];
        }
      }
      for (int d = 0; d < 3; d++)
        Ci->L[3 * jks + d] += L[d];
    }
  }
}

void L2L(Cell *Cj)
{
  complex_t Ynm[P * P], YnmTheta[P * P];
  for (Cell *Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++) {
    for (int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int j = 0; j < P; j++) {
      for (int k = 0; k <= j; k++) {
        int       jks = j * (j + 1) / 2 + k;
        complex_t L[3]{0, 0, 0};
        for (int n = j; n < P; n++) {
          for (int m = j + k - n; m < 0; m++) {
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nms  = n * (n + 1) / 2 - m;
            for (int d = 0; d < 3; d++)
              L[d] += std::conj(Cj->L[3 * nms + d]) * Ynm[jnkm] * real_t(oddOrEven(k));
          }
          for (int m = 0; m <= n; m++) {
            if (n - j >= abs(m - k)) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nms  = n * (n + 1) / 2 + m;
              for (int d = 0; d < 3; d++)
                L[d] += Cj->L[3 * nms + d] * Ynm[jnkm] * real_t(oddOrEven((m - k) * (m < k)));
            }
          }
        }
        for (int d = 0; d < 3; d++)
          Ci->L[3 * jks + d] += L[d];
      }
    }
  }
}
      
void L2P(Cell *Ci)
{
  multicomplex Ynm[P * P], YnmTheta[P * P];
  const double COMPLEX_STEP = 1e-32;
  for (Body *B = Ci->BODY; B != Ci->BODY + Ci->NBODY; B++) {
    real_t dadt = 0.0;
    for (int H_ind = 0; H_ind < 3; H_ind++) {
      complex_t dX[3];
      for (int d = 0; d < 3; d++)
        dX[d] = B->X[d] - Ci->X[d];
      dX[H_ind] += cplx(0, COMPLEX_STEP);
      complex_t spherical1[3];
      complex_t spherical2[3];
      complex_t spherical3[3];
      complex_t cartesian[3];
      for (int d = 0; d < 3; d++) {
        spherical1[d] = 0;
        spherical2[d] = 0;
        spherical3[d] = 0;
        cartesian[d]  = 0;
      }
      cplx r, theta, phi;
      cart2sph(dX, r, theta, phi);
      evalMultipole(r, theta, phi, Ynm, YnmTheta);
      for (int n = 0; n < P; n++) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        spherical1[0] += (product(init_from_C1(Ci->L[3 * nms + 0]), Ynm[nm])).A / r * n;
        spherical1[1] += (product(init_from_C1(Ci->L[3 * nms + 0]), YnmTheta[nm])).A;
        spherical2[0] += (product(init_from_C1(Ci->L[3 * nms + 1]), Ynm[nm])).A / r * n;
        spherical2[1] += (product(init_from_C1(Ci->L[3 * nms + 1]), YnmTheta[nm])).A;
        spherical3[0] += (product(init_from_C1(Ci->L[3 * nms + 2]), Ynm[nm])).A / r * n;
        spherical3[1] += (product(init_from_C1(Ci->L[3 * nms + 2]), YnmTheta[nm])).A;
        for (int m = 1; m <= n; m++) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          spherical1[0] += 2 * (product(init_from_C1(Ci->L[3 * nms + 0]), Ynm[nm])).A / r * n;
          spherical1[1] += 2 * (product(init_from_C1(Ci->L[3 * nms + 0]), YnmTheta[nm])).A;
          spherical1[2] += 2 * (product(product(init_from_C1(Ci->L[3 * nms + 0]), Ynm[nm]), init_from_C1(I))).A * m;
          spherical2[0] += 2 * (product(init_from_C1(Ci->L[3 * nms + 1]), Ynm[nm])).A / r * n;
          spherical2[1] += 2 * (product(init_from_C1(Ci->L[3 * nms + 1]), YnmTheta[nm])).A;
          spherical2[2] += 2 * (product(product(init_from_C1(Ci->L[3 * nms + 1]), Ynm[nm]), init_from_C1(I))).A * m;
          spherical3[0] += 2 * (product(init_from_C1(Ci->L[3 * nms + 2]), Ynm[nm])).A / r * n;
          spherical3[1] += 2 * (product(init_from_C1(Ci->L[3 * nms + 2]), YnmTheta[nm])).A;
          spherical3[2] += 2 * (product(product(init_from_C1(Ci->L[3 * nms + 2]), Ynm[nm]), init_from_C1(I))).A * m;
        }
      }
      for (int ind = 0; ind < 3; ind++)
        cartesian[ind] = 0;
      sph2cart(r, theta, phi, spherical1, cartesian);

      if (H_ind == 0) {
        B->velocity[2] -= 0.25 / M_PI * real_t(cartesian[1]);
        B->velocity[1] += 0.25 / M_PI * real_t(cartesian[2]);
        B->dadt[1] += 0.25 / M_PI * B->alpha[0] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt -= B->alpha[0] * imag(cartesian[1]) / COMPLEX_STEP;
        // B->dadt[2] -= 0.25 / M_PI * B->alpha[0] * imag(cartesian[1]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx1[3 * 0 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 1) {
        B->dadt[1] += 0.25 / M_PI * B->alpha[1] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt -= B->alpha[1] * imag(cartesian[1]) / COMPLEX_STEP;
        // B->dadt[2] -= 0.25 / M_PI * B->alpha[1] * imag(cartesian[1]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx2[3 * 0 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 2) {
        B->dadt[1] += 0.25 / M_PI * B->alpha[2] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt -= B->alpha[2] * imag(cartesian[1]) / COMPLEX_STEP;
        // B->dadt[2] -= 0.25 / M_PI * B->alpha[2] * imag(cartesian[1]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx3[3 * 0 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      }

      for (int ind = 0; ind < 3; ind++)
        cartesian[ind] = 0;
      sph2cart(r, theta, phi, spherical2, cartesian);

      if (H_ind == 0) {
        B->velocity[2] += 0.25 / M_PI * real_t(cartesian[0]);
        B->velocity[0] -= 0.25 / M_PI * real_t(cartesian[2]);
        B->dadt[0] -= 0.25 / M_PI * B->alpha[0] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt += B->alpha[0] * imag(cartesian[0]) / COMPLEX_STEP;
        // B->dadt[2] += 0.25 / M_PI * B->alpha[0] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx1[3 * 1 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 1) {
        B->dadt[0] -= 0.25 / M_PI * B->alpha[1] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt += B->alpha[1] * imag(cartesian[0]) / COMPLEX_STEP;
        // B->dadt[2] += 0.25 / M_PI * B->alpha[1] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx2[3 * 1 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 2) {
        B->dadt[0] -= 0.25 / M_PI * B->alpha[2] * imag(cartesian[2]) / COMPLEX_STEP;
        dadt += B->alpha[2] * imag(cartesian[0]) / COMPLEX_STEP;
        // B->dadt[2] += 0.25 / M_PI * B->alpha[2] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx3[3 * 1 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      }

      for (int ind = 0; ind < 3; ind++)
        cartesian[ind] = 0;
      sph2cart(r, theta, phi, spherical3, cartesian);

      if (H_ind == 0) {
        B->velocity[1] -= 0.25 / M_PI * real_t(cartesian[0]);
        B->velocity[0] += 0.25 / M_PI * real_t(cartesian[1]);
        B->dadt[0] += 0.25 / M_PI * B->alpha[0] * imag(cartesian[1]) / COMPLEX_STEP;
        B->dadt[1] -= 0.25 / M_PI * B->alpha[0] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx1[3 * 2 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 1) {
        B->dadt[0] += 0.25 / M_PI * B->alpha[1] * imag(cartesian[1]) / COMPLEX_STEP;
        B->dadt[1] -= 0.25 / M_PI * B->alpha[1] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx2[3 * 2 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      } else if (H_ind == 2) {
        B->dadt[0] += 0.25 / M_PI * B->alpha[2] * imag(cartesian[1]) / COMPLEX_STEP;
        B->dadt[1] -= 0.25 / M_PI * B->alpha[2] * imag(cartesian[0]) / COMPLEX_STEP;
        for (int ind = 0; ind < 3; ind++) {
          // B->dJdx3[3 * 2 + ind] += imag(cartesian[ind]) / COMPLEX_STEP;
        }
      }
    }
    B->dadt[2] += 0.25 / M_PI * dadt;
  }
}






} // namespace exafmm
#endif