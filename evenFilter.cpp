// Angular low-pass filter for the even-harmonic coupling.
//
// The even-harmonic source is a transverse dipole (p*W + i*gam/(H*ks)*dW, see
// solve_Sc_3D), so its spatial spectrum grows with k.  The orbit-averaged model
// keeps the on-axis Bessel factors at every angle, so nothing suppresses that
// growth at large angle and the radiated power grows with grid resolution
// (static SASE, dx = 20/10/5 um: h2 x1.8, h4 x2.3; the power inside ~1 mrad
// converges).  The filter removes the angles where the model is not valid:
//
//    F(k) = exp[ -(theta/theta_c)^16 ],   theta = |k| / (H*ks)
//
// F is real and even in k, i.e. a symmetric (self-adjoint) convolution, so
// applying the same F to the deposited source and to the field the particles
// read keeps push and deposit exact adjoints (energy exchange is unchanged).
// The field itself is never filtered in place, otherwise F would compound
// every step.  Odd harmonics are not filtered.
//
// Only a window around the beam is transformed: the source is zero outside the
// particles, and the push only needs the filtered field at the particles.  The
// window is the particle bounding box plus a margin of MARGIN_KC/k_c (k_c for the
// lowest even harmonic), rounded up to a power-of-two FFT with zero padding, so
// the cost follows the beam size, not the grid, and prime grid sizes (nx+1 = 251)
// no longer matter.
//
// Input ([Domain]):  even_filter = ON|OFF          (default ON)
//                    even_filter_angle = <urad>     (default below)
// Default cutoff:    gamR*theta_c = 0.3*sqrt(1 + (1+ue^2)*K0^2/2)
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <fftw3.h>
#include <mpi.h>
#include "mesh.h"

namespace {
   const double MARGIN_KC = 20.0;   // window margin in units of 1/k_c
   int marginX = 0, marginY = 0;    // [cells]
   std::vector<char> filtered;      // per harmonic index

   struct FFT {
      std::vector<cplx> buf;
      fftw_plan fwd = nullptr, bwd = nullptr;
   };
   std::map<std::pair<int,int>, FFT> ffts;                   // (nfx, nfy)
   std::map<std::tuple<int,int,int>, std::vector<double>> filts;   // (h, nfx, nfy)

   int pow2(int n) { int m = 1; while (m < n) m <<= 1; return m; }

   FFT &getFFT(int nfx, int nfy)
   {
      FFT &f = ffts[{nfx, nfy}];
      if (!f.fwd) {
         f.buf.assign(static_cast<size_t>(nfx)*nfy, cplx(0.0, 0.0));
         fftw_complex *b = reinterpret_cast<fftw_complex*>(f.buf.data());
         f.fwd = fftw_plan_dft_2d(nfy, nfx, b, b, FFTW_FORWARD,  FFTW_MEASURE);
         f.bwd = fftw_plan_dft_2d(nfy, nfx, b, b, FFTW_BACKWARD, FFTW_MEASURE);
      }
      return f;
   }

   const std::vector<double> &getFilter(const Domain &D, int h, int nfx, int nfy)
   {
      std::vector<double> &F = filts[std::make_tuple(h, nfx, nfy)];
      if (F.empty()) {
         const double kc = D.harmony[h]*D.ks*D.evenFilterTheta;
         const double norm = 1.0/(static_cast<double>(nfx)*nfy);
         F.resize(static_cast<size_t>(nfx)*nfy);
         for (int j = 0; j < nfy; ++j) {
            int jj = (j <= nfy/2) ? j : j - nfy;
            double ky = 2.0*M_PI*jj/(nfy*D.dy);
            for (int i = 0; i < nfx; ++i) {
               int ii = (i <= nfx/2) ? i : i - nfx;
               double kx = 2.0*M_PI*ii/(nfx*D.dx);
               double q = std::sqrt(kx*kx + ky*ky)/kc;
               F[static_cast<size_t>(j)*nfx + i] = std::exp(-std::pow(q, 16))*norm;
            }
         }
      }
      return F;
   }

   // Filter the window w of slice src (full nx*ny, x fastest) into dst (full size;
   // only the window is written).  src and dst may be the same array.
   void filterWindow(const Domain &D, int h, const cplx *src, cplx *dst, const EvenBox &w)
   {
      const int nx = D.nx;
      const int lx = w.i1 - w.i0, ly = w.j1 - w.j0;
      const int nfx = pow2(lx), nfy = pow2(ly);
      FFT &f = getFFT(nfx, nfy);
      std::fill(f.buf.begin(), f.buf.end(), cplx(0.0, 0.0));
      for (int j = 0; j < ly; ++j)
         std::copy(src + static_cast<size_t>(w.j0 + j)*nx + w.i0,
                   src + static_cast<size_t>(w.j0 + j)*nx + w.i1,
                   f.buf.begin() + static_cast<size_t>(j)*nfx);
      fftw_execute(f.fwd);
      const std::vector<double> &F = getFilter(D, h, nfx, nfy);
      for (size_t n = 0; n < f.buf.size(); ++n) f.buf[n] *= F[n];
      fftw_execute(f.bwd);
      for (int j = 0; j < ly; ++j)
         std::copy(f.buf.begin() + static_cast<size_t>(j)*nfx,
                   f.buf.begin() + static_cast<size_t>(j)*nfx + lx,
                   dst + static_cast<size_t>(w.j0 + j)*nx + w.i0);
   }
}

void setupEvenFilter(Domain &D)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   filtered.assign(D.numHarmony, 0);
   if (D.dimension != 3 || !D.evenFilterON) {
      if (myrank == 0 && D.dimension == 3)
         std::cout << "even-harmonic angular filter: OFF" << std::endl;
      return;
   }
   if (D.evenFilterTheta <= 0.0)
      D.evenFilterTheta = 0.3*std::sqrt(1.0 + 0.5*(1.0 + D.ue*D.ue)*D.K0*D.K0)/D.gamR;

   int hmin = 0;
   for (int h = 0; h < D.numHarmony; ++h)
      if (D.harmony[h] % 2 == 0) {
         filtered[h] = 1;
         if (hmin == 0 || D.harmony[h] < hmin) hmin = D.harmony[h];
      }
   if (hmin == 0) return;
   const double kc = hmin*D.ks*D.evenFilterTheta;
   marginX = static_cast<int>(std::ceil(MARGIN_KC/(kc*D.dx)));
   marginY = static_cast<int>(std::ceil(MARGIN_KC/(kc*D.dy)));

   if (myrank == 0)
      std::cout << "even-harmonic angular filter: ON, theta_c = " << D.evenFilterTheta*1e6
                << " urad, window margin = " << marginX << " x " << marginY << " cells"
                << " (grid Nyquist for H=2: " << M_PI/(std::max(D.dx, D.dy)*2.0*D.ks)*1e6
                << " urad)" << std::endl;
}

bool evenFiltered(const Domain &D, int h)
{
   return D.dimension == 3 && h < static_cast<int>(filtered.size()) && filtered[h];
}

// Particle bounding box of one slice (all species) plus the filter margin.
// An empty box (i0 >= i1) means no particles: nothing to filter.
EvenBox evenFilterBox(const Domain &D, int sliceI)
{
   int i0 = D.nx, i1 = -1, j0 = D.ny, j1 = -1;
   for (int s = 0; s < D.nSpecies; ++s) {
      const ptclHead *hd = D.particle[sliceI].head[s];
      if (hd == nullptr) continue;
      for (const ptclList *p = hd->pt; p != nullptr; p = p->next)
         for (size_t n = 0; n < p->x.size(); ++n) {
            int i = static_cast<int>(std::floor((p->x[n] - D.minX)/D.dx));
            int j = static_cast<int>(std::floor((p->y[n] - D.minY)/D.dy));
            i0 = std::min(i0, i);  i1 = std::max(i1, i + 1);
            j0 = std::min(j0, j);  j1 = std::max(j1, j + 1);
         }
   }
   EvenBox w;
   if (i1 < i0) { w.i0 = w.i1 = w.j0 = w.j1 = 0; return w; }
   w.i0 = std::max(0, i0 - marginX);  w.i1 = std::min(D.nx, i1 + 1 + marginX);
   w.j0 = std::max(0, j0 - marginY);  w.j1 = std::min(D.ny, j1 + 1 + marginY);
   if (w.i1 <= w.i0 || w.j1 <= w.j0) w.i0 = w.i1 = w.j0 = w.j1 = 0;
   return w;
}

// Source: filter the window of one slice in place.
void applyEvenFilter(const Domain &D, int h, cplx *slice, const EvenBox &w)
{
   if (!evenFiltered(D, h) || w.i1 <= w.i0) return;
   filterWindow(D, h, slice, slice, w);
}

// Push: write the filtered window of src into dst (dst is full slice size;
// only the window is written, and only the window is ever read).
void applyEvenFilterCopy(const Domain &D, int h, const cplx *src, cplx *dst, const EvenBox &w)
{
   if (!evenFiltered(D, h) || w.i1 <= w.i0) return;
   filterWindow(D, h, src, dst, w);
}
