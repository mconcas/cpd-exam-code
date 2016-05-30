#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Utilities.h"

__KERNEL__
void Trackleter(
    const cluster*  clu0,
    const int*      lut0,
    const cluster*  clu1,
    const int*      lut1,
    const float     dr,
    tracklet* trkl) {

#pragma omp parallel for
  for (int iP = 0; iP < kNphi; ++iP) {
    int idx_trkl = n_tracklets(lut0,lut1,iP);
    for (int iC0 = lut0[iP * kNz]; iC0 < lut0[(iP + 1) * kNz]; ++iC0) {
      const cluster& cl0 = clu0[iC0];
      for (int j = iP - 1; j <= iP + 1; ++j) {
        const int idx = j & (kNphi - 1);
        for (int iC1 = lut1[idx * kNz]; iC1 < lut1[(idx + 1) * kNz]; ++iC1) {
          const cluster& cl1 = clu1[iC1];
          const float dz = cl1.fZ - cl0.fZ;
          trkl[idx_trkl].i0 = iC0;
          trkl[idx_trkl].i1 = iC1;
          trkl[idx_trkl].dzdr = dz / dr;
          trkl[idx_trkl].magic = 0xdeadbeef;
          ++idx_trkl;
        }
      }
    }
  }
}
