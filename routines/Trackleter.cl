#include "definitions.h"

int GetNumberOfClustersPhiZ(__global int*  lut, int iPhi, int iZ) {
  iPhi &= (kNphi - 1);
  return lut[iPhi * kNz + iZ + 1] - lut[iPhi * kNz + iZ];
};

__kernel void Trackleter(
        __global float* x0,
        __global float* y0,
        __global float* z0,
        __global int*   lut0,
        __global float* x1,
        __global float* y1,
        __global float* z1,
        __global int*   lut1,
        __global int*   coarseLUT,
        __global int*   id0,
        __global int*   id1,
        __global float* phi,
        __global float* dzdr,
        const float     dr) {

    /// Group ID: needed to search the LUT
    const int group_size = get_local_size(0);

    /// Local ID
    const int local_id = get_local_id(0);

    for (int iPhi0 = get_group_id(0); iPhi0 < kNphi; iPhi0 += get_num_groups(0)) {
      /// Loop over clusters
      for (int iZ0 = local_id; iZ0 < kNz; iZ0 += group_size) {
        const int cls0 = GetNumberOfClustersPhiZ(lut0,iPhi0,iZ0);
        int offset = coarseLUT[iPhi0 * kNz + iZ0];
        for (int iPhi1 = iPhi0 - 1 ; iPhi1 <= iPhi0 + 1; ++iPhi1) {
          const int iP1 = iPhi1 & (kNphi - 1); /// Adjust index
          for (int iZ1 = 0; iZ1 < kNz; ++iZ1) {
            const int cls1 = GetNumberOfClustersPhiZ(lut1,iP1,iZ1);
            for (int iC0 = lut0[iPhi0 * kNz + iZ0]; iC0 < lut0[iPhi0 * kNz + iZ0 + 1]; ++iC0) {
              /// Loop over upper layer Phi granularity
              const int iC0_norm = (iC0 - lut0[iPhi0 * kNz + iZ0]);
              const float x_0 = x0[iC0];
              const float y_0 = y0[iC0];
              const float z_0 = z0[iC0];
              int idx_trkl = (iC0_norm * cls1) + offset;
              for (int iC1 = lut1[iP1 * kNz + iZ1] ; iC1 < lut1[iP1 * kNz + iZ1 + 1]; ++iC1) {
                const float x_1 = x1[iC1];
                const float y_1 = y1[iC1];
                const float z_1 = z1[iC1];
                id0[idx_trkl]  = iC0;
                id1[idx_trkl]  = iC1;
                phi[idx_trkl] = atan2(y_1 - y_0,x_1 - x_0) + M_PI_F;
                dzdr[idx_trkl] = (z_1 - z_0) / dr;
                idx_trkl++;
              }
            }
            offset += cls1 * cls0;
          }
        }
      }
    }
}
