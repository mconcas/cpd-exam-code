#include "definitions.h"

#ifndef _OPENCL
inline int get_group_id(int) { return __GID; }
inline int get_local_id(int) { return __LID; }
#endif

int get_nclusters(__global int* lut, int iPhi) {
    iPhi &= (kNphi - 1);
    return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};
int n_tracklets(__global int* lut0, __global int* lut1, int nphi) {
    int n = 0;
    for (int i = 0; i < nphi; ++i)
        n += get_nclusters(lut0,i) * (get_nclusters(lut1,i + 1) + get_nclusters(lut1,i) + get_nclusters(lut1,i - 1));
    return n;
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
        __global int*   id0,
        __global int*   id1,
        __global float* phi,
        __global float* dzdr,
        const float     dr) {

    /// Group ID: needed to search the LUT
    const int group_size = get_local_size(0);

    /// Local ID
    const int local_id = get_local_id(0);

    for (int iteration = get_group_id(0); iteration < kNphi; iteration += get_num_groups(0)) {
        const int first_cluster0 = lut0[iteration * kNz];
        const int last_cluster0  = lut0[(iteration + 1) * kNz];

        const int first_trkl = n_tracklets(lut0,lut1,iteration);
        const int n_trkl = get_nclusters(lut1,iteration + 1) + get_nclusters(lut1,iteration) + get_nclusters(lut1,iteration - 1);

        __local float x_0;
        __local float y_0;
        __local float z_0;
        __local float x_1;
        __local float y_1;
        __local float z_1;

        for (int iC0 = first_cluster0 + local_id; iC0 < last_cluster0; iC0 += group_size) {
            x_0 = x0[iC0];
            y_0 = y0[iC0];
            z_0 = z0[iC0];
            int idx_trkl = first_trkl + (iC0 - first_cluster0) * n_trkl;
            for (int j = iteration - 1; j <= iteration + 1; ++j) {
                const int idx = j & (kNphi - 1);
                for (int iC1 = lut1[idx * kNz]; iC1 < lut1[(idx + 1) * kNz]; ++iC1) {
                    x_1 = x1[iC1];
                    y_1 = y1[iC1];
                    z_1 = z1[iC1];
                    id0[idx_trkl]  = iC0;
                    id1[idx_trkl]  = iC1;
                    dzdr[idx_trkl] = (z_1 - z_0) / dr;
                    phi[idx_trkl] = atan2(y_1 - y_0,x_1 - x_0) + M_PI_F;
                    idx_trkl++;
                }
            }
        }
    }
}
