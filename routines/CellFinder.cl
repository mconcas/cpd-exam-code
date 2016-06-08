#include "definitions.h"

int GetNumberOfClustersPhi(__global int* restrict lut, int iPhi) {
    iPhi &= (kNphi - 1);
    return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int GetNumberOfClustersPhiZ(__global int* restrict lut, int iPhi, int iZ) {
    iPhi &= (kNphi - 1);
    return lut[iPhi * kNz + iZ + 1] - lut[iPhi * kNz + iZ];
};

int GetNumberOfClustersBin(__global int* restrict lut, int idx) {
    return lut[idx + 1] - lut[idx];
};

int FirstTrackletForTheBinOnLayer1(__global int* restrict lut0, __global int* restrict lut1, int phi0, int z0, int phi1, int z1) {
    int n = 0;
    const int numberOfClustersInBin0 = GetNumberOfClustersBin(lut0, (phi0 * kNz) + z0);
    for (int iPhi0 = phi0 - 1; iPhi0 < phi1; ++iPhi0)
        n += numberOfClustersInBin0 * GetNumberOfClustersPhi(lut1, iPhi0);
    phi1 &= (kNphi - 1);
    return n + (numberOfClustersInBin0 * (lut1[phi1 * kNz + z1] - lut1[phi1 * kNz]));
};

__kernel void CellFinder(
        __global int*   restrict id1_0,
        __global float* restrict tphi0,
        __global float* restrict dzdr0,
        __global int*   restrict id0_1,
        __global float* restrict tphi1,
        __global float* restrict dzdr1,
        __global int*   restrict lut0,
        __global int*   restrict lut1,
        __global int*   restrict neigh0_1,
        __global int*   restrict neigh1_0,
        __global int*   restrict coarseLUT0,
        __global int*   restrict coarseLUT1
        ) {

    /// Group ID: needed to search the LUT
    const int group_id = get_group_id(0);
    const int group_size = get_local_size(0);

    /// Local ID
    const int local_id = get_local_id(0);

    for (int bin0 = 0; bin0 < kNz * kNphi; bin0++) {
        const int phi0 = (bin0) / kNz;
        const int z0 = (bin0) % kNz;
        const int t0 = coarseLUT0[bin0];
        const int ncls0 = GetNumberOfClustersPhiZ(lut0, phi0, z0);

        for (int bin1 = group_id * group_size + local_id; bin1 < 3 * kNz; bin1 += get_num_groups(0) * group_size) {
            const int phi1 = (bin1 / kNz) + phi0 - 1;
            const int z1 = (bin1 % kNz);
            const int t01 = t0 + FirstTrackletForTheBinOnLayer1(lut0,lut1,phi0,z0,phi1,z1);
            const int t01_next = t01 + ncls0 * GetNumberOfClustersPhiZ(lut1,phi1,z1);

            const int t1 = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1];
            const int t1_next = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1 + 1];
            for (int iT01 = t01; iT01 < t01_next; iT01++){
                for (int iT1 = t1; iT1 < t1_next; iT1++) {
                    const bool flag = (id0_1[iT01] == id1_0[iT1]);// &&
                        //(fabs(dzdr0[iT01] - dzdr1[iT1]) + fabs(tphi0[iT01] - tphi1[iT1]) < kDiagonalTol);
                    // (fabs(dzdr0[iT01] - dzdr1[iT1]) < kDzDrTol) &&
                    // (fabs(tphi0[iT01] - tphi1[iT1]) < kDphiTol);
                    if (flag) {
                        neigh0_1[iT01] = iT1;
                        neigh1_0[iT1]  = iT01;
                    }
                }
            }
        }
    }
}

