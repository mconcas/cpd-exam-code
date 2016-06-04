#include "definitions.h"

int get_nclusters_phi(__global int* restrict lut, int iPhi) {
    iPhi &= (kNphi - 1);
    return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int get_nclusters_phi_z(__global int* restrict lut, int iPhi, int iZ) {
    iPhi &= (kNphi - 1);
    return lut[iPhi * kNz + iZ + 1] - lut[iPhi * kNz + iZ];
};

int n_tracklets_local(__global int* restrict lut0, __global int* restrict lut1, int nphi) {
    int n = 0;
    for (int i = 0; i < nphi; ++i)
        n += get_nclusters_phi(lut0,i) * (get_nclusters_phi(lut1,i + 1) + get_nclusters_phi(lut1,i) + get_nclusters_phi(lut1,i - 1));
    return n;
};

int n_tracklet_phi_z(__global int* restrict lut0, __global int* restrict lut1, int nphi0, int nz0) {
    int n = n_tracklets_local(lut0,lut1,nphi0);

    int mult = (get_nclusters_phi(lut1,nphi0 + 1) + get_nclusters_phi(lut1,nphi0) + get_nclusters_phi(lut1,nphi0 - 1));
    for (int iZ0 = 0; iZ0 < nz0; ++iZ0) {
        n += get_nclusters_phi_z(lut0,nphi0,iZ0) * mult;
    }
    return n;
};

__kernel void CellFinder(
        __global int*   restrict id1_0,
        __global float* restrict phi0,
        __global float* restrict dzdr0,
        __global int*   restrict id0_1,
        __global float* restrict phi1,
        __global float* restrict dzdr1,
        __global int*   restrict lut0,
        __global int*   restrict lut1,
        __global int*   restrict lut2,
        __global int*   restrict neigh0,
        __global int*   restrict neigh1
        ) {

    /// Group ID: needed to search the LUT
    const int group_id = get_group_id(0);
    const int group_size = get_local_size(0);

    /// Local ID
    const int local_id = get_local_id(0);

    for (int iteration = group_id * group_size; iteration < kNz * kNphi; iteration += group_size * get_num_groups(0)) {
        const int current_phi = iteration / kNz;
        const int current_z = iteration % kNz;
        const int next_phi = (iteration + 1) / kNz;
        const int next_z = (iteration + 1) % kNz;

        const int first_tracklet = n_tracklet_phi_z(lut1,lut2,current_phi,current_z);
        const int last_tracklet = n_tracklet_phi_z(lut1,lut2,next_phi,next_z);

        for (int iT1 = first_tracklet; iT1 < last_tracklet; ++iT1) {
            const float dzdr_1 = dzdr1[iT1];
            const float phi_1 = phi1[iT1];
            const int id_1  = id1_0[iT1];

            for (int iPhi = current_phi - 1; iPhi <= current_phi + 1; ++iPhi) {
                const int iPhiN = iPhi & (kNphi - 1);

                int clusters_until_this_bin = 0;
                for (int iC = current_phi + 1 - iPhi; iC > 0; iC--)
                    clusters_until_this_bin += get_nclusters_phi(lut1,current_phi - iC);
                clusters_until_this_bin += lut1[(iPhiN * kNz) + current_z] - lut1[iPhiN * kNz];

                for (int iZ = 0; iZ < kNz; ++iZ) {
                    const int first_assoc = n_tracklet_phi_z(lut0,lut1,iPhiN,iZ) + get_nclusters_phi_z(lut0,iPhiN,iZ) * clusters_until_this_bin;
                    for (int iT0 = first_assoc; iT0 < first_assoc + get_nclusters_phi_z(lut0,iPhiN,iZ) * get_nclusters_phi_z(lut1,current_phi,current_z); ++iT0) {
                        const bool flag = (id_1 == id0_1[iT0]) && (fabs(dzdr0[iT0]-dzdr_1) < kDzDrTol) && (fabs(phi_1 -phi0[iT0]) < kDphiTol || (fabs(phi_1 -phi0[iT0]) - 2 * M_PI_F) < kDphiTol);
                        neigh0[iT0] = flag ? iT1 : -1;
                        neigh1[iT1] = flag ? iT0 : -1;
                    }
                }
            }
        }
    }
}

