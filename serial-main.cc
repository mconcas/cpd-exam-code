#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include <string>
#include "Event.h"
#include "VertexCandidate.h"
#include "definitions.h"
#include <omp.h>

using std::vector;
using std::begin;
using std::end;
using std::cout;
using std::endl;
using std::array;

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;
constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};
constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

template<typename T> T clamp(T n, T lower, T upper) {
  return std::max(lower, std::min(n, upper));
}

int get_nclusters_phi(int* lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int get_nclusters_phi_z(int* lut, int iPhi, int iZ) {
  iPhi &= (kNphi - 1);
  return lut[iPhi * kNz + iZ + 1] - lut[iPhi * kNz + iZ];
};

int n_tracklets_phi(int* lut0, int* lut1, int nphi) {
  int n = 0;
  for (int i = 0; i < nphi; ++i)
    n += get_nclusters_phi(lut0,i) * (get_nclusters_phi(lut1,i + 1) + get_nclusters_phi(lut1,i) + get_nclusters_phi(lut1,i - 1));
  return n;
};

int n_tracklet_phi_z(int* lut0, int* lut1, int nphi0, int nz0) {
  int n = n_tracklets_phi(lut0,lut1,nphi0);

  int mult = (get_nclusters_phi(lut1,nphi0 + 1) + get_nclusters_phi(lut1,nphi0) + get_nclusters_phi(lut1,nphi0 - 1));
  for (int iZ0 = 0; iZ0 < nz0; ++iZ0) {
    n += get_nclusters_phi_z(lut0,nphi0,iZ0) * mult;
  }
  return n;
};

int GetNumberOfClustersPhi( int*  lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int GetNumberOfClustersPhiZ( int*  lut, int iPhi, int iZ) {
  iPhi &= (kNphi - 1);
  return lut[iPhi * kNz + iZ + 1] - lut[iPhi * kNz + iZ];
};

int GetNumberOfClustersBin( int*  lut, int idx) {
  return lut[idx + 1] - lut[idx];
};

int FirstTrackletForTheBinOnLayer0( int*  lut0,  int*  lut1, int binidx) {
  int n = 0;
  for (int i = 0; i < binidx; ++i) {
    const int iPhi = binidx / kNphi;
    n += GetNumberOfClustersBin(lut0, i) * (GetNumberOfClustersPhi(lut1,iPhi - 1) + \
        GetNumberOfClustersPhi(lut1, iPhi) + GetNumberOfClustersPhi(lut1, iPhi + 1));
  }
  return n;
};

void PopulateCoarseLUT(int* coarseLUT,  int* lut0,  int* lut1) {
  coarseLUT[0] = 0;
  for (int iPhi = 0; iPhi < kNphi; ++iPhi) {
    const int mult = GetNumberOfClustersPhi(lut1, iPhi - 1) + \
                     GetNumberOfClustersPhi(lut1, iPhi) + \
                     GetNumberOfClustersPhi(lut1, iPhi + 1);
    for (int bin = kNz * iPhi; bin < (iPhi + 1) * kNz; ++bin) {
      coarseLUT[bin + 1] = coarseLUT[bin] + mult * GetNumberOfClustersBin(lut0, bin);
    }
  }
};

int FirstTrackletForTheBinOnLayer1( int*  lut0,  int*  lut1, int phi0, int z0, int phi1, int z1) {
  int n = 0;
  const int numberOfClustersInBin0 = GetNumberOfClustersBin(lut0, (phi0 * kNz) + z0);
  for (int iPhi0 = phi0 - 1; iPhi0 < phi1; ++iPhi0)
    n += numberOfClustersInBin0 * GetNumberOfClustersPhi(lut1, iPhi0);
  phi1 &= (kNphi - 1);
  return n + (numberOfClustersInBin0 * (lut1[phi1 * kNz + z1] - lut1[phi1 * kNz]));
};

int main(int argc, char** argv) {
  if( argv[1] == NULL ) {
    std::cerr<<"Please, provide a data file."<<std::endl;
    exit(EXIT_FAILURE);
  }

  vector<Event> events( load_data(argv[1]) );


  auto index = [&](const float phi, float z, int l) {
    return int(phi * kInvDphi) * kNz + int((z + kZ[l]) * kInvDz[l]);
  };


  /// Loop over the vector of events
  //for ( Event& e : events ) {
  Event &e = events[0];
  std::array<vector<int>, 7> LUT;
  vector<float> vX[7];
  vector<float> vY[7];
  vector<float> vZ[7];
  vector<float> vPhi[7];
  vector<int>   vMcl[7];


  int tot_tracklets = 0;

  /// Loop over layers
  for (int iL = 0; iL < 7; ++iL ) {
    vector<int>& tLUT = LUT[iL];
    vector<float>& x = vX[iL];
    vector<float>& y = vY[iL];
    vector<float>& z = vZ[iL];
    vector<float>& phi = vPhi[iL];
    vector<int>&   mcl = vMcl[iL];

    x = e.GetLayer(iL).x;
    y = e.GetLayer(iL).y;
    z = e.GetLayer(iL).z;
    phi = e.GetLayer(iL).phi;
    mcl = e.GetLayer(iL).mcl;
    const int size = x.size();

    /// Use an array of indexes to sort 4 arrays
    vector<int> idx(size);
    for (int iC = 0; iC < size; ++iC) idx[iC] = iC;
    std::sort(begin(idx),end(idx), [&](const int& i, const int& j) {
        return index(phi[i],z[i],iL) < index(phi[j],z[j],iL); });

    for (int iC = 0; iC < size; ++iC) {
      x[iC] = e.GetLayer(iL).x[idx[iC]];
      y[iC] = e.GetLayer(iL).y[idx[iC]];
      z[iC] = e.GetLayer(iL).z[idx[iC]];
      phi[iC] = e.GetLayer(iL).phi[idx[iC]];
      mcl[iC] = e.GetLayer(iL).mcl[idx[iC]];
    }

    /// Fill the lookup-table
    for (int iC = 0; iC < size; ++iC) {

      while (index(phi[iC],z[iC],iL) > tLUT.size()) {
        tLUT.push_back(iC);
      }
    }
    while (tLUT.size() <= kNz * kNphi ) tLUT.push_back(size);  // Fix LUT size
  }


  vector<int>   vtId0[6];
  vector<int>   vtId1[6];
  vector<float> vtdzdr[6];
  vector<float> vtphi[6];
  vector<int> vtmc[6];
  vector<int> vcid0[6];
  vector<int> vcid1[6];
  for (int iL = 0; iL < 6; ++iL) {
    /// Create data structures to save tracklets
    int ntrkls = numTracklets(LUT[iL].data(), LUT[iL+1].data(), kNphi);
    vtId0[iL].resize(ntrkls);
    vtId1[iL].resize(ntrkls);
    vtdzdr[iL].resize(ntrkls);
    vtphi[iL].resize(ntrkls);
    vcid0[iL].resize(ntrkls);
    vcid1[iL].resize(ntrkls);
    vtmc[iL].resize(ntrkls);
  }

  using std::chrono::high_resolution_clock;
  using std::chrono::microseconds;
  auto t0 = high_resolution_clock::now();

  /// Loop over the layers
  for (int iL = 0; iL < 6; ++iL) {
    /// Loop over the Phi granularity
    const float dr = kRadii[iL + 1] - kRadii[iL];
    int m = 0;
    for (int iPhi0 = 0; iPhi0 < kNphi; ++iPhi0) {
      /// Loop over clusters
      int idx_trkl = n_tracklets_phi(LUT[iL].data(),LUT[iL+1].data(),iPhi0);
      for (int iC0 = LUT[iL][iPhi0 * kNz]; iC0 < LUT[iL][(iPhi0 + 1) * kNz]; ++iC0) {
        /// Loop over upper layer Phi granularity
        const float& x_0 = vX[iL][iC0];
        const float& y_0 = vY[iL][iC0];
        const float& z_0 = vZ[iL][iC0];
        for (int iPhi1 = iPhi0 - 1 ; iPhi1 <= iPhi0 + 1; ++iPhi1) {
          /// Adjust index
          const int iP1 = iPhi1 & (kNphi - 1);
          for (int iC1 = LUT[iL + 1][iP1 * kNz] ; iC1 < LUT[iL + 1][(iP1 + 1) * kNz]; ++iC1) {
            const float& x_1 = vX[iL + 1][iC1];
            const float& y_1 = vY[iL + 1][iC1];
            const float& z_1 = vZ[iL + 1][iC1];
            vtId0[iL][idx_trkl]  = iC0;
            vtId1[iL][idx_trkl]  = iC1;
            vtdzdr[iL][idx_trkl] = (z_1 - z_0) / dr;
            vtphi[iL][idx_trkl] = atan2(y_1 - y_0,x_1 - x_0) + kPi;
            vtmc[iL][idx_trkl] = vMcl[iL][iC0] == vMcl[iL + 1][iC1] ? vMcl[iL][iC0] : -1;
            if(vMcl[iL][iC0] == vMcl[iL + 1][iC1]) m++;
            idx_trkl++;
          }
        }
      }
    }
    cout << "\tGood tracklets: " << m << endl;
  }

  array<array<int, kNphi * kNz +1>, 6> tLUT;
  for (int iL = 0; iL < 6; ++iL) {
    PopulateCoarseLUT(tLUT[iL].data(),LUT[iL].data(),LUT[iL + 1].data());
    if (vtphi[iL].size() != tLUT[iL][kNz*kNphi])
      cout << "Error for tracklets " << iL << ": " << vtphi[iL].size() << "\t" << tLUT[iL][kNz*kNphi] << endl;
  }

  /*for (int i = 0; i < kNphi * kNz + 1; ++i) {
    if (!(i%kNz)) cout << endl;
    cout << tLUT[1][i] << "\t";
  }
  cout << endl;
  */

  for (int iL = 0; iL < 5; ++iL) {
    cout << "????? LAYER ?????" << iL << endl;
    int*   id0_1 = vtId1[iL].data();
    float* tphi0 = vtphi[iL].data();
    float* dzdr0 = vtdzdr[iL].data();
    int*   id1_0 = vtId0[iL + 1].data();
    float* tphi1 = vtphi[iL + 1].data();
    float* dzdr1 = vtdzdr[iL + 1].data();
    int*   lut0 = LUT[iL].data();
    int*   lut1 = LUT[iL + 1].data();
    int*   lut2 = LUT[iL + 2].data();
    int*   neigh0_1 = vcid1[iL].data();
    int*   neigh1_0 = vcid0[iL + 1].data();
    int* coarseLUT0 = tLUT[iL].data();
    int* coarseLUT1 = tLUT[iL + 1].data();

    int good = 0;

    for (int bin0 = 0; bin0 < kNz * kNphi; bin0++) {
      const int phi0 = (bin0) / kNz;
      const int z0 = (bin0) % kNz;
      const int phi0_next = (bin0 + 1) / kNz;
      const int z0_next = (bin0 + 1) % kNz;
      const int t0 = coarseLUT0[bin0];
      const int ncls0 = GetNumberOfClustersPhiZ(lut0, phi0, z0);
      for (int bin1 = 0; bin1 < 3 * kNz; bin1++) {
        const int phi1 = (bin1 / kNz) + phi0 - 1;
        const int z1 = (bin1 % kNz);
        const int t01 = t0 + FirstTrackletForTheBinOnLayer1(lut0,lut1,phi0,z0,phi1,z1);
        const int t01_next = t01 + ncls0 * GetNumberOfClustersPhiZ(lut1,phi1,z1);
        const int t1 = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1];
        const int t1_next = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1 + 1];
        for (int iT01 = t01; iT01 < t01_next; ++iT01) {
          int mcLabel = -1;
          if (vMcl[iL][vtId0[iL][iT01]] == vMcl[iL+1][vtId1[iL][iT01]]) {
            good++;
            mcLabel = vMcl[iL][vtId0[iL][iT01]];
          }
          for (int iT1 = t1; iT1 < t1_next; ++iT1) {
            const bool flag = (id0_1[iT01] == id1_0[iT1]);
            //                  (fabs(dzdr0[iT01] - dzdr1[iT1]) < kDzDrTol) && \
            //                  (fabs(tphi0[iT01] - tphi1[iT1]) < kDphiTol || fabs(tphi0[iT01] - tphi1[iT1]) - kTwoPi < kDphiTol);
            neigh0_1[iT01] = flag ? iT1  : -1;
            neigh1_0[iT1]  = flag ? iT01 : -1;
          }
        }
      }
    }
    cout << "GOOOODDE " << good << endl;
  }

  auto t1 = high_resolution_clock::now();
  microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
  cout<<" Event: "<<e.GetId()<<" - the vertexing ran in "<<total_ms.count()<<" microseconds"<<endl;

  for (int iL = 0; iL < 6; ++iL) {
    int good = 0,fake=0;
    for (size_t iT = 0; iT < vtId0[iL].size(); ++iT) {
      if (vMcl[iL][vtId0[iL][iT]] == vMcl[iL + 1][vtId1[iL][iT]]) good++;
      else fake++;
    }
    cout << "\n\tLayer " << iL << ": fakes " << double(fake) / vtId0[iL].size();
    cout << ", goods: " << double(good) / vtId0[iL].size() << endl;
  }

  int good = 0;
  int fake = 0;
  for (int iT = 0; iT < vtId0[1].size(); ++iT) {
    if (vcid0[1][iT] >= 0 && vcid1[1][iT] >= 0) {
      int idx = vcid0[1][iT];
      if (vMcl[0][vtId0[0][idx]] == vMcl[1][vtId1[0][idx]]) good++;
      else fake++;
    }
  }
  cout << "\n\tValidated tracklets: fakes " << fake;
  cout << ", goods: " << good<< endl;

  /// Vertex Finding and comparison
  float vtx[3];
  // computeVertex(vtx);
  return 0;
}

void computeVertex(int* id0, int* id1, int lenIds,  // Trusted cluster id on layer 0,1
    int* lut0, int* lut1,               // LUTs layer0, layer1
    float* x0, float* y0, float* z0,    // Clusters layer0
    float* x1, float* y1, float* z1,    // Clusters layer1
    float* final_vertex                 // Vertex array
    )
{
  int threads = 1;

  VertexCandidate vtxcand;
#pragma omp parallel
  {
    threads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int n = 0;
    VertexCandidate candidate;
#pragma omp for
    for ( int id = 0; id < lenIds; ++id ) {
      /// Reconstruct line from data
      Line l;
      l.x[0] = x0[lut0[id0[id]]];
      l.x[1] = y0[lut0[id0[id]]];
      l.x[2] = z0[lut0[id0[id]]];
      l.c[0] = x1[lut1[id1[id]]] - x0[lut0[id0[id]]];
      l.c[1] = y1[lut1[id1[id]]] - y0[lut0[id0[id]]];
      l.c[2] = z1[lut1[id1[id]]] - z0[lut0[id0[id]]];
      candidate.Add(l);
    }
#pragma omp critical
    {
      vtxcand.Add(candidate);
    }
  }
  vtxcand.ComputeClusterCentroid();
  vtxcand.GetVertex(final_vertex);
}


