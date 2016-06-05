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

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;
constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};
constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

int get_nclusters(int* lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};

int n_tracklets(int* lut0, int* lut1, int nphi) {
  int n = 0;
  for (int i = 0; i < nphi; ++i)
    n += get_nclusters(lut0,i) * (get_nclusters(lut1,i + 1) + get_nclusters(lut1,i) + get_nclusters(lut1,i - 1));
  return n;
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
  }

  using std::chrono::high_resolution_clock;
  using std::chrono::microseconds;
  auto t0 = high_resolution_clock::now();

  /// Loop over the layers
  for (int iL = 0; iL < 6; ++iL) {
    /// Loop over the Phi granularity
    const float dr = kRadii[iL + 1] - kRadii[iL];
    for (int iPhi0 = 0; iPhi0 < kNphi; ++iPhi0) {
      /// Loop over clusters
      int idx_trkl = n_tracklets(LUT[iL].data(),LUT[iL+1].data(),iPhi0);
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
            idx_trkl++;
          }
        }
      }
    }
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

 /* int good = 0;
  int fake = 0;
  for (int iT = 0; iT < vtId0[1].size(); ++iT) {
    if (vcid0[1][iT] >= 0 && vcid1[1][iT] >= 0) {
      int idx = vcid0[1][iT];
      if (vMcl[0][vtId0[0][idx]] == vMcl[1][vtId1[0][idx]]) good++;
      else fake++;
    }
  }
  cout << "\n\tValidated tracklets: fakes " << fake;
  cout << ", goods: " << good<< endl;*/

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


