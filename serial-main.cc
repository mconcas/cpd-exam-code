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
#include <fstream>
using std::vector;
using std::begin;
using std::end;
using std::cout;
using std::endl;
using std::array;
using std::ios;

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;
constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};
constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

void ComputeVertex(int lenIds,  // Trusted cluster id on layer 0,1
    float* x0, float* y0, float* z0,    // Clusters layer0
    float* x1, float* y1, float* z1,    // Clusters layer1
    float* final_vertex                 // Vertex array
    );

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
    for (int bin = kNz * iPhi; bin < (iPhi + 1) * kNz; ++bin)
      coarseLUT[bin + 1] = coarseLUT[bin] + mult * GetNumberOfClustersBin(lut0, bin);
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

int GoodTracklet(int *mc0, int *mc1, int *cls0, int *cls1, int iD) {
  if (mc0[cls0[iD]] == mc1[cls1[iD]] && mc0[cls0[iD]] > -1)
    return mc0[cls0[iD]];
  else
    return -1;
}

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

      while (index(phi[iC],z[iC],iL) > int(tLUT.size())) {
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
    vcid0[iL].resize(ntrkls, -1);
    vcid1[iL].resize(ntrkls, -1);
  }

  array<array<int, kNphi * kNz +1>, 6> tLUT;
  for (int iL = 0; iL < 6; ++iL) {
    PopulateCoarseLUT(tLUT[iL].data(),LUT[iL].data(),LUT[iL + 1].data());
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
      for (int iZ0 = 0; iZ0 < kNz; ++iZ0) {
        const int cls0 = GetNumberOfClustersPhiZ(LUT[iL].data(),iPhi0,iZ0);
        int offset = tLUT[iL][iPhi0 * kNz + iZ0];
        for (int iPhi1 = iPhi0 - 1 ; iPhi1 <= iPhi0 + 1; ++iPhi1) {
          const int iP1 = iPhi1 & (kNphi - 1); /// Adjust index
          for (int iZ1 = 0; iZ1 < kNz; ++iZ1) {
            const int cls1 = GetNumberOfClustersPhiZ(LUT[iL + 1].data(),iP1,iZ1);
            for (int iC0 = LUT[iL][iPhi0 * kNz + iZ0]; iC0 < LUT[iL][iPhi0 * kNz + iZ0 + 1]; ++iC0) {
              /// Loop over upper layer Phi granularity
              const int iC0_norm = (iC0 - LUT[iL][iPhi0 * kNz + iZ0]);
              const float& x_0 = vX[iL][iC0];
              const float& y_0 = vY[iL][iC0];
              const float& z_0 = vZ[iL][iC0];
              int idx_trkl = (iC0_norm * cls1) + offset;
              for (int iC1 = LUT[iL + 1][iP1 * kNz + iZ1] ; iC1 < LUT[iL + 1][iP1 * kNz + iZ1 + 1]; ++iC1) {
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
            offset += cls1 * cls0;
          }
        }
      }
    }
  }

#ifdef CHECK_OUTPUT
  std::ofstream myfile("/tmp/serial_trackleter.txt",ios::out);
  for (int iL = 0; iL < 6; ++iL) {
    for (size_t iT = 0; iT < vtId0[iL].size(); ++iT)
      myfile << vtId0[iL][iT] << "\t" << vtId1[iL][iT] << "\n";
    myfile << endl;
  }
  myfile.close();
#endif

  for (int iL = 0; iL < 5; ++iL) {
    int*   id0_1 = vtId1[iL].data();
    float* tphi0 = vtphi[iL].data();
    float* dzdr0 = vtdzdr[iL].data();
    int*   id1_0 = vtId0[iL + 1].data();
    float* tphi1 = vtphi[iL + 1].data();
    float* dzdr1 = vtdzdr[iL + 1].data();
    int*   lut0 = LUT[iL].data();
    int*   lut1 = LUT[iL + 1].data();
    int*   neigh0_1 = vcid1[iL].data();
    int*   neigh1_0 = vcid0[iL + 1].data();
    int* coarseLUT0 = tLUT[iL].data();
    int* coarseLUT1 = tLUT[iL + 1].data();

    for (int bin0 = 0; bin0 < kNz * kNphi; bin0++) {
      const int phi0 = (bin0) / kNz;
      const int z0 = (bin0) % kNz;
      const int t0 = coarseLUT0[bin0];
      const int ncls0 = GetNumberOfClustersPhiZ(lut0, phi0, z0);

#pragma omp parallel for // uncomment to parallelise
      for (int bin1 = 0; bin1 < 3 * kNz; bin1++) {
        const int phi1 = (bin1 / kNz) + phi0 - 1;
        const int z1 = (bin1 % kNz);
        const int t01 = t0 + FirstTrackletForTheBinOnLayer1(lut0,lut1,phi0,z0,phi1,z1);
        const int t01_next = t01 + ncls0 * GetNumberOfClustersPhiZ(lut1,phi1,z1);

        const int t1 = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1];
        const int t1_next = coarseLUT1[(phi1 & (kNphi - 1)) * kNz + z1 + 1];
        for (int iT01 = t01; iT01 < t01_next; ++iT01) {
          for (int iT1 = t1; iT1 < t1_next; ++iT1) {
            const bool flag = (id0_1[iT01] == id1_0[iT1]) &&
                              (fabs(dzdr0[iT01] - dzdr1[iT1]) + fabs(tphi0[iT01] - tphi1[iT1]) < kDiagonalTol);
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

  auto t1 = high_resolution_clock::now();
  microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
  cout<<" Event: "<<e.GetId()<<" - the vertexing ran in "<<total_ms.count()<<" microseconds"<<endl;

#ifdef CHECK_OUTPUT
  std::ofstream myfileCF("/tmp/serial_cellfinder.txt",ios::out);
  for (int iL = 0; iL < 6; ++iL) {
    for (size_t iT = 0; iT < vcid0[iL].size(); ++iT)
      myfileCF << vcid0[iL][iT] << "\t" << vcid1[iL][iT] << "\n";
    myfileCF << endl;
  }
  myfileCF.close();
#endif

  for (int iL = 0; iL < 6; ++iL) {
    int good = 0,fake=0;
    for (size_t iT = 0; iT < vtId0[iL].size(); ++iT) {
      if (vMcl[iL][vtId0[iL][iT]] == vMcl[iL + 1][vtId1[iL][iT]]) good++;
      else fake++;
    }
    cout << "\n\tLayer " << iL << ": fakes " << double(fake) / vtId0[iL].size();
    cout << ", goods: " << double(good) / vtId0[iL].size() << endl;
  }

  vector<float> trusted_x[2];
  vector<float> trusted_y[2];
  vector<float> trusted_z[2];
  int good = 0;
  for (size_t iT = 0; iT < vtId0[0].size(); ++iT) {
    int assoc_id = vcid1[0][iT];
    if (assoc_id > -1) {
      if (vcid0[1][assoc_id] > -1 && vcid1[1][assoc_id] > -1) {
        int assoc_id1 = vcid1[1][assoc_id];
        if (vcid0[2][assoc_id1] > -1 && vcid1[2][assoc_id1] > -1) {
          if (GoodTracklet(vMcl[0].data(), vMcl[1].data(), vtId0[0].data(), vtId1[0].data(), iT) >= 0) good++;
          trusted_x[0].push_back(vX[0][vtId0[0][iT]]);
          trusted_y[0].push_back(vY[0][vtId0[0][iT]]);
          trusted_z[0].push_back(vZ[0][vtId0[0][iT]]);
          trusted_x[1].push_back(vX[1][vtId1[0][iT]]);
          trusted_y[1].push_back(vY[1][vtId1[0][iT]]);
          trusted_z[1].push_back(vZ[1][vtId1[0][iT]]);
        }
      }
    }
  }
  cout << "Good: " << good << endl;
  /// Vertex Finding and comparison
  auto mcV = e.GetVertex();
  float vtx[3];
  ComputeVertex(trusted_x[0].size(),trusted_x[0].data(),trusted_y[0].data(),trusted_z[0].data(),
      trusted_x[1].data(),trusted_y[1].data(),trusted_z[1].data(),vtx);
  cout << vtx[0] - mcV[0] << "\t" << vtx[1] - mcV[1] << "\t" << vtx[2] - mcV[2] << "\t" << trusted_x[0].size() << endl;
  return 0;
}

void ComputeVertex(int lenIds,  // Trusted cluster id on layer 0,1
    float* x0, float* y0, float* z0,    // Clusters layer0
    float* x1, float* y1, float* z1,    // Clusters layer1
    float* final_vertex                 // Vertex array
    )
{

  vector<Line> lines(lenIds);
  for ( int id = 0; id < lenIds; ++id ) {
    /// Reconstruct line from data
    lines[id].x[0] = x0[id];
    lines[id].x[1] = y0[id];
    lines[id].x[2] = z0[id];
    lines[id].c[0] = x1[id] - x0[id];
    lines[id].c[1] = y1[id] - y0[id];
    lines[id].c[2] = z1[id] - z0[id];
  }

  //TODO: More elaborate strategy to reject fakes at the vertexing level.

  VertexCandidate vtxcand;
  for (auto& l : lines)
    vtxcand.Add(l);
  vtxcand.ComputeClusterCentroid();
  vtxcand.GetVertex(final_vertex);
}


