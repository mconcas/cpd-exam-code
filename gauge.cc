#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include "Event.h"
#include "definitions.h"
#include "util.hpp"

int __GID = 0;
int __LID = 0;

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

    int tot_tracklets = 0;

    /// Loop over layers
    for (int iL = 0; iL < 7; ++iL ) {
      vector<int>& tLUT = LUT[iL];
      vector<float>& x = vX[iL];
      vector<float>& y = vY[iL];
      vector<float>& z = vZ[iL];
      vector<float>& phi = vPhi[iL];

      x = e.GetLayer(iL).x;
      y = e.GetLayer(iL).y;
      z = e.GetLayer(iL).z;
      phi = e.GetLayer(iL).phi;
      const int size = x.size();

      /// Use an array of indexes to sort 4 arrays
      vector<int> idx(size);
      for (int iC = 0; iC < size; ++iC) idx[iC] = iC;
      std::sort(begin(idx),end(idx), [&](const int& i, const int& j) {
          return index(phi[i],z[i],iL) < index(phi[j],z[j],iL); }
          );

      for (int iC = 0; iC < size; ++iC) {
        x[iC] = e.GetLayer(iL).x[idx[iC]];
        y[iC] = e.GetLayer(iL).y[idx[iC]];
        z[iC] = e.GetLayer(iL).z[idx[iC]];
        phi[iC] = e.GetLayer(iL).phi[idx[iC]];
      }

      /// Fill the lookup-table
      for (int iC = 0; iC < size; ++iC) {
        while (index(phi[iC],z[iC],iL) > tLUT.size())
          tLUT.push_back(iC);
      }
      while (tLUT.size() <= kNz * kNphi ) tLUT.push_back(size);  // Fix LUT size
      cout<<"\tSize of layer "<<iL<<": "<<size<<"\t Sizeof tLUT: "<<tLUT.size()<<endl;
      cout<<"\tAvg clust/bin: "<<x.size()/(kNphi*kNz)<<endl;
    }
    for (int iLut = 0; iLut < 6; ++iLut) {
      int n_avg = 0;
      for (int iPhi = 0; iPhi < kNphi; ++iPhi) {
        for (int iZ = 0; iZ < kNz; ++iZ) {
          int first_z1 = firstTracklet(LUT[iLut].data(), LUT[iLut + 1].data(), iPhi, iZ+1);
          int first_z0 = firstTracklet(LUT[iLut].data(), LUT[iLut + 1].data(), iPhi, iZ);
          n_avg+= (first_z1 - first_z0);
        }
      }
      cout<<"\tAvg tracklets/bin between layer: "<<iLut<<" and "<<iLut+1<<" are: "<<n_avg/(kNz*kNphi)<<endl;
    }
  return 0;
}
