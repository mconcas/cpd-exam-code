#include <iostream>
#include <vector>
#include <algorithm>
#include "Event.h"
#include "MeanRadius.h"
#include "Definitions.h"
#include "Utilities.h"

using std::vector;
using std::cout;
using std::endl;

#define SIZE (1024)

int main(int argc, char** argv) {

  parse_args(argc, argv);
  vector<Event> events( load_data(argv[1]) );

  Event &e = events[0];

  for ( Event& e : events ) {
    cluster* layer_clusters[7] = {nullptr};
    vector<float> phiv[7];
    float* phi[7] = {nullptr};
    float radius[7] = {0.f};
    int   size[7] = {0};

    auto index = [&](const cluster& i, int l) {
      return i.fP * kInvDphi * kNz + (i.fZ + kZ[l]) * kInvDz[l];
    };

    auto int_index = [&](const cluster& i, int l) {
      return int(i.fP * kInvDphi) * kNz + int((i.fZ + kZ[l]) * kInvDz[l]);
    };

    vector<int> LUT[7];
//#pragma omp parallel for
    for (int iL = 0; iL < 7; ++iL) {
      vector<int> &tLUT = LUT[iL];
      tLUT.push_back(0);
      // tLUT = vector<int>(kNz*kNphi, int_index())
      size[iL] = e.GetClustersFromLayer(iL).size();

      /* Cluster sort */
      std::sort(e.GetClustersFromLayer(iL).begin(), e.GetClustersFromLayer(iL).end(), [&](const cluster& i, const cluster& j) {
            return index(i,iL) < index(j,iL); });

      /* Lookup table fill */
      for (int iC = 0; iC < size[iL]; ++iC) {
        while (int_index(e.GetClustersFromLayer(iL)[iC],iL) > tLUT.size())
          tLUT.push_back(iC);
      }
      tLUT.push_back(size[iL]);
      while (tLUT.size() < kNz*kNphi ) tLUT.push_back(tLUT[tLUT.size() - 1]);  // Fix LUT size 
      DumpLUT(tLUT, kNz);                                                                          
      // cout<<"Last LUT index: "<<tLUT[tLUT.size()-1]<<"\t"<<"Layer length: "<<size[iL]<<endl;
      /* Prepare data for offloading */
      layer_clusters[iL] = e.GetClustersFromLayer(iL).data();
      // cout << tLUT.size() << "\t" << int_index(e.GetClustersFromLayer(iL)[size[iL]-1],iL)<< endl;
    }

    for (int iL = 0; iL < 7; ++iL) {
#pragma offload_transfer target(mic:0) in(layer_clusters[iL] : length(size[iL]) ALLOC RETAIN) signal(layer_clusters[iL])
#pragma offload target(mic:0) nocopy(layer_clusters[iL] : length(size[iL]) REUSE RETAIN) wait(layer_clusters[iL])
      {
        for (int iC = 0; iC < size[iL]; ++iC) {
          radius[iL] += sqrt(layer_clusters[iL][iC].fX * layer_clusters[iL][iC].fX + layer_clusters[iL][iC].fY * layer_clusters[iL][iC].fY);
        }
      }
    }
  }

  return 0;
}
