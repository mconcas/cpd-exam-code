#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "Event.h"
#include "MeanRadius.h"
#include "Trackleter.h"
#include "Utilities.h"

using std::vector;
using std::cout;
using std::endl;

#define SIZE (1024)

int main(int argc, char** argv) {

  parse_args(argc, argv);
  vector<Event> events( load_data(argv[1]) );

  Event &e = events[0];

  auto index = [&](const cluster& i, int l) {
    return int(i.fP * kInvDphi) * kNz + int((i.fZ + kZ[l]) * kInvDz[l]);
  };
  for ( Event& e : events ) {
    //int   size[7] = {0};
    cluster* layer_clusters[7] = {nullptr};
    float* phi[7] = {nullptr};
    float radius[7] = {0.f};
    int   size[7] = {0};
    vector<float> phiv[7];

    vector<int> vLUT[7];
    int*        LUT[7] = {nullptr};
//#pragma omp parallel for
    for (int iL = 0; iL < 7; ++iL) {

      vector<int> &tLUT = vLUT[iL];
      tLUT.reserve(kNz * kNphi + 1);
      tLUT.push_back(0);
      LUT[iL] = tLUT.data();

      size[iL] = e.GetClustersFromLayer(iL).size();

      /* Cluster sort */
      std::sort(e.GetClustersFromLayer(iL).begin(), e.GetClustersFromLayer(iL).end(), [&](const cluster& i, const cluster& j) {
          return index(i,iL) < index(j,iL); });

      /* Lookup table fill */
      for (int iC = 0; iC < size[iL]; ++iC) {
        while (index(e.GetClustersFromLayer(iL)[iC],iL) > tLUT.size())
          tLUT.push_back(iC);
      }
      while (tLUT.size() <= kNz*kNphi ) tLUT.push_back(size[iL]);  // Fix LUT size
      /* Prepare data for offloading */
      layer_clusters[iL] = e.GetClustersFromLayer(iL).data();

//#pragma offload_transfer target(mic:0) in(layer_clusters[iL] : length(size[iL]) ALLOC RETAIN) in(LUT[iL] : length(kNz * kNphi) ALLOC RETAIN) signal(layer_clusters[iL])
    }


    tracklet* trkl[6];
    int       n_trkl[6] = {0};
    for (int iL = 0; iL < 6; ++iL) {
      int* lut0 = LUT[iL];
      int* lut1 = LUT[iL+1];
      cluster* clu0 = layer_clusters[iL];
      cluster* clu1 = layer_clusters[iL+1];
      n_trkl[iL] = n_tracklets(LUT[iL],LUT[iL+1]);
      trkl[iL] = (tracklet*)malloc(n_trkl[iL] * sizeof(tracklet));
      cout << n_trkl[iL] << "\t" << sizeof(tracklet) << "\t" << n_trkl[iL] * sizeof(tracklet) / (1024 * 1024) << endl;
#pragma offload target(mic:0) \
      in(trkl[iL] : length(n_trkl[iL]) ALLOC RETAIN) \
      in(clu0 : length(size[iL]) ALLOC RETAIN) \
      in(clu1 : length(size[iL+1]) ALLOC RETAIN) \
      in(lut0 : length(kNz * kNphi) ALLOC RETAIN) \
      in(lut1 : length(kNz * kNphi) ALLOC RETAIN)
      //wait(layer_clusters[iL],layer_clusters[iL+1])*/
      {
        Trackleter(clu0,lut0,clu1,lut1,kRadii[iL+1] - kRadii[iL],trkl[iL]);
      }
    }
  }

  return 0;
}
