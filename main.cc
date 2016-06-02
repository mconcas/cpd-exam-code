#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include "Event.h"
#include "definitions.h"
#include "cl.hpp"
#include "util.hpp"

using std::vector;
using std::begin;
using std::end;
using std::cout;
using std::endl;

inline int get_nclusters(const std::vector<int> lut, int iPhi) {
  iPhi &= (kNphi - 1);
  return lut[(iPhi + 1) * kNz] - lut[iPhi * kNz];
};


inline int n_tracklets(const std::vector<int> lut0, const std::vector<int> lut1, int nphi = kNphi) {
  int n = 0;
  for (int i = 0; i < nphi; ++i)
    n += get_nclusters(lut0,i) * (get_nclusters(lut1,i + 1) + get_nclusters(lut1,i) + get_nclusters(lut1,i - 1));
  return n;
};

#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_ACCELERATOR
#endif

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;
constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};
constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

char* err_code(cl_int);

int main(int argc, char** argv) {

  if( argv[1] == NULL ) {
    std::cerr<<"Please, provide a data file."<<std::endl;
    exit(EXIT_FAILURE);
  }

  vector<Event> events( load_data(argv[1]) );

  auto index = [&](const float phi, float z, int l) {
    return int(phi * kInvDphi) * kNz + int((z + kZ[l]) * kInvDz[l]);
  };

  /// OpenCL initialisation
  cl::Context context(DEVICE);
  auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
  cl::Program program_trackleter(context, util::loadProgram("routines/Trackleter.cl"));
  try {
    program_trackleter.build(devices,"-Iutils");
  } catch(cl::Error err) {
    std::cout << "Build Status: " << program_trackleter.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]);
    std::cout << std::endl;
    std::cout << "Build Options:\t" << program_trackleter.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]);
    std::cout<< std::endl;
    std::cout << "Build Log:\t " << program_trackleter.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    std::cout << std::endl;
  }
  cl::CommandQueue queue(context);

  auto Trackleter = cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,
       cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer, float>(program_trackleter, "Trackleter");

  cl::Buffer d_x[7];      // device memory used for the input cluster vector
  cl::Buffer d_y[7];      // device memory used for the input cluster vector
  cl::Buffer d_z[7];      // device memory used for the input cluster vector
  cl::Buffer d_LUT[7];      // device memory used for the input cluster vector
  cl::Buffer d_radii[7];  // device memory used for the input cluster vector

  cl::Buffer d_tid0[6];   // device memory used for the input/output tracklet id on layer 0
  cl::Buffer d_tid1[6];   // device memory used for the input/output tracklet id on layer 1
  cl::Buffer d_tdzdr[6];  // device memory used for the input/output dz/dr
  cl::Buffer d_tphi[6];    // device memory used for the input cluster vector
  cl::Buffer d_cid0[5];   // device memory used for the input/output cell id on layer 0
  cl::Buffer d_cid1[5];   // device memory used for the input/output cell id on layer 0

  /// Events loop
  //for ( Event& e : events ) {
  Event &e = events[0];

    std::array<vector<int>, 7> LUT;

    /// Cluster sort
    try {
      for (int iL = 0; iL < 7; ++iL ) {

        vector<int>& tLUT = LUT[iL];

        auto x = e.GetLayer(iL).x;
        auto y = e.GetLayer(iL).y;
        auto z = e.GetLayer(iL).z;
        auto phi = e.GetLayer(iL).phi;
        const int size = x.size();

        vector<int> idx(size);
        for (int iC = 0; iC < size; ++iC) idx[iC] = iC;
        std::sort(begin(idx),end(idx), [&](const int& i, const int& j) {
            return index(phi[i],z[i],iL) < index(phi[j],z[j],iL); });

        for (int iC = 0; iC < size; ++iC) {
          x[iC] = e.GetLayer(iL).x[idx[iC]];
          y[iC] = e.GetLayer(iL).y[idx[iC]];
          z[iC] = e.GetLayer(iL).z[idx[iC]];
          phi[iC] = e.GetLayer(iL).phi[idx[iC]];
        }
        /// Lookup table fill
        for (int iC = 0; iC < size; ++iC) {

          while (index(phi[iC],z[iC],iL) > tLUT.size()) {
            tLUT.push_back(iC);
          }
        }
        while (tLUT.size() <= kNz * kNphi ) tLUT.push_back(size);  // Fix LUT size
        int mean = 0;
        for (int i = 0; i < tLUT.size()-1; ++i) {
          mean += (tLUT[i+1] - tLUT[i]);
        }

        // cout<<"Size tLUT: "<<tLUT.size()<<" layer: "<<iL<<" Event: "<<e.GetId()<<endl;
        // cout<<"Size of layer "<<iL<<" is: "<<size<<endl;

        d_x[iL] = cl::Buffer(context, begin(x), end(x), true);
        d_y[iL] = cl::Buffer(context, begin(y), end(y), true);
        d_z[iL] = cl::Buffer(context, begin(z), end(z), true);
        d_LUT[iL] = cl::Buffer(context, begin(tLUT), end(tLUT), true);

      }

      //for (int iL = 0; iL < 6; iL ++) {
      int iL = 0;
        int ntrkls = n_tracklets(LUT[iL], LUT[iL+1]);
        vector<int> vtId0(ntrkls);
        vector<int> vtId1(ntrkls);
        vector<float> vtdzdr(ntrkls);
        vector<float> vtphi(ntrkls);

        d_tid0[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        d_tid1[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        d_tdzdr[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        d_tphi[iL] = cl::Buffer(context, begin(vtphi), end(vtphi), true);

        using std::chrono::high_resolution_clock;
        using std::chrono::microseconds;
        auto t0 = high_resolution_clock::now();

        Trackleter(cl::EnqueueArgs(queue,cl::NDRange(kNphi * kGroupSize),cl::NDRange(kGroupSize)),
            d_x[iL],
            d_y[iL],
            d_z[iL],
            d_LUT[iL],
            d_x[iL + 1],
            d_y[iL + 1],
            d_z[iL + 1],
            d_LUT[iL + 1],
            d_tid0[iL],
            d_tid1[iL],
            d_tphi[iL],
            d_tdzdr[iL],
            kRadii[iL + 1]-kRadii[iL]);

        queue.finish();

        auto t1 = high_resolution_clock::now();
        microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
        printf("\nThe kernels ran in %lli microseconds\n", total_ms.count());
      //}

    } catch (cl::Error err) {
      std::cout << "Exception\n";
      std::cerr << "ERROR: " << err.what() << "(" << err_code(err.err()) << ")" << std::endl;
    }

    /*for(int iL = 0; iL < 6; ++iL) {
      cout<<"\t combinatorial beteen layer: "<<iL<<" and layer: "<<iL+1<<" -> "<< \
      n_tracklets(LUT[iL], LUT[iL+1])<<endl;
      }*/

  //}

  return 0;
}
