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
  cl::Program program = cl::Program(context, util::loadProgram("routines/MeanRadius.cl"));
  try {
    program.build(devices,"-Iutils");
    cl::CommandQueue queue(context);
  } catch(cl::Error err) {
    std::cout << "Build Status: " << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]);
    std::cout << std::endl;
    std::cout << "Build Options:\t" << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]);
    std::cout<< std::endl;
    std::cout << "Build Log:\t " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    std::cout << std::endl;
  }
  auto MeanRadius = cl::make_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int>(program, "MeanRadius");
  cl::Buffer d_x[7];      // device memory used for the input cluster vector
  cl::Buffer d_y[7];      // device memory used for the input cluster vector
  cl::Buffer d_z[7];      // device memory used for the input cluster vector
  cl::Buffer d_phi[7];      // device memory used for the input cluster vector
  cl::Buffer d_radii[7];         // device memory used for the input cluster vector
  cl::CommandQueue queue(context);

  for ( Event& e : events ) {
    vector<int> tLUT;
    /// Cluster sort
    int iL = 0;
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
      while (index(phi[iC],z[iC],iL) > tLUT.size())
        tLUT.push_back(iC);
    }
    while (tLUT.size() <= kNz * kNphi ) tLUT.push_back(size);  // Fix LUT size

    try {

      d_x[iL] = cl::Buffer(context, begin(x), end(x), true);
      d_y[iL] = cl::Buffer(context, begin(y), end(y), true);
      d_z[iL] = cl::Buffer(context, begin(z), end(z), true);
      d_phi[iL] = cl::Buffer(context, begin(phi), end(phi), true);
      d_radii[iL]    = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * size);

      using std::chrono::high_resolution_clock;
      using std::chrono::microseconds;
      auto t0 = high_resolution_clock::now();

      MeanRadius(cl::EnqueueArgs(queue,cl::NDRange(1024)),
          d_x[iL],
          d_y[iL],
          d_radii[iL],
          size);

      queue.finish();

      auto t1 = high_resolution_clock::now();
      microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
      printf("\nThe kernels ran in %lli microseconds\n", total_ms.count());
    } catch (cl::Error err) {
      std::cout << "Exception\n";
      std::cerr << "ERROR: " << err.what() << "(" << err_code(err.err()) << ")" << std::endl;
    }
  }


  return 0;
}

