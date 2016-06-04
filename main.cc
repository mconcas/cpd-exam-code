#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include <string>
#include <omp.h>
#include "Event.h"
#include "VertexCandidate.h"
#include "definitions.h"
#include "util.hpp"

#ifndef _OPENCL
  #define VERSION "Serial"
  #include "Trackleter.cl"
  int __GID = 0;
  int __LID = 0;
#else
  #include "cl.hpp"
  #define VERSION "OpenCL"
  #ifndef DEVICE
  //#define DEVICE CL_DEVICE_TYPE_CPU
  #define DEVICE CL_DEVICE_TYPE_ACCELERATOR
  #endif
#endif

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

#ifdef _OPENCL
  /// OpenCL initialisation
  char* err_code(cl_int);

  /// Set context with a DEVICE
  cl::Context context(DEVICE);
  auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
  std::string s;
  devices[0].getInfo(CL_DEVICE_NAME, &s);
  std::cout << "\n\t\t" << VERSION << " vertexer running on: " << s << std::endl << std::endl;

  /// Create the Program, load and compile kernel
  cl::Program program_trackleter(context, util::loadProgram("routines/Trackleter.cl"));
  cl::Program program_cellfinder(context, util::loadProgram("routines/CellFinder.cl"));
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
  try {
    program_cellfinder.build(devices,"-Iutils");
  } catch(cl::Error err) {
    std::cout << "Build Status: " << program_cellfinder.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]);
    std::cout << std::endl;
    std::cout << "Build Options:\t" << program_cellfinder.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]);
    std::cout<< std::endl;
    std::cout << "Build Log:\t " << program_cellfinder.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    std::cout << std::endl;
  }

  // Enqueue the context
  cl::CommandQueue queue(context);

  /// Create kernel function
  auto Trackleter = cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,
       cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer, float>(program_trackleter, "Trackleter");

  cl::Buffer d_x[7];
  cl::Buffer d_y[7];
  cl::Buffer d_z[7];
  cl::Buffer d_LUT[7];
  cl::Buffer d_radii[7];
  cl::Buffer d_tid0[6];
  cl::Buffer d_tid1[6];
  cl::Buffer d_tdzdr[6];
  cl::Buffer d_tphi[6];
  cl::Buffer d_cid1[5];
#endif

  /// Loop over the vector of events
  //for ( Event& e : events ) {
    Event &e = events[0];
    std::array<vector<int>, 7> LUT;
    vector<float> vX[7];
    vector<float> vY[7];
    vector<float> vZ[7];
    vector<float> vPhi[7];

#ifdef _OPENCL
    try {
#endif

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
            return index(phi[i],z[i],iL) < index(phi[j],z[j],iL); });

        for (int iC = 0; iC < size; ++iC) {
          x[iC] = e.GetLayer(iL).x[idx[iC]];
          y[iC] = e.GetLayer(iL).y[idx[iC]];
          z[iC] = e.GetLayer(iL).z[idx[iC]];
          phi[iC] = e.GetLayer(iL).phi[idx[iC]];
        }

        /// Fill the lookup-table
        for (int iC = 0; iC < size; ++iC) {

          while (index(phi[iC],z[iC],iL) > tLUT.size()) {
            tLUT.push_back(iC);
          }
        }
        while (tLUT.size() <= kNz * kNphi ) tLUT.push_back(size);  // Fix LUT size

#ifdef _OPENCL
        d_x[iL] = cl::Buffer(context, begin(x), end(x), true);
        d_y[iL] = cl::Buffer(context, begin(y), end(y), true);
        d_z[iL] = cl::Buffer(context, begin(z), end(z), true);
        d_LUT[iL] = cl::Buffer(context, begin(tLUT), end(tLUT), true);
#endif
      }

      for (int iL = 0; iL < 6; iL ++) {
      // int iL = 0;

        /// Create data structures to save tracklets
        int ntrkls = numTracklets(LUT[iL].data(), LUT[iL+1].data(), kNphi);

        vector<int> vtId0(ntrkls);
        vector<int> vtId1(ntrkls);
        vector<float> vtdzdr(ntrkls);
        vector<float> vtphi(ntrkls);

#ifdef _OPENCL
        // d_tid0[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        // d_tid1[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        // d_tdzdr[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);
        // d_tphi[iL] = cl::Buffer(context, CL_MEM_ALLOC_HOST_PTR, sizeof(float) * ntrkls);

        d_tid0[iL] = cl::Buffer(context, begin(vtId0), end(vtId0), true);
        d_tid1[iL] = cl::Buffer(context, begin(vtId1), end(vtId1), true);
        d_tdzdr[iL] = cl::Buffer(context, begin(vtdzdr), end(vtdzdr), true);
        d_tphi[iL] = cl::Buffer(context, begin(vtphi), end(vtphi), true);

        using std::chrono::high_resolution_clock;
        using std::chrono::microseconds;
        auto t0 = high_resolution_clock::now();

        Trackleter(cl::EnqueueArgs(queue, cl::NDRange(kNphi * kGroupSize),
            cl::NDRange(kGroupSize)),
            d_x[iL],
            d_y[iL],
            d_z[iL],
            d_LUT[iL],
            d_x[iL+1],
            d_y[iL+1],
            d_z[iL+1],
            d_LUT[iL+1],
            d_tid0[iL],
            d_tid1[iL],
            d_tphi[iL],
            d_tdzdr[iL],
            kRadii[iL+1]-kRadii[iL]);

        queue.finish();

        auto t1 = high_resolution_clock::now();
        microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
        // printf("The kernels ran in %lli microseconds\n", total_ms.count());
        cout<<" Event: "<<e.GetId()<<" - the kernel ran in "<<total_ms.count()<<" microseconds"<<endl;
      }
    } catch (cl::Error err) {
      std::cout << "Exception\n";
      std::cerr << "ERROR: " << err.what() << "(" << err_code(err.err()) << ")" << std::endl;
    }
#endif
#ifndef _OPENCL
    using std::chrono::high_resolution_clock;
    using std::chrono::microseconds;
    auto t0 = high_resolution_clock::now();

    for (__GID = 0; __GID < kNphi; ++__GID) {
      for (__LID = 0; __LID < kGroupSize; ++__LID) {
        Trackleter( vX[iL].data(), vY[iL].data(), vZ[iL].data(), LUT[iL].data(), 
          vX[iL+1].data(), vY[iL+1].data(), vZ[iL+1].data(), LUT[iL+1].data(), 
          vtId0.data(), vtId1.data(), vtphi.data(), vtdzdr.data(),
          kRadii[iL+1]-kRadii[iL] );
      }
    }
    auto t1 = high_resolution_clock::now();

    microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
    cout<<" Event: "<<e.GetId()<<" - the kernel ran in "<<total_ms.count()<<" microseconds"<<endl;
#endif
#ifndef _OPENCL
  }
#endif

  /// Vertex Finding and comparison
  
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


