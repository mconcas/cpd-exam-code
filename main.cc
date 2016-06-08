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
#include "cl.hpp"
#include <fstream>

using std::vector;
using std::begin;
using std::end;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;

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
  cl_device_type DEVICE = CL_DEV_CPU;

  if( argv[2] != NULL ) {
    if ( atoi(argv[2]) == 1 ) {
      DEVICE = CL_DEV_CPU;
    } else if ( atoi(argv[2]) == 2 ) 
      DEVICE = CL_DEV_ACC;
  } else {
    cout<<"\n\tEmpty or invalid architecture specified, defaulting to CPU version. "<<endl;
    cout<<"\tRun: "<<argv[0]<<" "<<argv[1]<<" [ 1 (CPU) | 2 (ACC) ] to change architecture"<<endl<<endl;
  }

  vector<Event> events( load_data(argv[1]) );


  auto index = [&](const float phi, float z, int l) {
    return int(phi * kInvDphi) * kNz + int((z + kZ[l]) * kInvDz[l]);
  };

  /// OpenCL initialisation
  char* err_code(cl_int);

  /// Set context with a DEVICE
  cl::Context context(DEVICE);
  auto devices = context.getInfo<CL_CONTEXT_DEVICES>();
  std::string s;
  devices[0].getInfo(CL_DEVICE_NAME, &s);
  std::cout << "\n\t\tOpenCl vertexer running on: " << s << std::endl << std::endl;

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
  auto CellFinder = cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,
       cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer>(program_cellfinder, "CellFinder");

  cl::Buffer d_x[7];
  cl::Buffer d_y[7];
  cl::Buffer d_z[7];
  cl::Buffer d_LUT[7];
  cl::Buffer d_tid0[6];
  cl::Buffer d_tid1[6];
  cl::Buffer d_tdzdr[6];
  cl::Buffer d_tphi[6];
  cl::Buffer d_cid0[6];
  cl::Buffer d_cid1[6];

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
    while (int(tLUT.size()) <= kNz * kNphi) tLUT.push_back(size);  // Fix LUT size

    d_x[iL] = cl::Buffer(context, begin(x), end(x), true);
    d_y[iL] = cl::Buffer(context, begin(y), end(y), true);
    d_z[iL] = cl::Buffer(context, begin(z), end(z), true);
    d_LUT[iL] = cl::Buffer(context, begin(tLUT), end(tLUT), true);
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
    vcid0[iL].resize(ntrkls,-1);
    vcid1[iL].resize(ntrkls,-1);
  }

  using std::chrono::high_resolution_clock;
  using std::chrono::microseconds;
  auto t0 = high_resolution_clock::now();

  try {
    for (int iL = 0; iL < 6; iL++) {
      d_tid0[iL]  = cl::Buffer(context, begin(vtId0[iL]),  end(vtId0[iL]),  false);
      d_tid1[iL]  = cl::Buffer(context, begin(vtId1[iL]),  end(vtId1[iL]),  false);
      d_tdzdr[iL] = cl::Buffer(context, begin(vtdzdr[iL]), end(vtdzdr[iL]), false);
      d_tphi[iL]  = cl::Buffer(context, begin(vtphi[iL]),  end(vtphi[iL]),  false);

      Trackleter(
          cl::EnqueueArgs(queue,cl::NDRange(kNphi * kGroupSize), cl::NDRange(kGroupSize)),
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

      queue.flush();

      d_cid0[iL] = cl::Buffer(context, begin(vcid0[iL]), end(vcid0[iL]), false);
      d_cid1[iL] = cl::Buffer(context, begin(vcid1[iL]), end(vcid1[iL]), false);
      //TODO: use non-blocking I/O
      cl::copy(queue,d_tid0[iL],begin(vtId0[iL]),end(vtId0[iL]));
      cl::copy(queue,d_tid1[iL],begin(vtId1[iL]),end(vtId1[iL]));
    }
  } catch (cl::Error err) {
    std::cout << "Exception during the Trackleter execution.\n";
    std::cerr << "ERROR: " << err.what() << "(" << err_code(err.err()) << ")" << std::endl;
  }

  std::ofstream myfile("/tmp/parallel.txt",ios::out);
  for (int iL = 0; iL < 6; ++iL) {
    for (size_t iT = 0; iT < vtId0[iL].size(); ++iT)
      myfile << vtId0[iL][iT] << "\t" << vtId1[iL][iT] << "\n";
    myfile << endl;
  }
  myfile.close();

  try {
    for (int iL = 0; iL < 5; ++iL) {
      CellFinder(
          cl::EnqueueArgs(queue,cl::NDRange(kNphi * kGroupSize), cl::NDRange(kGroupSize)),
          d_tid1[iL],
          d_tphi[iL],
          d_tdzdr[iL],
          d_tid0[iL + 1],
          d_tphi[iL + 1],
          d_tdzdr[iL + 1],
          d_LUT[iL],
          d_LUT[iL + 1],
          d_LUT[iL + 2],
          d_cid1[iL],
          d_cid0[iL + 1]);
      cl::copy(queue,d_cid0[iL + 1],begin(vcid0[iL + 1]),end(vcid0[iL + 1]));
      cl::copy(queue,d_cid1[iL],begin(vcid1[iL]),end(vcid1[iL]));
    }
  } catch (cl::Error err) {
    std::cout << "Exception during the CellFinder execution.\n";
    std::cerr << "ERROR: " << err.what() << "(" << err_code(err.err()) << ")" << std::endl;
  }
  queue.finish();

  auto t1 = high_resolution_clock::now();
  microseconds total_ms = std::chrono::duration_cast<microseconds>(t1 - t0);
  cout<<" Event: "<<e.GetId()<<" - the vertexing ran in "<<total_ms.count()<<" microseconds"<<endl;

  /// Vertex Finding and comparison
  float vtx[3];
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
  for (size_t iT = 0; iT < vtId0[1].size(); ++iT) {
    if (vcid0[1][iT] >= 0 && vcid1[1][iT] >= 0) {
      int idx = vcid0[1][iT];
      if (vMcl[0][vtId0[0][idx]] == vMcl[1][vtId1[0][idx]]) good++;
      else fake++;
    }
  }
  cout << "\n\tValidated tracklets: fakes " << fake;
  cout << ", goods: " << good <<endl;
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

  VertexCandidate vtxcand;
#pragma omp parallel
  {
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
