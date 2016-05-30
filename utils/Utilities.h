#ifndef _DEFINITIONS_
#define _DEFINITIONS_

# ifdef __INTEL_COMPILER
#  define __KERNEL__ __declspec(target (mic))
#  define ALLOC alloc_if(1)
#  define FREE free_if(1)
#  define RETAIN free_if(0)
#  define REUSE alloc_if(0)
# else
#  define __KERNEL__
# endif

#include <vector>
using std::vector;

constexpr float kPi = 3.14159265359f;
constexpr float kTwoPi = 2.f * 3.14159265359f;

__KERNEL__ constexpr int kNz = 16;
__KERNEL__ constexpr int kNphi = 256;

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;

constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};

__KERNEL__ constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

#pragma offload_attribute(push,target(mic))
struct cluster {
  float fX;
  float fY;
  float fZ;
  float fP;
};

struct tracklet {
  int i0,i1;
  float dzdr;
  int magic;
};
#pragma offload_attribute(pop)

__KERNEL__ int get_nclusters(const int* lut, int iPhi);
__KERNEL__ int n_tracklets(const int* lut0, const int* lut1, int nphi = kNphi);

void parse_args(int argc, char** argv);

//template<typename T> void print_elm(T t, const int& width) { cout<< left << setw(width) << t; }
//void DumpLUT(const vector<int>& LUT, const int size_x) {

#endif
