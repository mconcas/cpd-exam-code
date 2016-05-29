#ifndef _DEFINITIONS_
#define _DEFINITIONS_

#ifdef __INTEL_COMPILER

#define __KERNEL__ __declspec(target (mic))
#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)

#else
#define __KERNEL__
#endif

constexpr float kPi = 3.14159265359f;
constexpr float kTwoPi = 2.f * 3.14159265359f;

constexpr int kNz = 16;
constexpr int kNphi = 64;

constexpr float kDphi = kTwoPi / kNphi;
constexpr float kInvDphi = kNphi / kTwoPi;

constexpr float kZ[7] = {16.333f,16.333f,16.333f,42.140f,42.140f,73.745f,73.745f};
constexpr float kInvDz[7] = {
  0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 16.333f,0.5 * kNz / 42.140f,
  0.5 * kNz / 42.140f,0.5 * kNz / 73.745f,0.5 * kNz / 73.745f};

constexpr float kRadii[7] = {2.34,3.15,3.93,19.6,24.55,34.39,39.34};

#endif
