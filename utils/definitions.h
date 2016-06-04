// #ifndef _OPENCL
// #define _OPENCL
// #endif
#ifndef UTILITIES_H
#define UTILITIES_H

#define kPi    (3.14159265359f)
#define kTwoPi (2.f * 3.14159265359f)

#define kNz   2
#define kNphi 256

#define kDzDrTol 0.02f
#define kDphiTol 0.1f

#define kGroupSize 32
#define kNgroups (kNz * kNphi / kGroupSize)

#ifndef _OPENCL
  #define __kernel
  #define __global
  #define __local
  #define M_PI_F (3.14159265359f)
  extern int __GID;
  extern int __LID;
  extern int __NGP;

  // inline int get_local_size(int) { return kGroupSize; }
#endif

#endif
