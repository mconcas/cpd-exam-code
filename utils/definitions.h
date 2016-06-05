#ifndef _OPENCL
  #define _OPENCL
  #ifndef DEVICE
    #define DEVICE CL_DEVICE_TYPE_CPU
    //#define DEVICE CL_DEVICE_TYPE_ACCELERATOR
   #endif
#endif

#ifndef UTILITIES_H
#define UTILITIES_H

#define kPi    (3.14159265359f)
#define kTwoPi (2.f * 3.14159265359f)

#define kNz   8
#define kNphi 256

#define kDzDrTol 0.02f
#define kDphiTol 0.1f

#define kGroupSize 1024
#define kNgroups (kNz * kNphi / kGroupSize)

#endif
