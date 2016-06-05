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

#ifndef _OPENCL
  #define __kernel
  #define __global
  #define __local
  #define M_PI_F (3.14159265359f)
  extern int __GID;
  extern int __LID;
  extern int __NGP;

  #define get_group_id(x) __GID
  #define get_local_id(x) __LID
  #define get_num_groups(x) kNgroups
  #define get_local_size(x) kGroupSize
#endif

#endif
