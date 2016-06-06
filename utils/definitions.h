#define CL_DEV_CPU CL_DEVICE_TYPE_CPU
#define CL_DEV_ACC CL_DEVICE_TYPE_ACCELERATOR

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
