#define CL_DEV_CPU CL_DEVICE_TYPE_CPU
#define CL_DEV_ACC CL_DEVICE_TYPE_GPU

#ifndef UTILITIES_H
#define UTILITIES_H

#define kPi    (3.14159265359f)
#define kTwoPi (2.f * 3.14159265359f)

#define kNz   16
#define kNphi 256

#define kDzDrTol 0.01f
#define kDphiTol 0.08f
#define kDiagonalTol 0.05f

#define kGroupSize 1024
#define kNgroups (kNz * kNphi / kGroupSize)

//#define CHECK_OUTPUT

#endif
