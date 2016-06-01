#include "definitions.h"
__kernel
void MeanRadius(
    __global float  * x,
    __global float  * y,
    __global float  * radii,
    __const int length) {
  const int idx = get_global_id(0);
  if (idx < length)
    radii[idx] = sqrt(x[idx] * x[idx] + y[idx] * y[idx]);
}
