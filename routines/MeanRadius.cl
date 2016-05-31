#include "definitions.h"
__kernel
void MeanRadius(
    __global struct cluster* buffer,
    __global float  * radii,
    __const int length) {
  const int idx = get_global_id(0);
  if (idx < length)
    radii[idx] = sqrt(buffer[idx].fX * buffer[idx].fX + buffer[idx].fY * buffer[idx].fY);
}
