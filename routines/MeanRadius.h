#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Definitions.h"

__KERNEL__
float MeanRadius(const float* x, const float* y, /*const float* phi,*/ const int size) {
  float radius = 0.f;
#pragma omp parallel for reduction(+:radius) 
  for (int iC = 0; iC < size; ++iC) {
  	
    radius += sqrt(x[iC] * x[iC] + y[iC] * y[iC]);
  }
  return radius / size;
}
