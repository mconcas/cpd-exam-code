#ifdef __INTEL_COMPILER
#define __KERNEL__ __declspec(target (mic))
#else
#define __KERNEL__
#endif
