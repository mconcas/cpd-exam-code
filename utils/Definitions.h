#ifdef __INTEL_COMPILER

#define __KERNEL__ __declspec(target (mic))
#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)

#else
#define __KERNEL__
#endif
