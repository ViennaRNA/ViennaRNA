/*
 * Adapted from "numcpus - Prints out the number of CPUs in your system"
 * http://csgsoft.doc.ic.ac.uk/numcpus/
 *
 * see also:
 * https://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine/150971#150971
 *
 */

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


int
num_proc_cores(int  *num_cores,
               int  *num_cores_conf)
{
  long        nprocs      = -1;
  long        nprocs_max  = -1;

#if defined(_WIN32)
  SYSTEM_INFO info;
  GetSystemInfo(&info);
#endif

#if defined(_SC_NPROCESSORS_ONLN)
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);

  if (nprocs < 1)
    return 0; /* failure */

  nprocs_max = sysconf(_SC_NPROCESSORS_CONF);

  if (nprocs_max < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

  return 1; /* success */

#elif defined(_SC_NPROC_ONLN)
  nprocs = sysconf(_SC_NPROC_ONLN);

  if (nprocs < 1)
    return 0; /* failure */

  nprocs_max = sysconf(_SC_NPROC_CONF);

  if (nprocs_max < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

  return 1; /* success */

#elif defined(_SC_CRAY_NCPU)
  nprocs = sysconf(_SC_CRAY_NCPU);

  if (nprocs < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs;

  return 1; /* success */

#elif defined(_WIN32)
  nprocs = info.dwNumberOfProcessors;

  nprocs_max = 0;

  while (info.dwActiveProcessorMask) {
    info.dwActiveProcessorMask  >>= 1;
    info.dwActiveProcessorMask  &= 0x7FFFFFFF;
    ++nprocs_max;
  }

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

  return 1; /* success */

#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__)
  int     mib[4];
  int     numCPU;
  size_t  len = sizeof(numCPU);

  /* set the mib for hw.ncpu */
  mib[0]  = CTL_HW;
  mib[1]  = HW_AVAILCPU; /* alternatively, try HW_NCPU; */

  /* get the number of CPUs from the system */
  sysctl(mib, 2, &numCPU, &len, NULL, 0);

  if (numCPU < 1) {
    mib[1] = HW_NCPU;
    sysctl(mib, 2, &numCPU, &len, NULL, 0);
    if (numCPU < 1)
      numCPU = 1;
  }

  *num_cores      = numCPU;
  *num_cores_conf = numCPU;

  return 1; /* success */

#else
  return 0;
#endif
}


int
max_user_threads(void)
{
  int threadm = -1;

#if defined(_WIN32)
  /* currently not implmeneted, so let's limit to some arbitrary value */
  threadm = 128;
#else
  threadm = (int)sysconf(_SC_CHILD_MAX);
#endif

  return threadm;
}
