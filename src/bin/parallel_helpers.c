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

#ifdef _WIN32
  long        curslotno;
  SYSTEM_INFO info;
  GetSystemInfo(&info);
#endif

#ifdef _SC_NPROCESSORS_ONLN
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);

  if (nprocs < 1)
    return 0; /* failure */

  nprocs_max = sysconf(_SC_NPROCESSORS_CONF);

  if (nprocs_max < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

  return 1; /* success */

#elif defined _SC_NPROC_ONLN
  nprocs = sysconf(_SC_NPROC_ONLN);

  if (nprocs < 1)
    return 0; /* failure */

  nprocs_max = sysconf(_SC_NPROC_CONF);

  if (nprocs_max < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

  return 1; /* success */

#elif defined _SC_CRAY_NCPU
  nprocs = sysconf(_SC_CRAY_NCPU);

  if (nprocs < 1)
    return 0; /* failure */

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs;

  return 1; /* success */

#elif defined _WIN32
  num_procs = info.dwNumberOfProcessors;

  nprocs_max = 0;

  while (info.dwActiveProcessorMask) {
    info.dwActiveProcessorMask  >>= 1;
    info.dwActiveProcessorMask  &= 0x7FFFFFFF;
    ++nprocs_max;
  }

  *num_cores      = (int)nprocs;
  *num_cores_conf = (int)nprocs_max;

#else
  return 0;
#endif
}
