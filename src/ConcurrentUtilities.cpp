#include "AmpGen/ConcurrentUtilities.h"

namespace AmpGen
{
    NamedParameter<size_t> getNumThreads()
    {
        static const unsigned int THREAD_MOD = 10;
        unsigned int maxThreads = std::thread::hardware_concurrency();
        unsigned int threadsToUse = maxThreads / THREAD_MOD;
        if (threadsToUse <= 0)
        {
            threadsToUse = 1;
        }

        NamedParameter<size_t> nThreads = NamedParameter<size_t>(
            "nThreads",
            threadsToUse,
            "Number of threads to use");
        INFO("Using " << nThreads.getVal() << " / " << maxThreads << " threads (DECREASE IN ConcurrentUtilities.cpp IF FREEZING)");

        return nThreads;
    }

#ifdef _OPENMP
    NamedParameter<unsigned int> getNCores()
    {
        static const unsigned int CORE_MOD = 10;
        unsigned int concurrentThreadsSupported = std::thread::hardware_concurrency();
        unsigned int threadsToUse = concurrentThreadsSupported / CORE_MOD;
        if (threadsToUse <= 0)
        {
            threadsToUse = 1;
        }

        NamedParameter<unsigned int> nCores = NamedParameter<unsigned int>(
            "nCores",
            threadsToUse,
            "Number of cores to use (OpenMP only)");

        INFO("Using OpenMP " << nCores.getVal() << " / " << concurrentThreadsSupported << " threads (DECREASE IN ConcurrentUtilities.cpp IF FREEZING)");

        return nCores;
    }
#endif

} // namespace AmpGen