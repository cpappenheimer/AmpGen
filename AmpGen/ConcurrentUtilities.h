#pragma once

#include "AmpGen/NamedParameter.h"

namespace AmpGen
{
    NamedParameter<size_t> getNumThreads();

#ifdef _OPENMP
    NamedParameter<unsigned int> getNCores();
#endif
} // namespace AmpGen
