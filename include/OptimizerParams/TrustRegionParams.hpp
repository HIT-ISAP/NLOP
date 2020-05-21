#ifndef TRUSTREGIONPARAMS_HPP
#define TRUSTREGIONPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::TrustRegionParams
/// @brief Abstract class for params for all trust region optimization
class TrustRegionParams
{
protected:
    virtual ~TrustRegionParams() {}
};
}

#endif
