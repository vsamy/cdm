#pragma once

#include "cdm/typedefs.hpp"

namespace cdm {

template <int Order>
LB66 DualMul(const coma::CMTM<Order>& lhs, const LB66& rhs)
{
    assert(lhs.nMat() == rhs.nBlocks());
    LB66 out(rhs.rows(), rhs.cols());
    
}

} // namespace cdm
