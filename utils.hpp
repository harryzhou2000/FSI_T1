#pragma once

#include "utype_basics.hpp"

namespace UTILS
{

    template <size dim>
    real RBFInversedMultiQuadric(UTBasic::utype<dim> x, real alpha)
    {
        real r = x.norm2();
        //std::cout << r << std::endl;
        return 1.0 / std::sqrt(UTBasic::sqr(r) + UTBasic::sqr(alpha) + 1e-36);
    }
}