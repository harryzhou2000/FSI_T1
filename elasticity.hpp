#pragma once
#include "defs.h"
#include "include/Eigen/Dense"

namespace SolidMaterial
{
    struct ElasSet
    {
        enum ElasInputType
        {
            E_nu = 0,
            E_G = 1,
            G_nu = 2
            //2-D conditions...
        };
        real E;         //Young's modulous
        real nu;        //poissonRatio
        real rho = 1.0; //mass density
        ElasSet(real I1 = 1, real I2 = 0.3, ElasInputType inputMode = E_nu)
        {
            switch (inputMode)
            {
            case E_G:
                E = I1;
                nu = 2 * I2 / E - 2;
                break;

            case G_nu:
                nu = I2;
                E = I1 * 2 * (1 + nu);
                break;

            default:
            case E_nu:
                E = I1;
                nu = I2;
                break;
            }
        }

        ElasSet(real I1, real I2, real Rho, ElasInputType inputMode = E_nu) : rho(Rho)
        {
            switch (inputMode)
            {
            case E_G:
                E = I1;
                nu = 2 * I2 / E - 2;
                break;

            case G_nu:
                nu = I2;
                E = I1 * 2 * (1 + nu);
                break;

            default:
            case E_nu:
                E = I1;
                nu = I2;
                break;
            }
        }

        Eigen::Matrix<real, 6, 6> constitutMatEFree() //constitutive Mat D/E, for normalized scale in linear problems
        {
            Eigen::Matrix<real, 6, 6> ret = Eigen::Matrix<real, 6, 6>::Constant(0.0);
            ret(0, 0) = ret(1, 1) = ret(2, 2) = 1.0;
            ret(3, 3) = ret(4, 4) = ret(5, 5) = (1.0 - 2.0 * nu) / ((1.0 - nu)); // be ware of gamma12 and e12
            ret(0, 1) = ret(0, 2) = ret(1, 2) = ret(1, 0) = ret(2, 0) = ret(2, 1) = nu / (1.0 - nu);
            ret *= (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
            return ret;
        }
    };

}
