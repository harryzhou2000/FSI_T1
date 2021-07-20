#include "roe.hpp"
using namespace RiemannSolver;

int main()

{
    real ULd[4] = {1, -1, 1, 4};
    real URd[4] = {2, -2, 1, 15};
    utype<4> UL;
    utype<4> UR;
    utype<4> F;
    UL.set(ULd);
    UR.set(URd);
    real vmax;
    int aim;
    RoeSet set;
    set.gamma = 1.4;

    real ULd3[5] = {1, -1, 1, 0, 4};
    real URd3[5] = {2, -2, 1, 0, 15};
    ut UL3;
    ut UR3;
    ut F3,FF3;
    UL3.set(ULd3);
    UR3.set(URd3);

    RoeFluxT<2>(UL, UR, set, F, vmax, aim, 1.0);
    RoeFluxT<3>(UL3, UR3, set, F3, vmax, aim, 1.0);
    RoeFlux2(UL3, UR3, set, FF3, vmax, aim, 1.0);

    auto ret1 = UL * TVDFunc<2>(UL,UR);
    auto ret2 = -UR * TVDFunc<2>(-UR, -UL);

    return 0;
}
