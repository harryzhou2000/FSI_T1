#pragma once
#include "roesg_common.hpp"
#include <cmath>
#include <cstdlib>
#include <float.h>

namespace RiemannSolver
{
    using namespace MeshBasic;
    using ut5 = utype<5>;
    using ut3 = utype<3>;
    struct RoeSet
    {
        real gamma = 1.4;
        real mu = 1.0;
        real lambda = 0.0;
        real thresEntropyFix = 0.2;
    };

    static inline void U2W(const ut &U, const RoeSet &set, ut &W)
    {
        W.v[0] = sqrt(U.v[0]);
        real rsussum = 0;
        for (int i = 0; i < DIM; i++)
        {
            W.v[i + 1] = U.v[i + 1] / W.v[0];
            rsussum += sqr(U.v[i + 1]);
        }

        W.v[FDIM - 1] = (U.v[FDIM - 1] + (U.v[FDIM - 1] - 0.5 * rsussum / U.v[0]) * (set.gamma - 1)) / sqrt(U.v[0]);
    }

    static inline void EulerFlux(const ut &U, const RoeSet &set, ut &F)
    {
        real qs = 0;
        for (int i = 0; i < DIM; i++)
            qs += sqr(U[i + 1]);
        real p = (U[FDIM - 1] - 0.5 * qs / U[0]) * (set.gamma - 1);
        F[0] = U[1];
        F[1] = sqr(U[1]) / U[0] + p;
        for (int i = 1; i < DIM; i++)
            F[i + 1] = U[i + 1] * U[1] / U[0];
        F[FDIM - 1] = U[1] / U[0] * (p + U[FDIM - 1]);
    }

    //Roe's Reimann solver
    //vmax is inverse of time
    void RoeFlux(const ut &UL, const ut &UR, const RoeSet &set, ut &F, real &vmax, int &aim, real distance)
    {
        ut WL, WR;
        U2W(UL, set, WL);
        U2W(UR, set, WR);
        ut Waverage = (WL + WR) * 0.5;
        real vmid[DIM];
        for (int i = 0; i < DIM; i++)
            vmid[i] = Waverage[i + 1] / Waverage[0];
        real Hmid = Waverage[FDIM - 1] / Waverage[0];
        real qsmid = 0;
        for (int i = 0; i < DIM; i++)
            qsmid += sqr(vmid[i]);
        real asmid = (set.gamma - 1) * (Hmid - 0.5 * qsmid);
        if (asmid < __FLT_MIN__)
        {
            aim = 1;
            asmid = 0;
        }
        real LambdaA = (vmid[0] - sqrt(asmid));
        real LambdaB = (vmid[0] + sqrt(asmid));
        //Harten-Yee::
        real deltaEntropy = 0.01 * (sqrt(qsmid) + sqrt(asmid));
        if (std::abs(LambdaA) < deltaEntropy)
            LambdaA = (sqr(LambdaA) + sqr(deltaEntropy)) / deltaEntropy * 0.5 * sign(LambdaA);
        if (std::abs(LambdaB) < deltaEntropy)
            LambdaB = (sqr(LambdaB) + sqr(deltaEntropy)) / deltaEntropy * 0.5 * sign(LambdaB);
        //
        ut incU = UR - UL;
        real vincv = 0;
        for (int i = 0; i < DIM; i++)
            vincv += incU[i + 1] * vmid[i];

        real a2 = (set.gamma - 1) / asmid * ((Hmid - sqr(vmid[0])) * incU[0] + vincv - incU[FDIM - 1]);
        real a1 = 0.5 * (incU[0] - a2 - (incU[1] - vmid[0] * incU[0]) / sqrt(asmid));
        real aend = incU[0] - (a1 + a2);
        real a3p[DIM - 1];
        for (int i = 1; i < DIM; i++)
            a3p[i - 1] = incU[i + 1] - vmid[i] * incU[0];

        ut fleft, fright;
        EulerFlux(UL, set, fleft);
        EulerFlux(UR, set, fright);

        ///get eigen vectors
        ut e1, e2, eend, e3p[DIM - 1];
        e1[0] = 1;
        e1[1] = LambdaA;
        for (int i = 1; i < DIM; i++)
            e1[i + 1] = vmid[i];
        e1[FDIM - 1] = Hmid - vmid[0] * sqrt(asmid);
        eend[0] = 1;
        eend[1] = LambdaB;
        for (int i = 1; i < DIM; i++)
            eend[i + 1] = vmid[i];
        eend[FDIM - 1] = Hmid + vmid[0] * sqrt(asmid);
        e2[0] = 1;
        e2[1] = vmid[0];
        for (int i = 1; i < DIM; i++)
            e2[i + 1] = vmid[i];
        e2[FDIM - 1] = 0.5 * qsmid;

        for (int ii = 1; ii < DIM; ii++)
        {
            e3p[ii - 1][0] = 0;
            e3p[ii - 1][1] = 0;
            e3p[ii - 1][FDIM - 1] = vmid[ii];
            for (int i = 1; i < DIM; i++)
                if (i == ii)
                    e3p[ii - 1][i + 1] = 1;
                else
                    e3p[ii - 1][i + 1] = 0;
        }

        ///sumsum
        for (int i = 0; i < FDIM; i++)
            F[i] = 0;
        F += (fleft + fright);

        if (asmid != 0)
        {
            F += e1 * (-a1 * std::abs(LambdaA)) + eend * (-aend * std::abs(LambdaB)) + e2 * (-a2 * std::abs(vmid[0]));
            for (int i = 1; i < DIM; i++)
                F += e3p[i - 1] * (-a3p[i - 1] * std::abs(vmid[0])); // vmid 0 but not vmid i!!!
        }

        ///vmax
        vmax = vmid[0] / distance;
        vmax = std::max(vmax, std::abs(LambdaA) / distance);
        vmax = std::max(vmax, std::abs(LambdaB) / distance);

        F *= 0.5;
    }

    void RoeFlux2(const ut5 &UL, const ut5 &UR, const RoeSet &set, ut5 &F, real &vmax, int &aim, real distance)
    {

        real rsvL = (sqr(UL[1]) + sqr(UL[2]) + sqr(UL[3])) / (UL[0]);
        real rsvR = (sqr(UR[1]) + sqr(UR[2]) + sqr(UR[3])) / (UR[0]);
        real pL = (UL[4] - 0.5 * rsvL) * (set.gamma - 1);
        real pR = (UR[4] - 0.5 * rsvR) * (set.gamma - 1);
        real HL = (UL[4] + pL) / UL[0];
        real HR = (UR[4] + pR) / UR[0];
        ut3 Umid(0.0);
        Umid[0] = (UL[1] / sqrt(UL[0]) + UR[1] / sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        Umid[1] = (UL[2] / sqrt(UL[0]) + UR[2] / sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        Umid[2] = (UL[3] / sqrt(UL[0]) + UR[3] / sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real Hmid = (HL * sqrt(UL[0]) + HR * sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real qsmid = sqr(Umid[0]) + sqr(Umid[1]) + sqr(Umid[2]);
        real asmid = (set.gamma - 1) * (Hmid - 0.5 * qsmid);

        if (asmid < __FLT_MIN__)
        {
            aim = 1;
            asmid = 0;
        }

        //getk
        ut5 K[5];
        K[0][0] = 1, K[0][1] = Umid[0] - sqrt(asmid), K[0][2] = Umid[1], K[0][3] = Umid[2], K[0][4] = Hmid - Umid[0] * sqrt(asmid);
        K[1][0] = 1, K[1][1] = Umid[0], K[1][2] = Umid[1], K[1][3] = Umid[2], K[1][4] = 0.5 * qsmid;
        K[2][0] = 0, K[2][1] = 0, K[2][2] = 1, K[2][3] = 0, K[2][4] = Umid[1];
        K[3][0] = 0, K[3][1] = 0, K[3][2] = 0, K[3][3] = 1, K[3][4] = Umid[2];
        K[4][0] = 1, K[4][1] = Umid[0] + sqrt(asmid), K[4][2] = Umid[1], K[4][3] = Umid[2], K[4][4] = Hmid + Umid[0] * sqrt(asmid);

        //geta
        ut5 incU = UR - UL;
        real a[5];
        a[2] = incU[2] - Umid[1] * incU[0];
        a[3] = incU[3] - Umid[2] * incU[0];
        real inc4 = incU[4] - a[2] * Umid[1] - a[3] * Umid[2];
        a[1] = (set.gamma - 1) / asmid * (incU[0] * (Hmid - sqr(Umid[0])) + Umid[0] * incU[1] - inc4);
        a[0] = 0.5 * (incU[0] * (Umid[0] + sqrt(asmid)) - incU[1] - sqrt(asmid) * a[1]) / sqrt(asmid);
        a[4] = incU[0] - (a[0] + a[1]);

        //lambdas & harten-yee
        ut5 Lams;
        Lams[0] = Umid[0] - sqrt(asmid), Lams[1] = Lams[2] = Lams[3] = Umid[0], Lams[4] = Umid[0] + sqrt(asmid);
        real deltaEntropy = 0.05 * (sqrt(qsmid) + sqrt(asmid));
        for (int i = 0; i < 5; i++)
            if (std::abs(Lams[i]) < deltaEntropy)
                Lams[i] = (sqr(Lams[i]) + sqr(deltaEntropy)) / deltaEntropy * 0.5 * sign(Lams[i]);

        //
        ut5 FL, FR;
        EulerFlux(UL, set, FL), EulerFlux(UR, set, FR);
        F = FL + FR;
        if (asmid != 0)
        {
            for (int i = 0; i < 5; i++)
                F -= K[i] * (a[i] * std::abs(Lams[i]));
        }
        vmax = std::max(std::abs(Lams[0]), std::abs(Lams[4])) / distance;
        F *= 0.5;
    }

    void NewtonianViscosFlux(const ut &UL, const ut &UR, const RoeSet &set, ut &F, real &vmax, real distance)
    {
        real dv1dx1 = (UR[1] - UL[1]) / distance;
        real dvidx1[DIM - 1];
        F[0] = 0;
        F[1] = set.mu * 2 * dv1dx1;
        F[FDIM - 1] = F[1] * 0.5 * (UL[1] + UR[1]);
        for (int i = 1; i < DIM; i++)
        {
            dvidx1[i - 1] = (UR[i + 1] - UL[i + 1]) / distance;
            F[i + 1] = set.mu * dvidx1[i - 1];
            F[FDIM - 1] += F[i + 1] * 0.5 * (UR[i + 1] + UL[i + 1]);
        }
        vmax = max(vmax, set.mu * 2 / (UL[0] + UR[0]) / sqr(distance));
    }

    ut5 getPressureBound(const ut &UIN, const ut &UR, const RoeSet &set)
    {

        real rsvR = (sqr(UR[1]) + sqr(UR[2]) + sqr(UR[3])) / (UR[0]);
        real pR = (UR[4] - 0.5 * rsvR) * (set.gamma - 1);
        real aB = sqrt(set.gamma * pR / UR[0]);
        real u = (UR[1] / sqrt(UR[0]) + UIN[1] / sqrt(UIN[0])) / (sqrt(UIN[0]) + sqrt(UR[0]));
        if (u - aB > 0) //0 to inner
        {
            return UIN;
        }
        else if (u > 0) //1 to inner ,P
        {
            ut5 ret = UIN;
            ret[4] = pR / (set.gamma - 1) + 0.5 * (sqr(UIN[1]) + sqr(UIN[2]) + sqr(UIN[3])) / (UIN[0]);
            return ret;
        }
        else if (u + aB > 0) // 4 to inner ,P rho v w
        {
            ut5 ret = UR;
            ret[4] += -sqr(UR[1]) / (UR[0]) * 0.5 + sqr(UIN[1]) / (UR[0]) * 0.5;
            ret[1] = UIN[1];
            return ret;
        }
        else // all to inner
        {
            return UR;
        }
    }

    ut5 UPhysToUconserv(const ut5 &UPhys, const RoeSet &set) //[rho u v w p] -> [rho rhou rhov rhow E]
    {
        ut5 ret(UPhys);
        ret[1] = UPhys[0] * UPhys[1], ret[3] = UPhys[0] * UPhys[3], ret[2] = UPhys[0] * UPhys[2];
        ret[4] = UPhys[0] * (0.5 * (sqr(UPhys[1]) + sqr(UPhys[2]) + sqr(UPhys[3])) + UPhys[4] / (UPhys[0] * (set.gamma - 1)));
        return ret;
    }

    template <unsigned int dim>
    static inline real VecUSQR(const utype<dim + 2> &U)
    {
        real ret = 0.0;
        for (int i = 0; i < dim; i++)
            ret += sqr(U[i + 1]);
        return ret;
    }

    template <unsigned int dim>
    static inline void EulerFlux(const utype<dim + 2> &U, const RoeSet &set, utype<dim + 2> &F)
    {
        real qs = 0;
        for (int i = 0; i < dim; i++)
            qs += sqr(U[i + 1]);
        real p = (U[dim + 1] - 0.5 * qs / U[0]) * (set.gamma - 1);
        F[0] = U[1];
        F[1] = sqr(U[1]) / U[0] + p;
        for (int i = 1; i < dim; i++)
            F[i + 1] = U[i + 1] * U[1] / U[0];
        F[dim + 1] = U[1] / U[0] * (p + U[dim + 1]);
    }

    template <unsigned int dim>
    static inline void EulerFlux_Lin(const utype<dim + 2> &U, const utype<dim + 2> &Uval, const RoeSet &set, utype<dim + 2> &F)
    {
        EulerFlux<dim>(U, set, F);
        //F0
        F[0] += Uval[1] - U[1];
        real QQ = 0.0;
        for (int i = 0; i < dim; i++)
            QQ += sqr(U[i + 1]);

        //F1
        F[1] += (Uval[0] - U[0]) * (QQ * (set.gamma - 1) - sqr(U[1]) * 2) / (2 * sqr[U[0]]) +
                (Uval[1] - U[1]) * (3 - set.gamma) * U[1] / U[0] +
                (Uval[dim + 1] - U[dim + 1]) * (set.gamma - 1);
        for (int i = 1; i < dim; i++)
            F[1] += (Uval[i + 1] - U[i + 1]) * (-(set.gamma - 1) * U[i + 1] / U[0]);

        //F2 3
        for (int j = 1; j < dim; j++)
        {
            F[j + 1] += (Uval[0] - U[0]) * (-U[1] * U[j + 1] / sqr(U[0])) +
                        (Uval[1] - U[1]) * (U[j + 1] / U[0]);
            F[j + 1] += (Uval[j + 1] - U[j + 1]) * (U[1] / U[0]);
        }

        //Fend
        F[dim + 1] += (Uval[0] - U[0]) * (QQ * (set.gamma - 1) - set.gamma * U[0] * U[dim + 1]) * U[1] / cube(U[0]) +
                      (Uval[1] - U[1]) * ((U[dim + 1] - (set.gamma - 1) * (QQ - 2 * U[0] * U[dim + 1]) / (2 * U[0])) / U[0] -
                                          sqr(U[1]) * (set.gamma - 1) / sqr(U[0])) +
                      (Uval[dim + 1] - U[dim + 1]) * set.gamma * U[1] / U[0];
        for (int i = 0; i < dim; i++)
            F[dim + 1] += (Uval[i + 1] - U[i + 1]) * (-(set.gamma - 1) * U[1] * U[i + 1] / sqr(U[0]));
    }

    template <unsigned int dim>
    void RoeFluxT(utype<dim + 2> &UL, utype<dim + 2> &UR, const RoeSet &set, utype<dim + 2> &F, real &vmax, int &aim, real distance)
    {
        bool retcent = false;
        if (UL[0] <= 0 || UR[0] <= 0)
        {
            std::cerr << "===RoeFlux===: RHO not Positive!!!" << std::endl;
            if (UL[0] <= 0)
            {
                UL[0] = 1e-10;
                for (int i = 0; i < dim; i++)
                    UL[i + 1] = 1e-20;
            }
            if (UR[0] <= 0)
            {
                UR[0] = 1e-10;
                for (int i = 0; i < dim; i++)
                    UR[i + 1] = 1e-20;
            }
            retcent = true;
        }
        real rsvL = VecUSQR<dim>(UL) / (UL[0]);
        real rsvR = VecUSQR<dim>(UR) / (UR[0]);
        real pL = (UL[dim + 1] - 0.5 * rsvL) * (set.gamma - 1);
        real pR = (UR[dim + 1] - 0.5 * rsvR) * (set.gamma - 1);
        real HL = (UL[dim + 1] + pL) / UL[0];
        real HR = (UR[dim + 1] + pR) / UR[0];
        utype<dim> Umid(0.0);
        for (int i = 0; i < dim; i++)
            Umid[i] = (UL[i + 1] / sqrt(UL[0]) + UR[i + 1] / sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real Hmid = (HL * sqrt(UL[0]) + HR * sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real qsmid = Umid.dot(Umid);
        real asmid = (set.gamma - 1) * (Hmid - 0.5 * qsmid);
        if (asmid < __FLT_MIN__)
        {
            aim = 1;
            asmid = std::abs(asmid);
        }

        //get eigen vectors
        utype<dim + 2> K[dim + 2];
        //k0
        K[0][0] = 1;
        K[0][1] = Umid[0] - sqrt(asmid);
        for (int i = 1; i < dim; i++)
            K[0][i + 1] = Umid[i];
        K[0][dim + 1] = Hmid - Umid[0] * sqrt(asmid);
        //k1
        K[1][0] = 1;
        for (int i = 0; i < dim; i++)
            K[1][i + 1] = Umid[i];
        K[1][dim + 1] = 0.5 * qsmid;
        //k2~kd
        for (int j = 1; j < dim; j++)
        {
            K[j + 1][0] = K[j + 1][1] = 0;
            for (int i = 1; i < dim; i++)
                K[j + 1][i + 1] = j == i ? 1 : 0;
            K[j + 1][dim + 1] = Umid[j];
        }
        //kd+1
        K[dim + 1][0] = 1;
        K[dim + 1][1] = Umid[0] + sqrt(asmid);
        for (int i = 1; i < dim; i++)
            K[dim + 1][i + 1] = Umid[i];
        K[dim + 1][dim + 1] = Hmid + Umid[0] * sqrt(asmid);

        //geta
        utype<dim + 2> incU = UR - UL;
        real a[dim + 2];
        for (int i = 1; i < dim; i++)
            a[i + 1] = incU[i + 1] - Umid[i] * incU[0];
        real inc4 = incU[dim + 1];
        for (int i = 1; i < dim; i++)
            inc4 -= a[i + 1] * Umid[i];

        a[1] = (set.gamma - 1) / asmid * (incU[0] * (Hmid - sqr(Umid[0])) + Umid[0] * incU[1] - inc4);
        a[0] = 0.5 * (incU[0] * (Umid[0] + sqrt(asmid)) - incU[1] - sqrt(asmid) * a[1]) / sqrt(asmid);
        a[dim + 1] = incU[0] - (a[0] + a[1]);

        //lambdas & harten-yee
        utype<dim + 2> Lams;
        Lams[0] = Umid[0] - sqrt(asmid);
        for (int i = 0; i < dim; i++)
            Lams[i + 1] = Umid[0];
        Lams[dim + 1] = Umid[0] + sqrt(asmid);
        real deltaEntropy = set.thresEntropyFix * (sqrt(qsmid) + sqrt(asmid));
        for (int i = 0; i < dim + 2; i++)
            if (std::abs(Lams[i]) < deltaEntropy)
                Lams[i] = (sqr(Lams[i]) + sqr(deltaEntropy)) / deltaEntropy * 0.5 * sign(Lams[i]);

        // sums
        utype<dim + 2> FL, FR;
        EulerFlux<dim>(UL, set, FL), EulerFlux<dim>(UR, set, FR);
        F = FL + FR;
        vmax = Lams[1] / distance;
        if (asmid != 0 && (!retcent))
        {
            for (int i = 0; i < dim + 2; i++)
                F -= K[i] * (a[i] * std::abs(Lams[i]));
            vmax = std::max(std::abs(Lams[0]), std::abs(Lams[dim + 1])) / distance;
        }

        F *= 0.5;
        for (int i = 0; i < dim + 2; i++)
            if (std::isnan(F[i]) || std::isinf(F[i]))
            {
                std::cerr << "===RoeFlux===: Return Value " << i << " Inf or Nan!!!" << std::endl;
                exit(1);
                return;
            }
    }

    // this RoeFlux is identical with RoeFluxT, but the convection uses UL,UR, convected uses ULval, URval,
    // which is linear mapping
    // for ULval-UL and URval-UR
    template <unsigned int dim>
    void RoeFluxT_Lin(utype<dim + 2> &UL, utype<dim + 2> &UR,
                      const utype<dim + 2> &ULval, const utype<dim + 2> &URval,
                      const RoeSet &set, utype<dim + 2> &F, real &vmax, int &aim, real distance)
    {
        bool retcent = false;
        if (UL[0] <= 0 || UR[0] <= 0)
        {
            std::cerr << "===RoeFluxLin===: RHO not Positive!!!" << std::endl;
            if (UL[0] <= 0)
            {
                UL[0] = 1e-10;
                for (int i = 0; i < dim; i++)
                    UL[i + 1] = 1e-20;
            }
            if (UR[0] <= 0)
            {
                UR[0] = 1e-10;
                for (int i = 0; i < dim; i++)
                    UR[i + 1] = 1e-20;
            }
            retcent = true;
        }
        real rsvL = VecUSQR<dim>(UL) / (UL[0]);
        real rsvR = VecUSQR<dim>(UR) / (UR[0]);
        real pL = (UL[dim + 1] - 0.5 * rsvL) * (set.gamma - 1);
        real pR = (UR[dim + 1] - 0.5 * rsvR) * (set.gamma - 1);
        real HL = (UL[dim + 1] + pL) / UL[0];
        real HR = (UR[dim + 1] + pR) / UR[0];
        utype<dim> Umid(0.0);
        for (int i = 0; i < dim; i++)
            Umid[i] = (UL[i + 1] / sqrt(UL[0]) + UR[i + 1] / sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real Hmid = (HL * sqrt(UL[0]) + HR * sqrt(UR[0])) / (sqrt(UL[0]) + sqrt(UR[0]));
        real qsmid = Umid.dot(Umid);
        real asmid = (set.gamma - 1) * (Hmid - 0.5 * qsmid);
        if (asmid < __FLT_MIN__)
        {
            aim = 1;
            asmid = 0;
        }

        //get eigen vectors
        utype<dim + 2> K[dim + 2];
        //k0
        K[0][0] = 1;
        K[0][1] = Umid[0] - sqrt(asmid);
        for (int i = 1; i < dim; i++)
            K[0][i + 1] = Umid[i];
        K[0][dim + 1] = Hmid - Umid[0] * sqrt(asmid);
        //k1
        K[1][0] = 1;
        for (int i = 0; i < dim; i++)
            K[1][i + 1] = Umid[i];
        K[1][dim + 1] = 0.5 * qsmid;
        //k2~kd
        for (int j = 1; j < dim; j++)
        {
            K[j + 1][0] = K[j + 1][1] = 0;
            for (int i = 1; i < dim; i++)
                K[j + 1][i + 1] = j == i ? 1 : 0;
            K[j + 1][dim + 1] = Umid[j];
        }
        //kd+1
        K[dim + 1][0] = 1;
        K[dim + 1][1] = Umid[0] + sqrt(asmid);
        for (int i = 1; i < dim; i++)
            K[dim + 1][i + 1] = Umid[i];
        K[dim + 1][dim + 1] = Hmid + Umid[0] * sqrt(asmid);

        //geta
        utype<dim + 2> incU = URval - ULval;
        real a[dim + 2];
        for (int i = 1; i < dim; i++)
            a[i + 1] = incU[i + 1] - Umid[i] * incU[0];
        real inc4 = incU[dim + 1];
        for (int i = 1; i < dim; i++)
            inc4 -= a[i + 1] * Umid[i];

        a[1] = (set.gamma - 1) / asmid * (incU[0] * (Hmid - sqr(Umid[0])) + Umid[0] * incU[1] - inc4);
        a[0] = 0.5 * (incU[0] * (Umid[0] + sqrt(asmid)) - incU[1] - sqrt(asmid) * a[1]) / sqrt(asmid);
        a[dim + 1] = incU[0] - (a[0] + a[1]);

        //lambdas & harten-yee
        utype<dim + 2> Lams;
        Lams[0] = Umid[0] - sqrt(asmid);
        for (int i = 0; i < dim; i++)
            Lams[i + 1] = Umid[0];
        Lams[dim + 1] = Umid[0] + sqrt(asmid);
        real deltaEntropy = set.thresEntropyFix * (sqrt(qsmid) + sqrt(asmid));
        for (int i = 0; i < dim + 2; i++)
            if (std::abs(Lams[i]) < deltaEntropy)
                Lams[i] = (sqr(Lams[i]) + sqr(deltaEntropy)) / deltaEntropy * 0.5 * sign(Lams[i]);

        // sums
        utype<dim + 2> FL, FR;
        EulerFlux_Lin<dim>(UL, ULval, set, FL), EulerFlux_Lin<dim>(UR, URval, set, FR);
        F = FL + FR;
        vmax = Lams[1] / distance;
        if (asmid != 0 && (!retcent))
        {
            for (int i = 0; i < dim + 2; i++)
                F -= K[i] * (a[i] * std::abs(Lams[i]));
            vmax = std::max(std::abs(Lams[0]), std::abs(Lams[dim + 1])) / distance;
        }

        F *= 0.5;
        for (int i = 0; i < dim + 2; i++)
            if (std::isnan(F[i]) || std::isinf(F[i]))
            {
                std::cerr << "===RoeFlux===: Return Value " << i << " Inf or Nan!!!" << std::endl;
                exit(1);
                return;
            }
    }
}