#pragma once
#include "roe.hpp"
#include <vector>
#include <iostream>
#include <omp.h>
#include <string>
#include <string.h>
////////////must be 2d
#include <utility>

namespace StructuralGrid
{
#define vector std::vector
    using namespace RiemannSolver;

    struct SGGrid
    {
        vector<utype<4>> d;
        vector<utype<4>> dtemp;
        vector<utype<4>> recval[4];
        vector<gridNeighbourSG<4>> N;
        vector<TDNodes<4, 4, 1>> pnodes;
        vector<TDNodesE<4>> pnodese;

        // boundary
        vector<SGGrid *> crossTarget; // if using outher SGGrids
        vector<size> crossLoc;        //node index of crossto location
        vector<size> crossEdge;       //edge index of crossto location
        vector<utype<4>> infinityData;
        size ld, mm, nn;

        // additional
        bool viscous, moving;
        size nwall;
        vector<real> wallload;

        vector<real> vmaxs;
        vector<real> vmaxsC;
        vector<int> aims;
        vector<utype<4>> dU0, dU1, dU2, dU3, U1, U2, U3;

        inline size ind(size i, size j) // suppose: i is direction of eta, j is xi
        {
            return i + j * ld;
        }

        //input by col-leading Whole boundary
        // x and y are ___ m+1 by n+1 ___
        /* 
           |   2   |
    i(eta) |3     1|
           |   0   |
               j(xi)
    */
        SGGrid(bool nviscous = false, bool nmoving = false) : viscous(nviscous), moving(nmoving)
        {
            ld = 0, mm = 0, nn = 0;
        }

        void build(size m, size n, vector<real> x, vector<real> y,
                   vector<ntype> bound, vector<SGGrid *> ncrossTarget, vector<size> ncrossLoc, vector<size> ncrossEdge,
                   vector<utype<4>> infdataN)
        {
            ld = mm = m, nn = n;
            size gridsz = m * n;
            size gridPeri = m * 2 + n * 2;
            infinityData = infdataN;
            d.resize(gridsz), dtemp.resize(gridsz), N.resize(gridsz), pnodes.resize(gridsz);
            recval[0].resize(gridsz), recval[1].resize(gridsz), recval[2].resize(gridsz), recval[3].resize(gridsz);
            pnodese.resize(gridsz);

            crossTarget = ncrossTarget;
            crossLoc = ncrossLoc;
            crossEdge = ncrossEdge;
            nwall = 0;
            for (size j = 0; j < n; j++)
            {
                for (size i = 0; i < m; i++)
                {
                    /// points
                    pnodes[ind(i, j)].d[0][0] = x[i + 0 + (m + 1) * (j + 0)], pnodes[ind(i, j)].d[0][1] = y[i + 0 + (m + 1) * (j + 0)];
                    pnodes[ind(i, j)].d[1][0] = x[i + 0 + (m + 1) * (j + 1)], pnodes[ind(i, j)].d[1][1] = y[i + 0 + (m + 1) * (j + 1)];
                    pnodes[ind(i, j)].d[2][0] = x[i + 1 + (m + 1) * (j + 1)], pnodes[ind(i, j)].d[2][1] = y[i + 1 + (m + 1) * (j + 1)];
                    pnodes[ind(i, j)].d[3][0] = x[i + 1 + (m + 1) * (j + 0)], pnodes[ind(i, j)].d[3][1] = y[i + 1 + (m + 1) * (j + 0)];
                    pnodes[ind(i, j)].getAreaBary();
                    pnodes[ind(i, j)].getNW();

                    /// neighbour & boundary
                    for (int iter = 0; iter < 4; iter++)
                        N[ind(i, j)]
                            .t[iter] = inner;
                    N[ind(i, j)].n[0] = ind(i - 1, j);
                    N[ind(i, j)].n[1] = ind(i, j + 1);
                    N[ind(i, j)].n[2] = ind(i + 1, j);
                    N[ind(i, j)].n[3] = ind(i, j - 1);

                    if (bound.size() == ncrossLoc.size())
                    {
                        if (i == 0)
                        {
                            N[ind(i, j)].t[0] = bound[j];
                            N[ind(i, j)].n[0] = j;
                            if (N[ind(i, j)].t[0] == wall)
                                nwall++;
                        }
                        if (j == n - 1)
                        {
                            N[ind(i, j)].t[1] = bound[n + i];
                            N[ind(i, j)].n[1] = n + i;
                            if (N[ind(i, j)].t[1] == wall)
                                nwall++;
                        }
                        if (i == m - 1)
                        {
                            N[ind(i, j)].t[2] = bound[m + n + j];
                            N[ind(i, j)].n[2] = m + n + j;
                            if (N[ind(i, j)].t[2] == wall)
                                nwall++;
                        }
                        if (j == 0)
                        {
                            N[ind(i, j)].t[3] = bound[m + n + n + i];
                            N[ind(i, j)].n[3] = m + n + n + i;
                            if (N[ind(i, j)].t[3] == wall)
                                nwall++;
                        }
                    }
                    else //size is 4
                    {
                        if (i == 0)
                        {
                            N[ind(i, j)].t[0] = bound[0];
                            N[ind(i, j)].n[0] = j;
                            if (bound[0] == wall)
                                nwall++;
                        }
                        if (j == n - 1)
                        {
                            N[ind(i, j)].t[1] = bound[1];
                            N[ind(i, j)].n[1] = n + i;
                            if (bound[1] == wall)
                                nwall++;
                        }
                        if (i == m - 1)
                        {
                            N[ind(i, j)].t[2] = bound[2];
                            N[ind(i, j)].n[2] = m + n + j;
                            if (bound[2] == wall)
                                nwall++;
                        }
                        if (j == 0)
                        {
                            N[ind(i, j)].t[3] = bound[3];
                            N[ind(i, j)].n[3] = m + n + n + i;
                            if (bound[3] == wall)
                                nwall++;
                        }
                    }
                }
            }
            wallload.resize(nwall * 2, 0);
            for (auto &x : crossTarget)
                if (!x)
                    x = this;
        }

        SGGrid(size m, size n, vector<real> x, vector<real> y,
               vector<ntype> bound, vector<SGGrid *> ncrossTarget, vector<size> ncrossLoc, vector<size> ncrossEdge,
               vector<utype<4>> infdataN)
        {
            viscous = moving = false;
            build(m, n, x, y, bound, ncrossTarget, ncrossLoc, ncrossEdge, infdataN);
        }

        void uploadGrid(vector<real> &x, vector<real> &y, vector<real> &dx, vector<real> &dy)
        {
            moving = true;
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size j = 0; j < nn; j++)
            {
                for (size i = 0; i < mm; i++)
                {
                    /// points
                    pnodes[ind(i, j)].d[0][0] = x[i + 0 + (mm + 1) * (j + 0)], pnodes[ind(i, j)].d[0][1] = y[i + 0 + (mm + 1) * (j + 0)];
                    pnodes[ind(i, j)].d[1][0] = x[i + 0 + (mm + 1) * (j + 1)], pnodes[ind(i, j)].d[1][1] = y[i + 0 + (mm + 1) * (j + 1)];
                    pnodes[ind(i, j)].d[2][0] = x[i + 1 + (mm + 1) * (j + 1)], pnodes[ind(i, j)].d[2][1] = y[i + 1 + (mm + 1) * (j + 1)];
                    pnodes[ind(i, j)].d[3][0] = x[i + 1 + (mm + 1) * (j + 0)], pnodes[ind(i, j)].d[3][1] = y[i + 1 + (mm + 1) * (j + 0)];
                    pnodes[ind(i, j)].getAreaBary();
                    pnodes[ind(i, j)].getNW();
                    pnodese[ind(i, j)].ddt[0][0] = dx[i + 0 + (mm + 1) * (j + 0)], pnodese[ind(i, j)].ddt[0][1] = dy[i + 0 + (mm + 1) * (j + 0)];
                    pnodese[ind(i, j)].ddt[1][0] = dx[i + 0 + (mm + 1) * (j + 1)], pnodese[ind(i, j)].ddt[1][1] = dy[i + 0 + (mm + 1) * (j + 1)];
                    pnodese[ind(i, j)].ddt[2][0] = dx[i + 1 + (mm + 1) * (j + 1)], pnodese[ind(i, j)].ddt[2][1] = dy[i + 1 + (mm + 1) * (j + 1)];
                    pnodese[ind(i, j)].ddt[3][0] = dx[i + 1 + (mm + 1) * (j + 0)], pnodese[ind(i, j)].ddt[3][1] = dy[i + 1 + (mm + 1) * (j + 0)];
                }
            }
        }

        void buildInitByFile(int m, int n, std::istream &meshin, std::istream &datain, vector<ntype> bound, vector<SGGrid *> ncrossTarget, vector<size> ncrossLoc, vector<size> ncrossEdge,
                             vector<utype<4>> infdataN)
        {
            vector<real> x((n + 1) * (m + 1));
            vector<real> y((n + 1) * (m + 1));
            for (size_t j = 0; j <= n; j++)
                for (size_t i = 0; i <= m; i++)
                {
                    meshin >> x[(m + 1) * j + i];
                    meshin >> y[(m + 1) * j + i];
                }
            real dump10;
            d.resize(m * n);
            char dbuf[2048];
            for (size i = 0; i < m * n; i++)
            {
                for (int dum = 0; dum < 10; dum++)
                    datain >> dump10;
                for (int read = 0; read < 4; read++)
                    datain >> d[i][read];
                datain.getline(dbuf, 2048);
            }
            build(m, n, x, y, bound, ncrossTarget, ncrossLoc, ncrossEdge, infdataN);
        }

        void buildByFile(int m, int n, std::istream &meshin, vector<ntype> bound, vector<SGGrid *> ncrossTarget, vector<size> ncrossLoc, vector<size> ncrossEdge,
                         vector<utype<4>> infdataN)
        {
            if (!meshin)
                std::cerr << "File Not Found"
                          << std::endl;
            vector<real> x((n + 1) * (m + 1));
            vector<real> y((n + 1) * (m + 1));
            for (size_t j = 0; j <= n; j++)
                for (size_t i = 0; i <= m; i++)
                {
                    meshin >> x[(m + 1) * j + i];
                    meshin >> y[(m + 1) * j + i];
                }
            build(m, n, x, y, bound, ncrossTarget, ncrossLoc, ncrossEdge, infdataN);
        }

        void SetDistance() //must be operated after all other grids are settled
        {
#ifdef OMP_ON
#pragma omp parallel for schedule(guided)
#endif
            for (size i = 0; i < d.size(); i++)
            {
                real xb = pnodes[i].bary[0], yb = pnodes[i].bary[1];
                real dis2Wall, dis2Wall1, xb1, yb1;
                for (int ii = 0; ii < 4; ii++)
                    switch (N[i].t[ii])
                    {
                    case inner:
                        xb1 = pnodes[N[i].n[ii]].bary[0], yb1 = pnodes[N[i].n[ii]].bary[1];
                        dis2Wall = pointToLine(xb, yb, pnodes[i].d[ii][0], pnodes[i].d[ii][1],
                                               pnodes[i].d[(ii + 1) % 4][0], pnodes[i].d[(ii + 1) % 4][1]);
                        dis2Wall1 = pointToLine(xb1, yb1, pnodes[N[i].n[ii]].d[ii][0], pnodes[N[i].n[ii]].d[ii][1],
                                                pnodes[N[i].n[ii]].d[(ii + 1) % 4][0], pnodes[N[i].n[ii]].d[(ii + 1) % 4][1]);
                        N[i].d[ii] = dis2Wall1 + dis2Wall;
                        N[i].dn[ii] = dis2(xb1, yb1, xb, yb);
                        N[i].dp[ii] = dis2Wall / N[i].d[ii];
                        break;
                    case infinity:
                        dis2Wall = pointToLine(xb, yb, pnodes[i].d[ii][0], pnodes[i].d[ii][1],
                                               pnodes[i].d[(ii + 1) % 4][0], pnodes[i].d[(ii + 1) % 4][1]);
                        N[i].dn[ii] = N[i].d[ii] = 2 * dis2Wall;
                        N[i].dp[ii] = 1. / 2;
                        break;
                    case wall:
                        dis2Wall = pointToLine(xb, yb, pnodes[i].d[ii][0], pnodes[i].d[ii][1],
                                               pnodes[i].d[(ii + 1) % 4][0], pnodes[i].d[(ii + 1) % 4][1]);
                        N[i].dn[ii] = N[i].d[ii] = 2 * dis2Wall;
                        N[i].dp[ii] = 1. / 2;
                        break;
                    case cross:
                        xb1 = crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].bary[0];
                        yb1 = crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].bary[1];
                        dis2Wall = pointToLine(xb, yb, pnodes[i].d[ii][0], pnodes[i].d[ii][1],
                                               pnodes[i].d[(ii + 1) % 4][0], pnodes[i].d[(ii + 1) % 4][1]);
                        dis2Wall1 = pointToLine(xb1, yb1,
                                                crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].d[crossEdge[N[i].n[ii]]][0], crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].d[crossEdge[N[i].n[ii]]][1],
                                                crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].d[(crossEdge[N[i].n[ii]] + 1) % 4][0], crossTarget[N[i].n[ii]]->pnodes[crossLoc[N[i].n[ii]]].d[(crossEdge[N[i].n[ii]] + 1) % 4][1]);
                        N[i].d[ii] = dis2(xb1, yb1, xb, yb);
                        N[i].dn[ii] = dis2Wall1 + dis2Wall;
                        N[i].dp[ii] = dis2Wall / N[i].d[ii];
                        break;
                    }
            }
        }

        void SetInitial(utype<4> Uuniform)
        {
            for (auto &x : d)
                x = Uuniform;
        }

        void Reconstruct(double tvdOverride) //tvd 2nd order reconstruct
        {
#ifdef OMP_ON
#pragma omp parallel for schedule(guided)
#endif
            for (size i = 0; i < d.size(); i++)
            {
                utype<4> deltaU[4];
                //utype<4> grad[2];
                utype<4> Uim;
                utype<4> Umax(d[i]), Umin(d[i]);
                for (size ii = 0; ii < 4; ii++)
                {
                    switch (N[i].t[ii])
                    {
                    case inner:
                        deltaU[ii] = (d[N[i].n[ii]] - d[i]);
                        Umax = Umax.max(d[N[i].n[ii]]);
                        Umin = Umin.min(d[N[i].n[ii]]);
                        break;

                    case infinity:
                        deltaU[ii] = (infinityData[N[i].n[ii]] - d[i]);
                        Umax = Umax.max(infinityData[N[i].n[ii]]);
                        Umin = Umin.min(infinityData[N[i].n[ii]]);
                        break;

                    case wall:
                        Uim = pnodes[i].getUOut(d[i], ii);
                        Uim[1] = -Uim[1];
                        if (viscous)
                            Uim[2] = -Uim[2];
                        Uim = pnodes[i].getUCent(Uim, ii);
                        deltaU[ii] = (Uim - d[i]);
                        Umax = Umax.max(Uim);
                        Umin = Umin.min(Uim);
                        break;

                    case cross:
                        deltaU[ii] = (crossTarget[N[i].n[ii]]->d[crossLoc[N[i].n[ii]]] - d[i]);
                        Umax = Umax.max(crossTarget[N[i].n[ii]]->d[crossLoc[N[i].n[ii]]]);
                        Umin = Umin.min(crossTarget[N[i].n[ii]]->d[crossLoc[N[i].n[ii]]]);
                        break;

                    default:
                        exit(-12);
                    }
                }
                //grad[0] = (deltaU[0] * TVDFunc<2>(deltaU[0], -deltaU[2]))
                recval[0][i] = (deltaU[0] * TVDFunc<2>(deltaU[0], -deltaU[2])) * (tvdOverride * N[i].dp[0]);
                recval[1][i] = (deltaU[1] * TVDFunc<2>(deltaU[1], -deltaU[3])) * (tvdOverride * N[i].dp[1]);
                recval[2][i] = (deltaU[2] * TVDFunc<2>(deltaU[2], -deltaU[0])) * (tvdOverride * N[i].dp[2]);
                recval[3][i] = (deltaU[3] * TVDFunc<2>(deltaU[3], -deltaU[1])) * (tvdOverride * N[i].dp[3]);
                auto recvalmaxD = utype<4>(0.0).max(recval[0][i]).max(recval[1][i]).max(recval[2][i]).max(recval[3][i]);
                auto recvalminD = utype<4>(0.0).min(recval[0][i]).min(recval[1][i]).min(recval[2][i]).min(recval[3][i]);
                auto UmaxD = Umax - d[i], UminD = Umin - d[i];
                UmaxD = (UmaxD / (recvalmaxD + 1e-3)).abs();
                UminD = (UminD / (recvalminD - 1e-3)).abs().min(UmaxD);
                real cc = std::min(1.0, UminD.minelem());
                recval[0][i] = (deltaU[0] * TVDFunc<2>(deltaU[0], -deltaU[2])) * (tvdOverride * cc * N[i].dp[0]) + d[i];
                recval[1][i] = (deltaU[1] * TVDFunc<2>(deltaU[1], -deltaU[3])) * (tvdOverride * cc * N[i].dp[1]) + d[i];
                recval[2][i] = (deltaU[2] * TVDFunc<2>(deltaU[2], -deltaU[0])) * (tvdOverride * cc * N[i].dp[2]) + d[i];
                recval[3][i] = (deltaU[3] * TVDFunc<2>(deltaU[3], -deltaU[1])) * (tvdOverride * cc * N[i].dp[3]) + d[i];
                // std::cout << cc << std::endl;
                for (int ii = 0; ii < 4; ii++)
                    recval[ii][i] = recval[ii][i].max(Umin).min(Umax);
            }
        }

        void GetIncrement(vector<utype<4>> &increment, double tvdOverride, const RoeSet &set, real &vmax, int &aim) // for explicit time stepping
        {
            increment.resize(d.size());
            Reconstruct(tvdOverride);
            vmaxs.assign(d.size(), 0.);
            aims.assign(d.size(), 0);
#ifdef OMP_ON
#pragma omp parallel for schedule(guided)
#endif
            for (size i = 0; i < d.size(); i++)
            {
                utype<4> fluxes[4];
                utype<4> vf;
                utype<4> Uthis[4];
                utype<4> Uthat[4];
                utype<4> Uwall[4];
                real vmaxa = 0;
                int aima = 0;
                size wallcount = 0;
                for (size ii = 0; ii < 4; ii++)
                {
                    switch (N[i].t[ii])
                    {
                    case inner:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(recval[(ii + 2) % 4][N[i].n[ii]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    case infinity:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(infinityData[N[i].n[ii]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    case wall:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii][1] = -Uthat[ii][1];
                        if (viscous)
                            Uthat[ii][2] = -Uthat[ii][2];
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        wallload[wallcount * 2 + 0] = fluxes[ii][1];
                        wallload[wallcount * 2 + 1] = fluxes[ii][2];
                        wallcount++;
                        break;

                    case cross:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(crossTarget[N[i].n[ii]]->recval[crossEdge[N[i].n[ii]]][crossLoc[N[i].n[ii]]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    default:
                        break;
                    }
                }
                // if (viscous)
                //     for (int ii = 0; ii < 4; ii++)
                //     {
                //         NewtonianViscosFlux(Uthis[ii], Uthat[ii], set, vf, vmax, N[i].d[ii]);
                //         fluxes[ii] -= pnodes[i].getUCent(vf, ii) * pnodes[i].w[ii];
                //     }
                increment[i] = -(fluxes[0] + fluxes[1] + fluxes[2] + fluxes[3]) * (1. / pnodes[i].area);
                vmaxs[i] = vmaxa;
                aims[i] = aima;
                // if (aims[i])
                //     aims[i] = aima;
            }
            for (size i = 0; i < d.size(); i++)
            {
                if (vmaxs[i] > vmax)
                    vmax = vmaxs[i];
                if (aims[i])
                    aim = 1;
            }
        }

        void GetIncrement_Lin(vector<utype<4>> &Uval, vector<utype<4>> &increment, double tvdOverride, const RoeSet &set, real &vmax, int &aim) // for explicit time stepping
        {
            increment.resize(d.size()); //TODO
            Reconstruct(tvdOverride);
            vmaxs.assign(d.size(), 0.);
            aims.assign(d.size(), 0);
#ifdef OMP_ON
#pragma omp parallel for schedule(guided)
#endif
            for (size i = 0; i < d.size(); i++)
            {
                utype<4> fluxes[4];
                utype<4> vf;
                utype<4> Uthis[4];
                utype<4> Uthat[4];
                utype<4> Uwall[4];
                real vmaxa = 0;
                int aima = 0;
                size wallcount = 0;
                for (size ii = 0; ii < 4; ii++)
                {
                    switch (N[i].t[ii])
                    {
                    case inner:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(recval[(ii + 2) % 4][N[i].n[ii]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    case infinity:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(infinityData[N[i].n[ii]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    case wall:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii][1] = -Uthat[ii][1];
                        if (viscous)
                            Uthat[ii][2] = -Uthat[ii][2];
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        wallload[wallcount * 2 + 0] = fluxes[ii][1];
                        wallload[wallcount * 2 + 1] = fluxes[ii][2];
                        wallcount++;
                        break;

                    case cross:
                        Uthis[ii] = pnodes[i].getUOut(recval[ii][i], ii);
                        Uthat[ii] = pnodes[i].getUOut(crossTarget[N[i].n[ii]]->recval[crossEdge[N[i].n[ii]]][crossLoc[N[i].n[ii]]], ii);
                        if (moving)
                        {
                            Uwall[ii][0] = Uwall[ii][3] = 0;
                            Uwall[ii][1] = 0.5 * (pnodese[i].ddt[ii][0] + pnodese[i].ddt[(ii + 1) % 4][0]);
                            Uwall[ii][2] = 0.5 * (pnodese[i].ddt[ii][1] + pnodese[i].ddt[(ii + 1) % 4][1]);
                            Uwall[ii] = pnodes[i].getUOut(Uwall[ii], ii);
                            Uthis[ii][1] -= Uwall[ii][1] * Uthis[ii][0];
                            Uthat[ii][1] -= Uwall[ii][1] * Uthat[ii][0];
                            Uthis[ii][2] -= Uwall[ii][2] * Uthis[ii][0];
                            Uthat[ii][2] -= Uwall[ii][2] * Uthat[ii][0];
                        }
                        RoeFluxT<2>(Uthis[ii], Uthat[ii], set, fluxes[ii], vmaxa, aima, N[i].dn[ii]);
                        if (moving)
                        {
                            fluxes[ii][1] += Uwall[ii][1] * fluxes[ii][0];
                            fluxes[ii][2] += Uwall[ii][2] * fluxes[ii][0];
                        }
                        fluxes[ii] = pnodes[i].getUCent(fluxes[ii], ii) * pnodes[i].w[ii];
                        break;

                    default:
                        break;
                    }
                }
                // if (viscous)
                //     for (int ii = 0; ii < 4; ii++)
                //     {
                //         NewtonianViscosFlux(Uthis[ii], Uthat[ii], set, vf, vmax, N[i].d[ii]);
                //         fluxes[ii] -= pnodes[i].getUCent(vf, ii) * pnodes[i].w[ii];
                //     }
                increment[i] = -(fluxes[0] + fluxes[1] + fluxes[2] + fluxes[3]) * (1. / pnodes[i].area);
                vmaxs[i] = vmaxa;
                aims[i] = aima;
                // if (aims[i])
                //     aims[i] = aima;
            }
            for (size i = 0; i < d.size(); i++)
            {
                if (vmaxs[i] > vmax)
                    vmax = vmaxs[i];
                if (aims[i])
                    aim = 1;
            }
        }

        std::tuple<double, utype<4>, int> StepRK4(double tvdOverride, const RoeSet &set, double CFL, double dtmax, double dtend, bool &conclude)
        {
            double vmax = 0.0;
            int aim = 0;
            dU0.resize(d.size()), dU1.resize(d.size()), dU2.resize(d.size()), dU3.resize(d.size()), U1.resize(d.size()), U2.resize(d.size()), U3.resize(d.size());
            GetIncrement(dU0, tvdOverride, set, vmax, aim);
            double dt = std::min(CFL / vmax, dtmax);
            if (dt >= dtend)
            {
                dt = dtend;
                conclude = true;
            }
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U1[i] = d[i] + dU0[i] * (dt * 0.5);
            std::swap(U1, d);
            GetIncrement(dU1, tvdOverride, set, vmax, aim);
            std::swap(U1, d);
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U2[i] = d[i] * (649.0 / 1600.0) - dU0[i] * (dt * 10890423.0 / 25193600.0) +
                        U1[i] * (951.0 / 1600.0) + dU1[i] * (dt * 5000.0 / 7873.0);
            std::swap(U2, d);
            GetIncrement(dU2, tvdOverride, set, vmax, aim);
            std::swap(U2, d);
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U3[i] = d[i] * (53989.0 / 2500000.0) - dU0[i] * (dt * 102261.0 / 5000000.0) +
                        U1[i] * (4806213.0 / 20000000.0) - dU1[i] * (dt * 5121.0 / 20000.0) +
                        U2[i] * (23619.0 / 32000.0) + dU2[i] * (dt * 7873.0 / 10000.0);
            std::swap(U3, d);
            GetIncrement(dU3, tvdOverride, set, vmax, aim);
            std::swap(U3, d);

            //No parallel
            utype<4> incmax(0.0);
            for (int i = 0; i < d.size(); i++)
            {
                utype<4> inc = d[i] * (-0.8) + dU0[i] * (dt * 0.1) +
                               U1[i] * (6127.0 / 30000.0) + dU1[i] * (dt * 1.0 / 6.0) +
                               U2[i] * (7873.0 / 30000.0) +
                               U3[i] * (1.0 / 3.0) + dU3[i] * (dt * 1.0 / 6.0);
                d[i] += inc, incmax = incmax.max(inc.abs());
            }

            if (aim)
            {
                std::cerr << "WARING: RK4 AIM!" << std::endl;
                //exit(-1);
            }
            return std::make_tuple(dt, incmax, aim);
        } //https://www.ece.uvic.ca/~bctill/papers/numacoust/Gottlieb_Shu_1998.pdf

        std::tuple<double, utype<4>, int> StepRK4Steady(double tvdOverride, const RoeSet &set, double CFL, double dtmax, double dtend, bool &conclude)
        {
            double vmax = 0.0;
            int aim = 0;
            dU0.resize(d.size()), dU1.resize(d.size()), dU2.resize(d.size()), dU3.resize(d.size()), U1.resize(d.size()), U2.resize(d.size()), U3.resize(d.size());
            GetIncrement(dU0, tvdOverride, set, vmax, aim);
            double dt = std::min(CFL / vmax, dtmax);
            std::swap(vmaxsC, vmaxs);
            if (dt >= dtend)
            {
                dt = dtend;
                conclude = true;
            }
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U1[i] = d[i] + dU0[i] * (CFL / vmaxsC[i] * 0.5);
            std::swap(U1, d);
            GetIncrement(dU1, tvdOverride, set, vmax, aim);
            std::swap(U1, d);
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U2[i] = d[i] * (649.0 / 1600.0) - dU0[i] * (CFL / vmaxsC[i] * 10890423.0 / 25193600.0) +
                        U1[i] * (951.0 / 1600.0) + dU1[i] * (CFL / vmaxsC[i] * 5000.0 / 7873.0);
            std::swap(U2, d);
            GetIncrement(dU2, tvdOverride, set, vmax, aim);
            std::swap(U2, d);
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U3[i] = d[i] * (53989.0 / 2500000.0) - dU0[i] * (CFL / vmaxsC[i] * 102261.0 / 5000000.0) +
                        U1[i] * (4806213.0 / 20000000.0) - dU1[i] * (CFL / vmaxsC[i] * 5121.0 / 20000.0) +
                        U2[i] * (23619.0 / 32000.0) + dU2[i] * (CFL / vmaxsC[i] * 7873.0 / 10000.0);
            std::swap(U3, d);
            GetIncrement(dU3, tvdOverride, set, vmax, aim);
            std::swap(U3, d);

            //No parallel
            utype<4> incmax(0.0);
            for (int i = 0; i < d.size(); i++)
            {
                utype<4> inc = d[i] * (-0.8) + dU0[i] * (CFL / vmaxsC[i] * 0.1) +
                               U1[i] * (6127.0 / 30000.0) + dU1[i] * (CFL / vmaxsC[i] * 1.0 / 6.0) +
                               U2[i] * (7873.0 / 30000.0) +
                               U3[i] * (1.0 / 3.0) + dU3[i] * (CFL / vmaxsC[i] * 1.0 / 6.0);
                d[i] += inc, incmax = incmax.max(inc.abs());
            }

            if (aim)
            {
                std::cerr << "WARING: RK4 AIM!" << std::endl;
                //exit(-1);
            }
            return std::make_tuple(dt, incmax, aim);
        } //https://www.ece.uvic.ca/~bctill/papers/numacoust/Gottlieb_Shu_1998.pdf

        std::tuple<double, utype<4>, int> StepEuler(double tvdOverride, const RoeSet &set, double CFL, double dtmax, double dtend, bool &conclude,
                                                    utype<4> absthres, int iitmax, double release)
        {
            double vmax = 0.0;
            int aim = 0;
            dU0.resize(d.size()), dU1.resize(d.size()), U1.resize(d.size()), U2.resize(d.size());
            GetIncrement(dU0, tvdOverride, set, vmax, aim);
            double dt = std::min(CFL / vmax, dtmax);
            if (dt >= dtend)
            {
                dt = dtend;
                conclude = true;
            }
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
                U2[i] = d[i] + dU0[i] * (release * dt * 1.0);
            std::swap(U2, d);
            dataFix(1e-6, set);
            std::swap(U2, d);
            auto inc0 = maxabs<4>(dU0);
            inc0 *= dt;
            std::cout << std::scientific << std::setw(15) << "\t\tEuler IIter = " << 0 << "  IncInit " << inc0[0] << "  " << inc0[1] << "  "
                      << inc0[2] << "  " << inc0[3] << "  " << std::endl;
            for (int iiter = 1; iiter <= iitmax; iiter++)
            {
                std::swap(U2, d); //TODO
                aim = 0;
                GetIncrement(dU0, tvdOverride, set, vmax, aim);
                std::swap(U2, d);
                utype<4> incU2(0.0);
                for (int i = 0; i < d.size(); i++)
                {
                    //auto inc = dU0[i] * (release / vmaxs[i]) + (d[i] - U2[i]) * (release / vmaxs[i] / dt);
                    //incU2 = incU2.max(inc.abs()), U2[i] += inc;
                    auto inc = dU0[i] * dt + d[i] - U2[i];
                    incU2 = incU2.max(inc.abs()), U2[i] += inc * release;
                }
                incU2 /= inc0;
                std::cout << std::scientific << std::setw(15) << "\t\tEuler IIter = " << iiter << "  IncRela " << incU2[0] << "  " << incU2[1] << "  "
                          << incU2[2] << "  " << incU2[3] << "  " << std::endl;
                std::swap(U2, d);
                dataFix(1e-6, set);
                std::swap(U2, d);
                bool end = true;
                for (int dim = 0; dim < 4; dim++)
                    if (incU2[dim] > absthres[dim])
                        end = false;
                if (end)
                    break;
            }

            utype<4> incmax(0.0); //TODO
            for (int i = 0; i < d.size(); i++)
            {
                incmax = incmax.max((d[i] - U2[i]).abs());
                d[i] = U2[i];
            }
            if (aim)
            {
                std::cerr << "===WARING===: EULER AIM!" << std::endl;
                //exit(-1);
                //conclude = true;
            }
            return std::make_tuple(dt, incmax, aim);
        }

        void dataIncrement(vector<utype<4>> &inc) //warning size of inc not checked
        {
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < d.size(); i++)
                d[i] += inc[i];
        }

        void dataFix(double fixratio, const RoeSet &set)
        {
            double maxP = 0, maxrho = 0;
            for (auto x : d)
            {
                double p = (x[4 - 1] - 0.5 * (sqr(x[1]) + sqr(x[2])) / x[0] * 0.5) * (set.gamma - 1);
                maxP = max(p, maxP), maxrho = max(x[0], maxrho);
            }
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (int i = 0; i < d.size(); i++)
            {
                auto &x = d[i];
                if (x[0] < fixratio * maxrho)
                {
                    x[0] = fixratio * maxrho;
                    x[1] = x[2] = 0;
                    double p = (x[4 - 1] - 0.5 * (sqr(x[1]) + sqr(x[2])) / x[0] * 0.5) * (set.gamma - 1);
                    x[3] = p / (set.gamma - 1);
                    //std::cout << "Fixed " << i << " Rho" << std::endl;
                    aims[i] |= 2;
                }
                double p = (x[4 - 1] - 0.5 * (sqr(x[1]) + sqr(x[2])) / x[0] * 0.5) * (set.gamma - 1);
                if (p < fixratio * maxP)
                {
                    aims[i] |= 4;
                    x[4 - 1] = 0.5 * (sqr(x[1]) + sqr(x[2])) / x[0] + maxP * fixratio * 2 / (set.gamma - 1);
                }
            }
        }

        utype<4> dataIncrementFrom(vector<utype<4>> &inc, vector<utype<4>> &dout) //warning size of inc not checked
        {
            /* #ifdef OMP_ON
#pragma omp parallel for
#endif */
            utype<4> ret;
            for (size i = 0; i < FDIM; i++)
                ret[i] = 0;
            for (size i = 0; i < d.size(); i++)
            {
                ret = ret.max((dout[i] + inc[i] - d[i]).abs());
                d[i] = dout[i] + inc[i];
            }
            return ret;
        }

        void xietaLaplacianFix(double fixratio)
        {
            dtemp = d;
            for (int j = 0; j < nn; j++)
                for (int i = 1; i < mm - 1; i++)
                {
                    for (int dd = 0; dd < FDIM; dd++)
                    {
                        real d0 = d[ind(i, j)][dd] - d[ind(i - 1, j)][dd];
                        real d1 = d[ind(i, j)][dd] - d[ind(i + 1, j)][dd];
                        real dm = (d[ind(i, j)][dd] * 2 + d[ind(i + 1, j)][dd] + d[ind(i - 1, j)][dd]) * 0.25;
                        if (1)
                        {
                            real rat = fixratio * abs(d1 * d0) / (sqr(dm) + abs(d1 * d0));
                            dtemp[ind(i, j)][dd] = dm * fixratio + d[ind(i, j)][dd] * (1 - fixratio);
                        }
                    }
                }
            d = dtemp;
            dtemp = d;
            for (int j = 1; j < nn - 1; j++)
                for (int i = 0; i < mm; i++)
                {
                    for (int dd = 0; dd < FDIM; dd++)
                    {
                        real d0 = d[ind(i, j)][dd] - d[ind(i, j - 1)][dd];
                        real d1 = d[ind(i, j)][dd] - d[ind(i, j - 1)][dd];
                        real dm = (d[ind(i, j)][dd] * 2 + d[ind(i, j - 1)][dd] + d[ind(i, j + 1)][dd]) * 0.25;
                        if (1)
                        {
                            real rat = fixratio * abs(d1 * d0) / (sqr(dm) + abs(d1 * d0));
                            dtemp[ind(i, j)][dd] = dm * fixratio + d[ind(i, j)][dd] * (1 - fixratio);
                        }
                    }
                }
            d = dtemp;
        }

        void fieldout(std::ostream &out)
        {
            for (size i = 0; i < d.size(); i++)
            {
                for (int ii = 0; ii < 4; ii++)
                    out << pnodes[i].d[ii][0] << '\t';
                for (int ii = 0; ii < 4; ii++)
                    out << pnodes[i].d[ii][1] << '\t';

                out << pnodes[i].bary[0] << '\t';
                out << pnodes[i].bary[1] << '\t';

                for (int ii = 0; ii < FDIM; ii++)
                    out << d[i][ii] << '\t';

                out << "\n";
            }
        }

        void wallLoadCircleOut(vector<real> &loadout) // only when the load target is really one circle
        {
            loadout.resize(wallload.size());
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 1; i < nwall; i++)
            {
                loadout[i * 2 + 0] = 0.5 * (wallload[(i - 1) * 2 + 0] + wallload[i * 2 + 0]);
                loadout[i * 2 + 1] = 0.5 * (wallload[(i - 1) * 2 + 1] + wallload[i * 2 + 1]);
            }
            loadout[0 * 2 + 0] = 0.5 * (wallload[(nwall - 1) * 2 + 0] + wallload[0 * 2 + 0]);
            loadout[0 * 2 + 1] = 0.5 * (wallload[(nwall - 1) * 2 + 1] + wallload[0 * 2 + 1]);
        }

        void fieldoutTecPlotDat(std::ostream &out)
        {
            out << "VARIABLES = \"X\", \"Y\", \"rho\", \"U\", \"V\", \"E\"\n";
            out << "ZONE N=" << (mm + 1) * (nn + 1) << " , E=" << mm * nn << ", ZONETYPE=FEQUADRILATERAL\n";
            out << "DATAPACKING = BLOCK\n";
            out << "VARLOCATION=([1-2]=NODAL, [3-6]=CELLCENTERED)\n";
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << pnodes[ind(i, j)].d[0][0] << '\n';
                    if (i == mm - 1)
                        out << pnodes[ind(i, j)].d[3][0] << '\n';
                }
            for (size i = 0; i < mm; i++)
            {
                out << pnodes[ind(i, nn - 1)].d[1][0] << '\n';
                if (i == mm - 1)
                    out << pnodes[ind(i, nn - 1)].d[2][0] << '\n';
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << pnodes[ind(i, j)].d[0][1] << '\n';
                    if (i == mm - 1)
                        out << pnodes[ind(i, j)].d[3][1] << '\n';
                }
            for (size i = 0; i < mm; i++)
            {
                out << pnodes[ind(i, nn - 1)].d[1][1] << '\n';
                if (i == mm - 1)
                    out << pnodes[ind(i, nn - 1)].d[2][1] << '\n';
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][1] / d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][2] / d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][3] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << i + 0 + (j + 0) * (mm + 1) + 1 << '\t'
                        << i + 0 + (j + 1) * (mm + 1) + 1 << '\t'
                        << i + 1 + (j + 1) * (mm + 1) + 1 << '\t'
                        << i + 1 + (j + 0) * (mm + 1) + 1 << '\n';
                }
        }

        void fieldoutTecPlotDatDebug(std::ostream &out)
        {
            out << "VARIABLES = \"X\", \"Y\", \"rho\", \"U\", \"V\", \"E\", \"AIM\"\n";
            out << "ZONE N=" << (mm + 1) * (nn + 1) << " , E=" << mm * nn << ", ZONETYPE=FEQUADRILATERAL\n";
            out << "DATAPACKING = BLOCK\n";
            out << "VARLOCATION=([1-2]=NODAL, [3-7]=CELLCENTERED)\n";
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << pnodes[ind(i, j)].d[0][0] << '\n';
                    if (i == mm - 1)
                        out << pnodes[ind(i, j)].d[3][0] << '\n';
                }
            for (size i = 0; i < mm; i++)
            {
                out << pnodes[ind(i, nn - 1)].d[1][0] << '\n';
                if (i == mm - 1)
                    out << pnodes[ind(i, nn - 1)].d[2][0] << '\n';
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << pnodes[ind(i, j)].d[0][1] << '\n';
                    if (i == mm - 1)
                        out << pnodes[ind(i, j)].d[3][1] << '\n';
                }
            for (size i = 0; i < mm; i++)
            {
                out << pnodes[ind(i, nn - 1)].d[1][1] << '\n';
                if (i == mm - 1)
                    out << pnodes[ind(i, nn - 1)].d[2][1] << '\n';
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][1] / d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][2] / d[ind(i, j)][0] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << d[ind(i, j)][3] << '\n';
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    out << aims[ind(i, j)] << '\n';

            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    out << i + 0 + (j + 0) * (mm + 1) + 1 << '\t'
                        << i + 0 + (j + 1) * (mm + 1) + 1 << '\t'
                        << i + 1 + (j + 1) * (mm + 1) + 1 << '\t'
                        << i + 1 + (j + 0) * (mm + 1) + 1 << '\n';
                }
        }

        void fieldintTecPlotDat(std::istream &in)
        {
            std::string buf;
            std::getline(in, buf);
            std::getline(in, buf);
            std::getline(in, buf);
            std::getline(in, buf);

            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    std::getline(in, buf);
                    if (i == mm - 1)
                        std::getline(in, buf);
                }
            for (size i = 0; i < mm; i++)
            {
                std::getline(in, buf);
                if (i == mm - 1)
                    std::getline(in, buf);
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    std::getline(in, buf);
                    if (i == mm - 1)
                        std::getline(in, buf);
                }
            for (size i = 0; i < mm; i++)
            {
                std::getline(in, buf);
                if (i == mm - 1)
                    std::getline(in, buf);
            }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    in >> d[ind(i, j)][0];
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    in >> d[ind(i, j)][1];
                    d[ind(i, j)][1] *= d[ind(i, j)][0];
                }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                {
                    in >> d[ind(i, j)][2];
                    d[ind(i, j)][2] *= d[ind(i, j)][0];
                }
            for (size j = 0; j < nn; j++)
                for (size i = 0; i < mm; i++)
                    in >> d[ind(i, j)][3];
        }
    };

    class SGSolver
    {

    public:
    };

#undef vector

}