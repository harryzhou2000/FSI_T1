#pragma once

#include "utype_basics.hpp"

namespace MeshBasic
{
    using namespace UTBasic;
    template <size nn>
    struct gridNeighbour
    {
        size n[nn];        //neighbour tetraindex (inner) or boundaryindex (notinner)
        size nthatidx[nn]; //index for otherside
        ntype t[nn];       //neighbour type
        utype<3> br[nn];   //bary relatively
        real d[nn];        //distance
        real dn[nn];       //distance normal
        real dp[nn];       //distance portion
    };

    template <size nn>
    struct gridNeighbourSG
    {
        size n[nn];
        ntype t[nn]; //neighbour type
        real d[nn];  //distance
        real dn[nn]; //distance normal
        real dp[nn]; //distance portion
    };

    template <size np, size nu, size ustart>
    struct TDNodes
    {
        real d[np][2]; //input data
        // cauculated data
        real bary[2];
        real nout[np][2];
        real w[np];
        real area;
        inline void getAreaBary()
        {
            area = 0;
            bary[0] = 0;
            bary[1] = 0;
            real x0 = d[0][0], y0 = d[0][1];
            for (int i = 2; i < np; i++)
            {
                real x1 = d[i - 1][0], y1 = d[i - 1][1];
                real x2 = d[i][0], y2 = d[i][1];
                real dx1 = x2 - x0, dy1 = y2 - y0, dx0 = x1 - x0, dy0 = y1 - y0;
                real areak = dx0 * dy1 - dx1 * dy0;
                area += areak;
                bary[0] += areak * (dx0 + dx1);
                bary[1] += areak * (dy0 + dy1);
            }
            bary[0] /= area * 3, bary[1] /= area * 3;
            area *= 0.5;
            bary[0] += x0;
            bary[1] += y0;
        }

        inline void getNW()
        {
            for (int i = 0; i < np; i++)
            {
                real x0 = d[i][0], y0 = d[i][1], x1 = d[(i + 1) % np][0], y1 = d[(i + 1) % np][1];
                real dx = x1 - x0, dy = y1 - y0;
                w[i] = sqrt(sqr(dx) + sqr(dy));
                nout[i][0] = dy / w[i];
                nout[i][1] = -dx / w[i];
            }
        }

        utype<nu> getUOut(const utype<nu> &U, size e)
        {
            utype<nu> nU;
            for (size i = 0; i < ustart; i++)
                nU[i] = U[i];
            for (size i = ustart + 2; i < nu; i++)
                nU[i] = U[i];
            nU[ustart + 0] = nout[e][0] * U[ustart + 0] + nout[e][1] * U[ustart + 1];
            nU[ustart + 1] = nout[e][0] * U[ustart + 1] - nout[e][1] * U[ustart + 0];
            return nU;
        }

        utype<nu> getUCent(const utype<nu> &UO, size e)
        {
            utype<nu> nU;
            for (size i = 0; i < ustart; i++)
                nU[i] = UO[i];
            for (size i = ustart + 2; i < nu; i++)
                nU[i] = UO[i];
            nU[ustart + 0] = nout[e][0] * UO[ustart + 0] - nout[e][1] * UO[ustart + 1];
            nU[ustart + 1] = nout[e][0] * UO[ustart + 1] + nout[e][1] * UO[ustart + 0];
            return nU;
        }
    };

    struct TetraNodes
    {
        size ivert[4];
        //+++++ now using face acquisition with std::map
        //size iface[4];
        //use face by default:
        /*
         3
        | \   
        |  2
        0----1
        0 -{0 2 1}
        1 -{0 1 3}
        2 -{2 0 3}
        3 -{1 2 3}
    */
        static const int fce2pidx[4][3]; //the mat above, to define face-vert
        utype<3> myPoints[4];            //locally stored coords;
        utype<3> nout[4];                //norm out(axis 1)
        utype<3> noutV[4];               //out normed axis 2
        utype<3> noutW[4];               //out normed axis 3
        real farea[4];                   //face area
        utype<3> bary;                   //barycenter
        utype<3> fbaryrel[4];            //facebarycenter
        real vol;                        //volume

        void __getVolume(const std::vector<utype<3>> &pointTable)
        {
            const utype<3> &p0 = pointTable[ivert[0]];
            const utype<3> &p1 = pointTable[ivert[1]];
            const utype<3> &p2 = pointTable[ivert[2]];
            const utype<3> &p3 = pointTable[ivert[3]];
            real nvol = (p3 - p0).dot(crossUT3(p1 - p0, p0 - p2));
            vol = std::abs(nvol) / 6.0;
            if (nvol > 0.0) // correct the left-right hand of the element
            {
                size mid = ivert[1];
                ivert[1] = ivert[2];
                ivert[2] = mid;
                //+++++ now using face acquisition with std::map
                // size fmid = iface[1];
                // iface[1] = iface[2];
                // iface[2] = mid;
            };
            if (vol <= __FLT_MIN__)
                exit(-2);
        }

        void getAreaBary(const std::vector<utype<3>> &pointTable)
        {
            __getVolume(pointTable);
            const utype<3> &p0 = myPoints[0] = pointTable[ivert[0]];
            const utype<3> &p1 = myPoints[1] = pointTable[ivert[1]];
            const utype<3> &p2 = myPoints[2] = pointTable[ivert[2]];
            const utype<3> &p3 = myPoints[3] = pointTable[ivert[3]];

            utype<3> l03 = p3 - p0;
            utype<3> l02 = p2 - p0;
            utype<3> l01 = p1 - p0;
            bary = p0 + (l03 + l02 + l01) * 0.25;

            nout[0] = crossUT3(l02, l01);
            nout[1] = crossUT3(l01, l03);
            nout[2] = crossUT3(l03, l02);
            nout[3] = crossUT3(p2 - p1, p3 - p2);
            farea[0] = std::sqrt(nout[0].dot(nout[0])) * 0.5;
            farea[1] = std::sqrt(nout[1].dot(nout[1])) * 0.5;
            farea[2] = std::sqrt(nout[2].dot(nout[2])) * 0.5;
            farea[3] = std::sqrt(nout[3].dot(nout[3])) * 0.5;
            nout[0] *= 0.5 / farea[0];
            nout[1] *= 0.5 / farea[1];
            nout[2] *= 0.5 / farea[2];
            nout[3] *= 0.5 / farea[3];
            noutV[0] = l02;
            noutV[1] = l01;
            noutV[2] = -l02;
            noutV[3] = p2 - p1;
            noutV[0].normalize(), noutV[1].normalize(), noutV[2].normalize(), noutV[3].normalize();
            noutW[0] = crossUT3(nout[0], noutV[0]);
            noutW[1] = crossUT3(nout[1], noutV[1]);
            noutW[2] = crossUT3(nout[2], noutV[2]);
            noutW[3] = crossUT3(nout[3], noutV[3]);
            fbaryrel[0] = triangleBary(p0, p2, p1) - bary;
            fbaryrel[1] = triangleBary(p0, p1, p3) - bary;
            fbaryrel[2] = triangleBary(p2, p0, p3) - bary;
            fbaryrel[3] = triangleBary(p1, p2, p3) - bary;
        }
        template <size nu, size ustart>
        utype<nu> getUOut(const utype<nu> &ucent, size nfce)
        {
            utype<3> ucentU;
            for (int i = 0; i < 3; i++)
                ucentU[i] = ucent[i + ustart];
            real magnitude = nout[nfce].dot(ucentU);
            real magnitudeV = noutV[nfce].dot(ucentU);
            real magnitudeW = noutW[nfce].dot(ucentU);
            utype<nu> ret(ucent);
            ret[ustart + 0] = magnitude;
            ret[ustart + 1] = magnitudeV;
            ret[ustart + 2] = magnitudeW;
            return ret;
        }
        template <size nu, size ustart>
        utype<nu> getUCent(const utype<nu> &uout, size nfce)
        {
            utype<3> ucentU = nout[nfce] * uout[ustart + 0] +
                              noutV[nfce] * uout[ustart + 1] +
                              noutW[nfce] * uout[ustart + 2];

            utype<nu> ret(uout);
            for (int i = 0; i < 3; i++)
                ret[i + ustart] = ucentU[i];
            return ret;
        }

#define eps 1e-40f
        //L & R are all positive
        template <size cnU>
        static inline utype<cnU> safePositiveDivide(const utype<cnU> &L, const utype<cnU> &R)
        {
            utype<cnU> ret;
            for (size i = 0; i < cnU; i++)
                ret[i] = (L[i] + eps) / (R[i] + eps);
            return ret;
        }
#undef eps

        //after getAreaBary, U[0~3]=uthat U[4]=uthis
        template <size cnU>
        void linearReconstruction(real steer, const gridNeighbour<4> &neighbour,
                                  const std::vector<utype<3>> &points,
                                  const utype<cnU> *U, //5-fce+cent
                                  mtype<cnU, 3> &grad, utype<cnU> &maxU, utype<cnU> &minU,
                                  utype<cnU> *Urec //4-fce
        )
        {
            maxU = U[4], minU = U[4];
            grad.set0();
            for (int i = 0; i < 4; i++)
            {
                maxU = maxU.max(U[i]);
                minU = minU.min(U[i]);
                for (int id = 0; id < 3; id++)
                {
                    utype<cnU> uface = U[i] * neighbour.dp[i] + U[4] * (1 - neighbour.dp[i]);
                    for (size iU = 0; iU < cnU; iU++)
                        grad(iU, id) += uface[iU] * (farea[i] * nout[i][id]);
                }
            }
            utype<cnU> maxUrel = maxU - U[4], minUrel = minU - U[4]; //relative to U[4]
            for (int id = 0; id < 3; id++)
                for (size iU = 0; iU < cnU; iU++)
                    grad(iU, id) /= vol;

            // barth limiter
            utype<cnU> maxUrec(0.0), minUrec(0.0); //relative to U[4]
            for (int inode = 0; inode < 4; inode++)
            {
                //utype<cnU> unode = grad * (points[ivert[inode]] - bary);
                utype<cnU> unode = grad * (myPoints[inode] - bary);
                maxUrec = maxUrec.max(unode);
                minUrec = minUrec.min(unode);
            }

            utype<cnU> compress(steer);
            compress = compress.min(safePositiveDivide<5>(maxUrel.abs(), maxUrec.abs()));
            compress = compress.min(safePositiveDivide<5>(minUrel.abs(), minUrec.abs()));
            for (size iU = 0; iU < cnU; iU++)
                for (int id = 0; id < 3; id++)
                    grad(iU, id) *= compress[iU];

            //get rec value
            for (int ifce = 0; ifce < 4; ifce++)
            {
                Urec[ifce] = U[4] + grad * fbaryrel[ifce];
                //Urec[ifce] = Urec[ifce].min(maxU).max(minU);
            }
        }
    };

    const int TetraNodes::fce2pidx[4][3] = {{0, 2, 1}, {0, 1, 3}, {2, 0, 3}, {1, 2, 3}};

    template <size np>
    struct TDNodesE
    {
        real ddt[np][2];
    };

    template <size nu>
    struct TetraMesh
    {
        //ntets
        std::vector<TetraNodes> tets;
        std::vector<gridNeighbour<4>> neighbours;
        //npoints
        std::vector<utype<3>> points;
        std::vector<std::pair<ntype, size>> pointntype; //boundtype + boundsetkey
        //nfaces
        std::map<std::tuple<size, size, size>, std::pair<size, size>> fcemap;
        //nbsts boundset key: valid if != 0
        std::vector<utype<nu>> boundset;

        real MinVol = 0.0;
        real MaxVol = 0.0;

        //use this to create unique tiplet key for each face
        static void FaceVertIdxAscend(size &v1, size &v2, size &v3)
        {
            if (v1 > v2)
                std::swap(v1, v2);
            if (v2 > v3)
                std::swap(v2, v3);
            if (v1 > v2)
                std::swap(v1, v2);
        }

        void LoadMesh(const std::string filename, std::ostream &logout = std::cout)
        {
            std::ifstream fin(filename);
            if (!fin)
            {
                logout << "File Open Faliure: " << filename << std::endl;
                exit(-1);
            }
            std::string name;
            size npoints, ntets, nbounds;
            //points
            fin >> name >> npoints;
            logout << "READ: " << name << std::endl;
            points.resize(npoints);
            for (size i = 0; i < npoints; i++)
                fin >> points[i][0] >> points[i][1] >> points[i][2];

            //tets
            fin >> name >> ntets;
            logout << "READ: " << name << std::endl;
            tets.resize(ntets);
            for (size i = 0; i < ntets; i++)
                fin >> tets[i].ivert[0] >> tets[i].ivert[1] >> tets[i].ivert[2] >> tets[i].ivert[3];

            //bounds
            boundset.resize(1);
            boundset.reserve(10);
            fin >> name >> nbounds;
            logout << "READ: " << name << std::endl;
            pointntype.resize(npoints); // not nbounds
            for (size i = 0; i < npoints; i++)
                pointntype[i].first = inner, pointntype[i].second = 0;
            for (size i = 0; i < nbounds; i++)
            {
                size tochange;
                fin >> tochange;
                int t;
                fin >> t >> pointntype[tochange].second; // for a boundary node, corresponding boundset key is set from file
                pointntype[tochange].first = ntype(t);
            }
            fin.close();

            //init neighbours
            MinVol = DBL_MAX;
            MaxVol = 0.0;
            neighbours.resize(ntets); //all tets' neighbour obj
            for (size i = 0; i < ntets; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    size p1idx = tets[i].ivert[TetraNodes::fce2pidx[j][0]];
                    size p2idx = tets[i].ivert[TetraNodes::fce2pidx[j][1]];
                    size p3idx = tets[i].ivert[TetraNodes::fce2pidx[j][2]];
                    FaceVertIdxAscend(p1idx, p2idx, p3idx);
                    //init neighbour types
                    // neighbours[i].t[j] = deriveFtype3(
                    //     pointntype[p1idx].first,
                    //     pointntype[p2idx].first,
                    //     pointntype[p3idx].first); //default values
                    //start pairing if inner
                    //if (neighbours[i].t[j] == inner)
                    //{
                    auto foundfce = fcemap.find(std::make_tuple(p1idx, p2idx, p3idx));
                    if (foundfce == fcemap.end())
                    {
                        auto &mappair = fcemap[std::make_tuple(p1idx, p2idx, p3idx)];
                        mappair.first = i, mappair.second = -1;
                    }
                    else
                        foundfce->second.second = i;
                    //}
                }
                tets[i].getAreaBary(points);
                MinVol = std::min(MinVol, tets[i].vol);
                MaxVol = std::max(MaxVol, tets[i].vol);
            }

            //Derive who's neighbour and which boundary
            for (size i = 0; i < ntets; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    size p1idx = tets[i].ivert[TetraNodes::fce2pidx[j][0]];
                    size p2idx = tets[i].ivert[TetraNodes::fce2pidx[j][1]];
                    size p3idx = tets[i].ivert[TetraNodes::fce2pidx[j][2]];
                    FaceVertIdxAscend(p1idx, p2idx, p3idx);
                    auto mappair = fcemap[std::make_tuple(p1idx, p2idx, p3idx)];
                    if (mappair.first != -1 && mappair.second != -1) //inner is a singled face
                    {
                        neighbours[i].t[j] == inner;
                        if (i == mappair.first)
                        {
                            //assert(mappair.second != -1);
                            neighbours[i].n[j] = mappair.second;
                        }
                        else
                        {
                            //assert(mappair.second == i);
                            neighbours[i].n[j] = mappair.first;
                        }

                        neighbours[i].br[j] = tets[neighbours[i].n[j]].bary - tets[i].bary; //caution if cycle boundary
                        neighbours[i].d[j] = neighbours[i].br[j].norm2();
                        neighbours[i].dn[j] = pointToFace(tets[i].bary,
                                                          points[p1idx], points[p2idx], points[p3idx]);
                        real that = pointToFace(tets[i].bary + neighbours[i].br[j],
                                                points[p1idx], points[p2idx], points[p3idx]);
                        neighbours[i].dp[j] = neighbours[i].dn[j] / (neighbours[i].dn[j] + that);
                        neighbours[i].dn[j] += that;
                    }
                    else
                    {
                        neighbours[i].t[j] = deriveFtype3(
                            pointntype[p1idx].first,
                            pointntype[p2idx].first,
                            pointntype[p3idx].first); //default values
                        neighbours[i].n[j] = deriveNbound3(pointntype[p1idx].first,
                                                           pointntype[p2idx].first,
                                                           pointntype[p3idx].first,
                                                           pointntype[p1idx].second,
                                                           pointntype[p2idx].second,
                                                           pointntype[p3idx].second);
                        neighbours[i].d[j] = neighbours[i].dn[j] =
                            2 * pointToFace(tets[i].bary,
                                            points[p1idx],
                                            points[p2idx],
                                            points[p3idx]);
                        neighbours[i].br[j] = tets[i].nout[j] * neighbours[i].dn[j];
                        neighbours[i].dp[j] = 0.5;
                    }
                }
            }
            logout << "Mesh Read; Min Vol: " << MinVol << " Max Vol: " << MaxVol << std::endl;
            for (size i = 0; i < ntets; i++)
                for (int j = 0; j < 4; j++)
                    for (int jj = 0; jj < 4; jj++)
                        if (neighbours[neighbours[i].n[j]].n[jj] == i)
                            neighbours[i].nthatidx[j] = jj;
        }

        bool checkAreInters(const std::vector<size> interPoints) const
        {
            for (auto i : interPoints)
                if (pointntype[i].first != inter)
                    return false;
            return true;
        }

        const std::vector<utype<3>> &showPoints() const
        {
            return points;
        }

        std::vector<utype<3>> &showPoints()
        {
            return points;
        }

        //must-do after you change points' coord
        void updateGeom()
        {
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < tets.size(); i++)
                tets[i].getAreaBary(points);
        }
    };

}
