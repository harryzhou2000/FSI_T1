#pragma once
#include "roe.hpp"
#include <iomanip>

namespace UGSover
{
    using namespace RiemannSolver;
    using ut3 = utype<3>;
    using ut5 = utype<5>;

    struct UrecTetra
    {
        ut5 u[4];
    };

    class RoeUGsolver : protected TetraMesh<5>
    {
    protected:
        std::vector<ut3> dpointsdt; //grid velocity

        std::vector<ut5> US;         //conservative field vec
        std::vector<UrecTetra> Urec; //reconstructed value (used for each Uinc)
        std::vector<int> aims;       //(used for each Uinc)
        std::vector<real> idtmaxs;   //(used for each Uinc)
        std::vector<ut5> uinc;       //(used for each timestep)
        std::vector<ut5> uincinc;    //(used for each timestep)
        std::vector<ut5> utemp;      //(used for each timestep)
        RoeSet set;
        bool ifuseInter = false;
        std::vector<ut3> interNodeForce;
        std::vector<int> forcecount;

    public:
        using TetraMesh::checkAreInters;
        void Load(const std::string &filename, const std::vector<utype<5>> &nboundset, std::ostream &logout = std::cout)
        {
            LoadMesh(filename, logout);
            boundset = nboundset;
            dpointsdt.assign(points.size(), ut3(0.0)); //set 0 for grid velocity
            US.resize(tets.size());
            Urec.resize(tets.size());
            aims.resize(tets.size()), idtmaxs.resize(tets.size());
            uincinc.resize(tets.size()), uinc.resize(tets.size()), utemp.resize(tets.size());
        }

        void setInter()
        {
            interNodeForce.resize(points.size());
            forcecount.resize(points.size());
        }

        void LoadTecDat(const std::string &meshfilename, const std::string &datfilename, std::vector<utype<5>> &nboundset, std::ostream &logout = std::cout)
        {
            LoadMesh(meshfilename, logout);
            boundset = nboundset;
            dpointsdt.assign(points.size(), ut3(0.0)); //set 0 for grid velocity
            US.resize(tets.size());
            Urec.resize(tets.size());
            aims.resize(tets.size()), idtmaxs.resize(tets.size());
            uincinc.resize(tets.size()), uinc.resize(tets.size()), utemp.resize(tets.size());
            fieldreadTecPlotDat(datfilename, logout);
        }

        void DerivedUdt(const std::vector<ut5> &Unow, std::vector<ut5> &Uinc, real steer, int &aim, real &idtmax)
        {
            // forall tet reconstruct values Urec
            //std::cout << "DBG\n";
            //std::cout << Unow;
            mtype<5, 3> grad;
            ut5 umax, umin;
#ifdef OMP_ON
#pragma omp parallel for private(grad, umax, umin)
#endif
            for (size i = 0; i < tets.size(); i++)
            {
                ut5 udata[5];
                udata[4] = Unow[i];
                for (size in = 0; in < 4; in++)
                {
                    switch (neighbours[i].t[in])
                    {
                    case inner:
                        udata[in] = Unow[neighbours[i].n[in]];
                        break;
                    case pressure:
                        udata[in] = boundset[neighbours[i].n[in]];
                        udata[in] = tets[i].getUCent<5, 1>(
                            getPressureBound(
                                tets[i].getUOut<5, 1>(udata[4], in),
                                tets[i].getUOut<5, 1>(udata[in], in), set),
                            in);
                        break;
                    case infinity:
                        udata[in] = boundset[neighbours[i].n[in]];
                        break;
                    case wall:
                    case inter:
                        udata[in] = tets[i].getUOut<5, 1>(udata[4], in);
                        udata[in][1] = -udata[in][1];
                        udata[in] = tets[i].getUCent<5, 1>(udata[in], in);
                        break;
                    default:
                        std::abort();
                        break;
                    }
                }
                tets[i].linearReconstruction<5>(steer, neighbours[i], points, udata, grad, umax, umin, Urec[i].u);
            }

            // forall tet calc face flux, derive aim & dtmax
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < tets.size(); i++)
            {
                aims[i] = 0;
                Uinc[i] = ut5(0.0);
                for (size in = 0; in < 4; in++)
                {
                    ut3 dUdt = (dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][0]]] +
                                dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][1]]] +
                                dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][2]]]) *
                               (1.0 / 3.0);
                    ut5 UL = Urec[i].u[in];
                    ut5 UR;
                    switch (neighbours[i].t[in])
                    {
                    case inner:
                        UR = Urec[neighbours[i].n[in]].u[neighbours[i].nthatidx[in]];
                        break;
                    case pressure:
                        UR = boundset[neighbours[i].n[in]];
                        UR = tets[i].getUCent<5, 1>(
                            getPressureBound(
                                tets[i].getUOut<5, 1>(UL, in),
                                tets[i].getUOut<5, 1>(UR, in), set),
                            in);
                        break;
                    case infinity:
                        UR = boundset[neighbours[i].n[in]];
                        break;
                    case wall:
                    case inter:
                        UR = tets[i].getUOut<5, 1>(UL, in);
                        UR[1] = -UR[1];
                        UR = tets[i].getUCent<5, 1>(UR, in);
                        break;
                    default:
                        std::abort();
                        break;
                    }
                    TransAddTo(UL, -dUdt), TransAddTo(UR, -dUdt); // frame change
                    UL = tets[i].getUOut<5, 1>(UL, in);
                    UR = tets[i].getUOut<5, 1>(UR, in);
                    real rhoUnmid = (UL[1] + UR[1]) * 0.5;
                    ut5 F;
                    RoeFlux2(UL, UR, set, F, idtmaxs[i], aims[i], neighbours[i].dn[in]);
                    ut5 Fcent = tets[i].getUCent<5, 1>(F, in);
                    Fcent[1] += rhoUnmid * dUdt[0], Fcent[2] += rhoUnmid * dUdt[1], Fcent[3] += rhoUnmid * dUdt[2];
                    Uinc[i] += Fcent * tets[i].farea[in];
                    for (int d = 0; d < 5; d++)
                        if (std::isnan(F[d]))
                        {
                            std::cerr << "NAN out of ROE" << std::endl;
                            exit(-5);
                        }
                }
                Uinc[i] *= (1.0 / tets[i].vol);
            }
            aim = 0;
            idtmax = 0;
            for (size i = 0; i < tets.size(); i++)
            {
                if (aims[i])
                    aim = 1;
                idtmax = std::max(idtmax, idtmaxs[i]);
            }
        }

        void DerivedUdtSteady(const std::vector<ut5> &Unow, std::vector<ut5> &Uinc, real steer, int &aim, real CFL)
        {
            // forall tet reconstruct values Urec

            mtype<5, 3> grad;
            ut5 umax, umin;
#ifdef OMP_ON
#pragma omp parallel for private(grad, umax, umin)
#endif
            for (size i = 0; i < tets.size(); i++)
            {
                ut5 udata[5];
                udata[4] = Unow[i];
                for (size in = 0; in < 4; in++)
                {
                    switch (neighbours[i].t[in])
                    {
                    case inner:
                        udata[in] = Unow[neighbours[i].n[in]];
                        break;
                    case pressure:
                        udata[in] = boundset[neighbours[i].n[in]];
                        udata[in] = tets[i].getUCent<5, 1>(
                            getPressureBound(
                                tets[i].getUOut<5, 1>(udata[4], in),
                                tets[i].getUOut<5, 1>(udata[in], in), set),
                            in);
                        break;
                    case infinity:
                        udata[in] = boundset[neighbours[i].n[in]];
                        break;
                    case wall:
                    case inter:
                        udata[in] = tets[i].getUOut<5, 1>(udata[4], in);
                        udata[in][1] = -udata[in][1];
                        udata[in] = tets[i].getUCent<5, 1>(udata[in], in);
                        break;
                    default:
                        std::abort();
                        break;
                    }
                }
                tets[i].linearReconstruction<5>(steer, neighbours[i], points, udata, grad, umax, umin, Urec[i].u);
            }

            // forall tet calc face flux, derive aim & dtmax
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < tets.size(); i++)
            {
                Uinc[i] = ut5(0.0);
                aims[i] = 0;
                for (size in = 0; in < 4; in++)
                {
                    ut3 dUdt = (dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][0]]] +
                                dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][1]]] +
                                dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][2]]]) *
                               (1.0 / 3.0);
                    ut5 UL = Urec[i].u[in];
                    //UL[1] -= dUdt[0], UL[2] -= dUdt[1], UL[3] -= dUdt[2];
                    ut5 UR;
                    switch (neighbours[i].t[in])
                    {
                    case inner:
                        UR = Urec[neighbours[i].n[in]].u[neighbours[i].nthatidx[in]];
                        break;
                    case pressure:
                        UR = boundset[neighbours[i].n[in]];
                        UR = tets[i].getUCent<5, 1>(
                            getPressureBound(
                                tets[i].getUOut<5, 1>(UL, in),
                                tets[i].getUOut<5, 1>(UR, in), set),
                            in);
                        break;
                    case infinity:
                        UR = boundset[neighbours[i].n[in]];
                        break;
                    case wall:
                    case inter:
                        UR = tets[i].getUOut<5, 1>(UL, in);
                        UR[1] = -UR[1];
                        UR = tets[i].getUCent<5, 1>(UR, in);
                        break;
                    default:
                        std::abort();
                        break;
                    }
                    //UR[1] -= dUdt[0], UR[2] -= dUdt[1], UR[3] -= dUdt[2];
                    TransAddTo(UL, -dUdt), TransAddTo(UR, -dUdt); // frame change
                    UL = tets[i].getUOut<5, 1>(UL, in);
                    UR = tets[i].getUOut<5, 1>(UR, in);
                    ut5 F;
                    RoeFlux2(UL, UR, set, F, idtmaxs[i], aims[i], neighbours[i].dn[in]);
                    Uinc[i] += tets[i].getUCent<5, 1>(F, in) * tets[i].farea[in];
                    for (int d = 0; d < 5; d++)
                        if (std::isnan(F[d]))
                        {
                            std::cerr << "NAN out of ROE" << std::endl;
                            exit(-5);
                        }
                }
                Uinc[i] *= (1.0 / tets[i].vol);
            }
            aim = 0;
            for (size i = 0; i < tets.size(); i++)
            {
                if (aims[i])
                    aim = 1;
                
            }
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for(size i = 0; i<tets.size(); i++)
                Uinc[i] *= CFL / idtmaxs[i];
        }

        real StepUEuler(double dtmax, double CFL, int iitmax, ut5 absresth, double release, double steer, std::ostream &logout = std::cout)
        {
            // get dUdt
            //std::cout << "re " << release << " CFL " << CFL << " resth " << absresth << " steer " << steer;
            if (release <= 0 || release > 1 || steer < 0 || steer > 1 || CFL >= 0)
            {
                std::cerr << "StepUEuler Innvalid Paramas!! "
                          << "Release " << release << " CFL " << CFL << " resth " << absresth << " steer " << steer
                          << std::endl;
            }
            real idt;
            int aim;
            real dt;
            DerivedUdt(US, uinc, steer, aim, idt);
            ut maxdif = maxabs<5>(uinc);
            ut maxdif0 = maxdif;
            dt = std::min(dtmax, CFL / idt);
            uinc *= dt;
            utemp = US;
            //utemp *= 1 / 1.0;
            utemp += uinc;
            //utemp *= 1.0;
            std::swap(uincinc, uinc);

            FIX(utemp);
            logout << "\t"
                   << "IIter = " << std::setw(4) << 0 << " AIM: " << aim << " Increment:   "
                   << std::scientific << std::setw(12)
                   << maxdif[0] << '\t' << maxdif[1] << '\t' << maxdif[2] << '\t' << maxdif[3] << '\t' << maxdif[4] << '\t' << "dt= " << dt << std::endl;
            for (int i = 1; i <= iitmax; i++)
            {
                DerivedUdt(utemp, uinc, steer, aim, idt);
                uincinc -= uinc;
                maxdif = maxabs<5>(uincinc);
                maxdif /= maxdif0;
                uincinc = uinc;
                logout << "\t"
                       << "IIter = " << std::setw(4) << i << " AIM: " << aim << " IncreChange: "
                       << std::scientific << std::setw(12)
                       << maxdif[0] << '\t' << maxdif[1] << '\t' << maxdif[2] << '\t' << maxdif[3] << '\t' << maxdif[4] << '\t' << std::endl;
                uinc *= dt;
                utemp *= (1 - release) / release;
                utemp += US;
                utemp += uinc;
                utemp *= release;
                FIX(utemp);

                int end = 1;
                for (int d = 0; d < 5; d++)
                    if (maxdif[d] > absresth[d])
                        end = false;
                if (end)
                    break;
            }

            std::swap(US, utemp);
            return dt;
            // increment if aim good
        }

        void StepUEulerSteady(double dtmax, double CFL, int iitmax, ut5 absresth, double release, double steer, std::ostream &logout = std::cout)
        {
            // get dUdt
            int aim;
            DerivedUdtSteady(US, uinc, 1.0, aim, CFL);
            ut maxdif = maxabs<5>(uinc);
            ut maxdif0 = maxdif;

            utemp = US;
            utemp += uinc;
            std::swap(uincinc, uinc);

            //FIX(utemp);
            logout << "\t"
                   << "STEADY IIter = " << std::setw(4) << 0 << " AIM: " << aim << " Increment:   "
                   << std::scientific << std::setw(12)
                   << maxdif[0] << '\t' << maxdif[1] << '\t' << maxdif[2] << '\t' << maxdif[3] << '\t' << maxdif[4] << std::endl;
            for (int i = 1; i <= iitmax; i++)
            {
                DerivedUdtSteady(utemp, uinc, steer, aim, CFL);
                uincinc -= uinc;
                maxdif = maxabs<5>(uincinc);
                maxdif /= maxdif0;
                uincinc = uinc;
                logout << "\t"
                       << "STEADY IIter = " << std::setw(4) << i << " AIM: " << aim << " IncreChange: "
                       << std::scientific << std::setw(12)
                       << maxdif[0] << '\t' << maxdif[1] << '\t' << maxdif[2] << '\t' << maxdif[3] << '\t' << maxdif[4] << '\t' << std::endl;
                utemp *= (1 - release) / release;
                utemp += US;
                utemp += uinc;
                utemp *= release;
                //FIX(utemp);

                int end = 1;
                for (int d = 0; d < 5; d++)
                    if (maxdif[d] > absresth[d])
                        end = false;
                if (end)
                    break;
            }

            std::swap(US, utemp);

            // increment if aim good
        }

        void FIX(std::vector<ut5> &US)
        {
            // for (int i = 0; i < US.size(); i++)
            // {
            // }
            ut maxU = maxabs<5>(US);

#ifdef OMP_ON
#pragma omp parallel for
#endif

            for (int i = 0; i < US.size(); i++)
            {
                if (US[i][0] < maxU[0] * 1e-4)
                {
                    US[i][0] = maxU[0] * 1e-4;
                    US[i][1] = US[i][2] = US[i][3] = 0;
                    //exit(3);
                }
                real rusqr = 0.5 * (US[i][1] * US[i][1] + US[i][2] * US[i][2] + US[i][3] * US[i][3]) / (US[i][0] + __FLT_MIN__);
                if (US[i][4] - rusqr < (maxU[4] - rusqr) * 1e-4)
                    US[i][4] = rusqr + (maxU[4] - rusqr) * 1e-4;
                //exit(2);
            }
        }

        void NANCheck(std::vector<ut5> &US)
        {
            for (int i = 0; i < US.size(); i++)
            {
                for (int d = 0; d < 5; d++)
                    if (std::isnan(US[i][d]))
                    {
                        std::cout << "=== ===NAN=== ===" << std::endl;
                        exit(-1);
                    }
            }
        }

        void RenewGrid(const std::vector<ut3> &npoints, const std::vector<ut3> &ndpointsdt)
        {
            points = npoints, dpointsdt = ndpointsdt;
        }

        ut3 calculateLoad(std::ostream &logout = std::cout)
        {
            interNodeForce.assign(interNodeForce.size(), ut3(0.0));
            forcecount.assign(forcecount.size(), 0);
            for (size i = 0; i < tets.size(); i++)
                for (size in = 0; in < 4; in++)
                    if (neighbours[i].t[in] == inter)
                    {
                        ut3 dUdt = (dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][0]]] +
                                    dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][1]]] +
                                    dpointsdt[tets[i].ivert[TetraNodes::fce2pidx[in][2]]]) *
                                   (1.0 / 3.0);
                        ut5 UL = Urec[i].u[in];
                        TransAddTo(UL, -dUdt); // frame change
                        UL = tets[i].getUOut<5, 1>(UL, in);
                        ut5 UR = UL;
                        UR[1] = -UL[1];
                        ut5 F;
                        RoeFlux2(UL, UR, set, F, idtmaxs[i], aims[i], neighbours[i].dn[in]);
                        F = tets[i].getUCent<5, 1>(F, in) * tets[i].farea[in];
                        ut3 FaceForce;
                        FaceForce[0] = F[1], FaceForce[1] = F[2], FaceForce[2] = F[3];
                        FaceForce *= (1.0 / 3.0);
                        for (int ip = 0; ip < 3; ip++)
                        {
                            size ipoint = tets[i].ivert[TetraNodes::fce2pidx[in][ip]];
                            interNodeForce[ipoint] += FaceForce;
                            forcecount[ipoint]++;
                        }
                    }
            ut3 TotalForce(0.0);
            for (size i = 0; i < points.size(); i++)
                if (pointntype[i].first == inter)
                    if (forcecount[i] > 0)
                    {
                        interNodeForce[i] /= real(forcecount[i]);
                        TotalForce += interNodeForce[i];
                    }
            logout << std::scientific << std::setprecision(8)
                   << "===Inter=== Got Fluid Inter Force, Total ["
                   << TotalForce[0] << ", " << TotalForce[1] << ", " << TotalForce[2] << "] "
                   << std::endl;
            return TotalForce;
        }

        void fieldprintTecPlotDat(const std::string &filename, std::ostream &logout = std::cout)
        {
            logout << "---FILE TO--- : " << filename << std::endl;
            std::ofstream out(filename);
            if (!out)
            {
                logout << "===<<= FILE OPEN FAILED =>>===" << std::endl;
                exit(-1);
            }
            out << "VARIABLES = \"X\", \"Y\", \"Z\", \"rho\", \"U\", \"V\", \"W\", \"E\"\n";
            out << "ZONE N=" << points.size() << " , E=" << tets.size() << ", ZONETYPE=FETETRAHEDRON\n";
            out << "DATAPACKING = BLOCK\n";
            out << "VARLOCATION=([1-3]=NODAL, [4-8]=CELLCENTERED)\n";
            out << std::scientific << std::setw(19) << std::setprecision(15);

            for (size it = 0; it < points.size(); it++)
                out << points[it][0] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][1] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][2] << '\n';

            for (size it = 0; it < tets.size(); it++)
                out << US[it][0] << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][1] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][2] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][3] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][4] << '\n';

            for (size it = 0; it < tets.size(); it++)
                out << tets[it].ivert[0] + 1 << '\t'
                    << tets[it].ivert[1] + 1 << '\t'
                    << tets[it].ivert[2] + 1 << '\t'
                    << tets[it].ivert[3] + 1 << '\n';

            out.close();
        }

        void fieldprintTecPlotDatDebug(const std::string &filename)
        {
            std::ofstream out(filename);
            out << "VARIABLES = \"X\", \"Y\", \"Z\", \"rho\", \"U\", \"V\", \"W\", \"E\", \"IS_EDGE\"\n";
            out << "ZONE N=" << points.size() << " , E=" << tets.size() << ", ZONETYPE=FETETRAHEDRON\n";
            out << "DATAPACKING = BLOCK\n";
            out << "VARLOCATION=([1-3]=NODAL, [4-8]=CELLCENTERED, [9]=NODAL)\n";

            for (size it = 0; it < points.size(); it++)
                out << points[it][0] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][1] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][2] << '\n';

            for (size it = 0; it < tets.size(); it++)
                out << US[it][0] << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][1] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][2] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][3] / (US[it][0] + 1e-100) << '\n';
            for (size it = 0; it < tets.size(); it++)
                out << US[it][4] << '\n';

            for (size it = 0; it < points.size(); it++)
                out << (int)(pointntype[it].first) << '\n';

            for (size it = 0; it < tets.size(); it++)
                out << tets[it].ivert[0] + 1 << '\t'
                    << tets[it].ivert[1] + 1 << '\t'
                    << tets[it].ivert[2] + 1 << '\t'
                    << tets[it].ivert[3] + 1 << '\n';

            out.close();
        }

        void fieldreadTecPlotDat(const std::string &filename, std::ostream &logout = std::cout)
        {
            logout << "---FILE FROM--- : " << filename << std::endl;
            std::ifstream in(filename);
            if (!in)
            {
                logout << "===<<= FILE OPEN FAILED =>>===" << std::endl;
                exit(-1);
            }

            char stall[512];
            in.getline(stall, 512);
            in.getline(stall, 512);
            in.getline(stall, 512);
            in.getline(stall, 512);
            for (size it = 0; it < 3 * points.size(); it++)
                in.getline(stall, 512);

            for (size it = 0; it < tets.size(); it++)
                in >> US[it][0];
            for (size it = 0; it < tets.size(); it++)
                in >> US[it][1];
            for (size it = 0; it < tets.size(); it++)
                in >> US[it][2];
            for (size it = 0; it < tets.size(); it++)
                in >> US[it][3];
            for (size it = 0; it < tets.size(); it++)
                in >> US[it][4];
            for (size it = 0; it < tets.size(); it++)
                US[it][1] *= US[it][0], US[it][2] *= US[it][0], US[it][3] *= US[it][0];

            in.close();
        }

        void InitByConst(const ut5 &initval)
        {
            US.assign(US.size(), initval);
        }

        void AlterByBox(const ut5 &initval,
                        real xmin, real xmax,
                        real ymin = -DBL_MAX, real ymax = DBL_MAX,
                        real zmin = -DBL_MAX, real zmax = DBL_MAX)
        {
            for (int i = 0; i < US.size(); i++)
            {
                ut3 lb, ub;
                lb[0] = xmin, lb[1] = ymin, lb[2] = zmin;
                ub[0] = xmax, ub[1] = ymax, ub[2] = zmax;
                if (inbox<3>(tets[i].bary, lb, ub))
                    US[i] = initval;
            }
        }

        void setRoeset(const RoeSet &nset)
        {
            set = nset;
        }

        using TetraMesh::showPoints;
        using TetraMesh::updateGeom;
        const std::vector<ut3> &showdPointsdt() const
        {
            return dpointsdt;
        }
        std::vector<ut3> &showdPointsdt()
        {
            return dpointsdt;
        }
        const std::vector<ut3> &showInterForces() const
        {
            return interNodeForce;
        }
    };

}