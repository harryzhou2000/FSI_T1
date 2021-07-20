#pragma once
#include "roeugsolver.hpp"
#include "femT4.hpp"
#include "utils.hpp"
#include <stdio.h>

namespace FSI
{
    using namespace FEM;
    using namespace UGSover;
    class UGFSISolver
    {
        FemT4Solver SolidSolver;
        RoeUGsolver FluidSolver;
        FemT4Solver FluidGridder;
        std::vector<size> Finter;
        std::vector<std::vector<size>> F2SIndex; //sizeof FInter * k1
        std::vector<size> Sinter;
        std::vector<std::vector<size>> S2FIndex; //sizeof FInter * k2

        std::vector<std::vector<real>> FluidGridMode;
        std::vector<std::vector<real>> modes;
        std::vector<real> eigs;
        std::vector<real> eigforce; //force upon eigs
        std::vector<real> eigx;
        std::vector<real> deigxdt;
        std::vector<real> solidLoad;

        real solidForceScale = 1;
        real solidTimeScale = 1;

        real rlAlpha = 1e-2, rlBeta = 1e-2;

        size Kset = 1;
    public:
        ut3 totalForceFace; //total conveyed force


        void SaveClass(const std::string &ClassFile,
                       std::ostream &logout = std::cout)
        {
            logout << "Saving Class To " << ClassFile << std::endl;
        }

        void Load(const std::string &SolidMeshFile, const std::string &FluidMeshFile,
                  const std::string &FluidInitVal, const std::string &InterKNNFile,
                  const std::vector<utype<3>> &SolidBndset, const std::vector<utype<5>> &FluidBndset,
                  std::ostream &logout = std::cout)
        {
            ElasSet setS(1, 0.3, ElasSet::E_nu), setS2(1, -0.5, ElasSet::E_nu);
            double GAM = 1.4;
            RoeSet setF;
            setF.gamma = GAM;

            logout << "===FSI=== Load Solid: " << std::endl;
            SolidSolver.Load(SolidMeshFile, SolidBndset, logout);
            SolidSolver.setElas(setS);
            solidLoad.resize(SolidSolver.showPoints().size() * 3);
            logout << "===FSI=== Load Fluid: " << std::endl;
            FluidSolver.Load(FluidMeshFile, FluidBndset, logout);
            FluidSolver.InitByConst(ut5(1.0));
            FluidSolver.setInter(), FluidSolver.setRoeset(setF),
                FluidSolver.fieldreadTecPlotDat(FluidInitVal, logout);
            FluidGridder.Load(FluidMeshFile, SolidBndset, logout); //SolidBndSet is a filler here
            FluidGridder.setElas(setS2);
            FluidGridder.setFluidGridder();

            //read the KNN file
            std::ifstream fin(InterKNNFile);
            if (!fin)
            {
                logout << "File Open Faliure: " << InterKNNFile << std::endl;
                exit(-1);
            }
            std::string buf1;
            fin >> buf1;
            logout << buf1 << std::endl;
            size FIntersize, kF;
            fin >> FIntersize;
            Finter.resize(FIntersize);
            F2SIndex.resize(FIntersize);
            fin >> buf1;
            logout << buf1 << std::endl;
            fin >> kF;
            for (size i = 0; i < FIntersize; i++)
            {
                fin >> Finter[i];
                F2SIndex[i].resize(kF);
                for (size k = 0; k < kF; k++)
                    fin >> F2SIndex[i][k];
            }

            fin >> buf1;
            logout << buf1 << std::endl;
            size SIntersize, kS;
            fin >> SIntersize;
            Sinter.resize(SIntersize);
            S2FIndex.resize(SIntersize);
            fin >> buf1;
            logout << buf1 << std::endl;
            fin >> kS;
            for (size i = 0; i < SIntersize; i++)
            {
                fin >> Sinter[i];
                S2FIndex[i].resize(kS);
                for (size k = 0; k < kS; k++)
                    fin >> S2FIndex[i][k];
            }
            fin.close();

            if ((!SolidSolver.checkAreInters(Sinter)) || (!FluidSolver.checkAreInters(Finter)))
            {
                logout << "===FSI found No Inter Points In KNN set=== ERROR" << std::endl;
                exit(-3);
            }
        }

        void Process(const std::string &ModeFileHead, const std::string &FluidGridModeFileHead,
                     std::ostream &logout = std::cout)
        {

            SolidSolver.BuildMKBound(logout);
            SolidSolver.SolveEigen(Kset, 1e-10, logout);
            SolidSolver.modeprintTecPlotDat(ModeFileHead, 0.0, logout);

            std::vector<ut3> curDisp(Finter.size());
            std::vector<ut3> modeByUT(SolidSolver.showPoints().size());
            modes = SolidSolver.GetModes();
            eigs = SolidSolver.GetEigs();
            eigforce.resize(eigs.size()), eigx.resize(eigs.size(), 0), deigxdt.resize(eigs.size(), 0); //0 init mode
            FluidGridMode.resize(modes.size());

            std::vector<size> interPos(FluidGridder.showPoints().size(), -2);
            for (size i = 0; i < Finter.size(); i++)
                interPos[Finter[i]] = i; //, std::cout << Finter[i] << std::endl;

            FluidGridder.BuildMKBound(logout);
            for (size k = 0; k < modes.size(); k++)
            {
                for (size i = 0; i < modeByUT.size(); i++)
                    modeByUT[i][0] = modes[k][i * 3 + 0],
                    modeByUT[i][1] = modes[k][i * 3 + 1],
                    modeByUT[i][2] = modes[k][i * 3 + 2];
                DistributeDisp(SolidSolver.showPoints(), FluidGridder.showPoints(),
                               modeByUT, Finter, F2SIndex, curDisp);
                FluidGridder.setInterDisplacement(interPos, curDisp);
                FluidGridder.SolveStatic(logout, FemT4Solver::SOR);
                //for (int i = 0; i < curDisp.size(); i++)
                //   std::cout << i << '\t' << curDisp[i][0] << '\t' << curDisp[i][1] << '\t' << curDisp[i][2] << std::endl;
                // for (int i = 0; i < FluidGridder.GetU().size(); i++)
                //     std::cout<< i <<'\t' << FluidGridder.GetU()[i] << std::endl;
                // for (int i = 0; i < FluidGridder.showPoints().size(); i++)
                //     std::cout << i <<'\t'<< FluidGridder.showPoints()[i][0] << std::endl;

                char nbuf[512];
                sprintf(nbuf, "_R_Mode_%d.dat", k);
                FluidGridder.getStrainStress(FluidGridder.GetU(), logout);
                FluidGridder.fieldprintTecPlotDat(FluidGridModeFileHead + nbuf, 0.0, FluidGridder.GetU(), logout);
                FluidGridMode[k] = (FluidGridder.GetU());
            }
        }

        void TestScale(std::ostream &logout = std::cout)
        {
            logout << "Testing Couple Scale" << std::endl;
            RoeUGsolver tempFsolver = FluidSolver;
            real rescontrol[] = {1, 1, 1, 1, 1};
            auto dtFluid = tempFsolver.StepUEuler(100, -0.1, 20, ut5(rescontrol) * 1e-3, 0.3, 1, logout);
            tempFsolver.calculateLoad(logout);
            DistributeForce(FluidSolver.showPoints(), SolidSolver.showPoints(),
                            tempFsolver.showInterForces(), Finter, F2SIndex, solidLoad);
            CalEigForces(1.0, logout);
            for (size k = 0; k < eigs.size(); k++)
                logout << std::scientific << std::setprecision(6)
                       << "Mode " << k << " Est_Force: " << eigforce[k] << " ModeMag: " << eigforce[k] / eigs[k]
                       << " MFreq: " << std::sqrt(eigs[k]) << std::endl;
        }

        real StepTimeEulerA(real dtmax, real CFL, size iitmax,
                            ut5 absresth, real release, real steer,
                            std::ostream &logout = std::cout)
        {
            real dt = FluidSolver.StepUEuler(dtmax, CFL, iitmax, absresth, release, steer, logout);
            totalForceFace = FluidSolver.calculateLoad(logout);
            DistributeForce(FluidSolver.showPoints(), SolidSolver.showPoints(),
                            FluidSolver.showInterForces(), Finter, F2SIndex, solidLoad);
            CalEigForces(solidForceScale, logout);
            real soliddt = std::abs(dt) * solidTimeScale;
            pushModeEuler(soliddt, logout);
            AssembleFluidGrid(logout);
            FluidSolver.updateGeom();
            logout << "===FSI=== Step, Primary Modes: " << std::scientific << std::setprecision(6);
            for (size i = 0; i < std::min(std::size_t(3), eigs.size()); i++)
                logout << i << " : " << eigx[i] << ", " << deigxdt[i] << "; ";
            logout << std::endl;
            return dt;
        }

        void CalEigForces(real scale, std::ostream &logout)
        {
            logout << "Mode Force: " << std::scientific << std::setprecision(6);
            for (size k = 0; k < eigs.size(); k++)
            {
                eigforce[k] = 0;
                for (size d = 0; d < solidLoad.size(); d++)
                    eigforce[k] += solidLoad[d] * modes[k][d] * scale; //using force scale
                logout << eigforce[k] << "  ";
            }
            logout << std::endl;
        }

        void AssembleFluidGrid(std::ostream &logout)
        {
            auto &targetp = FluidSolver.showPoints();
            auto &targetdpdt = FluidSolver.showdPointsdt();
            targetp = FluidGridder.showPoints(); //0 pos
            targetdpdt.assign(targetdpdt.size(), utype<3>(0.0));
            for (size k = 0; k < eigs.size(); k++)
#ifdef OMP_ON
#pragma omp parallel for
#endif
                for (size p = 0; p < targetp.size(); p++)
                {
                    targetp[p][0] += eigx[k] * FluidGridMode[k][p * 3 + 0];
                    targetp[p][1] += eigx[k] * FluidGridMode[k][p * 3 + 1];
                    targetp[p][2] += eigx[k] * FluidGridMode[k][p * 3 + 2];
                    targetdpdt[p][0] += deigxdt[k] * FluidGridMode[k][p * 3 + 0] * solidTimeScale;
                    targetdpdt[p][1] += deigxdt[k] * FluidGridMode[k][p * 3 + 1] * solidTimeScale;
                    targetdpdt[p][2] += deigxdt[k] * FluidGridMode[k][p * 3 + 2] * solidTimeScale;
                }
        }

        //For current versions of the interpolaions below, scales are not very safe,
        //currenly designed for only values around 1

        //pull to fluid disp(partial) from solid disp(total)
        static void DistributeDisp(const std::vector<ut3> &fromp, const std::vector<ut3> &top,
                                   const std::vector<ut3> &fromDisp, const std::vector<size> &toPinter,
                                   const std::vector<std::vector<size>> &KNN, std::vector<ut3> &toDispSub)
        {
            toDispSub.assign(toDispSub.size(), ut3(0.0));
            for (size i = 0; i < toDispSub.size(); i++)
            {
                real allw = 0;
                for (size k = 0; k < KNN[i].size(); k++)
                {
                    real w = UTILS::RBFInversedMultiQuadric<3>(fromp[KNN[i][k]] - top[toPinter[i]], 0);
                    toDispSub[i] += fromDisp[KNN[i][k]] * w;
                    //::cout << w << ", " << fromDisp[KNN[i][k]][0] << std::endl;
                    allw += w;
                }
                toDispSub[i] /= allw;
                //std::cout << toDispSub[i][0] << std::endl;
            }
        }

        //push to solid force(all) from fluid force(all)
        static void DistributeForce(const std::vector<ut3> &fromp, const std::vector<ut3> &top,
                                    const std::vector<ut3> &fromForce, const std::vector<size> &fromPinter,
                                    const std::vector<std::vector<size>> &KNN, std::vector<real> &toForce)
        {
            //toForce.assign(toForce.size(), ut3(0.0));
            toForce.assign(toForce.size(), 0.0);
            std::vector<real> weights;
            for (size i = 0; i < fromPinter.size(); i++)
            {
                weights.resize(KNN[i].size());
                real wtotal = 0;
                for (size k = 0; k < KNN[i].size(); k++)
                    wtotal += (weights[k] = UTILS::RBFInversedMultiQuadric<3>(fromp[fromPinter[i]] - top[KNN[i][k]], 0));
                for (size k = 0; k < KNN[i].size(); k++)
                {
                    auto inc = fromForce[fromPinter[i]] * (weights[k] / wtotal);
                    toForce[KNN[i][k] * 3 + 0] += inc[0], toForce[KNN[i][k] * 3 + 1] += inc[1], toForce[KNN[i][k] * 3 + 2] += inc[2];
                }
            }
        }

        void pushModeEuler(real dt, std::ostream &logout = std::cout)
        {
            for (size k = 0; k < eigs.size(); k++)
            {
                real c = eigs[k] * rlAlpha + rlBeta - 1.0 * dt / std::sqrt(eigs[k]) * pi * 0.125; //latter term is 2nd order fix
                real div = 1 + sqr(dt) * eigs[k] + dt * c;
                real eigxnew = (sqr(dt) * eigforce[k] + dt * deigxdt[k] + eigx[k] + c * dt * eigx[k]) / div;
                real deigxdtnew = (dt * eigforce[k] + deigxdt[k] - dt * eigs[k] * eigx[k]) / div;
                eigx[k] = eigxnew, deigxdt[k] = deigxdtnew;
            }
        }

        void setSolidScale(real fscale, real tscale)
        {
            solidForceScale = fscale, solidTimeScale = tscale;
        }

        void setAlphaBeta(real alpha, real beta)
        {
            rlAlpha = alpha, rlBeta = beta;
        }

        void statprintTecPlotDat(const std::string &FluidFileName, const std::string &SolidFileNameU,
                                 const std::string &SolidFileNamedUdt, real magnifier,
                                 std::ostream &logout = std::cout)
        {
            FluidSolver.fieldprintTecPlotDat(FluidFileName, logout);
            SolidSolver.modeSum(eigx, SolidSolver.GetU(), SolidFileNameU, magnifier, logout);
            SolidSolver.modeSum(deigxdt, SolidSolver.GetU(), SolidFileNamedUdt, 0.0, logout);
        }

        //x dxdt f x dxdt f...
        void spitEigenDofs(std::ostream &out)
        {
            out << std::scientific << std::setprecision(16);
            for (real k = 0; k < eigs.size(); k++)
                out << eigx[k] << "  " << deigxdt[k] << "  " << eigforce[k] << "  ";
        }

        RoeUGsolver &FSolver()
        {
            return FluidSolver;
        }

        void setModeOrder(size nk)
        {
            Kset = nk;
        }

        void SaveProcessResult(const std::string &processFile, std::ostream &logout = std::cout)
        {
            logout << "Process To File: " << processFile << std::endl;
            std::ofstream out(processFile);
            if (!out)
            {
                logout << "Process File not opened: " << processFile << std::endl;
                exit(-4);
            }
            out << eigs << eigx << deigxdt << modes << FluidGridMode;
        }

        void LoadProcessResult(const std::string &processFile, std::ostream &logout = std::cout)
        {
            logout << "Process From File: " << processFile << std::endl;
            std::ifstream in(processFile);
            if (!in)
            {
                logout << "Process File not opened: " << processFile << std::endl;
                exit(-6);
            }
            in >> eigs >> eigx >> deigxdt >> modes >> FluidGridMode;
            SolidSolver.GetModes() = modes;
            SolidSolver.GetEigs() = eigs;
            Kset = eigs.size();
            eigforce.resize(Kset);
        }
    };
}
