#include <vector>
#include <cstdio>
#include "fsi.hpp"

using namespace FSI;

void test();
void testWING();
void testTUBE();
void testRoe();
void testTRI12();
void runTRI();
void runTRI2();
void runTRI2F();
void testTri();
void testFemBEAM();

int main()
{
    testTRI12();
    return 0;
}

void test()
{
    //202:v2 = 36 rho2 = 1.4/20
    //206:v2 = 36 rho2 = 1.4/
    double U = 0.5, U2 = 000.1 / (1.4 / 20) / (3.14 * 0.4 * 0.4 * 0.25) / 316;
    double RHO = 1.4, RHO2 = 1.4 * 2;
    double GAM = 1.4;
    double PINF = 1, PINF2 = 200;
    double uinf[] = {RHO, 0.0, 0.0, -RHO * U, RHO * (U * U * 0.5 + PINF / (RHO * (GAM - 1)))};
    double uinf2[] = {RHO2, 0.0, 0.0, -RHO2 * U2, RHO2 * (U2 * U2 * 0.5 + PINF2 / (RHO2 * (GAM - 1)))};
    double rescontrol[] = {1e-2, 1e-2, 1e-2, 1e-2, 1e-2};
    int see = 200;
    double CFL = -2.0e-1;
    int itmax = 1000000;
    int iitmax = 20;
    // int see = 100;
    // double CFL = -2.0e1 ;
    // int itmax = 1000000;
    // int iitmax = 200;
    RoeSet set;
    set.gamma = GAM;
    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    boundset[1] = ut5(uinf), boundset[2] = ut5(uinf2);
    boundset[3] = ut5(uinf2);
    //solver.Load("MeshCAD.txt", boundset);
    solver.Load("NZL_MeshCADLRD.txt", boundset);
    //solver.Load("BLK_MeshCADC.txt", boundset);

    solver.setRoeset(set);
    solver.InitByConst(ut5(uinf));

    //std::cout << "Mesh OUTPUT\n";

    solver.fieldreadTecPlotDat("dout/TEST_211_STEP_60800.dat");

    solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");
    for (int i = 1; i <= itmax; i++)
    {
        solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5, 0.5, 1);
        std::cout << "===STEP: " << i << " DONE===" << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_212_STEP_%d.dat", i);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void testTUBE()
{
    double GAM = 1.4;
    RoeSet set;
    set.gamma = GAM;
    double vl = 0, vr = 0, rhol = 1, rhor = 1, PL = 1, pl = 10, pr = 1;
    double uL[] = {rhol, vl, 0, 0, pl};
    double uR[] = {rhor, vr, 0, 0, pr};
    ut5 UL = UPhysToUconserv(ut5(uL), set), UR = UPhysToUconserv(ut5(uR), set);
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 100;
    int itmax = 1000;
    int iitmax = 20;
    double CFL = -1e-1;

    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    //1:inf1 2:inter 3: inf2
    boundset[1] = ut5(UL), boundset[2] = ut5(UL);
    boundset[3] = ut5(UR);
    solver.Load("Meshes/MeshCAD_TUBE_T1.txt", boundset);
    solver.setRoeset(set);
    solver.InitByConst(ut5(UR));
    solver.AlterByBox(UL, -1, 0.5);

    int istart = 0;
    double t = 0;
    solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");
    for (int i = 1; i <= itmax; i++)
    {
        //solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5e-2, 0.5, 1);
        auto dt = solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 0.0);
        t += dt;
        std::cout << "===STEP: " << i << " DONE=== current t = " << std::scientific << std::setw(12) << t << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_TUBE_STEP_A_%d.dat", i + istart);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void testWING()
{
    //9: MA.8 alpha = 0 SW
    //7: MA.8 alpha = 5 SW
    //11: MA.9 alpha = 5 SW_M
    //22: MA2alpha = 0 SW_M
    double alphaWING = 0.0 / 180.0 * acos(-1);
    double U = 2, U2 = 1;
    double RHO = 1.4, RHO2 = 1.4;
    double GAM = 1.4;
    double PINF = 1, PINF2 = 1;
    double uinf[] = {RHO, 0.0, RHO * sin(alphaWING) * U, -RHO * cos(alphaWING) * U, RHO * (U * U * 0.5 + PINF / (RHO * (GAM - 1)))};
    double uinf2[] = {RHO2, 0.0, 0.0, -RHO2 * U2, RHO2 * (U2 * U2 * 0.5 + PINF2 / (RHO2 * (GAM - 1)))};
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 400;
    double CFL = -1e-1;
    int itmax = 2000;
    int iitmax = 20;
    // int see = 100;
    // double CFL = -2.0e1 ;
    // int itmax = 1000000;
    // int iitmax = 200;
    RoeSet set;
    set.gamma = GAM;
    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    boundset[1] = ut5(uinf), boundset[2] = ut5(uinf2);
    boundset[3] = ut5(uinf2);
    solver.Load("MeshCAD_SW_M.txt", boundset);

    solver.setRoeset(set);
    solver.InitByConst(ut5(uinf2));

    int istart = 2200;
    double t = 0;
    solver.fieldreadTecPlotDat("dout/TEST_210100_STEP_A_2200.dat");

    //solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");

    for (int i = 1; i <= itmax; i++)
    {
        //solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5e-2, 0.5, 1);
        auto dt = solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 1);
        t += dt;
        std::cout << "===STEP: " << i << " DONE=== current t = " << std::setprecision(5) << std::ios::scientific << t << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_210100_STEP_A_%d.dat", i + istart);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void testTRI12()
{
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 1000;
    int itmax = 100000;
    int iitmax = 100;
    double CFL = -1e-1;
    double dtmax = 1e-1, release = 0.5, steer = 1;

    UGFSISolver solver;
    solver.setModeOrder(10);
    double alpha = 5.0 / 180.0 * acos(-1);

    double GAM = 1.4;
    RoeSet set;
    set.gamma = GAM;
    double vl = 2, vr = 0, rhol = 1.4, rhor = 1, pl = 1, pr = 1;
    double uL[] = {rhol, 0, vl * sin(alpha), -vl * cos(alpha), pl};
    double uR[] = {rhor, 0, vr * sin(alpha), -vr * cos(alpha), pr};
    ut5 UL = UPhysToUconserv(ut5(uL), set), UR = UPhysToUconserv(ut5(uR), set);
    std::vector<ut5> boundset(4);
    //1:inf1 2:inter 3: inf2
    boundset[1] = ut5(UL), boundset[2] = ut5(UL);
    boundset[3] = ut5(UR);
    std::vector<utype<3>> bset(4);
    bset[0].set0(), bset[1].set0(), bset[2][0] = bset[2][2] = 0.0, bset[2][1] = 1.0, bset[3].set0();

    solver.Load("Meshes/MeshCSD_Tri_T1.txt", "Meshes/MeshCAD_TRI_K2.txt",
                "dout/TEST_TRI2_A5_STEP_A_2000.dat", "Meshes/FSIInter_TRI_K12.txt",
                bset, boundset);
    //solver.Process("dout_structural/MeshCSD_Tri_T1_Out", "dout_structural/MeshCAD_TRI_K2_Out");
    //solver.SaveProcessResult("dout_mode/TRI_K12_Process_10.txt");
    solver.LoadProcessResult("dout_mode/TRI_K12_Process_10.txt");

    //solver.TestScale();
    solver.setSolidScale(1.0 / 20e3, 1000.0);
    solver.setAlphaBeta(0.001, 0.0);

    std::string SolidUResultFile("dout_structural/MeshCSD_Tri_T1_UtOut");
    std::string SolidVResultFile("dout_structural/MeshCSD_Tri_T1_VtOut");
    std::string FluidResultFile("dout/MeshCAD_TRI_K2_tOut");

    std::ofstream ForceTimeSeqOut("dout_mode/TRI_K12_Force_tOut.txt");
    std::ofstream ModeTimeSeqOut("dout_mode/TRI_K12_Mode_tOut.txt");

    int istart = 0;
    double t = 0;

    for (int iter = 1; iter <= itmax; iter++)
    {

        ut5 absresth = ut5(5e-2);
        real dt;
        dt = solver.StepTimeEulerA(1e-1, CFL, iitmax, ut5(rescontrol) * 4e-3, 0.5, 1);
        t += dt;
        ModeTimeSeqOut << std::setw(10) << std::abs(t) << std::setw(20) << "  ";
        solver.spitEigenDofs(ModeTimeSeqOut);
        ModeTimeSeqOut << std::endl;
        ForceTimeSeqOut << std::abs(t) << solver.totalForceFace << std::endl;
        std::cout << "===FSI JOB=== Step Finished, iter = " << iter << std::endl;
        if (iter % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "_STEP_%d.dat", iter + istart);
            solver.statprintTecPlotDat(FluidResultFile + nbuf,
                                       SolidUResultFile + nbuf,
                                       SolidVResultFile + nbuf,
                                       0.0);
        }
    }
    ModeTimeSeqOut.close();
    ForceTimeSeqOut.close();
}

void runTRI()
{
    double alpha = 0 / 180 * acos(-1);

    double GAM = 1.4;
    RoeSet set;
    set.gamma = GAM;
    double vl = 2, vr = 0, rhol = 1.4, rhor = 1, pl = 1, pr = 1;
    double uL[] = {rhol, 0, vl * sin(alpha), -vl * cos(alpha), pl};
    double uR[] = {rhor, 0, vr * sin(alpha), -vr * cos(alpha), pr};
    ut5 UL = UPhysToUconserv(ut5(uL), set), UR = UPhysToUconserv(ut5(uR), set);
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 1000;
    int itmax = 10000;
    int iitmax = 100;
    double CFL = -2e-1;

    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    //1:inf1 2:inter 3: inf2
    boundset[1] = ut5(UL), boundset[2] = ut5(UL);
    boundset[3] = ut5(UR);
    solver.Load("Meshes/MeshCAD_TRI_K1.txt", boundset);
    solver.setRoeset(set);
    solver.InitByConst(ut5(UL));
    //solver.AlterByBox(UL, -1, 0.5);

    int istart = 3800;
    double t = -1.121693e+00;
    solver.fieldreadTecPlotDat("dout/TEST_TRI_STEP_A_3800.dat");
    //solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");
    for (int i = 1; i <= itmax; i++)
    {
        //solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5e-2, 0.5, 1);
        auto dt = solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 0.0);
        t += dt;
        std::cout << "===STEP: " << i + istart << " DONE=== current t = " << std::scientific << std::setw(12) << t << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_TRI_STEP_A_%d.dat", i + istart);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void runTRI2()
{
    double alpha = 5.0 / 180.0 * acos(-1);

    double GAM = 1.4;
    RoeSet set;
    set.gamma = GAM;
    double vl = 2, vr = 0, rhol = 1.4, rhor = 1, pl = 1, pr = 1;
    double uL[] = {rhol, 0, vl * sin(alpha), -vl * cos(alpha), pl};
    double uR[] = {rhor, 0, vr * sin(alpha), -vr * cos(alpha), pr};
    ut5 UL = UPhysToUconserv(ut5(uL), set), UR = UPhysToUconserv(ut5(uR), set);
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 1000;
    int itmax = 10000;
    int iitmax = 100;
    double CFL = -2e-1;

    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    //1:inf1 2:inter 3: inf2
    boundset[1] = ut5(UL), boundset[2] = ut5(UL);
    boundset[3] = ut5(UR);
    solver.Load("Meshes/MeshCAD_TRI_K2.txt", boundset);
    solver.setRoeset(set);
    solver.InitByConst(ut5(UL));
    //solver.AlterByBox(UL, -1, 0.5);

    // int istart = 1800;
    // double t = -5.313698e-01;
    // solver.fieldreadTecPlotDat("dout/TEST_TRI_STEP_A_1800.dat");
    double t = 0;
    int istart = 0;

    //solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");
    for (int i = 1; i <= itmax; i++)
    {
        //solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5e-2, 0.5, 1);
        solver.StepUEulerSteady(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 1);
        double dt = 0;
        //auto dt = solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 0.0);
        t += dt;
        std::cout << "===STEP: " << i + istart << " DONE=== current t = " << std::scientific << std::setw(12) << t << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_TRI2_A5_STEP_A_%d.dat", i + istart);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void runTRI2F()
{
    double alpha = 5.0 / 180.0 * acos(-1);

    double GAM = 1.4;
    RoeSet set;
    set.gamma = GAM;
    double vl = 2, vr = 0, rhol = 1.4, rhor = 1, pl = 1, pr = 1;
    double uL[] = {rhol, 0, vl * sin(alpha), -vl * cos(alpha), pl};
    double uR[] = {rhor, 0, vr * sin(alpha), -vr * cos(alpha), pr};
    ut5 UL = UPhysToUconserv(ut5(uL), set), UR = UPhysToUconserv(ut5(uR), set);
    double rescontrol[] = {1, 1, 1, 1, 1};
    int see = 1000;
    int itmax = 10000;
    int iitmax = 100;
    double CFL = -2e-1;

    RoeUGsolver solver;
    std::vector<ut5> boundset(4);
    //1:inf1 2:inter 3: inf2
    boundset[1] = ut5(UL), boundset[2] = ut5(UL);
    boundset[3] = ut5(UR);
    solver.Load("Meshes/MeshCAD_TRI_K2.txt", boundset);
    solver.setRoeset(set);
    solver.InitByConst(ut5(UL));
    //solver.AlterByBox(UL, -1, 0.5);

    int istart = 2000;
    double t = 0;
    solver.fieldreadTecPlotDat("dout/TEST_TRI2_A5_STEP_A_2000.dat");
    solver.setInter();

    //solver.fieldprintTecPlotDatDebug("./MeshCADtecdebug.dat");
    for (int i = 1; i <= itmax; i++)
    {
        solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 5e-2, 0.5, 1);

        // solver.StepUEulerSteady(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 1);
        solver.calculateLoad();
        double dt = 0;
        //auto dt = solver.StepUEuler(1e-1, CFL, iitmax, ut(rescontrol) * 1e-2, 0.5, 0.0);
        t += dt;
        std::cout << "===STEP: " << i + istart << " DONE=== current t = " << std::scientific << std::setw(12) << t << std::endl;
        if (i % see == 0)
        {
            char nbuf[512];
            sprintf(nbuf, "dout/TEST_TRI2_A5_STEP_A_%d.dat", i + istart);
            solver.fieldprintTecPlotDat(nbuf);
        }
    }
}

void testFemBEAM()
{
    ElasSet set(1, 0.3, ElasSet::E_nu);
    FemT4Solver solver;
    std::vector<utype<3>> bset(4);
    bset[0].set0(), bset[1].set0(), bset[2][0] = bset[2][2] = 0.0, bset[2][1] = 1.0, bset[3].set0();
    solver.Load("Meshes/MeshCSD_BEAM_T1.txt", bset);
    solver.setElas(set);
    solver.BuildMKBound();
    solver.SolveStatic();
    solver.getStrainStress(solver.GetU());
    solver.fieldprintTecPlotDat("dout_structural/MeshCSD_Out.dat", 0.0, solver.GetU());
    return;
}

void testTri()
{
    ElasSet set(1, 0.3, ElasSet::E_nu);
    FemT4Solver solver;
    std::vector<utype<3>> bset(4);
    bset[0].set0(), bset[1].set0(), bset[2][0] = bset[2][2] = 0.0, bset[2][1] = 1.0, bset[3].set0();
    solver.Load("Meshes/MeshCSD_Tri_T1.txt", bset);
    solver.setElas(set);
    solver.BuildMKBound();
    solver.SolveEigen(10);
    solver.modeprintTecPlotDat("dout_structural/MeshCSD_Out", 0.0);
    return;
}