#include <vector>
#include <cstdio>
#include "fsi.hpp"

using namespace FSI;

void testTRI12();


int main()
{
    testTRI12();
    return 0;
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




