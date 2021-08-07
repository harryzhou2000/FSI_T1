#include "roesgsolver.hpp"

using namespace StructuralGrid;
void coupleMain();
void coupleMainSteady();
void coupleMainSCMesh();
void shockMain();
void shockMain2();
void shockMainCol();
void stepMain();

int main()
{
    //coupleMainSCMesh();
    //coupleMain();
    coupleMainSteady();
    //shockMain();
    //shockMain2();
    //shockMainCol();
    //stepMain();
    return 0;
}

void coupleMain() //0012
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 1.25;
    double vB = 0.8 * sqrt(1.4);
    rhoInf = 1.0;
    std::string casename = "SG_0012_T2_A125_M08_";

    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.5;
    double fixRatio = 0.0001;
    int see = 5000;
    double ttvd = 0.0001; //20
    double tend = 400;
    bool vis = false, mov = false;
    int m = 400, n = 200;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/Mesh0012_IVA1.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);
    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();
    int istart = 0;

    // //READ breakpoint
    // std::ifstream fin("dout/" + casename + "data_0_AT14000.dat");
    // if (!fin)
    //     std::cerr << "Infile bad" << std::endl;
    // grid1.fieldintTecPlotDat(fin);
    // fin.close();
    // istart = 14000;
    // //READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        auto RK4ret = grid1.StepRK4(tvd, set, cfl, 100.0, tend - t, conclude);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void coupleMainSteady() //0012
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 1.25;
    double vB = 0.8 * sqrt(1.4);
    rhoInf = 1.0;
    std::string casename = "SG_0012_T2_A125_M08_Steady_";

    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.3;
    double fixRatio = 0.0001;
    int see = 5000;
    double ttvd = 0.0001; //20
    double tend = 400;
    bool vis = false, mov = false;
    int m = 400, n = 200;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/Mesh0012_IVA1.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);
    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();
    int istart = 0;

    // //READ breakpoint
    // std::ifstream fin("dout/" + casename + "data_0_AT14000.dat");
    // if (!fin)
    //     std::cerr << "Infile bad" << std::endl;
    // grid1.fieldintTecPlotDat(fin);
    // fin.close();
    // istart = 14000;
    // //READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        auto RK4ret = grid1.StepRK4Steady(tvd, set, cfl, 100.0, tend - t, conclude);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void coupleMainSCMesh()
{
    double rhoInf = 1.4;
    double pInf = 1;
    double alpha = pi / 180.0 * 3.5;
    double vB = 0.725;
    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 0.5;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.2;
    double fixRatio = 0.0001;
    int see = 1000;
    double ttvd = 0.00001; //20
    double tend = 200;
    bool vis = false, mov = false;
    std::string casename = "SG_20714_T2_Steady";
    int m = 600, n = 200;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/Mesh20714_IVA1_600x200.txt");

    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);
    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    std::vector<utype<4>> increment0;
    int istart = 0;

    // std::ifstream fin("dout/" + casename + "data_0_AT90000.dat");
    // if (!fin)
    //     std::cerr << "Infile bad" << std::endl;
    // grid1.fieldintTecPlotDat(fin);
    // fin.close();
    // istart = 90000;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        auto RK4ret = grid1.StepRK4Steady(tvd, set, cfl, 100.0, tend - t, conclude);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i + istart, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i + istart);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void shockMain()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 2.05;
    double vB = 2;
    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1.0;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.1;
    double fixRatio = 0.0001;
    int see = 5000;
    double ttvd = 3; //20
    double tend = 200;
    bool vis = false, mov = false;
    std::string casename = "SG_0012_T1_SHOCK";
    int m = 400, n = 200;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/Mesh0012_IVA1.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);
    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);
        double tvdo = tvd;
        double vmax;
        int aim;
        int throttle = 0;
        do
        {
            vmax = 0;
            aim = 0;
            if (tvd != tvdo)
                printf("\twarning:::TVD THROTTLING\n");
            grid1.GetIncrement(increment0, tvd, set, vmax, aim);
            tvd *= tvdDecrease;
            throttle++;
            if (throttle >= 128)
            {
                printf("BAD TERMINATION: TVD THROTTLE TOP\n");
                exit(0);
            }
            throttle++;
        } while (aim && false);
        double dt = cfl / vmax;
        bool conclude = false;
        if (t + dt >= tend)
        {
            dt = tend - t;
            conclude = true;
        }
        increment0 *= dt;
        t += dt;
        grid1.dataIncrement(increment0);
        auto incmax = maxabs<4>(increment0);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i, tvd / tvdDecrease, dt, t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void shockMain2()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 90;
    double vB = 0.775;
    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1.0;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.2;
    double fixRatio = 0.0001;
    int see = 5000;
    double ttvd = 0.0001; //20
    double tend = 200;
    bool vis = false, mov = false;
    std::string casename = "SG_0012_T1_SHOCK90";
    int m = 400, n = 200;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/Mesh0012_IVA1.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);
    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);
        double tvdo = tvd;
        double vmax;
        int aim;
        int throttle = 0;
        do
        {
            vmax = 0;
            aim = 0;
            if (tvd != tvdo)
                printf("\twarning:::TVD THROTTLING\n");
            grid1.GetIncrement(increment0, tvd, set, vmax, aim);
            tvd *= tvdDecrease;
            throttle++;
            if (throttle >= 128)
            {
                printf("BAD TERMINATION: TVD THROTTLE TOP\n");
                exit(0);
            }
            throttle++;
        } while (aim && false);
        double dt = cfl / vmax;
        bool conclude = false;
        if (t + dt >= tend)
        {
            dt = tend - t;
            conclude = true;
        }
        increment0 *= dt;
        t += dt;
        grid1.dataIncrement(increment0);
        auto incmax = maxabs<4>(increment0);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i, tvd / tvdDecrease, dt, t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void shockMainColOld()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 0;
    double vB = 0.8;
    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1.0;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.2;
    double fixRatio = 0.0001;
    int see = 100;
    double ttvd = 0.0001; //20
    double tend = 200;
    bool vis = false, mov = false;
    std::string casename = "SG_Col_T1_SHOCK90_RKTest";
    int m = 250, n = 250;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/MeshCol_250x250.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);

    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    int istart = 0;
    // READ breakpoint
    // std::ifstream fin("dout/" + casename + "data_0_AT44000.dat");
    // if (!fin)
    //     std::cerr << "Infile bad" << std::endl;
    // grid1.fieldintTecPlotDat(fin);
    // fin.close();
    // istart = 44000;
    // READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);
        double tvdo = tvd;
        double vmax;
        int aim;
        int throttle = 0;
        do
        {
            vmax = 0;
            aim = 0;
            if (tvd != tvdo)
                printf("\twarning:::TVD THROTTLING\n");
            grid1.GetIncrement(increment0, tvd, set, vmax, aim);
            tvd *= tvdDecrease;
            throttle++;
            if (throttle >= 128)
            {
                printf("BAD TERMINATION: TVD THROTTLE TOP\n");
                exit(0);
            }
            throttle++;
        } while (aim && false);
        double dt = cfl / vmax;
        bool conclude = false;
        if (t + dt >= tend)
        {
            dt = tend - t;
            conclude = true;
        }
        increment0 *= dt;
        t += dt;
        grid1.dataIncrement(increment0);
        auto incmax = maxabs<4>(increment0);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd / tvdDecrease, dt, t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i + istart, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i + istart);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void shockMainCol()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 0;
    double vB = 3.0;
    int itmax = 1000000;
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1.0;
    double tvdDecrease = 0.5;
    double tstart = 0;
    double cfl = 0.2;
    double fixRatio = 0.0001;
    int see = 100;
    double ttvd = 0.0001; //20
    double tend = 2000;
    bool vis = false, mov = false;
    std::string casename = "SG_Col_T1_SHOCK90_M3Test";
    int m = 250, n = 250;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;

    ///bound
    std::vector<ntype> bound(4);
    bound[0] = cross, bound[1] = infinity, bound[2] = cross, bound[3] = wall;
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k] = m * k + m - 1;
        crossEdge[k] = 2;
    }
    for (size_t k = 0; k < n; k++)
    {
        crossLoc[k + n + m] = m * k;
        crossEdge[k + n + m] = 0;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    std::ifstream MeshIN("SGMeshes/MeshCol_250x250.txt");
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);

    grid1.SetInitial(but);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    int istart = 0;
    //READ breakpoint
    std::ifstream fin("dout/" + casename + "data_0_AT14000.dat");
    if (!fin)
        std::cerr << "Infile bad" << std::endl;
    grid1.fieldintTecPlotDat(fin);
    fin.close();
    istart = 14000;
    //READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        auto RK4ret = grid1.StepRK4(tvd, set, cfl, 100.0, tend - t, conclude);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        //grid1.dataFix(fixRatio, set);
        //grid1.xietaLaplacianFix(0.05);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i + istart, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i + istart);
            // std::ofstream fout("dout/" + casename + "data_0_AT" + std::string(stepname) + ".txt");
            // if (!fout)
            // {
            //     printf("Ofile faliure!!!!\n");
            //     exit(-1);
            // }
            // grid1.fieldout(fout);
            // fout.close();

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDat(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void stepMainSmall()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 0;
    double vB = 3.0, vin = 1.0;
    int itmax = 1000000;
    int iitmax = 20;
    double absth[] = {1, 1, 1, 1, 1};
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1;
    double tvdDecrease = 0.2;
    double tstart = 0;
    double cfl = 0.2;
    set.thresEntropyFix = 0.2;
    double fixRatio = 0.0001;
    int see = 1000;
    double ttvd = 0.00001; //20
    double tend = 2000;
    bool vis = false, mov = false;
    std::string casename = "SG_Step_B";
    int m = 100, n = 200, nstep = 50;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    double eIn = pInf / (set.gamma - 1) + 0.5 * vin * vin * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;
    utype<4> bin;
    bin[0] = rhoInf, bin[1] = cos(alpha) * rhoInf * vin, bin[2] = sin(alpha) * rhoInf * vin, bin[3] = eIn;
    ///bound
    std::vector<ntype> bound(peri, wall);
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < m; k++)
    {
        bound[k + n] = infinity;
    }
    for (size_t k = 0; k < nstep; k++)
    {
        bound[k + m + n] = infinity;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    // std::ifstream MeshIN("SGMeshes/MeshStep_200x400x100.txt");//T1
    std::ifstream MeshIN("SGMeshes/MeshStep_100x200x50.txt"); //T2
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);

    grid1.SetInitial(bin);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    int istart = 0;
    // //READ breakpoint
    std::ifstream fin("dout/" + casename + "data_0_AT5000.dat");
    if (!fin)
        std::cerr << "Infile bad" << std::endl;
    grid1.fieldintTecPlotDat(fin);
    fin.close();
    istart = 5000;
    // //READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        //auto RK4ret = grid1.StepRK4(tvd, set, cfl, 100.0, tend - t, conclude);
        auto RK4ret = grid1.StepEuler(tvd, set, cfl, 100.0, tend - t, conclude, utype<4>(absth) * 1e-2, 10, 0.5);
        //grid1.dataFix(0.0, set);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i + istart, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i + istart);

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDatDebug(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}

void stepMain()
{
    double rhoInf = 1.4;
    double pInf = 1;
    //double alpha = pi / 180.0 * 2.05;
    double alpha = pi / 180.0 * 0;
    double vB = 3.0, vin = 3.0;
    int itmax = 1000000;
    int iitmax = 20;
    double absth[] = {1, 1, 1, 1, 1};
    RoeSet set;
    set.gamma = 1.4;
    double tvdOverride = 1.0;
    double tvdDecrease = 0.2;
    double tstart = 0;
    double cfl = 0.2;
    set.thresEntropyFix = 0.2;
    double fixRatio = 0.0001;
    int see = 1000;
    double ttvd = 0.00001; //20
    double tend = 2000;
    bool vis = false, mov = false;
    std::string casename = "SG_Step_C";
    int m = 200, n = 550, nstep = 50;
    int peri = n * 2 + m * 2;
    double eInf = pInf / (set.gamma - 1) + 0.5 * vB * vB * rhoInf;
    double eIn = pInf / (set.gamma - 1) + 0.5 * vin * vin * rhoInf;
    utype<4> but;
    but[0] = rhoInf, but[1] = cos(alpha) * rhoInf * vB, but[2] = sin(alpha) * rhoInf * vB, but[3] = eInf;
    utype<4> bin;
    bin[0] = rhoInf, bin[1] = cos(alpha) * rhoInf * vin, bin[2] = sin(alpha) * rhoInf * vin, bin[3] = eIn;
    ///bound
    std::vector<ntype> bound(peri, wall);
    std::vector<utype<4>> infinityData(peri, but);
    std::vector<SGGrid *> crossTarget(peri, nullptr);
    std::vector<size> crossLoc(peri, -1);
    std::vector<size> crossEdge(peri, -1);
    for (size_t k = 0; k < m; k++)
    {
        bound[k + n] = infinity;
    }
    for (size_t k = 0; k < nstep; k++)
    {
        bound[k + m + n] = infinity;
    }

    //zero construct
    SGGrid grid1(vis, mov);
    //grid1.build(m, n, x, y, bound, crossTarget, crossLoc, crossEdge, infinityData);
    //std::ifstream MeshIN("SGMeshes/MeshStep_200x400x100.txt"); //A
    //std::ifstream MeshIN("SGMeshes/MeshStep_100x200x50.txt"); //B
    std::ifstream MeshIN("SGMeshes/MeshStep_200x550x50.txt"); //C
    grid1.buildByFile(m, n, MeshIN, bound, crossTarget, crossLoc, crossEdge, infinityData);

    grid1.SetInitial(bin);
    grid1.SetDistance();

    //write zero
    std::ofstream fout("./dout/" + casename + "data_00.dat");
    grid1.fieldoutTecPlotDat(fout);
    fout.close();

    int istart = 0;
    // //READ breakpoint
    std::ifstream fin("dout/" + casename + "data_0_AT24000.dat");
    if (!fin)
        std::cerr << "Infile bad" << std::endl;
    grid1.fieldintTecPlotDat(fin);
    fin.close();
    istart = 24000;
    // //READ breakpoint

    std::vector<utype<4>> increment0;

    double t = tstart;

    for (int i = 1; i <= itmax; i++)
    {
        clock_t s = clock();
        double tvd = tvdOverride * std::min(t / ttvd, 1.);

        bool conclude = false;
        //auto RK4ret = grid1.StepRK4(tvd, set, cfl, 100.0, tend - t, conclude);
        auto RK4ret = grid1.StepEuler(tvd, set, cfl, 100.0, tend - t, conclude, utype<4>(absth) * 1e-2, 10, 0.5);
        //grid1.dataFix(0.0, set);
        t += std::get<0>(RK4ret);
        auto incmax = std::get<1>(RK4ret);
        printf("Iter = %d, tvd = %e, dt = %e, t = %e cpu time = %g",
               i + istart, tvd, std::get<0>(RK4ret), t,
               (double)(clock() - s) / CLOCKS_PER_SEC);
        printf("\t incs: %10e %10e %10e %10e\n", incmax[0], incmax[1], incmax[2], incmax[3]);

        fflush(stdout);

        if (i % see == 0 || conclude)
        {
            std::cout << "Data Out" << std::endl;
            char buf[1];
            std::cin.clear();
            //std::cin.getline(buf, 1);

            char stepname[128];
            sprintf(stepname, "%08d_t_%12e", i + istart, t);
            char stepnameC[128];
            sprintf(stepnameC, "%d", i + istart);

            std::ofstream fout1("dout/" + casename + "data_0_AT" + std::string(stepnameC) + ".dat");
            if (!fout1)
            {
                printf("Ofile faliure!!!!\n");
                exit(-1);
            }
            grid1.fieldoutTecPlotDatDebug(fout1);
            fout1.close();
        }
        if (conclude)
        {
            printf("Iter = %d, Concluded\n", i);
            break;
        }
    }
}