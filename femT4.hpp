#include "roesg_common.hpp"
#include "elasticity.hpp"
#include "sparseMat.hpp"
#include <stdio.h>

namespace FEM
{
    using namespace MeshBasic;
    using namespace SolidMaterial;
    using namespace SparseMat;
    class FemT4Solver : private TetraMesh<3>
    {
        typedef Eigen::Matrix<real, 6, 1> ev6;
        std::vector<real> loadVec;
        std::vector<real> U;
        std::vector<real> dUdt;             //for intetgal transient
        std::vector<std::vector<real>> phi; //modes
        std::vector<real> lams;             //mode stiffness (or omega^2)
        SparseMatGeR<real> M;               //mass
        SparseMatGeR<real> K;               //stiffness
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>> *pKsolver = nullptr;
        Eigen::SparseMatrix<real> EK;
        std::vector<ev6> strainElem;

        std::vector<size> nodeContributeNum;
        std::vector<ev6> strainNode;
        std::vector<ev6> stressNode;
        std::vector<real> VMNode;

        ElasSet set;
        Eigen::Matrix<real, 6, 6> Dmat; //for linear constitution
        size tDOF;                      //total dof including fixed

        std::vector<size> bInterPos;
        std::vector<utype<3>> bInterData;
        std::vector<int> fixedDOF;

        bool fluidRegion = false;

    public:
        void SaveClass(const std::string &ClassFile,
                       std::ostream &logout = std::cout) //TODO ...
        {
        }

        using TetraMesh::checkAreInters;
        FemT4Solver() {}
        ~FemT4Solver()
        {
            if (pKsolver)
                delete pKsolver;
        }

        void Load(const std::string &filename, const std::vector<utype<3>> &nboundset, std::ostream &logout = std::cout)
        {
            LoadMesh(filename, logout);
            boundset = nboundset;
            tDOF = points.size() * 3;
        }

        void BuildMKBound(std::ostream &logout = std::cout)
        {
            logout << "Building Stiffness And M" << std::endl;
            M.reserve(tets.size(), 12);
            K.reserve(tets.size(), 144); // 48 = total num of NZ elems of elementary Matrix
            loadVec.assign(tDOF, 0);

            Dmat = set.constitutMatEFree();

            //set dN_i/dxii_i, where N1 = xii1, N2 = xii2, N3 = xii3, N4 = 1 - xii1 - xii2 - xii3
            Eigen::Matrix<real, 3, 4>
                dNjdxiii = Eigen::Matrix<real, 3, 4>::Zero();
            dNjdxiii(0, 0) = 1, dNjdxiii(1, 1) = 1, dNjdxiii(2, 2) = 1, dNjdxiii(0, 3) = dNjdxiii(1, 3) = dNjdxiii(2, 3) = -1;
            //only one gaussian point(at centre)

            Eigen::Matrix<real, 3, 4> dNjdxi;   // cauculation: J^-1 * dNjdxiii
            Eigen::Matrix<real, 3, 3> dxjdxiii; //'Jacobian'
            //Eigen::Matrix<real, 3, 3> dxiijdxi; // inv J
            Eigen::Matrix<real, 6, 12> Bmat; // epsilon = B * a, a = [x1y1z1x2y2z2...]^T, epsilon = [e11 e22 e33 e23 e31 e12]^T
            Eigen::Matrix<real, 6, 3> BmatSlice;
            Eigen::Matrix<real, 4, 3> xpi; //x_pi, p is node idx, i is dim_idx
            Eigen::Matrix<real, 12, 12> Kelem;
            Eigen::Matrix<real, 12, 12> Melem;
            BmatSlice.setZero(); //only need setting zero once
            for (size it = 0; it < tets.size(); it++)
            {
                size ps[4] = {tets[it].ivert[0], tets[it].ivert[1], tets[it].ivert[2], tets[it].ivert[3]};
                for (int ii = 0; ii < 4; ii++)
                    for (int jj = 0; jj < 3; jj++)
                        xpi(ii, jj) = points[ps[ii]][jj];
                dxjdxiii = dNjdxiii * xpi;
                dNjdxi = dxjdxiii.colPivHouseholderQr().solve(dNjdxiii); // J^-1 *
                for (int iv = 0; iv < 4; iv++)
                {
                    //BmatSlice.setZero();
                    BmatSlice(0, 0) = dNjdxi(0, iv), BmatSlice(1, 1) = dNjdxi(1, iv), BmatSlice(2, 2) = dNjdxi(2, iv);
                    BmatSlice(3, 1) = 0.5 * dNjdxi(2, iv), BmatSlice(3, 2) = 0.5 * dNjdxi(1, iv);
                    BmatSlice(4, 2) = 0.5 * dNjdxi(0, iv), BmatSlice(4, 0) = 0.5 * dNjdxi(2, iv);
                    BmatSlice(5, 0) = 0.5 * dNjdxi(1, iv), BmatSlice(5, 1) = 0.5 * dNjdxi(0, iv); // geometry relations here
                    Bmat.block(0, 3 * iv, 6, 3) = BmatSlice;
                }
                Kelem = (Bmat.transpose() * Dmat * Bmat) * tets[it].vol; //integration at only one point
                if (fluidRegion)
                    Kelem *= 1.0 / std::pow(tets[it].vol / MinVol, 2.0) * std::exp(-0.0 * (tets[it].vol - MinVol) / (MaxVol - MinVol));
                Melem = Eigen::Matrix<real, 12, 12>::Identity() * tets[it].vol * 0.25; //the M is density-free
#ifdef _DEBUG
                logout << Bmat << std::endl
                       << std::endl
                       << Kelem << std::endl;
#endif
                for (int ii = 0; ii < 12; ii++)
                    for (int jj = 0; jj < 12; jj++)
                    {
                        size iworld = localDOF2indDOF(ii, ps);
                        size jworld = localDOF2indDOF(jj, ps);
                        if (std::abs(Kelem(ii, jj)) > DBL_MIN)
                            K.addTo(iworld, jworld, Kelem(ii, jj));
                        if (std::abs(Melem(ii, jj)) > DBL_MIN)
                            M.addTo(iworld, jworld, Melem(ii, jj));
                    }

                // face forces
                for (int ifce = 0; ifce < 4; ifce++)
                {
                    if (neighbours[it].t[ifce] == pressuresurf)
                    {
                        utype<3> loadV = boundset[neighbours[it].n[ifce]];
                        loadV *= tets[it].farea[ifce] * (1.0 / 3.0);
                        for (int ifcend = 0; ifcend < 3; ifcend++) //for every face node add force (equally distributed for linear elems)
                        {
                            loadVec[tets[it].ivert[TetraNodes::fce2pidx[ifce][ifcend]] * 3 + 0] += loadV[0];
                            loadVec[tets[it].ivert[TetraNodes::fce2pidx[ifce][ifcend]] * 3 + 1] += loadV[1];
                            loadVec[tets[it].ivert[TetraNodes::fce2pidx[ifce][ifcend]] * 3 + 2] += loadV[2];
                        }
                    }
                }
            }
#ifdef _DEBUG2
            real sum = 0;
            for (int i = 0; i < loadVec.size(); i++)
                sum += loadVec[i];
            logout << std::setprecision(10) << sum << std::endl
                   << std::endl;
#endif

            K.reduce(), M.reduce();
            K.sort();

            //deal with the displacements, large real method
            std::vector<int> diagNZ(tDOF, 0); //if diag is positive
            for (int inz = 0; inz < K.seeNZ(); inz++)
                if (K.triplet(inz).i == K.triplet(inz).j) // is maj
                {
                    size idof = K.triplet(inz).i;
                    auto local = indDof2indNode(idof);
                    if ((K.triplet(inz).v) >= DBL_MIN)
                        diagNZ[K.triplet(inz).i] = 1;
                    if (pointntype[local.first].first == fixed)
                    {
                        utype<3> displacement = boundset[pointntype[local.first].second];
                        K.triplet(inz).v = DBL_LARGE, loadVec[idof] = DBL_LARGE * displacement[local.second];
                    }
                }
            if (fluidRegion)
            {
                logout << "Bound Set For Fluid" << std::endl;
                for (int inz = 0; inz < K.seeNZ(); inz++)
                    if (K.triplet(inz).i == K.triplet(inz).j) // is maj
                    {
                        size idof = K.triplet(inz).i;
                        auto local = indDof2indNode(idof);
                        if (pointntype[local.first].first != inner)
                            K.triplet(inz).v = DBL_LARGE, loadVec[idof] = 0; // all set to 0 for fluid boundary
                    }
            }

            // check the diags
            for (int idof = 0; idof < tDOF; idof++)
                if (!diagNZ[idof])
                {
                    logout << "+++ Stiff Mat Got ZERO or NEGATIVE Diag, idof = "
                           << idof << " inode = " << idof / 3 << std::endl;
                    exit(-5);
                }
        }

        //interpos(size = np) is size of points, contains -1(notinter) or 0~size of InterData-1
        void setInterDisplacement(const std::vector<size> &InterPos, const std::vector<utype<3>> &InterData)
        {
            bInterPos = InterPos;
            bInterData = InterData;

            for (int inz = 0; inz < K.seeNZ(); inz++)
                if (K.triplet(inz).i == K.triplet(inz).j) // is maj
                {
                    size idof = K.triplet(inz).i;
                    auto local = indDof2indNode(idof);
                    if (pointntype[local.first].first == inter && InterPos[idof / 3] >= 0)
                    {
                        utype<3> displacement = InterData[InterPos[idof / 3]];
                        K.triplet(inz).v = DBL_LARGE, loadVec[idof] = DBL_LARGE * displacement[local.second];
                        //std::cout << "Inter Disp " << idof / 3 << '\t' << idof % 3 << '\t' << displacement[local.second]
                        //          << '\t' << points[idof / 3][0] << '\t' << points[idof / 3][1] << '\t' << points[idof / 3][2] << '\t' << std::endl;
                    }
                }
        }

        void setFluidGridder()
        {
            fluidRegion = true;
        }

        void getStrainStress(const std::vector<real> &UActual, std::ostream &logout = std::cout)
        {
            strainElem.resize(tets.size());
            strainNode.resize(points.size());
#ifndef OMP_ON
            for (auto &e : strainNode)
                e.setZero();
#else
#pragma omp parallel
            for (int it = 0; it < strainNode.size(); it++)
                strainNode[it].setZero();
#endif
            stressNode.resize(points.size());
            VMNode.resize(points.size());
            nodeContributeNum.resize(points.size());

            //set dN_i/dxii_i, where N1 = xii1, N2 = xii2, N3 = xii3, N4 = 1 - xii1 - xii2 - xii3
            Eigen::Matrix<real, 3, 4> dNjdxiii = Eigen::Matrix<real, 3, 4>::Zero();
            dNjdxiii(0, 0) = 1, dNjdxiii(1, 1) = 1, dNjdxiii(2, 2) = 1, dNjdxiii(0, 3) = dNjdxiii(1, 3) = dNjdxiii(2, 3) = -1;
            //only one gaussian point(at centre)

            Eigen::Matrix<real, 3, 4> dNjdxi;   // cauculation: J^-1 * dNjdxiii
            Eigen::Matrix<real, 3, 3> dxjdxiii; //'Jacobian'
            //Eigen::Matrix<real, 3, 3> dxiijdxi; // inv J
            Eigen::Matrix<real, 6, 12> Bmat; // epsilon = B * a, a = [x1y1z1x2y2z2...]^T, epsilon = [e11 e22 e33 e23 e31 e12]^T
            Eigen::Matrix<real, 6, 3> BmatSlice;
            Eigen::Matrix<real, 4, 3> xpi; //x_pi, p is node idx, i is dim_idx
            BmatSlice.setZero();           //only need setting zero once
            Eigen::Matrix<real, 12, 1> localU;
            for (size it = 0; it < tets.size(); it++)
            {
                size ps[4] = {tets[it].ivert[0], tets[it].ivert[1], tets[it].ivert[2], tets[it].ivert[3]};
                for (int ii = 0; ii < 4; ii++)
                    for (int jj = 0; jj < 3; jj++)
                        xpi(ii, jj) = points[ps[ii]][jj];
                dxjdxiii = dNjdxiii * xpi;
                dNjdxi = dxjdxiii.colPivHouseholderQr().solve(dNjdxiii);
                for (int iv = 0; iv < 4; iv++)
                {
                    //BmatSlice.setZero();
                    BmatSlice(0, 0) = dNjdxi(0, iv), BmatSlice(1, 1) = dNjdxi(1, iv), BmatSlice(2, 2) = dNjdxi(2, iv);
                    BmatSlice(3, 1) = 0.5 * dNjdxi(2, iv), BmatSlice(3, 2) = 0.5 * dNjdxi(1, iv);
                    BmatSlice(4, 2) = 0.5 * dNjdxi(0, iv), BmatSlice(4, 0) = 0.5 * dNjdxi(2, iv);
                    BmatSlice(5, 0) = 0.5 * dNjdxi(1, iv), BmatSlice(5, 1) = 0.5 * dNjdxi(0, iv); // geometry relations here
                    Bmat.block(0, 3 * iv, 6, 3) = BmatSlice;
                    localU(iv * 3 + 0) = UActual[ps[iv] * 3 + 0];
                    localU(iv * 3 + 1) = UActual[ps[iv] * 3 + 1];
                    localU(iv * 3 + 2) = UActual[ps[iv] * 3 + 2];
                }
                strainElem[it] = Bmat * localU;
                for (int iv = 0; iv < 4; iv++)
                {
                    strainNode[ps[iv]] += strainElem[it];
                    nodeContributeNum[ps[iv]]++;
                }
            }

            for (size ip = 0; ip < points.size(); ip++)
            {
                strainNode[ip] /= nodeContributeNum[ip];
                stressNode[ip] = Dmat * strainNode[ip];
                ev6 &s = stressNode[ip];
                VMNode[ip] = std::sqrt((sqr(s(0) - s(1)) + sqr(s(1) - s(2)) + sqr(s(2) - s(0)) +
                                        6 * (sqr(s(3)) + sqr(s(4)) + sqr(s(5)))) /
                                       2.0); //Von-Mises( = sqrt(3 * J2))
            }
        }

        void modeprintTecPlotDat(const std::string &filenameHead, real magnifier,
                                 std::ostream &logout = std::cout)
        {
            for (size k = 0; k < lams.size(); k++)
            {
                char nbuf[512];
                sprintf(nbuf, "_Mode_%d_%.6e.dat", k, lams[k]);
                getStrainStress(phi[k], logout);
                fieldprintTecPlotDat(filenameHead + nbuf, magnifier, phi[k], logout);
            }
        }

        //from modes assemble disp and wirte to file
        void modeSum(const std::vector<real> &modeMs, std::vector<real> Usum,
                     const std::string &filename, real magnifier,
                     std::ostream &logout = std::cout)
        {
            Usum.assign(tDOF, 0);
            for (size k = 0; k < lams.size(); k++)
            {
                if (phi[k].size() != tDOF)
                {
                    std::cerr << "fem modeSUM inputmode bad size" << std::endl;
                    exit(-5);
                }
                for (size d = 0; d < tDOF; d++)
                    Usum[d] += phi[k][d] * modeMs[k];
            }
            getStrainStress(Usum, logout);
            fieldprintTecPlotDat(filename, magnifier, Usum, logout);
        }

        //muse see to stress before go
        void fieldprintTecPlotDat(const std::string &filename, real magnifier,
                                  const std::vector<real> &Uactual, std::ostream &logout = std::cout) const
        {

            logout << "---FILE TO--- : " << filename << std::endl;
            std::ofstream out(filename);
            if (!out)
            {
                logout << "===<<= FILE OPEN FAILED =>>===" << std::endl;
                exit(-1);
            }
            //logout << "Probe1" << std::endl;
            out << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"VM\", \"BT\"\n";
            out << "ZONE N=" << points.size() << " , E=" << tets.size() << ", ZONETYPE=FETETRAHEDRON\n";
            //logout << "Probe2" << std::endl;
            out << "DATAPACKING = BLOCK\n";
            out << "VARLOCATION=([1-8]=NODAL)\n";
            out << std::scientific << std::setw(19) << std::setprecision(15);

            if (Uactual.size() < tDOF)
            {
                logout << "=== ERROR: Missing Data === Uactual size = " << Uactual.size() << std::endl;
                exit(-1);
            }
            for (size it = 0; it < points.size(); it++)
                out << points[it][0] + magnifier * Uactual[it * 3 + 0] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][1] + magnifier * Uactual[it * 3 + 1] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << points[it][2] + magnifier * Uactual[it * 3 + 2] << '\n';

            for (size it = 0; it < points.size(); it++)
                out << Uactual[it * 3 + 0] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << Uactual[it * 3 + 1] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << Uactual[it * 3 + 2] << '\n';
            for (size it = 0; it < points.size(); it++)
                out << VMNode[it] << '\n';

            for (size it = 0; it < points.size(); it++)
                out << (int)(pointntype[it].first) << '\n';

            for (size it = 0; it < tets.size(); it++)
                out << tets[it].ivert[0] + 1 << '\t'
                    << tets[it].ivert[1] + 1 << '\t'
                    << tets[it].ivert[2] + 1 << '\t'
                    << tets[it].ivert[3] + 1 << '\n';

            out.close();
        }

        size indNode2indDOF(size indNode, size idim)
        {
            return indNode * 3 + idim;
        }

        size localDOF2indDOF(size localDOF, size *iNods)
        {
            return iNods[(localDOF / 3)] * 3 + (localDOF % 3);
        }

        //first = indNode, second = idim
        std::pair<size, size> indDof2indNode(size indDof)
        {
            return std::make_pair(indDof / 3, indDof % 3);
        }

        void setElas(const ElasSet &newElas)
        {
            set = newElas;
        }
        enum StaticSolverType
        {
            LDLT = 0,
            PCG = 1,
            SOR = 2
        };
        void SolveStatic(std::ostream &logout = std::cout, StaticSolverType sovid = LDLT)
        {
            U.resize(loadVec.size());
            dUdt.resize(loadVec.size());
            switch (sovid)
            {
            case LDLT:
                SolveSparseLDL(K, loadVec, U, logout);
                break;
            case PCG:
                SolveSparseCG(K, loadVec, U, logout);
                break;
            case SOR:
                fixedDOF.assign(tDOF, 0);
                for (int i = 0; i < tDOF; i++)
                    if (pointntype[i / 3].first == fixed || (fluidRegion && (pointntype[i / 3].first != inner)))
                        fixedDOF[i] = 1;
                SolveMatSOR(K, loadVec, U, fixedDOF, 1e-10, 1.5, 400);
                break;
            default:
                break;
            }
        }

        void SolveStaticPre(std::ostream &logout = std::cout)
        {
            if (pKsolver)
                delete pKsolver;
            EK.resize(tDOF, tDOF);
            pKsolver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>>(EK);
            SolveSparseLDL_Step1(K, EK, *pKsolver, logout);
        }

        void
        SolveStaticCal(std::ostream &logout = std::cout)
        {

            SolveSparseLDL_Step2(*pKsolver, loadVec, U, logout);
        }

        void SolveEigen(size Ksolve = 10, real tol = 1e-10, std::ostream &logout = std::cout)
        {
            SpEigenIPowerSet set;
            set.tol = tol;
            SolveEigenInvLDL<real>(Ksolve, tDOF, K, M, phi, lams, logout, set);
        }

        std::vector<real> &GetU()
        {
            return U;
        }

        std::vector<std::vector<real>> &GetModes()
        {
            return phi;
        }

        std::vector<real> &GetEigs()
        {
            return lams;
        }

        using TetraMesh::showPoints;
    };

}
