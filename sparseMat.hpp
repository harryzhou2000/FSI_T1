#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <math.h>
#include <cmath>
#include <time.h>
#include <omp.h>
#include <utility>
#include "include/Eigen/Sparse"
#include "include/Eigen/Dense"
#include "defs.h"

#define vector std::vector
#define string std::string
#define ostream std::ostream
#define ifstream std::ifstream
#define stringstream std::stringstream

typedef double DOUBLE;

namespace SparseMat
{

    template <class R>
    struct MatTriplet // Row Major Seq
    {
        size i;
        size j;
        R v;
        MatTriplet(size ii, size jj, R vv = 0) : i(ii), j(jj), v(vv) {}
        bool operator>(const MatTriplet &Ri) const
        {
            return i == Ri.i ? (j > Ri.j) : i > Ri.i;
        }
        bool operator<(const MatTriplet &Ri) const
        {
            return i == Ri.i ? (j < Ri.j) : i < Ri.i;
        }
        bool operator==(const MatTriplet &Ri) const
        {
            return i == Ri.i && j == Ri.j;
        }
    };

    template <class R>
    class SparseMatGeR
    {
    private:
        vector<MatTriplet<R>> dat;
        bool inorder = false;
        bool hasrowstarts = false;
        size szi, szj;
        vector<size> RowStarts;
        bool reduced;

    public:
        MatTriplet<R> &triplet(size ind)
        {
            return dat[ind];
        }

        void clear()
        {
            dat.clear();
            inorder = false;
            hasrowstarts = false;
            szi = szj = -1;
            reduced = false;
        }

        void reserve(size nrow, size rowLoad)
        {
            dat.reserve(rowLoad * nrow);
        }

        auto find(size i, size j)
        {
            reduce();
            auto ans = std::lower_bound(dat.begin(), dat.end(), MatTriplet<R>(i, j));
            return ans;
        }

        inline void addTo(size i, size j, R const &value)
        {
            inorder = false;
            hasrowstarts = false;
            reduced = false;
            szi = -1, szj = -1;
            dat.push_back(MatTriplet<R>(i, j, value));
        }

        inline void addToFast(size i, size j, R const &value)
        {
            dat.push_back(MatTriplet<R>(i, j, value));
        }

        void outputEigen(vector<Eigen::Triplet<DOUBLE>> &o, const vector<size> &imap)
        {
            o.resize(dat.size());
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < dat.size(); i++)
            {
                o[i] = Eigen::Triplet<DOUBLE>(imap[dat[i].i], imap[dat[i].j], dat[i].v);
            }
        }

        void outputEigen(vector<Eigen::Triplet<DOUBLE>> &o) const
        {
            o.resize(dat.size());
#ifdef OMP_ON
#pragma omp parallel for
#endif
            for (size i = 0; i < dat.size(); i++)
            {
                o[i] = Eigen::Triplet<DOUBLE>(dat[i].i, dat[i].j, dat[i].v);
            }
        }

        void transverseTo(SparseMatGeR<R> &target)
        {
            if (&target - this)
            {
                target.clear();
                target.dat = dat;
            }
            else
            {
                inorder = false;
                hasrowstarts = false;
            }
#ifdef OMP_ON
#pragma omp parallel for
            for (size iter = 0; iter < dat.size(); iter++)
                std::swap(target.dat[iter].i, target.dat[iter].j);
#else
            for (auto &p : target.dat)
                std::swap(p.i, p.j);
#endif
        }

        void MatDenseVec(const vector<R> &b, vector<R> &target)
        {
            size szp = target.size();
            target.clear();
            target.resize(szp, 0);
#ifdef OMP_ON
            buildRowStarts();
#pragma omp parallel for //schedule(guided)//num_threads(2)
            for (size iterRow = 0; iterRow < RowStarts.size() - 1; iterRow++)
                for (size iter = RowStarts[iterRow]; iter < RowStarts[iterRow + 1]; iter++)
                {
                    if (std::abs(dat[iter].v) < DBL_LARGE * 0.95)
                        target[dat[iter].i] += dat[iter].v * b[dat[iter].j];
                }
#else
            for (auto const &p : dat)
                if (p.i < target.size() && p.j < b.size())
                    target[p.i] += p.v * b[p.j];
#endif
        }

        void MatDenseVec(const vector<R> &b, vector<R> &target, const vector<int> ignore)
        {
            size szp = target.size();
            target.clear();
            target.resize(szp, 0);
#ifdef OMP_ON
            buildRowStarts();
#pragma omp parallel for //schedule(guided)//num_threads(2)
            for (size iterRow = 0; iterRow < RowStarts.size() - 1; iterRow++)
                for (size iter = RowStarts[iterRow]; iter < RowStarts[iterRow + 1]; iter++)
                {
                    if (!ignore[dat[iter].i])
                        target[dat[iter].i] += dat[iter].v * b[dat[iter].j];
                }
#else
            for (auto const &p : dat)
                if (p.i < target.size() && p.j < b.size())
                    target[p.i] += p.v * b[p.j];
#endif
        }

        const R NormInf() const
        {
            R ret = 0;
            for (auto i : dat)
                ret = std::max(std::abs(i.v), ret);
            return ret;
        }

        void sort()
        {
            if (!inorder)
                std::sort(dat.begin(), dat.end());
            inorder = true;
        }

        void getsize()
        {
            sort();
            if (dat.size() == 0)
            {
                szi = szj = -1;
                return;
            }
            auto &nsz = dat[dat.size() - 1];
            szi = nsz.i;
            szj = nsz.j;
        }

        void reduce()
        {
            if (reduced)
                return;
            sort();
            vector<MatTriplet<R>> datn;
            datn.reserve(dat.size());

            for (const auto &p : dat)
            {
                if (datn.size() > 0 && p == datn[datn.size() - 1])
                    datn[datn.size() - 1].v += p.v;
                else
                    datn.push_back(p);
            }
            datn.reserve(datn.size());
            dat = std::move(datn);
            reduced = true;
        }

        void Print(string filename)
        {
            std::ofstream file(filename);
            for (const auto &p : dat)
            {
                file << p.i << '\t' << p.j << '\t' << p.v << '\n';
            }
        }

        void Print(ostream &out)
        {
            for (const auto &p : dat)
            {
                out << p.i << '\t' << p.j << '\t' << p.v << '\n';
            }
        }

        void normalizeRowInf(vector<R> b) //to be implemented
        {
        }

        std::pair<size, size> seeSize()
        {
            return std::make_pair(szi + 1, szj + 1);
        }

        size seeNZ()
        {
            return dat.size();
        }

        void buildRowStarts()
        {
            if (dat.empty() || hasrowstarts)
                return;
            sort();
            size nRow = 0;
            size curRowi = dat[0].i;
            for (const auto &elem : dat)
            {
                if (curRowi - elem.i)
                {
                    nRow++;
                    curRowi = elem.i;
                }
            }
            RowStarts.resize(nRow + 2);
            nRow = 0;
            curRowi = dat[0].i;
            RowStarts[nRow] = 0;
            for (size plc = 0; plc < dat.size(); plc++)
            {
                if (curRowi - dat[plc].i)
                {
                    nRow++;
                    curRowi = dat[plc].i;
                    RowStarts[nRow] = plc;
                }
            }
            RowStarts[nRow + 1] = dat.size();
            hasrowstarts = true;
        }

        void SORStep(vector<R> const &b, vector<R> const &x, vector<R> &xnew, R omega, size ns)
        {
            xnew.clear();
            xnew.resize(ns, 0);
            sort();
            getsize();

            size it(0), jt(0);
            for (; it < dat.size();)
            {
                size curi = dat[it].i;
                xnew[curi] = (1 - omega) * x[curi];
                R aii = 0;
                R gs1 = b[curi];
                for (jt = 0; it + jt < dat.size() && dat[it + jt].i == curi; jt++)
                {
                    size curj = dat[it + jt].j;
                    if (curj < curi)
                    {
                        gs1 -= dat[it + jt].v * xnew[curj];
                    }
                    else if (curj == curi)
                    {
                        aii = dat[it + jt].v;
                    }
                    else if (curj > curi)
                    {
                        gs1 -= dat[it + jt].v * x[curj];
                    }
                }
                it += jt;
                gs1 *= omega / aii;
                xnew[curi] += gs1;
            }
        }
    };

#undef string
#undef vector
#undef ostream
#undef ifstream

    template <class R>
    R normInf(std::vector<R> a, std::vector<R> b, const std::vector<int> &ignore)
    {
        R result = 0;
        size len = std::min(a.size(), b.size());
        for (size i = 0; i < len; i++)
            if (!ignore[i])
            {
                R diff = a[i] - b[i];
                diff = diff > 0 ? diff : -diff;
                result = diff > result ? diff : result;
            }
        return result;
    }

    template <class R>
    void SolveMatSOR(SparseMatGeR<R> &A, std::vector<R> &b, std::vector<R> &sol, std::vector<int> &ignore,
                     real res_threshold = 1e-5, DOUBLE omega = 1., size itmax = 500)
    {
        SparseMatGeR<DOUBLE> A1(A);
        std::vector<DOUBLE> b1(b);
        A1.normalizeRowInf(b1);
        sol.clear();
        sol.resize(b.size(), 0.);
        size ns = b.size();
        std::vector<DOUBLE> soln(sol);
        std::vector<DOUBLE> Asol(sol);
        A1.MatDenseVec(sol, Asol, ignore);
        DOUBLE res0 = normInf<DOUBLE>(b1, Asol, ignore);
        DOUBLE restarget = res0 * res_threshold;
        std::cout << "ResTarget = " << restarget << std::endl;

        for (int iter = 1; iter <= itmax; iter++)
        {
            A1.SORStep(b1, sol, soln, omega, ns);
            sol = soln;
            if (iter % 100 == 0)
            {
                A1.MatDenseVec(sol, Asol);
                DOUBLE res = normInf<DOUBLE>(b1, Asol, ignore);
                //logout << "R^2 = " << rdotr << endl;
                std::cout << "ResInf = " << res << std::endl;
                if (res <= restarget)
                    break;
            }
        }

        //for (auto &e : sol)
        //    e *= bRatio;
    }

    //input should be symmetric
    template <class R>
    void SolveSparseLDL(const SparseMatGeR<R> &Mat, const std::vector<R> &Vec, std::vector<R> &Sol, std::ostream &logout)
    {
#ifdef OMP_ON
        Eigen::setNbThreads(128);
#endif
        std::vector<Eigen::Triplet<R>> mids;
        Mat.outputEigen(mids);
        Eigen::SparseMatrix<R> EMat(Vec.size(), Vec.size());
        EMat.setFromTriplets(mids.begin(), mids.end());

        Eigen::Matrix<R, -1, 1> EVec(Vec.size());
        for (int i = 0; i < Vec.size(); i++)
            EVec(i) = Vec[i];
        Eigen::Matrix<R, -1, 1> ESol(Vec.size());
        Eigen::SimplicialLDLT<decltype(EMat)> solver;
        logout << "LDL Decomposition..."
               << std::endl;
        solver.compute(EMat);
        if (solver.info() != Eigen::Success)
        {
            logout << "LDL Faliure!!!"
                   << std::endl;
            exit(0);
        }
        logout << "LDL Solving..."
               << std::endl;
        ESol = solver.solve(EVec);
        for (int i = 0; i < Vec.size(); i++)
            Sol[i] = ESol(i);
#ifdef _DEBUG
        logout << ESol << std::endl;
#endif
    }

    template <class R>
    void SolveSparseCG(const SparseMatGeR<R> &Mat, const std::vector<R> &Vec, std::vector<R> &Sol,
                       std::ostream &logout, real tol = 1e-50)
    {
#ifdef OMP_ON
        Eigen::setNbThreads(128 );
#endif
        std::vector<Eigen::Triplet<R>> mids;
        Mat.outputEigen(mids);
        Eigen::SparseMatrix<R> EMat(Vec.size(), Vec.size());
        EMat.setFromTriplets(mids.begin(), mids.end());

        Eigen::Matrix<R, -1, 1> EVec(Vec.size());
        for (int i = 0; i < Vec.size(); i++)
            EVec(i) = Vec[i];
        Eigen::Matrix<R, -1, 1> ESol(Vec.size());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<R>, Eigen::Lower | Eigen::Upper, Eigen::LeastSquareDiagonalPreconditioner<R>> cg;

        cg.setTolerance(tol);
        cg.setMaxIterations(500);
        logout << "CG Computing..." << std::endl;
        cg.compute(EMat);
        logout << "CG Solving..." << std::endl;
        ESol = cg.solve(EVec);
        logout << "CG Done #iter = " << cg.iterations() << " Res = " << cg.error() << std::endl;
        for (int i = 0; i < Vec.size(); i++)
            Sol[i] = ESol(i);
    }

    template <class R>
    void SolveSparseLDL_Step1(const SparseMatGeR<R> &Mat,
                              Eigen::SparseMatrix<R> &EMat, Eigen::SimplicialLDLT<Eigen::SparseMatrix<R>> &solver,
                              std::ostream &logout)
    {
        std::vector<Eigen::Triplet<R>> mids;
        Mat.outputEigen(mids);
        EMat.setFromTriplets(mids.begin(), mids.end());
        logout << "LDL Decomposition...\n";
        solver.compute(EMat);
        if (solver.info() != Eigen::Success)
        {
            logout << "LDL Faliure!!!\n";
            exit(0);
        }
    }

    template <class R>
    void SolveSparseLDL_Step2(Eigen::SimplicialLDLT<Eigen::SparseMatrix<R>> &solver, const std::vector<R> &Vec, std::vector<R> &Sol, std::ostream &logout)
    {

        Eigen::Matrix<R, -1, 1> EVec(Vec.size());
        for (int i = 0; i < Vec.size(); i++)
            EVec(i) = Vec[i];

        logout << "LDL Solving...\n";
        Eigen::Matrix<R, -1, 1> ESol = solver.solve(EVec);
        for (int i = 0; i < Vec.size(); i++)
            Sol[i] = ESol(i);
    }

    namespace PE //inherited from planarElasticity
    {
        using Msize = size;
        using ostream = std::ostream;
        void projectToNormalSubspace(const Eigen::SparseMatrix<DOUBLE> &M,
                                     const Eigen::MatrixXd &V, Eigen::VectorXd &v, Msize i)
        {
            Eigen::VectorXd v1 = v;
            //std::cout << "Vproj\n"
            //          << v << std::endl;
            for (int a = 0; a < i; a++)
                v -= (V.col(a).transpose() * M * v1) * V.col(a);
        }

        const DOUBLE tolfloat = 1e-15;
        void generateFromSubspace(const Eigen::SparseMatrix<DOUBLE> &M, Eigen::MatrixXd &V,
                                  Eigen::VectorXd &v, Msize i, ostream &logout)
        {
            //v.setConstant(1.);
            v.setRandom();
            v(0) = 0.;
            for (Msize iter = 0; iter < v.size(); iter++)
            {
                projectToNormalSubspace(M, V, v, i);
                DOUBLE vsize = sqrt(v.transpose() * M * v);
                if (vsize > v.size() * tolfloat)
                {
                    v /= vsize;
                    return;
                }
                if (iter == v.size() - 1)
                {
                    logout << "\tIP: Subspace generation ERROR!!" << std::endl;
                    exit(0);
                }
                v.setConstant(1.);
                v(iter + 1) = 0;
            }
        }

        void inversedPowerIteration(Eigen::SimplicialLDLT<Eigen::SparseMatrix<DOUBLE>> &Asolver,
                                    const Eigen::SparseMatrix<DOUBLE> &A, const Eigen::SparseMatrix<DOUBLE> &M,
                                    const Eigen::MatrixXd &V, Eigen::VectorXd &v, Msize i,
                                    DOUBLE tol, Msize maxit, Msize itcheck, ostream &logout)
        {
            for (Msize iter = 1; iter <= maxit; iter++)
            {
                v = M * v;
                Eigen::VectorXd v1 = Asolver.solve(v);
                auto pres = (A * v1 - v).norm();
                v = v1;
                projectToNormalSubspace(M, V, v, i);
                v /= sqrt(v.transpose() * M * v);

                if (iter % itcheck == 0)
                {
                    auto res = A * v - M * v * ((v.transpose() * A * v) / (v.transpose() * M * v));
                    //std::cout << "VMDAV " << (v.transpose() * MDA * v) << std::endl;
                    DOUBLE resnorm = res.lpNorm<Eigen::Infinity>();
                    logout << "\t--IP iter " << iter << ", res = " << resnorm << std::endl;
                    if (resnorm <= (tol * (1 << i / 2)))
                    {
                        logout << "\t--IP done " << iter << ", number " << i << " of the eigens" << std::endl;
                        break;
                    }
                }
            }
        }
    }

    struct SpEigenIPowerSet
    {
        real tol = 1e-10;
        real reltol = 1e-2; //currently unused
        size maxit = 120;
        size itcheck = 8;
        SpEigenIPowerSet(){};
    };

    //input should be symmetric  and (semi for KMat) positive-definite
    template <class R>
    void SolveEigenInvLDL(size K, size n, const SparseMatGeR<R> &KMat, const SparseMatGeR<R> &MMat,
                          std::vector<std::vector<R>> &Phi, std::vector<R> &Eig,
                          std::ostream &logout = std::cout, SpEigenIPowerSet set = SpEigenIPowerSet())
    {
#ifdef OMP_ON
        Eigen::setNbThreads(64);
#endif
        Phi.resize(K, std::vector<R>(n)), Eig.resize(K);
        std::vector<Eigen::Triplet<R>> mids;
        //Get MMat
        MMat.outputEigen(mids);
        Eigen::SparseMatrix<R> EMMat(n, n);
        EMMat.setFromTriplets(mids.begin(), mids.end());
        //Get KMat
        KMat.outputEigen(mids);
        Eigen::SparseMatrix<R> EKMat(n, n);
        EKMat.setFromTriplets(mids.begin(), mids.end());
        //modification
        real MRatio = MMat.NormInf();
        EMMat /= MRatio;
        Eigen::SparseMatrix<R> EKMatT = EKMat.transpose();
        EKMat = (EKMat + EKMatT) * 0.5;
        EKMatT.setIdentity();
        EKMat = EKMat + EKMatT * 1e-9; //deflection
        //K LDL
        Eigen::SimplicialLDLT<decltype(EKMat)> solver;
        logout << "===Eigen IPower=== A Mat LDL Decomposition..." << std::endl;
        solver.compute(EKMat);
        if (solver.info() != Eigen::Success)
        {
            logout << "LDL Faliure!!!\n";
            exit(0);
        }
        logout << "===Eigen IPower=== A Mat LDL Decomposition Done" << std::endl;

        Eigen::VectorXd vv(n);
        Eigen::MatrixXd VV(n, K);
        vv.setConstant(1.);
        VV.setConstant(0.);
#ifdef _DEBUG
        logout << EMMat << std::endl
               << std::endl
               << EKMat << std::endl;

#endif

        for (size i = 0; i < K; i++)
        {
            PE::generateFromSubspace(EMMat, VV, vv, i, logout);
            //std::cout << vv << std::endl;
            PE::inversedPowerIteration(solver, EKMat, EMMat, VV, vv, i, set.tol, set.maxit, set.itcheck, logout);
            VV.col(i) = vv;
            Eig[i] = (vv.transpose() * EKMat * vv);
            Eig[i] /= (vv.transpose() * EMMat * vv);
            Eig[i] = (Eig[i] - 1e-9) / MRatio;
            for (size p = 0; p < n; p++)
                Phi[i][p] = vv(p) / std::sqrt(MRatio);
        }
    }

}
