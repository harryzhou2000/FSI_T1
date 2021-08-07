#pragma once
#include <vector>
#include <iostream>
#include <omp.h>
#include <string>
#include <string.h>
#include <iomanip>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "defs.h"

namespace UTBasic
{

    constexpr real sqr(real x)
    {
        return x * x;
    }
    constexpr real cube(real x)
    {
        return x * x * x;
    }
    constexpr real sign(real x)
    {
        return x > 0 ? 1 : -1;
    }
    constexpr real max(real x, real y)
    {
        return x > y ? x : y;
    }
    inline real dis2(real x0, real y0, real x1, real y1)
    {
        return sqrt(sqr(x0 - x1) + sqr(y0 - y1));
    }

    template <int dim>
    struct utype
    {
        real v[dim];
        utype()
        {
            for (int i = 0; i < dim; i++)
                v[i] = 0;
        }
        utype(real x)
        {
            for (int i = 0; i < dim; i++)
                v[i] = x;
        }
        utype(const real x[dim])
        {
            for (int i = 0; i < dim; i++)
                v[i] = x[i];
        }
        utype(const utype<dim> &R)
        {
            for (int i = 0; i < dim; i++)
                v[i] = R.v[i];
        }
        real &operator[](int i) { return v[i]; }

        real operator[](int i) const { return v[i]; }

        utype<dim> operator+(const utype<dim> &R) const //add R
        {
            utype<dim> ret;
            for (int i = 0; i < dim; i++)
                ret.v[i] = R.v[i] + v[i];
            return ret;
        }
        utype<dim> operator-(const utype<dim> &R) const //minus R
        {
            utype<dim> ret;
            for (int i = 0; i < dim; i++)
                ret.v[i] = v[i] - R.v[i];
            return ret;
        }
        utype<dim> operator*(const real R) const //times scalar
        {
            utype<dim> ret;
            for (int i = 0; i < dim; i++)
                ret.v[i] = v[i] * R;
            return ret;
        }
        utype<dim> operator*(const utype<dim> &R) const //times R
        {
            utype<dim> ret;
            for (int i = 0; i < dim; i++)
                ret.v[i] = v[i] * R[i];
            return ret;
        }

        void operator+=(const utype<dim> &R) //add by R
        {
            for (int i = 0; i < dim; i++)
                v[i] += R.v[i];
        }

        void operator-=(const utype<dim> &R) //minus by R
        {
            for (int i = 0; i < dim; i++)
                v[i] -= R.v[i];
        }

        void operator/=(const utype<dim> &R) //divide by R
        {
            for (int i = 0; i < dim; i++)
                v[i] /= R.v[i];
        }

        utype<dim> operator/(const utype<dim> &R) //divide by R
        {
            utype<dim> ret(*this);
            for (int i = 0; i < dim; i++)
                ret[i] /= R.v[i];
            return ret;
        }

        void operator*=(real R) //times scalar
        {
            for (int i = 0; i < dim; i++)
                v[i] *= R;
        }

        void set(const real *a) //set uniform
        {
            for (int i = 0; i < dim; i++)
                v[i] = a[i];
        }

        void set0()
        {
            for (int i = 0; i < dim; i++)
                v[i] = 0;
        }

        utype<dim> operator-() //negative of
        {
            utype<dim> ret;
            for (int i = 0; i < dim; i++)
                ret.v[i] = -v[i];
            return ret;
        }

        utype<dim> abs() const
        {
            utype<dim> ret(*this);
            for (int i = 0; i < dim; i++)
                ret[i] = std::abs(ret[i]);
            return ret;
        }

        utype<dim> max(const utype<dim> &R) const //max by component
        {
            utype<dim> ret(*this);
            for (int i = 0; i < dim; i++)
                ret[i] = std::max(v[i], R[i]);
            return ret;
        }

        utype<dim> min(const utype<dim> &R) const //min by component
        {
            utype<dim> ret(*this);
            for (int i = 0; i < dim; i++)
                ret[i] = std::min(v[i], R[i]);
            return ret;
        }

        real minelem()
        {
            real ret = v[0];
            for (int i = 1; i < dim; i++)
                ret = std::min(v[i], ret);
            return ret;
        }

        real dot(const utype<dim> &R) const
        {
            real ans = 0.;
            for (int i = 0; i < dim; i++)
                ans += R[i] * v[i];
            return ans;
        }

        real norm2() const
        {
            real ans = 0.;
            for (int i = 0; i < dim; i++)
                ans += v[i] * v[i];
            return std::sqrt(ans);
        }

        void normalize()
        {
            real dn2 = 1.0 / norm2();
            for (int i = 0; i < dim; i++)
                v[i] *= dn2;
        }
    };

    utype<3> crossUT3(const utype<3> &L, const utype<3> &R)
    {
        utype<3> ret;
        ret[0] = L[1] * R[2] - L[2] * R[1];
        ret[1] = L[2] * R[0] - L[0] * R[2];
        ret[2] = L[0] * R[1] - L[1] * R[0];
        return ret;
    }

    template <size col, size row>
    struct mtype
    {
        real v[col][row];

        mtype()
        {
            for (size i = 0; i < col; i++)
                for (size j = 0; j < row; j++)
                    v[i][j] = 0;
        }

        mtype(const mtype<col, row> &R)
        {
            for (size i = 0; i < col; i++)
                for (size j = 0; j < row; j++)
                    v[i][j] = R(i, j);
        }

        void set0()
        {
            for (size i = 0; i < col; i++)
                for (size j = 0; j < row; j++)
                    v[i][j] = 0;
        }

        real &operator()(const size m, const size n)
        {
            return v[m][n];
        }

        utype<col> operator*(const utype<row> &r)
        {
            utype<col> ret;
            for (size i = 0; i < col; i++)
                for (size j = 0; j < row; j++)
                    ret[i] += r[j] * v[i][j];
            return ret;
        }
    };
    mtype<2, 2> a;

#define DIM 3
#define FDIM 5
    typedef utype<FDIM> ut;

    template <unsigned dim>
    utype<dim + 2> TVDFunc(const utype<dim + 2> &dl, const utype<dim + 2> &dr)
    {
        utype<dim + 2> ret;
        for (int i = 0; i < dim + 2; i++)
            if ((dl[i] <= 0 && dr[i] >= 0) || (dl[i] >= 0 && dr[i] <= 0))
                ret[i] = 0;
            else
                ret[i] = 2 * dr[i] / (dl[i] + dr[i]);
        return ret;
    }

    void operator*=(std::vector<ut> &a, real R)
    {
#ifndef OMP_ON
        for (auto &x : a)
            x *= R;
#endif
#ifdef OMP_ON
#pragma omp parallel for
        for (size i = 0; i < a.size(); i++)
            a[i] *= R;
#endif
    }

    void operator*=(std::vector<utype<4>> &a, real R)
    {
#ifndef OMP_ON
        for (auto &x : a)
            x *= R;
#endif
#ifdef OMP_ON
#pragma omp parallel for
        for (size i = 0; i < a.size(); i++)
            a[i] *= R;
#endif
    }

    void operator+=(std::vector<ut> &a, const std::vector<ut> &R) //warning, R size not checked
    {
#ifdef OMP_ON
#pragma omp parallel for
#endif
        for (size i = 0; i < a.size(); i++)
            a[i] += R[i];
    }

    void operator+=(std::vector<utype<4>> &a, const std::vector<utype<4>> &R) //warning, R size not checked
    {
#ifdef OMP_ON
#pragma omp parallel for
#endif
        for (size i = 0; i < a.size(); i++)
            a[i] += R[i];
    }

    void operator-=(std::vector<ut> &a, const std::vector<ut> &R) //warning, R size not checked
    {
#ifdef OMP_ON
#pragma omp parallel for
#endif
        for (size i = 0; i < a.size(); i++)
            a[i] -= R[i];
    }

    void operator-=(std::vector<utype<4>> &a, const std::vector<utype<4>> &R) //warning, R size not checked
    {
#ifdef OMP_ON
#pragma omp parallel for
#endif
        for (size i = 0; i < a.size(); i++)
            a[i] -= R[i];
    }

    template <unsigned dim>
    utype<dim> maxabs(const std::vector<utype<dim>> &a)
    {
        utype<dim> ret(0.0);
        for (size i = 0; i < a.size(); i++)
            ret = ret.max(a[i].abs());
        return ret;
    }

    real triangleArea2nd(utype<3> a, utype<3> b, utype<3> c)
    {
        return crossUT3(b - a, c - a).norm2();
    }

    utype<3> triangleBary(utype<3> a, utype<3> b, utype<3> c)
    {
        return a + ((b - a) + (c - a)) * (1.0 / 3.0);
    }

    real tetraVolume6th(utype<3> a, utype<3> b, utype<3> c, utype<3> d)
    {
        return std::abs(crossUT3(b - a, c - a).dot(d - a));
    }

    real pointToLine(real x, real y, real x0, real y0, real x1, real y1)
    {
        real area = fabs((x0 - x) * (y1 - y) - (x1 - x) * (y0 - y));
        return area / sqrt(sqr(x0 - x1) + sqr(y0 - y1));
    }

    real pointToFace(utype<3> p, utype<3> a, utype<3> b, utype<3> c)
    {
        return tetraVolume6th(p, a, b, c) / triangleArea2nd(a, b, c);
    }

    enum __ntype__
    {
        inner = 0,
        infinity = 1,        // CFD
        wall = 2,            // CFD
        pressure = 3,        // CFD
        freesurf = 2001,     // CSD
        fixed = 2002,        // CSD
        pressuresurf = 2003, // CSD
        cross = 9003,        // currently unused
        inter = 66,          // interaction for FSI. wall for CFD and force/free for CSD
        unknown = 10203
    };

#ifdef INF_FIRST
    __ntype__ deriveFtype3(__ntype__ t1, __ntype__ t2, __ntype__ t3)
    {
        if (t1 && t2 && t3)
        {
            if ((t1 == infinity) && (t2 == infinity) && (t3 == infinity))
                return infinity;
            if ((t1 == pressure) && (t2 == pressure) && (t3 == pressure))
                return pressure;
            if ((t1 == fixed) && (t2 == fixed) && (t3 == fixed))
                return fixed;
            if ((t1 == wall) || (t2 == wall) || (t3 == wall))
                return wall;
            if ((t1 == pressuresurf) || (t2 == pressuresurf) || (t3 == pressuresurf))
                return pressuresurf; // TODO should discuss here
            if ((t1 == inter) || (t2 == inter) || (t3 == inter))
                return inter;
            if ((t1 == freesurf) || (t2 == freesurf) || (t3 == freesurf))
                return freesurf;

            if ((t1 == cross) || (t2 == cross) || (t3 == cross))
                return cross;
            return unknown;
        }
        return inner;
    }
#else
    __ntype__ deriveFtype3(__ntype__ t1, __ntype__ t2, __ntype__ t3)
    {
        if (t1 && t2 && t3)
        {
            if ((t1 == infinity) || (t2 == infinity) || (t3 == infinity))
                return infinity;
            if ((t1 == pressure) || (t2 == pressure) || (t3 == pressure))
                return pressure;
            if ((t1 == fixed) && (t2 == fixed) && (t3 == fixed))
                return fixed;
            if ((t1 == wall) || (t2 == wall) || (t3 == wall))
                return wall;
            if ((t1 == pressuresurf) || (t2 == pressuresurf) || (t3 == pressuresurf))
                return pressuresurf; // TODO should discuss here
            if ((t1 == inter) || (t2 == inter) || (t3 == inter))
                return inter;
            if ((t1 == freesurf) || (t2 == freesurf) || (t3 == freesurf))
                return freesurf;

            if ((t1 == cross) || (t2 == cross) || (t3 == cross))
                return cross;
            return unknown;
        }
        return inner;
    }
#endif

#ifdef INF_FIRST
    size deriveNbound3(__ntype__ t1, __ntype__ t2, __ntype__ t3, size n1, size n2, size n3) //when all are not inner
    {
        if ((t1 == infinity) && (t2 == infinity) && (t3 == infinity))
            return n1;
        if ((t1 == pressure) && (t2 == pressure) && (t3 == pressure))
            return n1;
        if (t1 == wall || t1 == inter || t1 == pressuresurf)
            return n1;
        if (t2 == wall || t2 == inter || t2 == pressuresurf)
            return n2;
        if (t3 == wall || t3 == inter || t3 == pressuresurf)
            return n3;
        return 0;
    }
    typedef __ntype__ ntype;
#else
    size deriveNbound3(__ntype__ t1, __ntype__ t2, __ntype__ t3, size n1, size n2, size n3) //when all are not inner
    {
        if ((t1 == wall) && (t2 == wall) && (t3 == wall))
            return n1;
        if ((t1 == inter) && (t2 == inter) && (t3 == inter))
            return n1;
        if (t1 == infinity || t1 == pressure || t1 == pressuresurf)
            return n1;
        if (t2 == infinity || t2 == pressure || t2 == pressuresurf)
            return n2;
        if (t3 == infinity || t3 == pressure || t3 == pressuresurf)
            return n3;
        return 0;
    }

#endif
    typedef __ntype__ ntype;

    template <size dim>
    bool inbox(const utype<dim> p, const utype<dim> lb, const utype<dim> ub)
    {
        for (int i = 0; i < dim; i++)
            if (p[i] < lb[i] || p[i] > ub[i])
                return false;
        return true;
    }

    void TransAddTo(utype<5> &U, const utype<3> &velo)
    {
        U[1] += velo[0] * U[0];
        U[2] += velo[1] * U[0];
        U[3] += velo[2] * U[0];
    }

    // WHY ???
    template <size d>
    std::ostream &operator<<(std::ostream &out, const utype<d> &u)
    {
        for (size i = 0; i < d; i++)
            out << u[i] << '\t';
        return out;
    }

    std::ostream &operator<<(std::ostream &out, const utype<5> &u)
    {
        for (size i = 0; i < 5; i++)
            out << u[i] << '\t';
        return out;
    }

    std::ostream &operator<<(std::ostream &out, const utype<3> &u)
    {
        for (size i = 0; i < 3; i++)
            out << u[i] << '\t';
        return out;
    }
    ///////

    template <size d>
    std::istream &operator>>(std::istream &in, utype<d> &u)
    {
        for (size i = 0; i < d; i++)
            in >> u[i];
        return in;
    }
    template <typename T>
    void SaveVector(const std::vector<T> &vec, std::ostream &out)
    {
        out << "\nVECTOR " << std::scientific << std::setprecision(16) << std::setw(20)
            << vec.size() << '\n';
        for (size i = 0; i < vec.size(); i++)
            out << vec[i] << "\n";
    }

    template <typename T>
    bool ReadVector(std::vector<T> &vec, std::istream &in)
    {
        std::string name;
        in >> name;
        if (name != "VECTOR")
        {
            std::cerr << "VECTOR READ FALILIURE\n"
                      << std::endl;
            return false;
        }
        size vsize;
        in >> vsize;
        vec.resize(vsize);
        for (size i = 0; i < vec.size(); i++)
            in >> vec[i];
        return true;
    }

    template <typename T>
    std::istream &operator>>(std::istream &in, std::vector<T> &vec)
    {
        ReadVector(vec, in);
        return in;
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec)
    {
        SaveVector(vec, out);
        return out;
    }

    const real pi = std::acos(-1.0);
}