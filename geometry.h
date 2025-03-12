#pragma once

#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>
#include <initializer_list>

const int SpaceDim = 3;
using real = double;

class Point3D
{
protected:
    real X[SpaceDim];
public:
    Point3D(real x, real y, real z)
    {
        X[0] = x; X[1] = y; X[2] = z;
    }
    Point3D(){memset(X,0,sizeof(X));}
    Point3D(std::initializer_list<real> a)
    {
        assert(a.size() == SpaceDim && "Point3D can be constructed only from 3 coordinates");
        auto _a = std::data(a);
        X[0] = _a[0]; X[1] = _a[1]; X[2] = _a[2];
    }
    Point3D(const real *d)
    {
        memcpy(X,d,sizeof(X));
    }
    real operator[](int i) const
    {
        assert(i >=0 && i < SpaceDim && "Index out of bounds");
        return X[i];
    }
    real& operator[](int i)
    {
        assert(i >=0 && i < SpaceDim && "Index out of bounds");
        return X[i];
    }
};

real distance2(const Point3D& a, const Point3D& b)
{
    real r = 0;
    real d_i = 0;
    for(int i = 0; i < SpaceDim; ++i)
    {
        d_i = b[i] - a[i];
        r += d_i * d_i;
    }
    return r;
}

real distance(const Point3D& a, const Point3D& b)
{
    return std::sqrt(distance2(a,b));
}

class Vector3D: public Point3D
{
public:
    Vector3D(const Point3D& x): Point3D(x) {}
    Vector3D(real x, real y, real z): Point3D(x,y,z) {}
    Vector3D(): Point3D() {}
    Vector3D(const real *d): Point3D(d) {}
    real dot(const Vector3D& v) const
    {
        const Vector3D &x(*this);
        real r = 0;
        for(int i = 0; i < SpaceDim; ++i)
        {
            r += x[i] * v[i];
        }
        return r;
    }
    Vector3D vector_product(const Vector3D& u) const
    {
        auto& v(*this);
        return Vector3D(v[1] * u[2] - v[2] * u[1], 
                      -(v[0] * u[2] - v[2] * u[0]), 
                        v[0] * u[1] - v[1] * u[0]);
    }
    real norm2() const
    {
        return dot(*this);
    }
    real norm() const
    {
        return std::sqrt(norm2());
    }
    Vector3D& operator *= (real a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] *= a;
        }
        return *this;
    }
    Vector3D& operator /= (real a)
    {
        assert(a != 0 && "Trying to divide by zero!");
        return (*this) *= 1.0 / a;
    }
    Vector3D& operator += (const Vector3D& a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] += a[i];
        }
        return *this;
    }
    Vector3D& operator -= (const Vector3D& a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] -= a[i];
        }
        return *this;
    }
    Vector3D operator - ()
    {
        Vector3D r(*this);
        r *= -1;
        return r;
    }
};

class Segment3D
{
protected:
    Vector3D X[2];
    Vector3D l;
    real S;
public:
    Segment3D(const Point3D& x0, const Point3D& x1)
    {
        X[0] = x0;
        X[1] = x1;
        l = X[1]; l -= X[0];
        S = l.norm();
        if(S == 0)
        {
            l = Vector3D(1,0,0);
        }
        else
        {
            l /= S;
        }
    }
    const Vector3D& begin()const{return X[0];}
    const Vector3D& end()const{return X[1];}
    const Vector3D& ort()const{return l;}
    real lenght() const {return S;}
    real projection(const Point3D& x)const
    {
        return (Vector3D(x)-=begin()).dot(ort());
    }
    int is_projection_out(real p)const
    {
        if(p < 0) return -1;
        else if(p < lenght()) return 0;
        else return 1;
    }
    bool parallel(const Segment3D& s, real eps = std::numeric_limits<real>::epsilon())const
    {
        real l1l2 = s.ort().dot(ort());
    }
};

real distance2_to_line(const Point3D& x, const Segment3D& l)
{
    Vector3D c = x; c -= l.begin();
    real b = c.dot(l.ort());
    return c.norm2() - b * b;
}

real distance_to_line(const Point3D& x, const Segment3D& l)
{
    return sqrt(distance2_to_line(x, l));
}

std::pair<std::pair<real,real>, bool> crest(
    const Vector3D x10, const Vector3D& l1, 
    const Vector3D x20, const Vector3D& l2, 
    real eps = 10 * std::numeric_limits<real>::epsilon())
{
    std::pair<std::pair<real,real>, bool> result;
    real l1l2 = l1.dot(l2);
    real d = 1.0 - l1l2 * l1l2;
    if(fabs(d) < eps)
    {
        result.second = false;
    }
    else
    {
        Vector3D dx = x20; dx -= x10;
        real l1dx = l1.dot(dx);
        real l2dx = l2.dot(dx);
        result.first.first =  (l1dx - l1l2 * l2dx) / d;
        result.first.second = (-l2dx + l1l2 * l1dx) / d;
        result.second = true;
    }
    return result;
}

inline real segments_distance(const Segment3D& S1, const Segment3D& S2)
{
    auto cr = crest(S1.begin(),S1.ort(),S2.begin(),S2.ort());
    if(cr.second) //Lines contaning the segments are crossing
    {
        real x1 = cr.first.first;
        real x2 = cr.first.second;
        if(x1 >=0 && x1<=S1.lenght() && x2 >=0 && x2 <= S2.lenght()) //The crest solution lies inside the segmrnts
        {
            //Return the distance btw parallel planes containing the segments
            auto &u = S1.ort();
            auto &v = S2.ort();
            auto w = u.vector_product(v);
            Vector3D dx(S2.begin()); dx-=S1.begin();
            real result = fabs(w.dot(dx) / w.norm());
            return result;
        }
        else //The crest solution lies outside the segments
        {
            if(x1 < 0)
            {
                if(x2 < 0)
                {
                    return distance(S1.begin(),S2.begin());
                }
                else if(x2 < S2.lenght())
                {
                    return distance_to_line(S1.begin(),S2);
                }
                else
                {
                    return distance(S1.begin(), S2.end());
                }
            }
            else if(x1 < S1.lenght())
            {
                if(x2 < 0)
                {
                    return distance_to_line(S2.begin(),S1);
                }
                else
                {
                    return distance_to_line(S2.end(),S1);
                }
            }
            else
            {
                if(x2 < 0)
                {
                    return distance(S1.end(),S2.begin());
                }
                else if(x2 < S2.lenght())
                {
                    return distance_to_line(S1.end(),S2);
                }
                else
                {
                    return distance(S1.end(),S2.end());
                }
            }

        }
    }
    else //Lines are parallel
    {
        int sb = S1.is_projection_out(S1.projection(S2.begin()));
        int se = S1.is_projection_out(S1.projection(S2.end()));
        if(sb * se <= 0)
        {
            //Calculate and return the distance between lines
            Vector3D c(S2.begin()); c-=S1.begin();
            real a = c.dot(S2.ort());
            return std::sqrt(c.norm2() - a * a);
        }
        else if(sb < 0)
        {
            return std::min(distance(S1.begin(),S2.begin()),distance(S1.begin(),S2.end()));
        }
        else
        {
            return std::min(distance(S1.end(),S2.begin()),distance(S1.end(),S2.end()));
        }
    }
}

std::ostream& operator << (std::ostream& s, const Point3D& x)
{
    s << "(" << x[0] << ", " << x[1] << ", " << x[2] << ")";
    return s;
}

std::ostream& operator << (std::ostream& s, const Segment3D& x)
{
    s << "{" << x.begin() << " -- " << x.end() << "}";
    return s;
}

