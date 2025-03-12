#pragma once

#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>
#include <initializer_list>
#include <iostream>

const int SpaceDim = 3;
using real = double;


//! Clacc Point 3D represents the point in 3D space
class Point3D
{
protected:
    real X[SpaceDim]; /// The array with coordinates (x,y,z)
public:
    /// @brief Constructs the Point3D object by given coordinates
    /// @param x X coordinate
    /// @param y Y coordinate
    /// @param z Z coordinate
    Point3D(real x, real y, real z)
    {
        X[0] = x; X[1] = y; X[2] = z;
    }
    /// @brief Default constructor constructs the null point
    Point3D(){memset(X,0,sizeof(X));}
    /// @brief Constructs the Point3D object by coordinates given in initializer list
    /// @param a Initializer list of size = 3 containing te coordinates.
    Point3D(std::initializer_list<real> a)
    {
        assert(a.size() == SpaceDim && "Point3D can be constructed only from 3 coordinates");
        auto _a = std::data(a);
        X[0] = _a[0]; X[1] = _a[1]; X[2] = _a[2];
    }
    /// @brief Constructs the Point3D object by coordinates given in data array pointed with d
    /// @param d Pointer at data array containing point coordinates
    Point3D(const real *d)
    {
        memcpy(X,d,sizeof(X));
    }
    /// @brief Extracts the i-th coordinate value
    /// @param i The number of the coordinate axis (0 is x, 1 is y, 2 is z)
    /// @return The desired coordinate value
    real operator[](int i) const
    {
        return X[i];
    }
    /// @brief Extracts the i-th coordinate value to be modified
    /// @param i The number of the coordinate axis (0 is x, 1 is y, 2 is z)
    /// @return The reference to desired coordinate value
    real& operator[](int i)
    {
        return X[i];
    }
};

/// @brief Calculates the squared euclidian distance between two given points
/// @param a The first point
/// @param b The second point
/// @return The squared euclidian distance between two given points value
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

/// @brief Calculates the euclidian distance between two given points
/// @param a The first point
/// @param b The second point
/// @return The euclidian distance between two given points value
real distance(const Point3D& a, const Point3D& b)
{
    return std::sqrt(distance2(a,b));
}

/// @brief Vector3D class represents the radius-vector of a point
class Vector3D: public Point3D
{
public:
    /// @brief Constructs the radius-vector of a given point
    /// @param x The point that the vector represents
    Vector3D(const Point3D& x): Point3D(x) {}
    /// @brief Constructs the Vector3D object with the given coordinates
    /// @param x X coordinate
    /// @param y Y coordinate
    /// @param z Z coordinate
    Vector3D(real x, real y, real z): Point3D(x,y,z) {}
    /// @brief Default constructor constructs the null vector
    Vector3D(): Point3D() {}
    /// @brief Calculates the inner product of this vector and the given one
    /// @param v The vector to calculate the inner product with
    /// @return The value of the inner product
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
    /// @brief Calculates the vector product of this vector and the given one
    /// @param v The vector to calculate the vector product with
    /// @return The vector representing the vector product
    Vector3D vector_product(const Vector3D& u) const
    {
        const Vector3D& v(*this);
        return Vector3D(v[1] * u[2] - v[2] * u[1], 
                      -(v[0] * u[2] - v[2] * u[0]), 
                        v[0] * u[1] - v[1] * u[0]);
    }
    /// @brief Calculates the dot of this with itself
    /// @return The squared norm of this vector
    real norm2() const
    {
        return dot(*this);
    }
    /// @brief Calculates the norm of this vector
    /// @return The norm of this vector
    real norm() const
    {
        return std::sqrt(norm2());
    }
    /// @brief Modifies this vector by multiplying each coordinate on the given number
    /// @param a The given number
    /// @return This vector modified by multiplication on the given number
    Vector3D& operator *= (real a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] *= a;
        }
        return *this;
    }
    /// @brief Modifies this vector by dividing each coordinate on the given number
    /// @param a The given number
    /// @return This vector modified by division on the given number
    Vector3D& operator /= (real a)
    {
        assert(a != 0 && "Trying to divide by zero!");
        return (*this) *= 1.0 / a;
    }
    /// @brief Modifies this vector by adding another vector
    /// @param a Added vector
    /// @return Modified vector
    Vector3D& operator += (const Vector3D& a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] += a[i];
        }
        return *this;
    }
    /// @brief Modifies this vector by substraction another vector
    /// @param a Vector to be substracted
    /// @return Modified vector
    Vector3D& operator -= (const Vector3D& a)
    {
        for(int i = 0; i < SpaceDim; ++i)
        {
            (*this)[i] -= a[i];
        }
        return *this;
    }
    /// @brief Modifies this vector by reversing
    /// @return Modified vector
    Vector3D operator - ()
    {
        Vector3D r(*this);
        r *= -1;
        return r;
    }
};

/// @brief Class segment represents the line segment between the given points
class Segment3D
{
protected:
    /// @brief A storage of pair of points
    Vector3D X[2];
    /// @brief Vector representing the direction of the line
    Vector3D l;
    /// @brief The segment lenght value
    real S;
public:
    /// @brief Constructs the segment by a pair of given points
    /// @param x0 The begin point
    /// @param x1 The end point
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
    /// @brief Gets the point at segment begin
    /// @return Constant referens of the point
    const Vector3D& begin()const{return X[0];}
    /// @brief Gets the point at segment end
    /// @return Constant referens of the point
    const Vector3D& end()const{return X[1];}
    /// @brief Gets the direction vector const reference
    /// @return Constant referens of the point
    const Vector3D& ort()const{return l;}
    /// @brief Returns the segment lenght value
    /// @return Segment lenght value
    real lenght() const {return S;}
    /// @brief Calculates the projection of the given point onto cegment line
    /// @param x The point that projection to be calculated
    /// @return The distance from the begin point to the projection along the segment line direction
    real projection(const Point3D& x)const
    {
        return (Vector3D(x)-=begin()).dot(ort());
    }
    /// @brief Check if the projection given by the distance from the segment beginning poing lies inside or outside the segment
    /// @param p The Given distance from the segment beginning pointing the point on the segment line
    /// @return -1 if p < 0; 0 if p is owned by segment; 1 if point lies to the right of the segment ending point
    int is_projection_out(real p)const
    {
        if(p < 0) return -1;
        else if(p < lenght()) return 0;
        else return 1;
    }
};

/// @brief Calculates the squared distance from the given point to the given line
/// @param x The given point
/// @param l The segment representinting the line
/// @return The squared distance value
real distance2_to_line(const Point3D& x, const Segment3D& l)
{
    Vector3D c = x; c -= l.begin();
    real b = c.dot(l.ort());
    return c.norm2() - b * b;
}

/// @brief Calculates the distance from the given point to the given line
/// @param x The given point
/// @param l The segment representinting the line
/// @return The distance value
real distance_to_line(const Point3D& x, const Segment3D& l)
{
    return sqrt(distance2_to_line(x, l));
}

/// @brief Calculates the points on the two lines with a minimal distance from each other
/// @param x10 The point on the first line beginning
/// @param l1 The Direction vector of the first line
/// @param x20 The point on the second line beginning
/// @param l2 The Direction vector of the second line
/// @param eps The tolerance for comparing the floating numbers
/// @return The Pair. The ferst member is a pair of the desired points coordinates in lical line coordinate system. The second member is the sign if the given lines are skew
std::pair<std::pair<real,real>, bool> crest(
    const Vector3D& x10, const Vector3D& l1, 
    const Vector3D& x20, const Vector3D& l2, 
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

/// @brief Calculates the distance between the two given 3D segmrnts
/// @param S1 The first segment
/// @param S2 The second segment
/// @return The distance value
/** 
 Algorithm description:
 The parametric equation of a line is expressed as
 p = p0 + t*l
 where p0 is a given point, l is a line direction, and t is a coordinate of the point p along the line
 For the line from point p0 to p1 it can be expressed as
 p = p0 + t * (p1 - p0) / |p1 - p0|
 So the segment of the line between p0 and p1 could be defined via t limiting from 0 to |p1 - p0|.
 Minimizing the distance between the two points owned by the two non-intersecting 3D lines p and q yields the solution 
 for the coordinates of the closest points:
         tp_m =  ((lp,q0-p0) - (lp,lq) * (lq,q0-p0)) / d
         tq_m = (-(lq,q0-p0) + (lp,lq) * (lp,q0-p0)) / d
         d = 1 - (lp,lq)^2
 So the pair of (tp_m, tq_m) may be located in several positions relative to the rectangle 
 on the parametric plane which could be defined for the segments' points coordinates on the lines.
 The distance between segments is calculated from that location as it is implemented in the source code.
 */
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

