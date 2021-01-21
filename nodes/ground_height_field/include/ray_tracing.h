/*
 *  Copyright (c) 2018, Nagoya University
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of Autoware nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef RAY_TRACING_H
#define RAY_TRACING_H

//STD
#include <cmath> 
#include <vector>
#include <iostream>
#include <limits>
#include <type_traits>

//Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

//Boost
#include <boost/utility/enable_if.hpp>


#if !defined(DEBUG)
//#define DEBUG 0
#endif

namespace ground_field {

/**
    \brief Implements a single 3D point
*/
template<typename DataType, int dims = 3>
struct Point {
    //the actual coordinates
    typedef Eigen::Matrix<DataType,dims,1> Vector_t;
    enum { NeedsToAlign = (sizeof(Vector_t)%16)==0 };    
    Vector_t data_;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
    /**
        \brief Default <0,0,0> point
    */
    Point() 
    {
        clear();
    }
    /**
        \brief Point at <x,y>
    */
    // template < class T = DataType,
    //            typename boost::enable_if_c< dims == 2, T >::type = 0>
    Point(DataType x, DataType y)
    {
        init(x,y);
        static_assert(dims == 2, "XY constructor only usable in 2D");
    }
    /**
        \brief Point at <x,y,z>
    */
    // template < class T = DataType,
    //            typename boost::enable_if_c< dims == 3, T >::type = 0>
    Point(DataType x, DataType y, DataType z)
    {
        init(x,y,z);
        static_assert(dims == 3, "XYZ constructor only usable in 3D");
    }
    /**
        \brief Copy constructor
        \param[in] other is the point to copy from
    */
    Point(const Point<DataType,dims>& other) 
    {
        data_ = other.data_;
    }
    /**
        \brief Copy constructor
        \param[in] other is the point to copy from
    */
    template <typename Y, int dimensions = dims>
    Point(const Point<Y,dimensions>& other) 
    {
        static_assert(dimensions == dims, "Point dimensions must be equal");    
        for (int i=0; i<data_.rows(); i++) {
            data_(i) = static_cast<DataType>(other.data_(i));
        }
    }
    
    /**
        \brief Copy operator
        \param[in] other is the point to copy from
    */
    Point<DataType,dims>& operator=(const Point<DataType,dims>& other)
    {
        // check for self-assignment
        if(&other != this) {
            data_ = other.data_;
        }
        return *this;
    }
    /**
        \brief Copy operator from other point class
        \param[in] other is the point to copy from
    */
    template <typename Y, int dimensions = dims>
    Point<DataType,dims>& operator=(const Point<Y,dimensions>& other)
    {
        for (int i=0; i<data_.rows(); i++) {
            data_(i) = static_cast<DataType>(other.data_(i));
        }
        return *this;
    }
    /**
        \brief Sum operator
        \param[in] other is the point to add to
    */
    Point<DataType,dims>& operator+=(const Point<DataType,dims>& other)
    {
        data_ += other.data_;
        return *this;
    }
    /**
        \brief Subtract operator
        \param[in] other is the point to subtract from
    */
    Point<DataType,dims>& operator-=(const Point<DataType,dims>& other)
    {
        data_ -= other.data_;
        return *this;
    }
    /**
        \brief Sum operator
        \param[in] other is the point to add to
    */
    friend Point<DataType,dims> operator+(Point<DataType,dims> lhs, const Point<DataType,dims> rhs)
    {
        lhs += rhs;
        return lhs;
    }
    /**
        \brief Subtract operator
        \param[in] other is the point to subtract from
    */
    friend Point<DataType,dims> operator-(Point<DataType,dims> lhs, const Point<DataType,dims> rhs)
    {
        lhs -= rhs;
        return lhs;
    }
    /**
        \brief Scalar multiplication and assigment
    */
    Point<DataType,dims>& operator*=(const DataType& s) {
        data_ *= s;
        return *this;
    }
    /**
        \brief Scalar multiplication
        \retval the product of the point p and scalar s
    */
    friend Point<DataType,dims> operator*(Point<DataType,dims> p, DataType s) 
    {
        p *= s;
        return p;
    }
    /**
        \brief Dot product
        \retval the product of this point and the other point
    */
    DataType operator*(const Point<DataType,dims>& other) const
    {
        return data_.dot(other.data_);
    }
    /**
        \brief Vector inner product and assigment
    */
    Point<DataType,dims>& operator%=(const Point<DataType,dims>& rhs) {
        data_.array() *= rhs.data_.array();
        return *this;
    }
    /**
        \brief Vector inner product
        \retval the .* (inner product) of this point and the other point
    */
    friend Point<DataType,dims> operator%(Point<DataType,dims> lhs, const Point<DataType,dims> rhs)
    {
        lhs %= rhs;
        return lhs;
    }
    /**
        \brief Scalar division and assigment
    */
    Point<DataType,dims>& operator/=(const DataType& s) {
        data_ /= s;
        return *this;
    }
    /**
        \brief Scalar division
        \retval the quotient of lhs point and scalar s
    */
    friend Point<DataType,dims> operator/(Point<DataType,dims> lhs, DataType s)
    {
        lhs /= s;
        return lhs;
    }
    /**
        \brief Vector inner division and assignment
        \retval the ./ (inner division) of this point and the other point
    */
    Point<DataType,dims>& operator/=(const Point<DataType,dims>& rhs) {
        data_.array() /= rhs.data_.array();
        return *this;
    }
    /**
        \brief Vector point division
        \retval the ./ (division) of lhs point and rhs point
    */
    friend Point<DataType,dims> operator/(Point<DataType,dims> lhs, const Point<DataType,dims> rhs) 
    {
        lhs /= rhs;
        return lhs;
    }

    /**
        \brief Vector cross product and assignment
        \retval the ./ (inner division) of this point and the other point
    */
    Point<DataType,dims>& operator*=(const Point<DataType,dims>& rhs) {
        data_ *= rhs.data_;
        return *this;
    }
    /**
        \brief Cross-product
        \retval the cross product of lhs point and rhs point
    */
    friend Point<DataType,dims> operator*(Point<DataType,dims> lhs, const Point<DataType,dims> rhs)
    {
        lhs *= rhs;
        return lhs;
    }
    /**
        \brief Equal comparison
        \retval true if both points are equal
    */
    friend bool operator==(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return (lhs.data_.isApprox(rhs.data_));
    }
    /**
        \brief Different comparison
        \retval true if both points are different
    */
    friend bool operator!=(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return !(lhs == rhs);
    }
    /**
        \brief Less-than comparison
        \retval true if this point is less than the other point
    */
    friend bool operator<(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return (lhs.data_.array() < rhs.data_.array()).any();
    }
    /**
        \brief Less-than or equal to comparison
        \retval true if this point is less than or equal to the other point
    */
    friend bool operator<=(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return !(lhs > rhs);
    }
    /**
        \brief Greater-than comparison
        \retval true if this point is greater than the other point
    */
    friend bool operator>(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return (rhs < lhs);
    }
    /**
        \brief Greater-than or equal to comparison
        \retval true if this point is greater than or equal to the other point
    */
    friend bool operator>=(const Point<DataType,dims>& lhs, const Point<DataType,dims>& rhs) 
    {
        return !(lhs < rhs);
    }
    /**
        \brief Sets the X coordinate of this point
        \param[in] x is the new value
    */
    inline void setX(const DataType x) 
    {
        data_[0] = x;
    }
    /**
        \brief Sets the Y coordinate of this point
        \param[in] y is the new value
    */
    inline void setY(const DataType y)
    {
        data_[1] = y;
    }
    /**
        \brief Sets the Z coordinate of this point
        \param[in] z is the new value
    */
    inline void setZ(const DataType z)
    {
        data_[2] = z;
    }
    /**
        \brief X accessor
        \retval the current value of x coordinate
    */
    DataType x() const
    {
        return data_[0];
    }
    /**
        \brief Y accessor
        \retval the current value of y coordinate
    */            
    DataType y() const
    {
        return data_[1];
    }
    /**
        \brief Z accessor
        \retval the current value of Z coordinate
    */
    DataType z() const
    {
        return data_[2];
    }
    /**
        \brief X accessor
        \retval the current value of x coordinate
    */
    DataType& x()
    {
        return data_[0];
    }
    /**
        \brief Y accessor
        \retval the current value of y coordinate
    */            
    DataType& y()
    {
        return data_[1];
    }
    /**
        \brief Z accessor
        \retval the current value of Z coordinate
    */
    DataType& z()
    {
        return data_[2];
    }
    /**
        \brief Initializer
    */
    void clear()
    {
        data_.setConstant(0);
    }
    void init(DataType x, DataType y)
    {
        data_ << x,y;
    }
    void init(DataType x, DataType y, DataType z)
    {
        data_ << x,y,z;
    }
    /**
        \brief Returns the length of this point
    */
    inline DataType length()
    {
        return static_cast<DataType>(data_.norm());
    }
    /**
        \brief Returns the maximum of two points p1 and p2
    */
    static const Point<DataType,dims>& max(const Point<DataType,dims>& p1, const Point<DataType,dims>& p2)
    {
        if (p1 > p2) {
            return p1;
        } else {
            return p2;
        }
    }
    /**
        \brief Returns the maximum value among the three axes of this point
    */
    DataType max() const 
    {
        return data_.array().max();
    }

    /**
        \brief Returns the minimum of two points p1 and p2
    */
    static const Point<DataType,dims>& min(const Point<DataType,dims>& p1, const Point<DataType,dims>& p2)
    {
        if (p1 < p2) {
            return p1;
        } else {
            return p2;
        }
    }
    /**
        \brief Returns the minimum value among the three axes of this point
    */
    DataType min() const 
    {
        return data_.array().min();
    }

    /**
        \brief Gets the normalized value of this point (unitary)
    */
    static const Point<DataType,dims>& norm(const Point<DataType,dims>& p)
    {
        Point<DataType,dims> q = p;
        q.data_.normalize();
        return q;
    }
    /**
        \brief Gets the normalized value of this point (unitary)
    */
    Point<DataType,dims>& norm() const
    {
        Point<DataType,dims> p = *this;
        p.data_.normalize();
        return p;
    }

    /**
        \brief Output stream operator
    */
    friend std::ostream& operator<<(std::ostream& out, const Point<DataType,dims>& p)
    {
    	//Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    	Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ",", ",", "", "", "[", "]");
        return out << p.data_.format(OctaveFmt);
    }
    
    /**
        \brief Returns a point which components are rounded to upper integer
        using ceil of the given point
    */
    static const Point<int,dims> ceil(const Point<DataType,dims>& p) 
    {
        Point<int,dims> q;
        q.data_.array() = p.data_.array().ceil().template cast < int >();
        return q;
    }
    /**
        \brief Returns a point which components are rounded to upper integer
        using ceil of this point
    */
    Point<int,dims> ceil() const
    {
        Point<int,dims> q;
        q.data_.array() = data_.array().ceil().template cast < int >();
        return q;
    }
    /**
        \brief Returns a point which components are the absolute value
        of the given point
    */
    static const Point<DataType,dims> abs(const Point<DataType,dims>& p)
    {
        Point<DataType,dims> q;
        q.data_.array() = p.data_.array().abs();
        return q;
    }
    /**
        \brief Returns a point which components are the absolute value
        of this point
    */
    Point<DataType,dims> abs() const
    {
        Point<DataType,dims> q;
        q.data_.array() = data_.array().abs();
        return q;
    }
    /**
        \brief Returns a point which components are fixed from INFINITE and NAN
        errors, using the given point
    */
    static Point<DataType,dims> fix(const Point<DataType,dims>& p)
    {
        Point<DataType,dims> q;
        q.data_.array() = p.data_.array().isFinite().select(p.data_, 0);
        return q;
    }
    /**
        \brief Fixes this point from INFINITE and NAN errors
    */
    void fix()
    {
        data_.array().isFinite().select(data_, 0);
    }
    /**
        \brief Returns a point which components are the inverse 
        of the given point
    */
    static Point<DataType,dims> inverse(const Point<DataType,dims>& p)
    {
        Point<DataType,dims> q;
        q.data_.array() = p.data_.array().inverse();
        return q;
    }
    /**
        \brief Returns a point which components are the inverse 
        of this point
    */
    Point<DataType,dims> inverse() const
    {
        Point<DataType,dims> q;
        q.data_.array() = data_.array().inverse();
        return q;
    }
    /**
        \brief Returns a point which components are the sign value
        (0 for positive, 1 for negative) of the given point
    */
    static Point<DataType,dims> sign(const Point<DataType,dims>& p)
    {
        Point<DataType,dims> q;
        q.data_.array() = (p.data_.array() < 0).select(Vector_t::Constant(dims,1),Vector_t::Constant(dims,0));
        return q;
    }
    /**
        \brief Returns a point which components are the sign value
        (0 for positive, 1 for negative) of this point
    */
    Point<DataType,dims> sign() const
    {
        Point<DataType,dims> q;
        q.data_.array() = (data_.array() < 0).select(Vector_t::Constant(dims,1),Vector_t::Constant(dims,0));
        return q;
    }
    /**
        \brief Returns a point which components are the sign value
        (1 for positive, -1 for negative, 0 for zero) of this point
    */
    Point<DataType,dims> sign2() const
    {
        Point<DataType,dims> q;
        q.data_.array() = data_.array().sign();
        return q;
    }
    /**
        \brief Typecasting for point elements
    */
    template <typename DataTypeOther, int dimensions = dims>
    Point<DataTypeOther,dimensions> cast() const
    {
        Point<DataTypeOther,dimensions> q;
        for (int i=0; i<q.data_.rows(); i++) {
            q.data_(i) = static_cast<DataTypeOther>(data_(i));
        }
        return q;
    }
    /**
        \brief Ascending sorting of point elements
    */
    Point<DataType,dims> sort() const
    {
        Point<DataType,dims> q;
        q.data_.array() = data_.array();
        std::sort(q.data_.data(), q.data_.data()+q.data_.size());
        return q;
    }
};

template<typename DataType>
using Point2D = Point<DataType,2>;

template<typename DataType>
using Point3D = Point<DataType,3>;


/**
    \brief Implements a 3D ray (line segment)

    Includes extra attributs to support A.Williams ray-box intersection algorithm
*/
template<typename DataType>
struct Ray {
    /**
        \brief Default empty ray
    */
    Ray()
    {
        origin_.clear();
        end_.clear();
        direction_.clear();
        inv_direction_.clear();
        sign_.clear();
    }
    /**
        \brief Ray using given points
        \param[in] p1 is the ray origin point
        \param[in] p2 is the ray end point
    */            
    Ray(Point3D<DataType>& p1, Point3D<DataType>& p2)
    {
        origin_ = p1;
        end_ = p2;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
    }
    /**
        \brief Ray using given points
        \param[in] p1 is the ray origin point
        \param[in] p2 is the ray end point
    */            
    Ray(Point3D<DataType> p1, Point3D<DataType> p2)
    {
        origin_ = p1;
        end_ = p2;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
    }
    /**
        \brief Copy constructor
        \param[in] other is the ray to copy from
    */
    Ray(const Ray<DataType>& other)
    {
        origin_ = other.origin_;
        end_ = other.end_;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
    }
    /**
        \brief Copy operator
        \param[in] other is the ray to copy from
    */
    Ray<DataType>& operator=(const Ray<DataType>& other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;
        origin_ = other.origin_;
        end_ = other.end_;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
        return *this;             
    }
    /**
        \brief Sets the ray origin
        \param[in] p is the new origin point
    */
    void setOrigin(Point3D<DataType> p)
    {
        origin_ = p;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
    }
    /**
        \brief Sets the ray end
        \param[in] p is the new end point
    */
    void setEnd(Point3D<DataType> p)
    {
        end_ = p;
        direction_ = end_ - origin_;
        inv_direction_ = direction_.inverse();
        sign_ = inv_direction_.sign();
    }
    /**
        \brief Gets the ray origin
        \retval this ray's current origin
    */
    inline Point3D<DataType> origin() const { return origin_; }
    /**
        \brief Gets the ray end
        \retval this ray's current end
    */
    inline Point3D<DataType> end() const { return end_; }
    /**
        \brief Gets the ray direction
        \retval this ray's current direction
    */
    inline Point3D<DataType> direction() const { return direction_; }
    /**
        \brief Gets the ray inverse direction
        \retval this ray's current inverse direction
    */
    inline Point3D<DataType> inv_direction() const { return inv_direction_; }

    /**
        \brief Gets the ray direction sign
        \retval this ray's current direction sign
    */
    inline Point3D<int> direction_sign() const { return sign_; }
    /**
        \brief Computes the Euclidean length of the ray
        \retval the length of the ray
    */
    inline DataType length() const 
    {
        Point3D<DataType> diff = end_ - origin_;
        return diff.length();
    }

    /**
        \brief Output stream operator
    */
    friend std::ostream& operator<<(std::ostream& out, const Ray<DataType>& ray)
    {
        return out << "origin " << ray.origin_ << ", end " << ray.end_;
    }

    Point3D<DataType> origin_;
    Point3D<DataType> end_;
    Point3D<DataType> direction_;
    Point3D<DataType> inv_direction_;
    Point3D<int> sign_;
};

/**
    \brief Defines a grid cell of a 2D grid map
*/
template<typename DataType>
class GridCell {
    /**
        \brief Default constructor
    */
    GridCell()
    {
        height_ = std::numeric_limits<DataType>::infinity();
        num_points_ = 0;
        end_point_ = false;
        mean_ = 0;
        variance_ = 0;
    }
    /**
        \brief Copy constructor
        \param[in] other is the grid cell to copy from
    */
    GridCell(const GridCell<DataType>& other) 
    {
        height_ = other.height_;
        num_points_ = other.num_points_;
        end_point_ = other.end_point_;
        mean_ = other.mean_;
        variance_ = other.variance_;
    }
    /**
        \brief Copy constructor
        \param[in] other is the grid cell to copy from
    */
    template <typename Y>
    GridCell(const GridCell<Y>& other) 
    {
        height_ = static_cast<DataType>(other.height_);
        num_points_ = other.num_points_;
        end_point_ = other.end_point_;
        mean_ = static_cast<DataType>(other.mean_);
        variance_ = static_cast<DataType>(other.variance_);
    }
    
    /**
        \brief Copy operator
        \param[in] other is the grid cell to copy from
    */
    GridCell<DataType>& operator=(const GridCell<DataType>& other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;
        height_ = other.height_;
        num_points_ = other.num_points_;
        end_point_ = other.end_point_;
        mean_ = other.mean_;
        variance_ = other.variance_;
        return *this;
    }
    /**
        \brief Gets the current height value
    */
    DataType height() const
    {
        return height_;
    }
    bool isEndPoint() const
    {
        return end_point_;
    }
    void setEndPoint(bool end_point)
    {
        end_point_ = end_point;
        // if (end_point_) {
        //
        // }
    }
    /**
        \brief Adds a height value to this cell
        \param[in] height is the new height value to add, if less than the current height value then replaces it
        \param[in] end_point tells whether this cell is a ray end-point (true) or a simple intersection cell (cell)
    */ 
    void add(DataType height, bool end_point=false)
    {
        if (height < height_) {
            height_ = height;
        }
        // if (end_point || end_point_) {
        //     end_point_ = true;
        //     num_points_++;
        // }
    }
    /**
        \brief Height assignment operator
        \param[in] height is the new height value to add, if less than the current height value then replaces it
    */
    GridCell<DataType>& operator=(const DataType height)
    {
        add(height, end_point_);
        return *this;
    }


    DataType height_;    //! this cell's height value (keeping lowest height)
    short int num_points_; //! number of points in this cell
    bool end_point_;  //! tells whether this cell is the ray's end-point (true) of just an intersecting cell (false)
    DataType mean_;  //! mean value for end-point cells
    DataType variance_; //! variance value for end-point cells
};

/**
   \brief Computation of ray intersection with several voxels on a grid
   (not a complete raytracer implementation)
   This implementation is based on:
          A. Williams et al., "An Efficient and Robust Ray–Box Intersection Algorithm"
          http://www.cs.utah.edu/~awilliam/box/box.pdf
   This class is selfcontained to reduce dependencies
 **/
template<typename DataType>
class RayTracing {
    public:
        typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> GridMap2D; 
        /**
            \brief Structure for data comparison 
        */
        struct CompareZ {
            /**
                \brief Compare two points, if X and Y values are equal return the point with minimum
                Z axis value, otherwise return p2 assuming it is the newest point
            */
            template<typename C>
            auto operator()(const Point3D<C>& p1, const Point3D<C>& p2)
            {
                if (p1.x() == p2.x() && p1.y() == p2.y()) {
                    return Point3D<C>(p1.x(), 
                                      p1.y(), 
                                      std::min(p1.z(), p2.z()));
                } else {

                    return p2;
                }
            }
        };

    public:
        /**
          \brief Default constructor
        */
        RayTracing();
        /**
          \brief Defines a grid with the given parameters
          \param grid_min_bounds defines the minimum bound coordinates of this grid
          \param grid_max_bounds defines the maximum bound coordinates of this grid
          \param cell_size_x defines the size of a cell (voxel) in X axis
          \param cell_size_y defines the size of a cell (voxel) in Y axis
          \param cell_size_z defines the size of a cell (voxel) in Z axis
        */
        RayTracing(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds, DataType cell_size_x, DataType cell_size_y, DataType cell_size_z);
        /**
          \brief Defines a grid with the given parameters
          \param grid_min_bounds defines the minimum bound coordinates of this grid
          \param grid_max_bounds defines the maximum bound coordinates of this grid
          \param cell_size defines the size of a cell (voxel) in each axis (X,Y,Z)
        */
        RayTracing(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds, Point3D<DataType> cell_size);
        /**
          \brief Default destructor
        */        
        ~RayTracing();

        /**
          \brief Sets the grid bounding limits
          \param grid_min_bounds defines the minimum bound coordinates of this grid
          \param grid_max_bounds defines the maximum bound coordinates of this grid  
        */
        void setGridBounds(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds);

        /**
          \brief Sets the cell (voxel) size
          \param cell_size defines the size of a cell (voxel) in each axis (X,Y,Z)
        */
        void setCellSize(Point3D<DataType> cell_size);

        /**
            \brief Adds a new ray to the list of rays for intersection computation
            \param ray is the new segment (ray)
        */
        void addRay(const Ray<DataType>& ray);

        /**
            \brief Gets the ray at the given index idx
            \param idx is the index (number) of the ray to extract
        */
        Ray<DataType>& getRay(const unsigned int idx);

        /**
            \brief Gets number of rays
        */
        unsigned int raysCount() const
        {
            return rays_.size();
        }

        /**
            \brief Remove all the rays and corresponding voxel data
        */
        void clearRays()
        {
            rays_.clear();
            //voxel_coordinates_.clear();
        }
        /**
            \brief Remove all voxel data
        */
        void clearVoxels()
        {
            for (auto v : voxel_coordinates_) {
                v.clear();
            }
        }

        /**
            \brief Computes the intersections of a ray with all voxels on the grid
            and return the list of voxels, empty if no intersection exists.
            \param idx is the index (number) of the ray to use for intersections
        */
        void computeRayIntersections(const unsigned int idx);
        /**
            \brief Overloaded version, assumes we want to use the first ray
        */
        void computeRayIntersections()
        {
            unsigned int idx = 0;
            computeRayIntersections(idx);
        }
        /**
            \brief Overloaded version, assumes we want to use the given ray 
            \param[in] ray is the segment (ray) to use for intersection computation
        */
        void computeRayIntersections(const Ray<DataType>& ray)
        {
            addRay(ray);
            unsigned int idx = raysCount() - 1;
            computeRayIntersections(idx);
        }

        /**
            \brief Computes the intersections of a ray with all voxels on the grid
            filtering cells with same <X,Y> coordinates keeping the lowest Z coordinate cell
            and return the list of voxels, empty if no intersection exists.
            \param[in] idx is the index (number) of the ray to use for intersections
        */
        void computeRayIntersectionsZFilter(const unsigned int idx);
        /**
            \brief Overloaded version, assumes we want to use the first ray
        */
        void computeRayIntersectionsZFilter()
        {
            unsigned int idx = 0;
            computeRayIntersectionsZFilter(idx);
        }
        /**
            \brief Overloaded version, assumes we want to use the given ray 
            \param[in] ray is the segment (ray) to use for intersection computation
        */
        void computeRayIntersectionsZFilter(const Ray<DataType>& ray)
        {
            addRay(ray);
            unsigned int idx = raysCount() - 1;
            computeRayIntersectionsZFilter(idx);
        }
        /**
            \brief Computes ray intersections and stores results on 2D grid map (using Point2D)
            \param[in] ray is the segment (ray) to use for intersection computation
        */
        void computeRayIntersectionsZFilter2DGrid(const Ray<DataType>& ray, DataType base_height = 0);
        /**
            \brief Overloaded version, assumes we want to use the given ray index
            \param[in] idx is the index (number) of the ray to use for intersections
            \param[in,optional] base_height is the basic height value to add to height computations
        */
        void computeRayIntersectionsZFilter2DGrid(const unsigned int idx, DataType base_height = 0)
        {
        	computeRayIntersectionsZFilter2DGrid(getRay(idx), base_height);
        }
        /**
			\brief Computes ray intersections and stores results on 2D grid map using Extra-Fast Line Algorithm
			\param[in] ray is the segment (ray) to use for intersection computation
			\param[in,optional] base_height is the basic height value to add to height computations
		*/
        void computeRayIntersectionsZFilter2DGridEFLA(const Ray<DataType>& ray, DataType base_height = 0,
													  DataType min_height = -1*std::numeric_limits<DataType>::infinity(),
													  DataType max_height = std::numeric_limits<DataType>::infinity());
        /**
			\brief Overloaded version, assumes we want to use the given ray index
			\param[in] idx is the index (number) of the ray to use for intersections
			\param[in,optional] base_height is the basic height value to add to height computations
			\param[in,optional] min_height is minimum height limit
			\param[in,optional] max_height is maximum height limit
		*/
        void computeRayIntersectionsZFilter2DGridEFLA(const unsigned int idx, DataType base_height = 0,
        											  DataType min_height = -1*std::numeric_limits<DataType>::infinity(),
													  DataType max_height = std::numeric_limits<DataType>::infinity())
        {
			computeRayIntersectionsZFilter2DGridEFLA(getRay(idx), base_height, min_height, max_height);
        }
        /**
			\brief Computes ray intersections and stores results on 2D grid map (using Point3D)
			\param[in] ray is the segment (ray) to use for intersection computation
			\param[in,optional] base_height is the basic height value to add to height computations
			\param[in,optional] min_height is minimum height limit
			\param[in,optional] max_height is maximum height limit
		*/
		void computeRayIntersectionsZFilter2DGridPoint3D(const Ray<DataType>& ray, DataType base_height = 0);
        /**
			\brief Overloaded version, assumes we want to use the given ray index
			\param[in] idx is the index (number) of the ray to use for intersections
			\param[in,optional] base_height is the basic height value to add to height computations
		*/
		void computeRayIntersectionsZFilter2DGridPoint3D(const unsigned int idx, DataType base_height = 0)
		{
			computeRayIntersectionsZFilter2DGridPoint3D(getRay(idx), base_height);
		}
		/**
		    \brief Gets the current grid map
		 */
        GridMap2D& getGridMap()
        {
            return grid_map_;
        }
		/**
		    \brief Clears the grid map
		 */
        void clearGridMap()
        {
            grid_map_.setConstant(std::numeric_limits<DataType>::infinity());
        }
		/**
		    \brief Creates a grid map of the given size
		 */
        void makeGridMap(unsigned int rows, unsigned int cols)
        {
            grid_map_.resize(rows,cols);
            grid_map_.setConstant(std::numeric_limits<DataType>::infinity());
        }
		/**
		    \brief Creates a grid map with size equal to the bounding box
		 */
        void makeGridMapSize()
        {
            grid_map_.resize(grid_size_.x(),grid_size_.y());
            grid_map_.setConstant(std::numeric_limits<DataType>::infinity());
        }

        /**
            \brief Access the list of voxels intersected by the ray
        */
        std::vector< Point3D<int> > getVoxels(const unsigned int idx) const 
        { 
            return voxel_coordinates_[idx]; 
        }
        /**
            \brief Access the list of voxels intersected by the ray
        */
        std::vector< Point3D<int> > getVoxels() const 
        { 
            unsigned int idx = 0;
            return voxel_coordinates_[idx]; 
        }

        /**
            \brief Get the grid size
        */
        Point3D<DataType> gridSize();

        /**
            \brief Get the cell size
        */
        inline Point3D<DataType> cellSize() const 
        { 
            return voxel_size_; 
        }

        /**
            \brief Filter voxels using functor
        */
        template<class Compare>
        void filterVoxels(Compare comp);

    protected:
        /**
            \brief Checks whether there is any intersection between ray and the defined grid-boundaries
            The ray will be expressed in parametric form, as o + td where o is the ray origin,
            d is the ray direction and t is a parameter.
            \param[in] r is the ray to test
            \param[in] t0 is the lower limit bound of t for the intersection test (example, t0=0)
            \param[in] t1 is the upper limit bound of t for the intersection test (example, t1=1)
            \param[out] tmin_val is the minimum limit of the intersection interval along the ray
            \param[out] tmax_val is the maximum limit of the intersection interval along the ray
            \retval true if there is intersection, false otherwise
        */
        bool intersectTest(const Ray<DataType>& r, float t0, float t1, float* tmin_val, float* tmax_val);
        /**
            \brief Overloaded version without intersection interval
        */
        bool intersectTest(const Ray<DataType>& r, float* tmin_val, float* tmax_val);

    private:
        /**
            \brief Fix the limits of the bounding box of the grid
        */
        void fixBoundingBoxLimits();

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:        
        Point3D<int> grid_size_;          //! number of cells in each axis
        Point3D<DataType> grid_bounds_[2];     //! start and end points of the grid
        Point3D<DataType> voxel_size_;         //! the dimensions of one voxel in each axis
        std::vector< Ray<DataType> > rays_;
        std::vector< std::vector< Point3D<int> > > voxel_coordinates_; //! the coordinates of the voxels intersected by the ray
        GridMap2D grid_map_;
};

//As this class is based on templates, HAVE to include the implementation here.
//(Awful)
#include "ray_tracing.tcc"

} //namespace

#endif //RAY_TRACING_H
