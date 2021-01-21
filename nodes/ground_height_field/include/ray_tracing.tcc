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

#include <omp.h>

//#define DEBUG

static const Point3D<int> pZero(0,0,0);
static const Point3D<int> pOnes(1,1,1);
static const Point2D<int> pZero2D(0,0);
static const Point2D<int> pOnes2D(1,1);

template<typename DataType>
RayTracing<DataType>::RayTracing():
    grid_size_(Point3D<int>()),
    voxel_size_(Point3D<DataType>())
{
    grid_bounds_[0].clear();
    grid_bounds_[1].clear();

//    voxel_coordinates_.clear();
}

template<typename DataType>
RayTracing<DataType>::RayTracing(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds, DataType cell_size_x, DataType cell_size_y, DataType cell_size_z):
    voxel_size_(Point3D<DataType>(cell_size_x, cell_size_y, cell_size_z))
{
    grid_bounds_[0] = grid_min_bounds;
    grid_bounds_[1] = grid_max_bounds;

//    voxel_coordinates_.clear();
    fixBoundingBoxLimits();
}

template<typename DataType>
RayTracing<DataType>::RayTracing(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds, Point3D<DataType> cell_size):
    voxel_size_(cell_size)
{
    setGridBounds(grid_min_bounds, grid_max_bounds);
}

template<typename DataType>
RayTracing<DataType>::~RayTracing()
{

}

template<typename DataType>
void RayTracing<DataType>::setGridBounds(Point3D<DataType> grid_min_bounds, Point3D<DataType> grid_max_bounds)
{
    grid_bounds_[0] = grid_min_bounds;
    grid_bounds_[1] = grid_max_bounds;

//    voxel_coordinates_.clear();
    fixBoundingBoxLimits();
}

template<typename DataType>
void RayTracing<DataType>::setCellSize(Point3D<DataType> cell_size)
{
    voxel_size_ = cell_size;
    fixBoundingBoxLimits();
}

template<typename DataType>
void RayTracing<DataType>::fixBoundingBoxLimits()
{
    //make sure limits are correctly ordered
    if (grid_bounds_[0].x() > grid_bounds_[1].x()) {
        DataType val = grid_bounds_[0].x();
        grid_bounds_[0].setX(grid_bounds_[1].x());
        grid_bounds_[1].setX(val);
    }
    if (grid_bounds_[0].y() > grid_bounds_[1].y()) {
        DataType val = grid_bounds_[0].y();
        grid_bounds_[0].setY(grid_bounds_[1].y());
        grid_bounds_[1].setY(val);
    }
    if (grid_bounds_[0].z() > grid_bounds_[1].z()) {
        DataType val = grid_bounds_[0].z();
        grid_bounds_[0].setZ(grid_bounds_[1].z());
        grid_bounds_[1].setZ(val);
    }

    grid_size_ = Point3D<int>((grid_bounds_[1] - grid_bounds_[0]) / voxel_size_);
}

template<typename DataType>
void RayTracing<DataType>::addRay(const Ray<DataType>& ray)
{
    rays_.push_back(ray);
//    voxel_coordinates_.push_back(std::vector< Point3D<int> >());
}

template<typename DataType>
Ray<DataType>& RayTracing<DataType>::getRay(const unsigned int idx)
{
    if (idx < rays_.size()) {
        return rays_.at(idx);
    } else {
        return rays_.back();
    }
}

template<typename DataType>
void RayTracing<DataType>::computeRayIntersections(const unsigned int idx)
{
    /*
    bool intersection = false;
    float tMin = 0, tMax = 0;
    Point<DataType> cell_start;
    Point<DataType> cell_end;
    Point<int> startIndex, endIndex;
    Point<int> step(1,1,1);
    Point<DataType> tDelta;
    Point<DataType> tMaxVal;
    Point<DataType> tDeltaZero;
    Point<DataType> tMaxValZero;
    Point<int> dirSign;
    const Point<int> pzero;
    Point<DataType> diff;
    Ray<DataType> ray;

    //Check if ray index is valid
    if (idx >= voxel_coordinates_.size()) {
        return;
    }

    //gets this ray
    ray = getRay(idx);

    //initializes the voxel list
    voxel_coordinates_[idx].clear();

    //check if the grid and/or the ray are valid
    if ((grid_size_ == Point<int>()) || ray.length() == 0) {
        return;
    }

#if defined(DEBUG)
    std::cout << "grid_size " << grid_size_ << std::endl;
    std::cout << "grid_bounds[0] " << grid_bounds_[0] << std::endl;
    std::cout << "grid_bounds[1] " << grid_bounds_[1] << std::endl;
    std::cout << "ray origin " << ray.origin() << std::endl;
    std::cout << "ray direction " << ray.direction() << std::endl;
#endif

    //check if the ray intersects any part of the grid
    intersection = intersectTest(ray, &tMin, &tMax);

    if (!intersection) {
        //no, nothing to return
        return;
    }

    //gets the cell start and end coordinates
    cell_start = ray.origin() + ray.direction()*tMin;
    cell_end = ray.origin() + ray.direction()*tMax;

#if defined(DEBUG)
    std::cout << "tMin " << tMin << std::endl;
    std::cout << "tMax " << tMax << std::endl;
    std::cout << "cell_start " << cell_start << std::endl;
    std::cout << "cell_end " << cell_end << std::endl;
    std::cout << "voxelSize" << voxel_size_ << std::endl;
#endif

    // Determine initial voxel coordinates and line traversal directions
    //NOTE: round could be better than ceil
    //consider using a rounding function F which minimizes
    // (cell_start.*voxel_size + grid_bounds_[0]) - start
    // and 
    // (cell_end.*voxel_size + grid_bounds_[0]) - end
    diff = (cell_start - grid_bounds_[0]);
    diff /= voxel_size_;
    startIndex = Point<int>::max(pzero,
                                 Point<int>::ceil(diff));
    diff = (cell_end - grid_bounds_[0]);
    diff /= voxel_size_;
    endIndex = Point<int>::max(pzero,
                            Point<int>::ceil(diff));

    dirSign = Point<int>(((ray.direction().x() > 0) ? 1 : ((ray.direction().x() < 0) ? -1 : 0)),
                         ((ray.direction().y() > 0) ? 1 : ((ray.direction().y() < 0) ? -1 : 0)),
                         ((ray.direction().z() > 0) ? 1 : ((ray.direction().z() < 0) ? -1 : 0))
                        );

    step = dirSign;
    tDeltaZero = Point<DataType>(tMax,tMax,tMax);
    tMaxValZero = Point<DataType>(tMax,tMax,tMax);
    Point<DataType> taux(tMax,tMax,tMax);
    tDelta = voxel_size_ / (Point<DataType>(dirSign) % ray.direction());
    //if direction > 0, use startIndex, 
    //otherwise (< 0), use startIndex - 1
    //to avoid branching we use the expression
    //startIndex + (dirSign - 1)/2 in integer arithmetic
    Point<int> auxP( startIndex + ((dirSign - Point<int>(1,1,1)) / 2) );
    tMaxVal = Point<DataType>(tMin, tMin, tMin) + 
                    (( grid_bounds_[0] + (( Point<DataType>(auxP) % voxel_size_) - cell_start) ) / ray.direction());

    if (!dirSign.x()) {
        tDelta.setX(tDeltaZero.x());
        tMaxVal.setX(tMaxValZero.x());
    }
    if (!dirSign.y()) {
        tDelta.setY(tDeltaZero.y());
        tMaxVal.setY(tMaxValZero.y());
    }
    if (!dirSign.z()) {
        tDelta.setZ(tDeltaZero.z());
        tMaxVal.setZ(tMaxValZero.z());
    }

#if defined(DEBUG)
    std::cout << "startIndex " << startIndex << std::endl;
    std::cout << "endIndex " << endIndex << std::endl;
    std::cout << "dirSign " << dirSign << std::endl;
    std::cout << "tDelta " << tDelta << std::endl;
    std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif
    //add first voxel grid coordinates
    voxel_coordinates_[idx].push_back(startIndex);

    tDelta.fix();
    tMaxVal.fix();
#if defined(DEBUG)
#if (0)
    std::cout << "tDelta fixed " << tDelta << std::endl;
    std::cout << "tMaxVal fixed " << tMaxVal << std::endl;
#endif
#endif

    // Step iteratively through the grid
    while (startIndex != endIndex) {
        Point<int> indxIncr;
        Point<DataType> tMaxIncr;
        if (tMaxVal.x() < tMaxVal.y()) {
            if (tMaxVal.x() < tMaxVal.z()) {
#if defined(DEBUG)
                std::cout << "(x<y && x<z)" << std::endl;
#endif
                indxIncr.setX(step.x());
                tMaxIncr.setX(tDelta.x());
            } else {
#if defined(DEBUG)
                std::cout << "(x<y && z<x)" << std::endl;
#endif
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        } else {
#if defined(DEBUG)            
            std::cout << "(y<x && y<z)" << std::endl;
#endif            
            if (tMaxVal.y() < tMaxVal.z()) {
                indxIncr.setY(step.y());
                tMaxIncr.setY(tDelta.y());
            } else {
#if defined(DEBUG)
                std::cout << "(y<x && z<y)" << std::endl;
#endif                
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        }
        startIndex += indxIncr;
        tMaxVal += tMaxIncr;
#if defined(DEBUG)
        std::cout << "indxIncr " << indxIncr << std::endl;
        std::cout << "tMaxIncr " << tMaxIncr << std::endl;
        std::cout << "startIndex " << startIndex << std::endl;
        std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif

        //add the voxel grid coordinates
        voxel_coordinates_[idx].push_back(startIndex);
    }
    */
}

template<typename DataType>
void RayTracing<DataType>::computeRayIntersectionsZFilter(const unsigned int idx)
{
/*
    bool intersection = false;
    float tMin = 0, tMax = 0;
    Point3D<DataType> cell_start;
    Point3D<DataType> cell_end;
    Point3D<int> startIndex, endIndex;
    Point3D<int> step(1,1,1);
    Point3D<DataType> tDelta;
    Point3D<DataType> tMaxVal;
    Point3D<DataType> tDeltaZero;
    Point3D<DataType> tMaxValZero;
    Point3D<int> dirSign;
    Point3D<DataType> diff;
    Point3D<int> pCeil;
    Ray<DataType> ray;
    Point3D<int> indxIncr;
    Point3D<DataType> tMaxIncr;

    //Check if ray index is valid
    if (idx >= voxel_coordinates_.size()) {
        return;
    }

    //gets this ray
    ray = getRay(idx);

    //initializes the voxel list
    voxel_coordinates_[idx].clear();
    
    //check if the grid and/or the ray are valid
    if ((grid_size_ == pZero) || ray.length() == 0) {
        return;
    }

    //check if the ray intersects any part of the grid
    intersection = intersectTest(ray, &tMin, &tMax);

#if defined(DEBUG)
    std::cout << "grid_size " << grid_size_ << std::endl;
    std::cout << "grid_bounds[0] " << grid_bounds_[0] << std::endl;
    std::cout << "grid_bounds[1] " << grid_bounds_[1] << std::endl;
    std::cout << "ray origin " << ray.origin() << std::endl;
    std::cout << "ray direction " << ray.direction() << std::endl;
    std::cout << "ray length " << ray.length() << std::endl;
    std::cout << "intersection " << intersection << std::endl;
#endif

    if (!intersection) {
        //no, nothing to return
        return;
    }

    //gets the cell start and end coordinates
    //cell_start = ray.origin() + ray.direction()*tMin;
    cell_start = ray.origin();
    cell_end = ray.origin() + ray.direction()*tMax;

#if defined(DEBUG)
    std::cout << "tMin " << tMin << std::endl;
    std::cout << "tMax " << tMax << std::endl;
    std::cout << "cell_start " << cell_start << std::endl;
    std::cout << "cell_end " << cell_end << std::endl;
    std::cout << "voxelSize" << voxel_size_ << std::endl;
#endif

	///**
	// * \note
    // * Determine initial voxel coordinates and line traversal directions
    // * NOTE: round could be better than ceil
    // * consider using a rounding function F which minimizes
    // * (cell_start.*voxel_size + grid_bounds_[0]) - start
    // * and
    // * (cell_end.*voxel_size + grid_bounds_[0]) - end
    // /
    diff = (cell_start - grid_bounds_[0]);
    diff /= voxel_size_;
    pCeil = diff.ceil();
    startIndex = Point3D<int>::max(pZero, pCeil);
    diff = (cell_end - grid_bounds_[0]);
    diff /= voxel_size_;
    pCeil = diff.ceil();
    endIndex = Point3D<int>::max(pZero, pCeil);

    dirSign = ray.direction().sign2();

    step = dirSign;
    tDeltaZero = Point3D<DataType>(tMax,tMax,tMax);
    tMaxValZero = Point3D<DataType>(tMax,tMax,tMax);
    Point3D<DataType> taux(tMax,tMax,tMax);
    tDelta = voxel_size_;
    tDelta /= (Point3D<DataType>(dirSign) % ray.direction());

    //if direction > 0, use startIndex, 
    //otherwise (< 0), use startIndex - 1
    //to avoid branching we use the expression
    //startIndex + (dirSign - 1)/2 in integer arithmetic
    Point3D<int> diffInt = (dirSign - pOnes);
    diffInt /= 2;
    diffInt += startIndex;
    Point3D<DataType> aux2 = diffInt.template cast<DataType>();
    aux2 %= voxel_size_;
    diff = aux2 - cell_start;
    aux2 = grid_bounds_[0] + diff;
    aux2 /= ray.direction();
    tMaxVal.init(tMin, tMin, tMin);
    tMaxVal += aux2;

    if (!dirSign.x()) {
        tDelta.setX(tDeltaZero.x());
        tMaxVal.setX(tMaxValZero.x());
    }
    if (!dirSign.y()) {
        tDelta.setY(tDeltaZero.y());
        tMaxVal.setY(tMaxValZero.y());
    }
    if (!dirSign.z()) {
        tDelta.setZ(tDeltaZero.z());
        tMaxVal.setZ(tMaxValZero.z());
    }

#if defined(DEBUG)
    std::cout << "startIndex " << startIndex << std::endl;
    std::cout << "endIndex " << endIndex << std::endl;
    std::cout << "dirSign " << dirSign << std::endl;
    std::cout << "tDelta " << tDelta << std::endl;
    std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif

    tDelta.fix();
    tMaxVal.fix();
#if defined(DEBUG)
#if (0)
    std::cout << "tDelta fixed " << tDelta << std::endl;
    std::cout << "tMaxVal fixed " << tMaxVal << std::endl;
#endif
#endif

    // //preallocate the vector of coordinates with a rough size estimation
    // Point<DataType> len = ray.direction() / tDelta;
    // len = len.abs();
    // len.sort();
    // if (!len.x()) len.x() = 1;
    // if (!len.y()) len.y() = 1;
    // if (!len.z()) len.z() = 1;
    // voxel_coordinates_[idx].resize((len.x()+len.y())*len.z(), pZero);
    // unsigned int count = 0;

    //add first voxel grid coordinates
    voxel_coordinates_[idx].push_back(startIndex);
    // voxel_coordinates_[idx].at(count++) = startIndex;

    // Step iteratively through the grid
    while (startIndex != endIndex) {
        indxIncr.clear();
        tMaxIncr.clear();
        if (tMaxVal.x() < tMaxVal.y()) {
            if (tMaxVal.x() < tMaxVal.z()) {
#if defined(DEBUG)
                std::cout << "(x<y && x<z)" << std::endl;
#endif
                indxIncr.setX(step.x());
                tMaxIncr.setX(tDelta.x());
            } else {
#if defined(DEBUG)
                std::cout << "(x<y && z<x)" << std::endl;
#endif
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        } else {
#if defined(DEBUG)            
            std::cout << "(y<x && y<z)" << std::endl;
#endif            
            if (tMaxVal.y() < tMaxVal.z()) {
                indxIncr.setY(step.y());
                tMaxIncr.setY(tDelta.y());
            } else {
#if defined(DEBUG)
                std::cout << "(y<x && z<y)" << std::endl;
#endif                
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        }
        startIndex += indxIncr;
        tMaxVal += tMaxIncr;
#if defined(DEBUG)
        std::cout << "indxIncr " << indxIncr << std::endl;
        std::cout << "tMaxIncr " << tMaxIncr << std::endl;
        std::cout << "startIndex " << startIndex << std::endl;
        std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif
        auto p = voxel_coordinates_[idx].back();
        //auto p = voxel_coordinates_[idx].at(count-1);
        if (p.x() == startIndex.x() &&
            p.y() == startIndex.y() &&
            p.z() > startIndex.z()) {
            //update the voxel coordinates
            voxel_coordinates_[idx].back().setZ(startIndex.z());
        } else if (p.x() == startIndex.x() &&
                   p.y() == startIndex.y() &&
                   p.z() <= startIndex.z()) {
            //do nothing
        } else {
            //add the new voxel grid coordinates
            voxel_coordinates_[idx].push_back(startIndex);
            //voxel_coordinates_[idx].at(count++) = startIndex;
        }
    }

    //voxel_coordinates_[idx].resize(count);
*/
}

template<typename DataType>
void RayTracing<DataType>::computeRayIntersectionsZFilter2DGrid(const Ray<DataType>& ray, DataType base_height)
{
    bool intersection = false;
    float tMin = 0, tMax = 0;
    Point2D<DataType> cell_start;
    Point2D<DataType> cell_end;
    Point2D<int> startIndex, endIndex;
    Point2D<int> step(1,1);
    Point2D<DataType> tDelta;
    Point2D<DataType> tMaxVal;
    Point2D<DataType> tDeltaZero;
    Point2D<DataType> tMaxValZero;
    Point2D<int> dirSign;
    Point2D<DataType> diff;
    Point2D<int> pCeil;
    Point2D<int> indxIncr;
    Point2D<DataType> tMaxIncr;
    Point2D<DataType> taux;
    Point2D<DataType> vox_size2D(voxel_size_.x(), voxel_size_.y());
    Point2D<DataType> gbounds2D(grid_bounds_[0].x(), grid_bounds_[0].y());
    Point2D<DataType> raydir2D;
    double tan_alpha = 0;

    //check if the grid and/or the ray are valid
    if ((grid_size_ == pZero) || ray.length() == 0) {
        return;
    }

    raydir2D = Point2D<DataType>(ray.direction().x(), ray.direction().y());
    //computes the tangent of the angle formed by this ray
    tan_alpha = tan(ray.direction().z() / raydir2D.length());

    //check if the ray intersects any part of the grid
    intersection = intersectTest(ray, &tMin, &tMax);

    if (!intersection) {
        //no, nothing to return
        return;
    }
    

    //gets the cell start and end coordinates
    //cell_start = ray.origin() + ray.direction()*tMin;
    cell_start = ray.origin();
    cell_end = ray.origin() + ray.direction()*tMax;

	/**
	 * \note
     * Determine initial voxel coordinates and line traversal directions
     * NOTE: round could be better than ceil
     * consider using a rounding function F which minimizes
     * (cell_start.*voxel_size + grid_bounds_[0]) - start
     * and
     * (cell_end.*voxel_size + grid_bounds_[0]) - end
     */
    diff = cell_start - gbounds2D;
    diff /= vox_size2D;
    pCeil = diff.ceil();
    startIndex = Point2D<int>::max(pZero2D, pCeil);
    diff = cell_end - gbounds2D;
    diff /= vox_size2D;
    pCeil = diff.ceil();
    endIndex = Point2D<int>::max(pZero2D, pCeil);

    dirSign = ray.direction().sign2();

    step = dirSign;
    tDeltaZero = Point2D<DataType>(tMax,tMax);
    tMaxValZero = tDeltaZero;
    taux = tDeltaZero;
    tDelta = vox_size2D;
    tDelta /= (Point2D<DataType>(dirSign) % raydir2D);

    /**
     * \note:
     * if direction > 0, use startIndex,
     * otherwise (< 0), use startIndex - 1
     * to avoid branching we use the expression
     * startIndex + (dirSign - 1)/2 in integer arithmetic
     */
    Point2D<int> diffInt = (dirSign - pOnes2D);
    diffInt /= 2;
    diffInt += startIndex;
    Point2D<DataType> aux2 = diffInt.template cast<DataType>();
    aux2 %= vox_size2D;
    diff = aux2 - cell_start;
    aux2 = gbounds2D + diff;
    aux2 /= raydir2D;
    tMaxVal.init(tMin, tMin);
    tMaxVal += aux2;

    if (!dirSign.x()) {
        tDelta.setX(tDeltaZero.x());
        tMaxVal.setX(tMaxValZero.x());
    }
    if (!dirSign.y()) {
        tDelta.setY(tDeltaZero.y());
        tMaxVal.setY(tMaxValZero.y());
    }
    //if (!dirSign.z()) {
    //    tDelta.setZ(tDeltaZero.z());
    //    tMaxVal.setZ(tMaxValZero.z());
    //}


    tDelta.fix();
    tMaxVal.fix();

    //if grid_map is empty, resize it to match the bounding box size
    if (grid_map_.size() == 0) {
        makeGridMapSize();
    }

    //add first voxel grid coordinates
    DataType height = startIndex.length()*tan_alpha + base_height;
    if (height < grid_map_(startIndex.x(), startIndex.y())) {
        grid_map_(startIndex.x(), startIndex.y()) = height;
    }

    // Step iteratively through the grid
    while (startIndex != endIndex) {
        indxIncr.clear();
        tMaxIncr.clear();
        if (tMaxVal.x() < tMaxVal.y()) {
            //if (tMaxVal.x() < tMaxVal.z()) {
                indxIncr.setX(step.x());
                tMaxIncr.setX(tDelta.x());
            //} else {
            //    indxIncr.setZ(step.z());
            //    tMaxIncr.setZ(tDelta.z());
            //}
        } else {
            //if (tMaxVal.y() < tMaxVal.z()) {
                indxIncr.setY(step.y());
                tMaxIncr.setY(tDelta.y());
            //} else {
            //    indxIncr.setZ(step.z());
            //    tMaxIncr.setZ(tDelta.z());
            //}
        }
        startIndex += indxIncr;
        tMaxVal += tMaxIncr;

        height = startIndex.length()*tan_alpha + base_height;
        if (height < grid_map_(startIndex.x(), startIndex.y())) {
            grid_map_(startIndex.x(), startIndex.y()) = height;
        }
        //std::cout << startIndex << std::endl;
    }
}

//template<typename DataType>
//void RayTracing<DataType>::computeRayIntersectionsZFilter2DGridEFLA(const unsigned int idx, DataType base_height, DataType min_height, DataType max_height)
//{
//    float tMin = 0, tMax = 0;
//    Point2D<DataType> cell_start;
//    Point2D<DataType> cell_end;
//    Point2D<int> startIndex, endIndex;
//    Point2D<int> step(1,1);
//    Point2D<DataType> tDelta;
//    Point2D<DataType> tMaxVal;
//    Point2D<DataType> tDeltaZero;
//    Point2D<DataType> tMaxValZero;
//    Point2D<int> dirSign;
//    Point2D<DataType> diff;
//    Point2D<int> pCeil;
//    Ray<DataType> ray;
//    Point2D<int> indxIncr;
//    Point2D<DataType> tMaxIncr;
//    Point2D<DataType> taux;
//    Point2D<DataType> vox_size2D(voxel_size_.x(), voxel_size_.y());
//    Point2D<DataType> gbounds2D(grid_bounds_[0].x(), grid_bounds_[0].y());
//    Point2D<DataType> raydir2D;
//    double tan_alpha = 0;
//	bool yLonger=false;
//	double incrementVal, endVal;
//	DataType height;
//
//    //gets this ray
//    ray = getRay(idx);
//    //check if the grid and/or the ray are valid
//    if ((grid_size_ == pZero) || ray.length() == 0) {
//        return;
//    }
//
//    raydir2D = Point2D<DataType>(ray.direction().x(), ray.direction().y());
//    //computes the tangent of the angle formed by this ray
//    tan_alpha = tan(ray.direction().z() / raydir2D.length());
//
//    //check if the ray intersects any part of the grid
//    if (!(intersectTest(ray, &tMin, &tMax))) {
//        //no, nothing to return
//        return;
//    }
//
//
//    //gets the cell start and end coordinates
//    //cell_start = ray.origin() + ray.direction()*tMin;
//    cell_start = ray.origin();
//    cell_end = ray.origin() + ray.direction()*tMax;
//
//	/**
//	 * \note
//     * Determine initial voxel coordinates and line traversal directions
//     * NOTE: round could be better than ceil
//     * consider using a rounding function F which minimizes
//     * (cell_start.*voxel_size + grid_bounds_[0]) - start
//     * and
//     * (cell_end.*voxel_size + grid_bounds_[0]) - end
//     */
//    diff = cell_start - gbounds2D;
//    diff /= vox_size2D;
//    pCeil = diff.ceil();
//    startIndex = Point2D<int>::max(pZero2D, pCeil);
//    diff = cell_end - gbounds2D;
//    diff /= vox_size2D;
//    pCeil = diff.ceil();
//    endIndex = Point2D<int>::max(pZero2D, pCeil);
//
//    dirSign = ray.direction().sign2();
//
//    step = dirSign;
//    tDeltaZero = Point2D<DataType>(tMax,tMax);
//    tMaxValZero = tDeltaZero;
//    taux = tDeltaZero;
//    tDelta = vox_size2D;
//    tDelta /= (Point2D<DataType>(dirSign) % raydir2D);
//    tDelta.fix();
//    //tDelta %= Point2D<DataType>(step.x(), step.y());
//
//    diff = cell_end - cell_start;
//	DataType shortLen=diff.y();
//	DataType longLen=diff.x();
//	if (abs(shortLen)>abs(longLen)) {
//		DataType swap=shortLen;
//		shortLen=longLen;
//		longLen=swap;
//		yLonger=true;
//	}
//
//	endVal=longLen;
//	if (longLen<0) {
//		incrementVal=-1;
//		longLen=-longLen;
//	} else incrementVal=1;
//
//	double decInc;
//	if (longLen==0) decInc=(double)shortLen;
//	else decInc=((double)shortLen/(double)longLen);
//	double j=0.0;
//	decInc *= vox_size2D.x();
//
//	Point2D<DataType> index = cell_start;
//	Point2D<int> indexInt = startIndex;
//	Point2D<DataType> p;
//
//	if (yLonger) {
//		incrementVal *= vox_size2D.y();
//		if (dirSign.y() > 0) {
//			for (double i=0; /*indexInt < endIndex*/i<endVal; i+=incrementVal) {
//				indexInt = startIndex + Point2D<int>(j/vox_size2D.x(), i/vox_size2D.y());
//				height = index.length()*tan_alpha + base_height;
//				if (height < min_height) {
//					height = min_height;
//				} else if (height < max_height) {
//					height = max_height;
//				}
//				if (height < grid_map_(indexInt.x(), indexInt.y())) {
//				    grid_map_(indexInt.x(), indexInt.y()) = height;
//				}
//				j+=decInc;
//				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
//			}
//		} else {
//			for (double i=0; /*indexInt > endIndex*/i>endVal; i+=incrementVal) {
//				indexInt = startIndex + Point2D<int>(j/vox_size2D.x(), i/vox_size2D.y());
//				height = index.length()*tan_alpha + base_height;
//				if (height < min_height) {
//					height = min_height;
//				} else if (height < max_height) {
//					height = max_height;
//				}
//				if (height < grid_map_(indexInt.x(), indexInt.y())) {
//				    grid_map_(indexInt.x(), indexInt.y()) = height;
//				}
//				j+=decInc;
//				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
//			}
//		}
//	} else {
//		incrementVal *= vox_size2D.x();
//		if (dirSign.x() > 0) {
//			for (double i=0; /*indexInt < endIndex*/i<endVal; i+=incrementVal) {
//				indexInt = startIndex + Point2D<int>(i/vox_size2D.x(), j/vox_size2D.y());
//				height = index.length()*tan_alpha + base_height;
//				if (height < min_height) {
//					height = min_height;
//				} else if (height < max_height) {
//					height = max_height;
//				}
//				if (height < grid_map_(indexInt.x(), indexInt.y())) {
//				    grid_map_(indexInt.x(), indexInt.y()) = height;
//				}
//				j+=decInc;
//				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
//			}
//		} else {
//			for (double i=0; /*indexInt > endIndex*/i>endVal; i+=incrementVal) {
//				indexInt = startIndex + Point2D<int>(i/vox_size2D.x(), j/vox_size2D.y());
//				height = index.length()*tan_alpha + base_height;
//				if (height < min_height) {
//					height = min_height;
//				} else if (height < max_height) {
//					height = max_height;
//				}
//				if (height < grid_map_(indexInt.x(), indexInt.y())) {
//				    grid_map_(indexInt.x(), indexInt.y()) = height;
//				}
//				j+=decInc;
//				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
//			}
//		}
//	}
//}

template<typename DataType>
void RayTracing<DataType>::computeRayIntersectionsZFilter2DGridEFLA(const Ray<DataType>& ray, DataType base_height, DataType min_height, DataType max_height)
{
    float tMin = 0, tMax = 0;
    Point2D<DataType> cell_start;
    Point2D<DataType> cell_end;
    Point2D<int> startIndex, endIndex;
    Point2D<int> step(1,1);
    Point2D<DataType> tDelta;
    Point2D<DataType> tMaxVal;
    Point2D<DataType> tDeltaZero;
    Point2D<DataType> tMaxValZero;
    Point2D<int> dirSign;
    Point2D<DataType> diff;
    Point2D<int> pCeil;
    Point2D<int> indxIncr;
    Point2D<DataType> tMaxIncr;
    Point2D<DataType> taux;
    Point2D<DataType> vox_size2D(voxel_size_.x(), voxel_size_.y());
    Point2D<DataType> gbounds2D(grid_bounds_[0].x(), grid_bounds_[0].y());
    Point2D<DataType> raydir2D;
    double tan_alpha = 0;
	bool yLonger=false;
	double incrementVal, endVal;
	DataType height;

    //check if the grid and/or the ray are valid
    if ((grid_size_ == pZero) || ray.length() == 0) {
        return;
    }

    raydir2D = Point2D<DataType>(ray.direction().x(), ray.direction().y());
    //computes the tangent of the angle formed by this ray
    tan_alpha = tan(ray.direction().z() / raydir2D.length());

    //check if the ray intersects any part of the grid
    if (!(intersectTest(ray, &tMin, &tMax))) {
        //no, nothing to return
        return;
    }


    //gets the cell start and end coordinates
    //cell_start = ray.origin() + ray.direction()*tMin;
    cell_start = ray.origin();
    cell_end = ray.origin() + ray.direction()*tMax;

	/**
	 * \note
     * Determine initial voxel coordinates and line traversal directions
     * NOTE: round could be better than ceil
     * consider using a rounding function F which minimizes
     * (cell_start.*voxel_size + grid_bounds_[0]) - start
     * and
     * (cell_end.*voxel_size + grid_bounds_[0]) - end
     */
    diff = cell_start - gbounds2D;
    diff /= vox_size2D;
    pCeil = diff.ceil();
    startIndex = Point2D<int>::max(pZero2D, pCeil);
    diff = cell_end - gbounds2D;
    diff /= vox_size2D;
    pCeil = diff.ceil();
    endIndex = Point2D<int>::max(pZero2D, pCeil);

    dirSign = ray.direction().sign2();

    step = dirSign;
    tDeltaZero = Point2D<DataType>(tMax,tMax);
    tMaxValZero = tDeltaZero;
    taux = tDeltaZero;
    tDelta = vox_size2D;
    tDelta /= (Point2D<DataType>(dirSign) % raydir2D);
    tDelta.fix();
    //tDelta %= Point2D<DataType>(step.x(), step.y());

    diff = cell_end - cell_start;
	DataType shortLen=diff.y();
	DataType longLen=diff.x();
	if (abs(shortLen)>abs(longLen)) {
		DataType swap=shortLen;
		shortLen=longLen;
		longLen=swap;
		yLonger=true;
	}

	endVal=longLen;
	if (longLen<0) {
		incrementVal=-1;
		longLen=-longLen;
	} else incrementVal=1;

	double decInc;
	if (longLen==0) decInc=(double)shortLen;
	else decInc=((double)shortLen/(double)longLen);
	double j=0.0;
	decInc *= vox_size2D.x();

	Point2D<DataType> index = cell_start;
	Point2D<int> indexInt = startIndex;
	Point2D<DataType> p;

	if (yLonger) {
		incrementVal *= vox_size2D.y();
		if (dirSign.y() > 0) {
			for (double i=0; /*indexInt < endIndex*/i<endVal; i+=incrementVal) {
				indexInt = startIndex + Point2D<int>(j/vox_size2D.x(), i/vox_size2D.y());
				height = index.length()*tan_alpha + base_height;
				if (height < min_height) {
					height = min_height;
				} else if (height < max_height) {
					height = max_height;
				}
				if (height < grid_map_(indexInt.x(), indexInt.y())) {
				    grid_map_(indexInt.x(), indexInt.y()) = height;
				}
				j+=decInc;
				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
			}
		} else {
			for (double i=0; /*indexInt > endIndex*/i>endVal; i+=incrementVal) {
				indexInt = startIndex + Point2D<int>(j/vox_size2D.x(), i/vox_size2D.y());
				height = index.length()*tan_alpha + base_height;
				if (height < min_height) {
					height = min_height;
				} else if (height < max_height) {
					height = max_height;
				}
				if (height < grid_map_(indexInt.x(), indexInt.y())) {
				    grid_map_(indexInt.x(), indexInt.y()) = height;
				}
				j+=decInc;
				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
			}
		}
	} else {
		incrementVal *= vox_size2D.x();
		if (dirSign.x() > 0) {
			for (double i=0; /*indexInt < endIndex*/i<endVal; i+=incrementVal) {
				indexInt = startIndex + Point2D<int>(i/vox_size2D.x(), j/vox_size2D.y());
				height = index.length()*tan_alpha + base_height;
				if (height < min_height) {
					height = min_height;
				} else if (height < max_height) {
					height = max_height;
				}
				if (height < grid_map_(indexInt.x(), indexInt.y())) {
				    grid_map_(indexInt.x(), indexInt.y()) = height;
				}
				j+=decInc;
				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
			}
		} else {
			for (double i=0; /*indexInt > endIndex*/i>endVal; i+=incrementVal) {
				indexInt = startIndex + Point2D<int>(i/vox_size2D.x(), j/vox_size2D.y());
				height = index.length()*tan_alpha + base_height;
				if (height < min_height) {
					height = min_height;
				} else if (height < max_height) {
					height = max_height;
				}
				if (height < grid_map_(indexInt.x(), indexInt.y())) {
				    grid_map_(indexInt.x(), indexInt.y()) = height;
				}
				j+=decInc;
				//std::cout << "(" << i << "," << j << ") -> " << index << ";" << indexInt << " => " << height << std::endl;
			}
		}
	}
}

template<typename DataType>
void RayTracing<DataType>::computeRayIntersectionsZFilter2DGridPoint3D(const Ray<DataType>& ray, DataType base_height)
{
    bool intersection = false;
    float tMin = 0, tMax = 0;
    Point3D<DataType> cell_start;
    Point3D<DataType> cell_end;
    Point3D<int> startIndex, endIndex;
    Point3D<int> step(1,1,1);
    Point3D<DataType> tDelta;
    Point3D<DataType> tMaxVal;
    Point3D<DataType> tDeltaZero;
    Point3D<DataType> tMaxValZero;
    Point3D<int> dirSign;
    Point3D<DataType> diff;
    Point3D<int> pCeil;
    Point3D<int> indxIncr;
    Point3D<DataType> tMaxIncr;


    //check if the grid and/or the ray are valid
    if ((grid_size_ == pZero) || ray.length() == 0) {
        return;
    }

    //check if the ray intersects any part of the grid
    intersection = intersectTest(ray, &tMin, &tMax);

#if defined(DEBUG)
    std::cout << "grid_size " << grid_size_ << std::endl;
    std::cout << "grid_bounds[0] " << grid_bounds_[0] << std::endl;
    std::cout << "grid_bounds[1] " << grid_bounds_[1] << std::endl;
    std::cout << "ray origin " << ray.origin() << std::endl;
    std::cout << "ray direction " << ray.direction() << std::endl;
    std::cout << "ray length " << ray.length() << std::endl;
    std::cout << "intersection " << intersection << std::endl;
#endif

    if (!intersection) {
        //no, nothing to return
        return;
    }


    //gets the cell start and end coordinates
    //cell_start = ray.origin() + ray.direction()*tMin;
    cell_start = ray.origin();
    cell_end = ray.origin() + ray.direction()*tMax;

#if defined(DEBUG)
    std::cout << "tMin " << tMin << std::endl;
    std::cout << "tMax " << tMax << std::endl;
    std::cout << "cell_start " << cell_start << std::endl;
    std::cout << "cell_end " << cell_end << std::endl;
    std::cout << "voxelSize" << voxel_size_ << std::endl;
#endif

	/**
	 * \note
     * Determine initial voxel coordinates and line traversal directions
     * NOTE: round could be better than ceil
     * consider using a rounding function F which minimizes
     * (cell_start.*voxel_size + grid_bounds_[0]) - start
     * and
     * (cell_end.*voxel_size + grid_bounds_[0]) - end
     */
    diff = (cell_start - grid_bounds_[0]);
    diff /= voxel_size_;
    pCeil = diff.ceil();
    startIndex = Point3D<int>::max(pZero, pCeil);
    diff = (cell_end - grid_bounds_[0]);
    diff /= voxel_size_;
    pCeil = diff.ceil();
    endIndex = Point3D<int>::max(pZero, pCeil);

    dirSign = ray.direction().sign2();

    step = dirSign;
    tDeltaZero = Point3D<DataType>(tMax,tMax,tMax);
    tMaxValZero = Point3D<DataType>(tMax,tMax,tMax);
    Point3D<DataType> taux(tMax,tMax,tMax);
    tDelta = voxel_size_;
    tDelta /= (Point3D<DataType>(dirSign) % ray.direction());

    /**
     * \note:
     * if direction > 0, use startIndex,
     * otherwise (< 0), use startIndex - 1
     * to avoid branching we use the expression
     * startIndex + (dirSign - 1)/2 in integer arithmetic
     */
    Point3D<int> diffInt = (dirSign - pOnes);
    diffInt /= 2;
    diffInt += startIndex;
    Point3D<DataType> aux2 = diffInt.template cast<DataType>();
    aux2 %= voxel_size_;
    diff = aux2 - cell_start;
    aux2 = grid_bounds_[0] + diff;
    aux2 /= ray.direction();
    tMaxVal.init(tMin, tMin, tMin);
    tMaxVal += aux2;

    if (!dirSign.x()) {
        tDelta.setX(tDeltaZero.x());
        tMaxVal.setX(tMaxValZero.x());
    }
    if (!dirSign.y()) {
        tDelta.setY(tDeltaZero.y());
        tMaxVal.setY(tMaxValZero.y());
    }
    if (!dirSign.z()) {
        tDelta.setZ(tDeltaZero.z());
        tMaxVal.setZ(tMaxValZero.z());
    }

#if defined(DEBUG)
    std::cout << "startIndex " << startIndex << std::endl;
    std::cout << "endIndex " << endIndex << std::endl;
    std::cout << "dirSign " << dirSign << std::endl;
    std::cout << "tDelta " << tDelta << std::endl;
    std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif

    tDelta.fix();
    tMaxVal.fix();
#if defined(DEBUG)
#if (0)
    std::cout << "tDelta fixed " << tDelta << std::endl;
    std::cout << "tMaxVal fixed " << tMaxVal << std::endl;
#endif
#endif

    //if grid_map is empty, resize it to match the bounding box size
    if (grid_map_.size() == 0) {
        makeGridMapSize();
    }

    //add first voxel grid coordinates
    if (startIndex.z() < grid_map_(startIndex.x(), startIndex.y())) {
        grid_map_(startIndex.x(), startIndex.y()) = startIndex.z() + base_height;
    }

    // Step iteratively through the grid
    while (startIndex != endIndex) {
        indxIncr.clear();
        tMaxIncr.clear();
        if (tMaxVal.x() < tMaxVal.y()) {
            if (tMaxVal.x() < tMaxVal.z()) {
#if defined(DEBUG)
                std::cout << "(x<y && x<z)" << std::endl;
#endif
                indxIncr.setX(step.x());
                tMaxIncr.setX(tDelta.x());
            } else {
#if defined(DEBUG)
                std::cout << "(x<y && z<x)" << std::endl;
#endif
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        } else {
#if defined(DEBUG)
            std::cout << "(y<x && y<z)" << std::endl;
#endif
            if (tMaxVal.y() < tMaxVal.z()) {
                indxIncr.setY(step.y());
                tMaxIncr.setY(tDelta.y());
            } else {
#if defined(DEBUG)
                std::cout << "(y<x && z<y)" << std::endl;
#endif
                indxIncr.setZ(step.z());
                tMaxIncr.setZ(tDelta.z());
            }
        }
        startIndex += indxIncr;
        tMaxVal += tMaxIncr;
#if defined(DEBUG)
        std::cout << "indxIncr " << indxIncr << std::endl;
        std::cout << "tMaxIncr " << tMaxIncr << std::endl;
        std::cout << "startIndex " << startIndex << std::endl;
        std::cout << "tMaxVal " << tMaxVal << std::endl;
#endif
        if (startIndex.z() < grid_map_(startIndex.x(), startIndex.y())) {
            grid_map_(startIndex.x(), startIndex.y()) = startIndex.z() + base_height;
        }
    }
}

template<typename DataType>
bool RayTracing<DataType>::intersectTest(const Ray<DataType>& r, float t0, float t1, float* tMin_val, float* tmax_val)
{
    Point3D<DataType> tmin;
    Point3D<DataType> tmax;

    tmin = Point3D<DataType>(grid_bounds_[r.direction_sign().x()].x(),
                           grid_bounds_[r.direction_sign().y()].y(),
                           grid_bounds_[r.direction_sign().z()].z());
    tmin -= r.origin();
    tmin %= r.inv_direction();
    tmax = Point3D<DataType>(grid_bounds_[1-r.direction_sign().x()].x(),
                           grid_bounds_[1-r.direction_sign().y()].y(),
                           grid_bounds_[1-r.direction_sign().z()].z());
    tmax -= r.origin();
    tmax %= r.inv_direction();

    if ((tmin.x() > tmax.y()) || (tmin.y() > tmax.x()) ||
        (tmin.x() > tmax.z()) || (tmin.z() > tmax.x())) {
        return false;
    }

    if ( tmin.y() > tmin.x() ) {
        tmin.x() = tmin.y();
    }
    if ( tmax.y() < tmax.x() ) {
        tmax.x() = tmax.y();
    }

    if ( tmin.z() > tmin.x() ) {
        tmin.x() = tmin.z();
    }
    if ( tmax.z() < tmax.x() ) {
        tmax.x() = tmax.z();
    }

    *tMin_val = tmin.x();
    *tmax_val = tmax.x();
    
    return ( (tmin.x() < t1) && (tmax.x() > t0) );
}

template<typename DataType>
bool RayTracing<DataType>::intersectTest(const Ray<DataType>& r, float* tMin_val, float* tmax_val)
{
    Point3D<DataType> tmin;
    Point3D<DataType> tmax;

    tmin = Point3D<DataType>(grid_bounds_[r.direction_sign().x()].x(),
                           grid_bounds_[r.direction_sign().y()].y(),
                           grid_bounds_[r.direction_sign().z()].z());
    tmin -= r.origin();
    tmin %= r.inv_direction();
    tmax = Point3D<DataType>(grid_bounds_[1-r.direction_sign().x()].x(),
                           grid_bounds_[1-r.direction_sign().y()].y(),
                           grid_bounds_[1-r.direction_sign().z()].z());
    tmax -= r.origin();
    tmax %= r.inv_direction();

    if ((tmin.x() > tmax.y()) || (tmin.y() > tmax.x()) ||
        (tmin.x() > tmax.z()) || (tmin.z() > tmax.x())) {
        return false;
    }

    if ( tmin.y() > tmin.x() ) {
        tmin.x() = tmin.y();
    }
    if ( tmax.y() < tmax.x() ) {
        tmax.x() = tmax.y();
    }

    if ( tmin.z() > tmin.x() ) {
        tmin.x() = tmin.z();
    }
    if ( tmax.z() < tmax.x() ) {
        tmax.x() = tmax.z();
    }

    *tMin_val = tmin.x();
    *tmax_val = tmax.x();
    
    return true;
}

template<typename DataType>
Point3D<DataType> RayTracing<DataType>::gridSize()
{
    return Point3D<DataType>::abs(grid_bounds_[1] - grid_bounds_[0]);
}

template<typename DataType>
template<class Compare>
void RayTracing<DataType>::filterVoxels(Compare comp)
{
    // std::vector< Point<int> > voxels_aux;

    // voxels_aux.push_back(voxel_coordinates_.front());
    // for (auto it=voxel_coordinates_.begin() + 1;
    //           it != voxel_coordinates_.end();
    //           it++) {
    //     //char method = ' ';
    //     auto p = voxels_aux.back();
    //     auto q = comp(p, *it);
    //     if (q == *it && it->z() < p.z()) {
    //         voxels_aux.back() = q;
    //         //method = '1';
    //     } else if (q != p) {
    //         voxels_aux.push_back(q);
    //         //method = '2';
    //     } else {
    //         //no need to update
    //         //method = '3';
    //     }
    //     ////uncomment these for debugging
    //     //std::cout << p.x() << "," << p.y() << "," << p.z() << " <-> ";
    //     //std::cout << it->x() << "," << it->y() << "," << it->z() << " result: ";
    //     //std::cout << q.x() << "," << q.y() << "," << q.z() << " method: " << method << std::endl;
    // }
    // voxel_coordinates_ = voxels_aux;
}
