#pragma once
#ifndef __UTILS_H__
#define __UTILS_H__

namespace distmesh {
    namespace utils {
        // select array elements based on mask，选择函数，用来筛选
        template <
            class type
        >
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> selectMaskedArrayElements(
            Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
            Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, 1> const> const mask) {
            Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(mask.count(), array.cols());

            int resultCount = 0;
            for (int row = 0; row < array.rows(); ++row) {
                if (mask(row)) {
                    result.row(resultCount) = array.row(row);
                    resultCount++;
                }
            }

            return result;
        }

        // select array elements based on indices
        template <
            class type
        >
        Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> selectIndexedArrayElements(
            Eigen::Ref<Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> const> const array,
            Eigen::Ref<Eigen::ArrayXi const> const indices) {
            Eigen::Array<type, Eigen::Dynamic, Eigen::Dynamic> result(indices.rows(), array.cols());

            for (int row = 0; row < indices.rows(); ++row) {
                result.row(row) = array.row(indices(row));
            }

            return result;
        }

        // calculate factorial recursively,阶乘
        inline unsigned factorial(unsigned const n) {
            if (n <= 1) {
                return 1;
            }
            else {
                return n * factorial(n - 1);
            }
        }

        // easy creation of n-dimensional bounding box
        Eigen::ArrayXXd boundingBox(unsigned const dimensions);

        // create initial points distribution
        Eigen::ArrayXXd createInitialPoints(Functional const& Fd,
            double const h0, Functional const& Fh,
            Eigen::Ref<Eigen::ArrayXXd const> const bbox,
            Eigen::Ref<Eigen::ArrayXXd const> const pfix);

        // create array with all unique combinations n over k
        Eigen::ArrayXXi nOverK(unsigned const n, unsigned const k);

        // get a unique list of all edges in given triangulation
        Eigen::ArrayXXi findUniqueEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation);

        // get indices of bars in triangulation
        Eigen::ArrayXXi getTriangulationEdgeIndices(Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
            Eigen::Ref<Eigen::ArrayXXi const> const edges);

        // determine boundary edges of given triangulation
        Eigen::ArrayXi boundEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
            Eigen::Ref<Eigen::ArrayXXi const> const edges = Eigen::ArrayXXi(),
            Eigen::Ref<Eigen::ArrayXXi const> const edgeIndices = Eigen::ArrayXXi());

        // fix orientation of edges located at the boundary
        Eigen::ArrayXXi fixBoundaryEdgeOrientation(Eigen::Ref<Eigen::ArrayXXd const> const nodes,
            Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
            Eigen::Ref<Eigen::ArrayXXi const> const edges,
            Eigen::Ref<Eigen::ArrayXXi const> const edgeIndices);

        // project points outside of domain back to boundary
        void projectPointsToBoundary(Functional const& Fd,
            double const h0, Eigen::Ref<Eigen::ArrayXXd> points);

        // check whether points lies inside or outside of polygon
        Eigen::ArrayXd pointsInsidePoly(
            Eigen::Ref<Eigen::ArrayXXd const> const points,
            Eigen::Ref<Eigen::ArrayXXd const> const polygon);
    }
}

#endif