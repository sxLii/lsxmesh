//核心算法
#include <vector>
#include <set>
#include <algorithm>

#include "distmesh.h"
#include "constants.h"
#include "triangulation.h"

// apply the distmesh algorithm
std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh::distmesh(
    Functional const& Fd, double const h0,
    Functional const& Fh, Eigen::Ref<Eigen::ArrayXXd const> const bbox,
    Eigen::Ref<Eigen::ArrayXXd const> const pfix) {
    // determine dimension of mesh
    unsigned const dimension = bbox.cols();

    // create initial distribution in bounding box
    Eigen::ArrayXXd points = utils::createInitialPoints(Fd,
        h0, Fh, bbox, pfix);

    // create initial triangulation
    Eigen::ArrayXXi triangulation = triangulation::delaunay(points);

    // create buffer to store old point locations to calculate
    // retriangulation and stop criterion
    Eigen::ArrayXXd retriangulationCriterionBuffer = Eigen::ArrayXXd::Constant(
        points.rows(), points.cols(), INFINITY);
    Eigen::ArrayXXd stopCriterionBuffer = Eigen::ArrayXXd::Zero(
        points.rows(), points.cols());

    // main distmesh loop
    Eigen::ArrayXXi edgeIndices;
    for (unsigned step = 0; step < constants::maxSteps; ++step) {
        // retriangulate if point movement is above threshold
        if ((points - retriangulationCriterionBuffer).square().rowwise().sum().sqrt().maxCoeff() >
            constants::ttol * h0) {
            // update triangulation
            triangulation = triangulation::delaunay(points);

            // reject triangles with circumcenter outside of the region
            Eigen::ArrayXXd circumcenter = Eigen::ArrayXXd::Zero(triangulation.rows(), dimension);
            for (int point = 0; point < triangulation.cols(); ++point) {
                circumcenter += utils::selectIndexedArrayElements<double>(
                    points, triangulation.col(point)) / triangulation.cols();
            }
            triangulation = utils::selectMaskedArrayElements<int>(triangulation,
                Fd(circumcenter) < -constants::geps * h0);

            // find unique edge indices，创建bar并使其unique化
            edgeIndices = utils::findUniqueEdges(triangulation);

            // store current points positions
            retriangulationCriterionBuffer = points;
        }

        // calculate edge vectors and their length
        auto const edgeVector = (utils::selectIndexedArrayElements<double>(points, edgeIndices.col(0)) -
            utils::selectIndexedArrayElements<double>(points, edgeIndices.col(1))).eval();
        auto const edgeLength = edgeVector.square().rowwise().sum().sqrt().eval();

        // evaluate elementSizeFunction at midpoints of edges
        auto const desiredElementSize = Fh(0.5 *
            (utils::selectIndexedArrayElements<double>(points, edgeIndices.col(0)) +
                utils::selectIndexedArrayElements<double>(points, edgeIndices.col(1)))).eval();

        // calculate desired edge length
        auto const desiredEdgeLength = (desiredElementSize * (1.0 + 0.4 / std::pow(2.0, dimension - 1)) *
            std::pow((edgeLength.pow(dimension).sum() / desiredElementSize.pow(dimension).sum()),
                1.0 / dimension)).eval();

        // calculate force vector for each edge
        auto const forceVector = (edgeVector.colwise() *
            ((desiredEdgeLength - edgeLength) / edgeLength).max(0.0)).eval();

        // store current points positions
        stopCriterionBuffer = points;

        // move points
        for (int edge = 0; edge < edgeIndices.rows(); ++edge) {
            if (edgeIndices(edge, 0) >= pfix.rows()) {
                points.row(edgeIndices(edge, 0)) += constants::deltaT * forceVector.row(edge);
            }
            if (edgeIndices(edge, 1) >= pfix.rows()) {
                points.row(edgeIndices(edge, 1)) -= constants::deltaT * forceVector.row(edge);
            }
        }

        // project points outside of domain to boundary
        utils::projectPointsToBoundary(Fd, h0, points);

        // stop, when maximum points movement is below threshold
        if ((points - stopCriterionBuffer).square().rowwise().sum().sqrt().maxCoeff() <
            constants::dptol * h0) {
            break;
        }
    }

    return std::make_tuple(points, triangulation);
}