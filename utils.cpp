#include <set>
#include <array>
#include <vector>

#include "distmesh.h"
#include "constants.h"

// easy creation of n-dimensional bounding box
Eigen::ArrayXXd distmesh::utils::boundingBox(unsigned const dimension) {
    Eigen::ArrayXXd box(2, dimension);
    box.row(0).fill(-1.0);
    box.row(1).fill(1.0);
    return box;
}

// create initial points distribution
Eigen::ArrayXXd distmesh::utils::createInitialPoints(
    Functional const& Fd, double const h0,
    Functional const& Fh, Eigen::Ref<Eigen::ArrayXXd const> const bbox,
    Eigen::Ref<Eigen::ArrayXXd const> const pfix) {
    // extract dimension of mesh
    unsigned const dimension = bbox.cols();

    // initially distribute points evenly in complete bounding box,最初在完整的边界框中均匀分布点
    Eigen::ArrayXi pointsPerDimension(dimension);
    for (int dim = 0; dim < dimension; ++dim) {
        pointsPerDimension(dim) = ceil((bbox(1, dim) - bbox(0, dim)) /
            (h0 * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0)));
    }//每一列维度，ceil返回最大整数，10.25->11

    Eigen::ArrayXXd points(pointsPerDimension.prod(), dimension);
    for (int point = 0; point < points.rows(); ++point)
        for (int dim = 0; dim < dimension; ++dim) {
            int const pointIndex = (point / std::max(pointsPerDimension.topRows(dim).prod(), 1)) %
                pointsPerDimension(dim);

            points(point, dim) = bbox(0, dim) + (double)pointIndex * h0 *
                (dim == 0 ? 1.0 : sqrt(3.0) / 2.0);

            if (dim > 0) {
                points(point, dim - 1) += pointIndex % 2 != 0 ? h0 / 2.0 : 0.0;
            }//x-> h0; y->sqrt3/2 h0; x/2->1/2h0
        }

    // reject points outside of region defined by distance function,删除大于Fd的点
    points = selectMaskedArrayElements<double>(points,
        Fd(points) < constants::geps * h0);

    // clear duplicate points，unique函数，删除重复点
    Eigen::Array<bool, Eigen::Dynamic, 1> isUniquePoint =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(points.rows(), true);
    for (int i = 0; i < pfix.rows(); ++i)
        for (int j = 0; j < points.rows(); ++j) {
            isUniquePoint(j) &= !(pfix.row(i) == points.row(j)).all();
        }
    points = selectMaskedArrayElements<double>(points, isUniquePoint);

    // calculate probability to keep points，r0=1/fh^2，按fh距离函数概率删点
    Eigen::ArrayXd probability = 1.0 / Fh(points).pow(dimension);
    probability /= probability.maxCoeff();

    // reject points with wrong probability
    points = selectMaskedArrayElements<double>(points,
        0.5 * (1.0 + Eigen::ArrayXd::Random(points.rows())) < probability);

    // combine fixed and variable points to one array
    if (pfix.rows()!=0) {
        Eigen::ArrayXXd finalPoints(points.rows() + pfix.rows(), dimension);
        finalPoints << pfix, points; //pfix不能为空
        return finalPoints;
    }
    
    return points;
}

// create array with all unique combinations n over k创建具有所有唯一组合 n 超过 k 的数组,即排列组合
Eigen::ArrayXXi distmesh::utils::nOverK(unsigned const n, unsigned const k) {
    // fill an array with all unique combinations n over k,
    // starting with 0, 1, 2, ...
    // n=3，k=2->{{0,1},{0,2},{1,2}}
    Eigen::ArrayXXi combinations = Eigen::ArrayXXi::Zero(factorial(n) /
        (factorial(k) * factorial(n - k)), k);//->3!/(2!*1!)
    combinations.row(0).setLinSpaced(k, 0, k - 1); //类似于matlab的linespace（start：step：end）

    for (int combination = 1; combination < combinations.rows(); ++combination) {
        combinations.row(combination) = combinations.row(combination - 1);
        for (int col = k - 1; col >= 0; --col) {
            combinations.block(combination, col, 1, k - col).row(0)
                .setLinSpaced(k - col, combinations(combination, col) + 1,
                    combinations(combination, col) + k - col);
            if (combinations(combination, k - 1) < (int)n) {
                break;
            }
        }
    }

    return combinations;
}

Eigen::ArrayXXi distmesh::utils::findUniqueEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation) {
    // find all unique combinations
    auto const combinations = nOverK(triangulation.cols(), 2);//3，2->3

    // find unique edges for all combinations
    // guarantee direction of edges with lower node index to higher index
    std::set<std::array<int, 2>> uniqueEdges;
    std::array<int, 2> edge = { {0, 0} };
    for (int combination = 0; combination < combinations.rows(); ++combination)//0,1,2
        for (int triangle = 0; triangle < triangulation.rows(); ++triangle) {
            edge[0] = triangulation(triangle, combinations(combination, 0));
            edge[1] = triangulation(triangle, combinations(combination, 1));

            edge = edge[1] < edge[0] ? std::array<int, 2>{edge[1], edge[0]} : edge;//排序，小的放前面

            uniqueEdges.insert(edge);
        }

    // copy set to eigen array
    Eigen::ArrayXXi edgeIndices(uniqueEdges.size(), 2);
    int index = 0;
    for (auto const& edge : uniqueEdges) {
        edgeIndices(index, 0) = edge[0];
        edgeIndices(index, 1) = edge[1];

        index++;
    }

    return edgeIndices;
}

Eigen::ArrayXXi distmesh::utils::getTriangulationEdgeIndices(
    Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
    Eigen::Ref<Eigen::ArrayXXi const> const edges) {
    // find indices for each edge of triangulation in edge index array
    Eigen::ArrayXXi edgeIndices(triangulation.rows(), triangulation.cols());
    for (int element = 0; element < triangulation.rows(); ++element)
        for (int node = 0; node < triangulation.cols(); ++node) {
            // create edge with direction from node with lower index
            // to node with higher index
            auto const edge = (Eigen::ArrayXi(2) << triangulation(element, node), triangulation(element, (node + 1) % triangulation.cols())).finished();

            // check if edge is in edges list, and get index
            int edgeIndex = 0;
            if (((edges.rowwise() - edge.transpose()).square().rowwise().sum().minCoeff(&edgeIndex) == 0) ||
                ((edges.rowwise() - edge.transpose().reverse()).square().rowwise().sum().minCoeff(&edgeIndex) == 0)) {
                edgeIndices(element, node) = edgeIndex;
            }
        }

    return edgeIndices;
}


// determine boundary edges of given triangulation
Eigen::ArrayXi distmesh::utils::boundEdges(
    Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
    Eigen::Ref<Eigen::ArrayXXi const> const _edges,
    Eigen::Ref<Eigen::ArrayXXi const> const _edgeIndices) {
    // create a new edge list, if none was given
    Eigen::ArrayXXi edges;
    if (_edges.rows() == 0) {
        edges = utils::findUniqueEdges(triangulation);
    }
    else {
        edges = _edges;
    }

    // get edge indices for each triangle in triangulation
    Eigen::ArrayXXi edgeIndices;
    if (_edgeIndices.rows() == 0) {
        edgeIndices = utils::getTriangulationEdgeIndices(triangulation, edges);
    }
    else {
        edgeIndices = _edgeIndices;
    }

    // find edges, which only appear once in triangulation
    std::set<int> uniqueEdges;
    std::vector<int> boundaryEdges;
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle)
        for (int edge = 0; edge < triangulation.cols(); ++edge) {
            auto const edgeIndex = edgeIndices(triangle, edge);

            // insert edge in set to get info about multiple appearance
            if (!std::get<1>(uniqueEdges.insert(edgeIndex))) {
                // find edge in vector and delete it
                auto const it = std::find(boundaryEdges.begin(), boundaryEdges.end(), edgeIndex);
                if (it != boundaryEdges.end()) {
                    boundaryEdges.erase(it);
                }
            }
            else {
                boundaryEdges.push_back(edgeIndex);
            }
        }

    // convert stl vector to eigen array
    Eigen::ArrayXi boundary(boundaryEdges.size());
    for (int edge = 0; edge < boundary.rows(); ++edge) {
        boundary(edge) = boundaryEdges[edge];
    }

    return boundary;
}

// fix orientation of edges located at the boundary
Eigen::ArrayXXi distmesh::utils::fixBoundaryEdgeOrientation(
    Eigen::Ref<Eigen::ArrayXXd const> const nodes,
    Eigen::Ref<Eigen::ArrayXXi const> const triangulation,
    Eigen::Ref<Eigen::ArrayXXi const> const _edges,
    Eigen::Ref<Eigen::ArrayXXi const> const edgeIndices) {
    Eigen::ArrayXXi edges = _edges;

    // for the 2-D case fix orientation of boundary edges
    if (nodes.cols() == 2) {
        auto const boundary = utils::boundEdges(triangulation, edges, edgeIndices);

        for (int edge = 0; edge < boundary.rows(); ++edge) {
            // find get index of element containing boundary edge
            int elementIndex = 0, edgeIndex = 0;
            (edgeIndices - boundary(edge)).square().minCoeff(&elementIndex, &edgeIndex);

            // get index of node not used in edge, but in the triangle
            int nodeIndex = 0;
            for (int node = 0; node < triangulation.cols(); ++node) {
                if ((triangulation(elementIndex, node) != edges(boundary(edge), 0)) &&
                    (triangulation(elementIndex, node) != edges(boundary(edge), 1))) {
                    nodeIndex = node;
                    break;
                }
            }

            // boundary edges with wrong orientation are marked with a negative sign
            auto const v1 = (nodes.row(edges(boundary(edge), 1)) - nodes.row(edges(boundary(edge), 0))).eval();
            auto const v2 = (nodes.row(triangulation(elementIndex, nodeIndex)) - nodes.row(edges(boundary(edge), 1))).eval();
            if (v1(0) * v2(1) - v1(1) * v2(0) < 0.0) {
                edges.row(boundary(edge)) = edges.row(boundary(edge)).reverse().eval();
            }
        }
    }

    return edges;
}

// project points outside of domain back to boundary
void distmesh::utils::projectPointsToBoundary(
    Functional const& Fd, double const h0,
    Eigen::Ref<Eigen::ArrayXXd> points) {
    Eigen::ArrayXd distance = Fd(points);

    // 找到在边界外的点的集合，bool类型，n*1的矩阵
    Eigen::Array<bool, Eigen::Dynamic, 1> outside = distance > 0.0;
    if (outside.any()) {//有一个为true，则为true
        // 计算梯度
        Eigen::ArrayXXd gradient(points.rows(), points.cols());
        Eigen::ArrayXXd deps = Eigen::ArrayXXd::Zero(points.rows(), points.cols());

        for (int dim = 0; dim < points.cols(); ++dim) {
            deps.col(dim).fill(constants::eps * h0);//第dim列用deps填充
            gradient.col(dim) = (Fd(points + deps) - distance) /(constants::eps * h0);
            deps.col(dim).fill(0.0);//更新deps
        }
        // project points back to boundary，即p(ix,:)=p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady]; 
        //int 范围+-2147483647
        /*points -= outside.replicate(points.rows(), points.cols()).select(
            gradient.colwise() * distance / gradient.square().rowwise().sum(), 0.0);*/
        for (int index = 0; index < points.rows(); ++index) { 
            if (outside(index, 0)) {
                points(index, 0) -= distance(index,0) * gradient(index, 0);//x
                points(index, 1) -= distance(index,0) * gradient(index, 1);//y
            }
        }
    }
}

// check whether points lies inside or outside of polygon
Eigen::ArrayXd distmesh::utils::pointsInsidePoly(
    Eigen::Ref<Eigen::ArrayXXd const> const points,
    Eigen::Ref<Eigen::ArrayXXd const> const polygon) {
    Eigen::ArrayXd inside = Eigen::ArrayXd::Zero(points.rows());

    for (int i = 0, j = polygon.rows() - 1; i < polygon.rows(); j = i++) {
        inside = (((points.col(1) < polygon(i, 1)) != (points.col(1) < polygon(j, 1))) &&
            (points.col(0) < (polygon(j, 0) - polygon(i, 0)) * (points.col(1) - polygon(i, 1)) /
                (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))).select(1.0 - inside, inside);
    }

    return inside;
}