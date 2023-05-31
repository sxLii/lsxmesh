#include <iostream>
#include "helper.h"

int main() {
    distmesh::helper::HighPrecisionTime time;

    Eigen::ArrayXXd polygon(10, 2);
    polygon << -0.4, -0.5, 0.4, -0.2, 0.4, -0.7,
        1.5, -0.4, 0.9, 0.1, 1.6, 0.8, 0.5, 0.5,
        0.2, 1.0, 0.1, 0.4, -0.7, 0.7;

    // create mesh
    Eigen::ArrayXXd points;
    Eigen::ArrayXXi elements;

    std::tie(points, elements) = distmesh::distmesh(
        distmesh::Fd::polygon(polygon),
        0.1, 1.0, distmesh::utils::boundingBox(2), polygon);

    // print mesh properties and elapsed time
    std::cout << "Created mesh with " << points.rows() << " points and " << elements.rows() <<
        " elements in " << time.elapsed() * 1e3 << " ms." << std::endl;

    //std::tie(points, elements) = distmesh::distmesh(
    //    distmesh::Fd::circular(1.0), 0.2);

    //// print mesh properties and elapsed time
    //std::cout << "Created mesh with " << points.rows() << " points and " << elements.rows() <<
    //    " elements in " << time.elapsed() * 1e3 << " ms." << std::endl;

    // save mesh to file
    distmesh::helper::savetxt<double>(points, "points.txt");
    distmesh::helper::savetxt<int>(elements, "triangulation.txt");

    // plot mesh using python
    //return system("python plot_mesh.py");
    system("pause");
    return 0;
}
