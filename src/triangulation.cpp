#include "distmesh.h"
#include "triangulation.h"

// CGAL
#include<CGAL/point_generators_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include<CGAL/Delaunay_triangulation_2.h>

using namespace CGAL;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
typedef Delaunay::Face_iterator Face_iterator;

Eigen::ArrayXXi distmesh::triangulation::delaunay(Eigen::Ref<Eigen::ArrayXXd const> const points) {
	//将points转换为vector<point2>类型
	std::vector<Point>pts;
	for (int i = 0; i < points.rows(); ++i) {
		double x=0.000, y=0.000;
		for (int j = 0; j < points.cols(); ++j) {
			if (j == 0) {
				x = points(i, j);
			}
			else {
				y= points(i, j);
			}
		}
		pts.push_back(Point(x, y));
	}
	
	//通过dt剖分
	Delaunay dt;
	dt.insert(pts.begin(), pts.end());

	//再将剖分后的点转换至triangulation中
	std::vector<int> tri;
	unsigned face_count = 0;
	for (Face_iterator fa = dt.faces_begin(); fa != dt.faces_end(); fa++) {
		face_count++;
		for (int i = 0; i < 3; i++) {
			Point p = fa->vertex(i)->point();
			tri.push_back(find(pts.begin(), pts.end(), p)- pts.begin());
			//std::cout << "(" << p.x() << "," << p.y() << ")";
		}
	}

	//返回的应该是points的索引，以0开始
	Eigen::ArrayXXi triangulation(face_count, points.cols() + 1);
	for (unsigned i = 0; i < face_count; i++) {
		for (unsigned j = 0; j < points.cols() + 1; j++) {
			triangulation(i, j) = tri[i * (points.cols() + 1) + j];
		}
	}
	return triangulation;

}