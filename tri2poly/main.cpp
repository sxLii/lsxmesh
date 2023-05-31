#include "quadrangulation.h"

#include "readtxt.h"
#include <iostream>


int main() {
	Eigen::ArrayXXd points= getTemplatePoints<double>("points.txt");
	Eigen::ArrayXXi triangulations = getTemplatePoints<int>("triangulation.txt");
	//std::cout <<"points:" << points << std::endl;
	//std::cout <<"elements:" << triangulations << std::endl;
	Eigen::ArrayXXi edge_mix;
	std::tie(points, edge_mix) = quadrangulation::tri2quad(points, triangulations, false);
	std::cout << "edge_mix:\n" << edge_mix<< "\nedge_mix.size():\n" << edge_mix.size() << std::endl;

	system("pause");
	return 0;
}