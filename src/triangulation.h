#pragma once
#ifndef __TRIANGULATION_H__
#define __TRIANGULATION_H__

#include "Eigen/Core"

namespace distmesh::triangulation {
    // create delsaunay triangulation from points array
    Eigen::ArrayXXi delaunay(Eigen::Ref<Eigen::ArrayXXd const> const points);
}

#endif