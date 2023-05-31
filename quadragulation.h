#pragma once
#ifndef __QUADRANGULATION_H__
#define __QUADRANGULATION_H__

#include "Eigen/Core"

namespace distmesh::quadrangulation {

    // change triangulation to quadrangulation
    Eigen::ArrayXXi tri2quad(Eigen::Ref<Eigen::ArrayXXd const> const points);
}

#endif