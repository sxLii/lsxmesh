#pragma once
#ifndef __DISTMESH_H__
#define __DISTMESH_H__

// standard c++ lib
#include <functional>
#include <tuple>
#include <iostream>

// Eigen lib for array handling
#include "Eigen/Core"

// libdistmesh includes
#include "functional.h"
#include "Fd.h"
#include "utils.h"

namespace distmesh {
    // apply the distmesh algorithm
    std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> distmesh(
        Functional const& Fd, double const h0,
        Functional const& Fh = 1.0,
        Eigen::Ref<Eigen::ArrayXXd const> const bbox = utils::boundingBox(2),
        Eigen::Ref<Eigen::ArrayXXd const> const pfix = Eigen::ArrayXXd());
}

#endif