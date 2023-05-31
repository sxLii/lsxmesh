#pragma once
#ifndef __FD_H__
#define __FD_H__

#include "functional.h"

namespace distmesh {
    namespace Fd {
        // creates distance function for a nd rectangular domain,创建矩形域的距离函数
        // Attention: Not a real distance function at the corners of domainm,注意：域的角点处的距离函数不是实数
        // you have to give the corners as fixed points to distmesh algorithm,必须将角点作为固定点提供给DistMesh算法
        Functional rectangular(Eigen::Ref<Eigen::ArrayXXd const> const rectangle);

        // creates the true distance function for a 2d rectangular domain,矩形
        Functional rectangle(Eigen::Ref<Eigen::ArrayXXd const> const rectangle);

        // creates distance function for elliptical domains，创建椭圆域的距离函数
        // Note: not a real distance function but a level function,注意：不是实数距离函数，而是水平函数
        // which is sufficient
        Functional elliptical(Eigen::Ref<Eigen::ArrayXd const> const radii = Eigen::ArrayXd(),
            Eigen::Ref<Eigen::ArrayXd const> const midpoint = Eigen::ArrayXd());

        // creates the true distance function for circular domains，圆
        Functional circular(double const radius = 1.0,
            Eigen::Ref<Eigen::ArrayXd const> const midpoint = Eigen::ArrayXd());

        // creates distance function for a 2d domain described by polygon
        // Attention: Not a real distance function at the corners of domainm
        // you have to give the corners as fixed points to distmesh algorithm，多边形
        Functional polygon(Eigen::Ref<Eigen::ArrayXXd const> const polygon);
    }
}

#endif