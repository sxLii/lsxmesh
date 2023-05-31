#pragma once
// 存放常数变量
#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <complex>
#include <limits>

namespace distmesh {
    namespace constants {
        // 当最大相对点移动低于阈值时，算法停止,pointsMovementThreshold
        static double const dptol = 1e-3;

        // 当最大相对点移动超过阈值时，更新三角测量,retriangulationThreshold
        static double const ttol = 1e-1;

        // 几何评估中的相对阈值,geometryEvaluationThreshold
        static double const geps = 1e-3;

        // 用欧拉法更新点位的时间步长
        static double const deltaT = 2e-1;

        // 数值微分的步长
        // EPSILON指的是浮点数可表示的最小值（the smallest value）,即matlab的eps,deltaX
        static double const eps = std::sqrt(std::numeric_limits<double>::epsilon());

        // algorithm will be terminated after the maximum number of iterations,
        // when no convergence can be achieved,算法将在最大迭代次数后终止
        static unsigned const maxSteps = 1000;
    }
}

#endif