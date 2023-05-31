#pragma once
#ifndef __QUADRANGULATION_H__
#define __QUADRANGULATION_H__

#include "Eigen/Core"
#include <iostream>

#include <functional>
#include <tuple>

#define PI (3.1415926535897932346f)
#define recombination_treshold_ (0.10000f)

namespace quadrangulation {
    // 用于第四步，创建一个结构体，用来比较Q-EDGE这一数据结构
    typedef struct qualityquad
    {
        double quality;
        int quad1,quad2,quad3,quad4;
        //对于向量元素是结构体的，可在结构体内部定义比较函数，下面按照id,length,width升序排序。  
        bool operator< (const qualityquad & a)  const
        {
            return quality < a.quality;
        }
    }QualityQuad;

    // 用于返回从小到大排列的edge
    template <typename type>
    std::array<type, 2> makeEdge(const type edge0, const type edge1) {
        std::array<type, 2> edge = { edge0,edge1 };
        edge = edge[1] < edge[0] ? std::array<type, 2>{edge[1], edge[0]} : edge;
        return edge;
    }

    //1. 发现edges，并区分为外部边 n*2、内部边 m*2、quad k*4
    std::tuple<Eigen::ArrayXXi, Eigen::ArrayXXi> createEdges(Eigen::Ref<Eigen::ArrayXXi const> const triangulation);

    //2. 将两个三角形合并为一个四边形，输出四边形节点quad
    Eigen::ArrayXXi merge2triangles(Eigen::Ref<Eigen::ArrayXXi const> const triangulation, Eigen::Ref<Eigen::ArrayXXi const> const edge_interior);

    //3. 计算每一个四边形的质量，导入quad，输出n*1矩阵Q存放质量，质量越大越好
    Eigen::ArrayXd element2quality(Eigen::Ref<Eigen::ArrayXXd const> const points, Eigen::Ref<Eigen::ArrayXXi const> const quadrangulation);
    
    //4. 删除公共边，导入edge数据和Q质量数据，优先删除质量小的edge，edge里减1(导出），delete_edge里减5（判断能不能删）
    std::tuple<Eigen::ArrayXXi, Eigen::ArrayXXi> deleteCommonEdge(Eigen::Ref<Eigen::ArrayXXi const> const quadrangulation, Eigen::Ref<Eigen::ArrayXd const> const quality);

    //5. 导出混合网格的edge
    Eigen::ArrayXXi mixEdge(Eigen::Ref<Eigen::ArrayXXi const> const edge_interior, Eigen::Ref<Eigen::ArrayXXi const> const edge_surface, Eigen::Ref<Eigen::ArrayXXi const> const edge_delete);

    //6. 找到剩下的三角形,导出剩余的tri
    Eigen::ArrayXXi remainTri(Eigen::Ref<Eigen::ArrayXXi const> const edge_mix);

    //7. 在剩下的三角形中添加重心点和edge的中点，将三角形分成三个四边形,导出points，quad
    std::tuple<Eigen::ArrayXXd, Eigen::ArrayXXi> refine(Eigen::Ref<Eigen::ArrayXXd const> const points, 
        Eigen::Ref<Eigen::ArrayXXi const> const new_tri, Eigen::Ref<Eigen::ArrayXXi const> const quad);
    
    //8. 对剩余节点进行光滑处理，每一个节点和其相邻节点进行求平均值的操作
    Eigen::ArrayXXd smooth(Eigen::Ref<Eigen::ArrayXXd const> const new_points, Eigen::Ref<Eigen::ArrayXXi const> const new_edge,Eigen::Ref<Eigen::ArrayXXi const> const edge_surface);
    
    //9. 导出绘图，混合网格即导出edge；如果是纯四边形网格，内置函数将quad转换为edge再导出
    std::tuple< Eigen::ArrayXXd, Eigen::ArrayXXi> tri2quad(Eigen::Ref<Eigen::ArrayXXd const> const points, Eigen::Ref<Eigen::ArrayXXi const> const triangulation, bool ismix);
}

#endif