#pragma once
#ifndef __FUNCTIONAL_H__
#define __FUNCTIONAL_H__

// macro for easier creation of distmesh lambda functions
#include <functional>
#include "Eigen/Core"

#define DISTMESH_FUNCTIONAL(function_body) \
    (distmesh::Functional([=](Eigen::Ref<Eigen::ArrayXXd const> const points) -> Eigen::ArrayXd \
        function_body))

namespace distmesh {
    // base class of all function expression for allowing easy function arithmetic，所有函数表达式的基类，允许进行简单的函数运算
    class Functional {
    public:
        // function type of Functional callable
        typedef std::function<Eigen::ArrayXd(Eigen::Ref<Eigen::ArrayXXd const> const)> function_t;

        // create class from function type
        Functional(function_t const& func) : function_(func) {}
        Functional(double const constant) :
            Functional(DISTMESH_FUNCTIONAL({
                return Eigen::ArrayXd::Constant(points.rows(),constant);
                })) {}

                // copy constructor
                Functional(Functional const& rhs) : function_(rhs.function()) {}
                Functional(Functional&& rhs) : function_(std::move(rhs.function())) {}

                // assignment operator
                Functional& operator=(Functional const& rhs);
                Functional& operator=(Functional&& rhs);

                // evaluate function by call
                Eigen::ArrayXd operator() (Eigen::Ref<Eigen::ArrayXXd const> const points) const;

                // basic arithmetic operations
                Functional operator+() const { return *this; }
                Functional operator-() const;
                Functional& operator+=(Functional const& rhs);
                Functional& operator+=(double const rhs);
                Functional& operator-=(Functional const& rhs);
                Functional& operator-=(double const rhs);
                Functional& operator*=(Functional const& rhs);
                Functional& operator*=(double const rhs);
                Functional& operator/=(Functional const& rhs);
                Functional& operator/=(double const rhs);
                friend Functional operator+(Functional const& lhs, Functional const& rhs);
                friend Functional operator+(Functional const& lhs, double const rhs);
                friend Functional operator+(double const lhs, Functional const& rhs);
                friend Functional operator-(Functional const& lhs, Functional const& rhs);
                friend Functional operator-(Functional const& lhs, double const rhs);
                friend Functional operator-(double const lhs, Functional const& rhs);
                friend Functional operator*(Functional const& lhs, Functional const& rhs);
                friend Functional operator*(Functional const& lhs, double const rhs);
                friend Functional operator*(double const lhs, Functional const& rhs);
                friend Functional operator/(Functional const& lhs, Functional const& rhs);
                friend Functional operator/(Functional const& lhs, double const rhs);
                friend Functional operator/(double const lhs, Functional const& rhs);

                // mathematical methods
                Functional min(Functional const& rhs) const;
                Functional max(Functional const& rhs) const;
                Functional abs() const;

                // geometric transform
                Functional shift(Eigen::Ref<Eigen::ArrayXd const> const offset) const;
                Functional rotate2D(double const angle) const;

                // accessors
                function_t& function() { return this->function_; }
                function_t const& function() const { return this->function_; }

    private:
        // stores std function
        function_t function_;
    };
}

#endif