#include "distmesh.h"

// assignment operator
distmesh::Functional& distmesh::Functional::operator=(
    Functional const& rhs) {
    this->function() = rhs.function();
    return *this;
}
distmesh::Functional& distmesh::Functional::operator=(
    Functional&& rhs) {
    this->function() = std::move(rhs.function());
    return *this;
}

Eigen::ArrayXd distmesh::Functional::operator()(
    Eigen::Ref<Eigen::ArrayXXd const> const points) const {
    return this->function()(points);
}

distmesh::Functional distmesh::Functional::operator-() const {
    auto const func = this->function();;
    return DISTMESH_FUNCTIONAL({
        return -func(points);
        });
}

distmesh::Functional& distmesh::Functional::operator+=(
    Functional const& rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs(points);
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator+=(
    double const rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) + rhs;
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    Functional const& rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs(points);
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator-=(
    double const rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) - rhs;
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    Functional const& rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs(points);
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator*=(
    double const rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) * rhs;
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    Functional const& rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs(points);
        });
    return *this;
}

distmesh::Functional& distmesh::Functional::operator/=(
    double const rhs) {
    auto const func = this->function();
    this->function() = DISTMESH_FUNCTIONAL({
        return func(points) / rhs;
        });
    return *this;
}

distmesh::Functional distmesh::operator+(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs(points);
        });
}

distmesh::Functional distmesh::operator+(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) + rhs;
        });
}

distmesh::Functional distmesh::operator+(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs + rhs(points);
        });
}

distmesh::Functional distmesh::operator-(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs(points);
        });
}

distmesh::Functional distmesh::operator-(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) - rhs;
        });
}

distmesh::Functional distmesh::operator-(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs - rhs(points);
        });
}

distmesh::Functional distmesh::operator*(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs(points);
        });
}

distmesh::Functional distmesh::operator*(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) * rhs;
        });
}

distmesh::Functional distmesh::operator*(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs * rhs(points);
        });
}

distmesh::Functional distmesh::operator/(
    Functional const& lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs(points);
        });
}

distmesh::Functional distmesh::operator/(
    Functional const& lhs, double const rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs(points) / rhs;
        });
}

distmesh::Functional distmesh::operator/(
    double const lhs, Functional const& rhs) {
    return DISTMESH_FUNCTIONAL({
        return lhs / rhs(points);
        });
}

distmesh::Functional distmesh::Functional::min(
    Functional const& rhs) const {
    auto const func = this->function();
    return DISTMESH_FUNCTIONAL({
        return func(points).min(rhs(points));
        });
}

distmesh::Functional distmesh::Functional::max(
    Functional const& rhs) const {
    auto const func = this->function();
    return DISTMESH_FUNCTIONAL({
        return func(points).max(rhs(points));
        });
}

distmesh::Functional distmesh::Functional::abs() const {
    auto const func = this->function();
    return DISTMESH_FUNCTIONAL({
        return func(points).abs();
        });
}

// geometric transform
distmesh::Functional distmesh::Functional::shift(Eigen::Ref<Eigen::ArrayXd const> const offset) const {
    auto const func = this->function();
    return DISTMESH_FUNCTIONAL({
        return func(points.rowwise() - offset.transpose());
        });
}

distmesh::Functional distmesh::Functional::rotate2D(double const angle) const {
    auto const func = this->function();
    return DISTMESH_FUNCTIONAL({
        Eigen::ArrayXXd transformedPoints(points.rows(), points.cols());
        transformedPoints.col(0) = points.col(0) * std::cos(angle) + points.col(1) * std::sin(angle);
        transformedPoints.col(1) = -points.col(0) * std::sin(angle) + points.col(1) * std::cos(angle);

        return func(transformedPoints);
        });
}