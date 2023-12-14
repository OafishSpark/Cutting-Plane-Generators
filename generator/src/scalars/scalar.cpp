#include "../../include/scalars/scalar.h"

#include <iomanip>
#include <sstream>

#if SCALAR == SCALAR_DOUBLE
Scalar str2scalar(std::string& s, std::size_t* pos) { return std::stod(s, pos); }
std::string ScalarToString(const Scalar& scalar, const int width, const int precision) {
    std::ostringstream ost;
    ost << std::setprecision(precision) << std::setw(width) << scalar;
    return ost.str();
}
#elif SCALAR == SCALAR_FLOAT
Scalar str2scalar(std::string& s, std::size_t* pos) { return std::stof(s, pos); }
std::string ScalarToString(const Scalar& scalar, const int width, const int precision) {
    std::ostringstream ost;
    ost << std::setprecision(precision) << std::setw(width) << scalar;
    return ost.str();
}
#elif SCALAR == SCALAR_RATIONAL
Scalar str2scalar(std::string& s, std::size_t* pos) {
    LOG_ERROR("str2scalar is yet not supported for Rational scalar\n");
    exit(1);
}
std::string ScalarToString(const Scalar& scalar, const int width, const int precision) {
    std::ostringstream ost;
    ost << std::setw(width)
        << std::to_string(scalar.numerator()) + "/" + std::to_string(scalar.denominator());
    return ost.str();
}
#endif
