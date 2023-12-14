#ifndef ESOLVER_RATIONAL_H
#define ESOLVER_RATIONAL_H

using std::abs;

#include "utils/logging.h"

template <typename I>
class Rational {
   public:
    // Constructors
    constexpr Rational() : numerator_(0), denominator_(1) {}

    constexpr explicit Rational(I n) : numerator_(n), denominator_(1) {}

    constexpr Rational(I n, I d) : numerator_(n), denominator_(d) {}

    // Normal copy constructors and assignment operators
    // Assignment from I
    Rational &operator=(I n);

    // Assign in place
    Rational &assign(I n, I d);

    // Representation
    constexpr I numerator() const { return numerator_; }

    constexpr I denominator() const { return denominator_; }

    // Two-arg arithmetic
    Rational operator+(const Rational &r) const { return Rational<I>(*this) += r; }

    Rational operator-(const Rational &r) const { return Rational<I>(*this) -= r; }

    Rational operator*(const Rational &r) const { return Rational<I>(*this) *= r; }

    Rational operator/(const Rational &r) const { return Rational<I>(*this) /= r; }

    // Arithmetic operators
    Rational &operator+=(const Rational &r) {
        // Optimized performance and less chance to overflow
        I r_num = r.numerator_;
        I r_den = r.denominator_;

        auto g = gcd(denominator_, r_den);
        denominator_ /= g;
        numerator_ = numerator_ * (r_den / g) + r_num * denominator_;
        g = gcd(numerator_, g);
        numerator_ /= g;
        denominator_ *= r_den / g;

        return *this;
    }

    Rational &operator-=(const Rational &r) {
        // Optimized performance and less chance to overflow
        I r_num = r.numerator_;
        I r_den = r.denominator_;

        auto g = gcd(denominator_, r_den);
        denominator_ /= g;
        numerator_ = numerator_ * (r_den / g) - r_num * denominator_;
        g = gcd(numerator_, g);
        numerator_ /= g;
        denominator_ *= r_den / g;

        return *this;
    }

    Rational &operator*=(const Rational &r) {
        this->numerator_ = this->numerator_ * r.numerator_;
        this->denominator_ = this->denominator_ * r.denominator_;
        this->simplify();
        return *this;
    }

    Rational &operator/=(const Rational &r) {
        if (r.numerator_ >= 0) {
            this->numerator_ = this->numerator_ * r.denominator_;
            this->denominator_ = this->denominator_ * r.numerator_;
        } else {
            this->numerator_ = -this->numerator_ * r.denominator_;
            this->denominator_ = -this->denominator_ * r.numerator_;
        }
        this->simplify();
        return *this;
    }

    // Arithmetic with integers
    template <typename II>
    Rational operator+(const II &i) const {
        return {this->numerator_ + this->denominator_ * i, this->denominator_};
    }

    template <typename II>
    Rational operator-(const II &i) const {
        return {this->numerator_ - this->denominator_ * i, this->denominator_};
    }

    template <typename II>
    Rational operator*(const II &i) const {
        Rational<I> result(this->numerator_ * i, this->denominator_);
        result.simplify();
        return result;
    }

    template <typename II>
    Rational operator/(const II &i) const {
        if (i >= 0) {
            Rational<I> result(this->numerator_, this->denominator_ * i);
            result.simplify();
            return result;
        } else {
            Rational<I> result(-this->numerator_, -this->denominator_ * i);
            result.simplify();
            return result;
        }
    }

    // Comparison operators
    bool operator<(const Rational &r) const {
        if (this->denominator_ == 0 and r.denominator_ == 0) {
            return this->numerator_ < r.numerator_;
        }
        return this->numerator_ * r.denominator_ < r.numerator_ * this->denominator_;
    }

    bool operator<=(const Rational &r) const { return !(*this > r); }

    bool operator>(const Rational &r) const {
        if (this->denominator_ == 0 and r.denominator_ == 0) {
            return this->numerator_ > r.numerator_;
        }
        return this->numerator_ * r.denominator_ > r.numerator_ * this->denominator_;
    }

    bool operator>=(const Rational &r) const { return !(*this < r); }

    constexpr bool operator==(const Rational &r) const {
        return (this->numerator_ == r.numerator_) and (this->denominator_ == r.denominator_);
    }

    constexpr bool operator!=(const Rational &r) const {
        return (this->numerator_ != r.numerator_) or (this->denominator_ != r.denominator_);
    }

   private:
    I numerator_;
    I denominator_;

    I gcd(I a, I b) {
        if (a < b) {
            std::swap(a, b);
        }
        while (b != 0) {
            a %= b;
            std::swap(a, b);
        }
        return a;
    }

    void simplify() {
        if (denominator_ == 0) {
            if (numerator_ > 0) {
                numerator_ = 1;
            }
            if (numerator_ < 0) {
                numerator_ = -1;
            }
            return;
        }

        auto a = gcd(abs(numerator_), denominator_);

        numerator_ /= a;
        denominator_ /= a;
    }
};

// Unary operators
template <typename I>
constexpr Rational<I> operator+(const Rational<I> &r);

template <typename I>
Rational<I> operator-(const Rational<I> &r) {
    return {-r.numerator(), r.denominator()};
}

// Reversed order operators
template <typename I, typename II>
inline Rational<I> operator+(II i, const Rational<I> &r) {
    return r + i;
}

template <typename I, typename II>
inline Rational<I> operator-(II i, const Rational<I> &r) {
    return -(r - i);
}

template <typename I, typename II>
inline Rational<I> operator*(II i, const Rational<I> &r) {
    return r * i;
}

template <typename I>
Rational<I> abs(const Rational<I> &r) {
    return {abs(r.numerator()), r.denominator()};
}

template <typename I>
Rational<I> floor(const Rational<I> &r) {
    if (r.numerator() >= 0) {
        return {r.numerator() / r.denominator(), 1};
    } else {
        return {r.numerator() / r.denominator() - 1, 1};
    }
}

template <typename I>
Rational<I> ceil(const Rational<I> &r) {
    // Todo: calculate explicitly
    return floor(Rational<I>(r.numerator() + r.denominator() - 1, r.denominator()));
}

template <typename I>
bool isnan(const Rational<I> &r) {
    return r.numerator() == 0 and r.denominator() == 0;
}

// Input and output
template <typename I>
std::istream &operator>>(std::istream &is, Rational<I> &r) {
    LOG_ERROR(">> operator is yet not supported for Rational scalar\n");
    exit(1);
}

template <typename I>
std::ostream &operator<<(std::ostream &os, const Rational<I> &r) {
    os << r.numerator() << "/" << r.denominator();
    return os;
}

#endif  // ESOLVER_RATIONAL_H
