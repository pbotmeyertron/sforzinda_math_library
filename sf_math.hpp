#pragma once

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Sforzinda Math Library                                                     */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* This library provides a set of mathematical functions and utilities that   */
/* can be used in graphics and game development.                              */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/

/* Standard Includes */
#include <random>
#include <assert.h>
#include <cmath>

/* Sforzinda Includes */
#include "sf_simd.hpp"
#include "sf_base.hpp"

/* SIMDE Includes */
#include "libs/simde/simde/x86/avx512.h"

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Mathematical Constants                                                     */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a number of mathematical constants.                                */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

namespace sf {
namespace math {

/* Machine epsilon */
template<typename T>
constexpr sf_inline T
EPSILON {
    static_cast<T>(std::numeric_limits<T>::epsilon())
};

/* A circle's circumference divided by its radius. */
template<typename T>
constexpr sf_inline T 
PI {
    static_cast<T>(3.141592653589793238462643383279502884)
};

/* The number of radians in a full circle. */
template<typename T>
constexpr sf_inline T
TAU {
    static_cast<T>(6.283185307179586476925286766559005768)
};

/* PI/2 */
template<typename T>
constexpr sf_inline T
PI_DIV_2 {
    static_cast<T>(1.570796326794896619231321691639751442)
};

/* Square root of 2 */
template<typename T>
constexpr sf_inline T
SQRT2 {
    static_cast<T>(1.414213562373095048801688724209698078)
};

/* e */
template<typename T>
constexpr sf_inline T 
E {
    static_cast<T>(2.718281828459045235360287471352662498)
};

/* log_2(e) */
template<typename T>
constexpr sf_inline T
LOG_2_E {
    static_cast<T>(1.442695040888963407359924681001892137)
};

/* log_10(e) */
template<typename T>
constexpr sf_inline T
LOG_10 {
    static_cast<T>(0.434294481903251827651128918916605082)
};

/* log_e(2) */
template<typename T>
constexpr sf_inline T
LN_2 {
    static_cast<T>(0.693147180559945309417232121458176568)
};

/* log_e(10) */
template<typename T>
constexpr sf_inline T
LN_10 {
    static_cast<T>(2.302585092994045684017991454684364208)
};

/* Euler-Mascheroni Constant */
template<typename T>
constexpr sf_inline T 
EULER_MASCHERONI {
    static_cast<T>(0.577215664901532860606512090082402431)
};

/* Golden ratio */
template<typename T>
constexpr sf_inline T 
PHI {
    static_cast<T>(1.618033988749894848204586834365638118)
};

/* Apery's Constant */
template<typename T>
constexpr sf_inline T 
APERY {
    static_cast<T>(1.202056903159594285399738161511449991)
};

/* Catalan's Constant */
template<typename T>
constexpr sf_inline T 
CATALAN {
    static_cast<T>(0.915965594177219015054603514932384110)
};

/* Khinchin's Constant */
template<typename T>
constexpr sf_inline T 
KHINCHIN {
    static_cast<T>(2.685452001065306445309714835481795693)
};

/* Feigenbaum's Constant */
template<typename T>
constexpr sf_inline T 
FEIGENBAUM {
    static_cast<T>(4.669201609102990671853203820466201617)
};

/* Landau's Constant */
template<typename T>
constexpr sf_inline T 
LANDAU {
    static_cast<T>(0.5)
};

/* Glaisher-Kinkelin Constant */
template<typename T>
constexpr sf_inline T 
GLAISHER_KINKELIN {
    static_cast<T>(1.282427129100622636875342568869791727)
};

/* Napier's Constant */
template<typename T>
constexpr sf_inline T 
NAPIER {
    static_cast<T>(2.685452001065306445309714835481795693)
};

/* Planck's Constant */
template<typename T>
constexpr sf_inline T 
PLANCK {
    static_cast<T>(6.626070040e-34)
};

/* Boltzmann's Constant */
template<typename T>
constexpr sf_inline T 
BOLTZMANN {
    static_cast<T>(1.38064852e-23)
};

/* Avogadro's Constant */
template<typename T>
constexpr sf_inline T 
AVOGADRO {
    static_cast<T>(6.022140857e23)
};

/* Faraday's Constant */
template<typename T>
constexpr sf_inline T 
FARADAY {
    static_cast<T>(96485.33289)
};

/* Gas Constant */
template<typename T>
constexpr sf_inline T 
GAS_CONSTANT {
    static_cast<T>(8.3144598)
};

/* Earth's Gravity */
template<typename T>
constexpr sf_inline T
EARTH_GRAVITY {
    static_cast<T>(9.80665)
};

/* Gravitational Constant */
template<typename T>
constexpr sf_inline T 
GRAVITATIONAL_CONSTANT {
    static_cast<T>(6.67408e-11)
};

/* Speed of Light */
template<typename T>
constexpr sf_inline T 
SPEED_OF_LIGHT {
    static_cast<T>(299792458)
};

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Elementary Functions                                                       */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a number of elementary functions, such as trigonometric,           */
/* hyperbolic, exponential, logarithmic, power, and rounding functions.       */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/*-----------------------------*/
/*-----------------------------*/
/* Scalar Elementary Functions */
/*-----------------------------*/
/*-----------------------------*/

/*-----------------*/
/* Basic Functions */
/*-----------------*/

/* Returns the absolute value of the given value. */
template<typename T>
constexpr sf_inline T
abs(const T& a) {
    return a < 0 ? -a : a;
}

/* Returns the sign of the given value. */
template<typename T>
constexpr sf_inline T
sign(const T& a) {
    return a < 0 ? -1 : 1;
}

/* Returns the minimum of two values. */
template<typename T>
constexpr sf_inline T
min(const T& a, const T& b) {
    return a < b ? a : b;
}

/* Returns the minimum of N values. */
template<typename T>
constexpr sf_inline T 
min(std::initializer_list<T> il) {
    T result = std::numeric_limits<T>::max();
    for (const T& value : il) {
        result = value < result ? value : result;
    }
    return result;
}

/* Returns the maximum of two values. */
template<typename T>
constexpr sf_inline T
max(const T& a, const T& b) {
    return a > b ? a : b;
}

/* Returns the maximum of N values. */
template<typename T>
constexpr sf_inline T
max(std::initializer_list<T> il) {
    T result = std::numeric_limits<T>::min();
    for (const T& value : il) {
        result = value > result ? value : result;
    }
    return result;
}

/* Rounds the given value to the nearest integer. */
template<typename T>
constexpr sf_inline T
round(const T& a) {
    if constexpr (std::is_integral<T>::value) {
        return a;
    } else {
        auto floor_value = static_cast<i64>(a);
        return static_cast<T>(floor_value + (a - floor_value >= T(0.5) ? 1 : 0));
    }
}

/* Rounds the given value to the nearest integer less than or equal to the 
 * given value. */
template<typename T>
constexpr sf_inline T
floor(const T& a) {
    if constexpr (std::is_integral<T>::value) {
        return a;
    } else {
        return static_cast<T>(static_cast<i64>(a));
    }
}

/* Rounds the given value to the nearest integer greater than or equal 
 * to the given value. */
template<typename T>
constexpr sf_inline T
ceil(const T& a) {
    if constexpr (std::is_integral<T>::value) {
        return a;
    } else {
        auto floor_value = static_cast<i64>(a);
        return static_cast<T>(floor_value + (a > floor_value ? 1 : 0));
    }
}

/* Rounds the given value towards zero. */
template<typename T>
constexpr sf_inline T
trunc(const T& a) {
    if constexpr (std::is_integral<T>::value) {
        return a;
    } else {
        return static_cast<T>(a >= T(0) ? 
               static_cast<i64>(a) : 
              -static_cast<i64>(-a));
    }
}

/* Returns the fractional part of the given value. */
template<typename T>
constexpr sf_inline T
fract(const T& a) {
    if constexpr (std::is_integral<T>::value) {
        return T(0);
    } else {
        return a - math::floor(a);
    }
}

/* Returns the mathematical modulo of a divided by b.
 * Note: This differs from C++'s % operator for negative numbers. */
template<typename T>
constexpr sf_inline T
mod(const T& a, const T& b) {
    if constexpr (std::is_integral<T>::value) {
        T m = a % b;
        return m >= T(0) ? m : m + math::abs(b);
    } else {
        return a - b * math::floor(a / b);
    }
}

/* Copies the sign of the second value to the first value. */
template<typename T>
constexpr sf_inline T
copysign(const T& a, const T& b) {
    return math::abs(a) * math::sign(b);
}

/*-----------------------*/
/* Exponential Functions */
/*-----------------------*/

/* Computes the exponential of the given value. */
template<typename T>
constexpr sf_inline T
exp(const T& a) {
    return std::exp(a);
}

/* Computes the natural logarithm of the given value. */
template<typename T>
constexpr sf_inline T
log(const T& a) {
    return std::log(a);
}

/* Computes the base 2 logarithm of the given value. */
template<typename T>
constexpr sf_inline T
log2(const T& a) {
    return std::log2(a);
}

/* Computes the base 10 logarithm of the given value. */
template<typename T>
constexpr sf_inline T
log10(const T& a) {
    return std::log10(a);
}

/* Computes the power of the given value. */
template<typename T>
constexpr sf_inline T 
pow(T base, i32 exponent) {
    return std::pow(base, exponent);
}

/* Computes the power of the given value. */
template<typename T>
constexpr sf_inline T 
pow(T base, T exponent) {
    return std::pow(base, exponent);
}

/* Computes the power of the given value. */
template<typename T, typename U>
constexpr sf_inline T
pow(const T& base, const U& exponent) {
    return math::pow(base, static_cast<T>(exponent));
}

/* Computes the inverse square root of the given value. */
template<typename T>
constexpr sf_inline T
rsqrt(const T& a) {
    if constexpr (std::is_same<T, f32>::value) {
        f32 x2 = a * 0.5f;
        f32 y  = a;
        i32 i  = *(i32*) & y;
        i      = 0x5f3759df - (i >> 1);
        y      = *(f32*) & i;
        y      = y * (1.5f - (x2 * y * y));
        return y;
    } else if constexpr (std::is_same<T, f64>::value) {
        f64 halfx = 0.5 * a;
        i64 i = *(i64*)&a;
        i = 0x5FE6ED2102DCBFDA - (i >> 1);
        f64 y = *(f64*)&i;
        f64 hyhy = halfx * y * y;
        y *= 1.50087895511633457 - hyhy;
        hyhy = halfx * y * y;
        y *= 1.50000057967625766 - hyhy;
        hyhy = halfx * y * y;
        y *= 1.5000000000002520 - hyhy;
        hyhy = halfx * y * y;
        y *= 1.5000000000000000 - hyhy;
        return y;
    } else {
        return static_cast<T>(1) / std::sqrt(a);
    }
}

template<typename T>
constexpr sf_inline T
sqrt(const T& x) {
    return std::sqrt(x);
}

/* Computes the cube root of the given value. */
template<typename T>
constexpr sf_inline T
cbrt(const T& a) {
    return std::cbrt(a);
}

/* Computes the reciprocal of the given value. */
template<typename T>
constexpr sf_inline T
rcp(const T& a) {
    return static_cast<T>(1) / a;
}

/*-------------------------*/
/* Trigonometric Functions */
/*-------------------------*/

/* Virtual lookup table for fast trigonometric functions. */
template<typename T, std::size_t SIN_BITS = 16>
class fast_trig {
private:

    constexpr sf_inline std::size_t SIN_MASK  = (1 << SIN_BITS) - 1;

    constexpr sf_inline std::size_t SIN_COUNT = SIN_MASK + 1;

    constexpr sf_inline T           radian_to_index = 
                                    static_cast<T>(SIN_COUNT) / math::TAU<T>;

    constexpr sf_inline T           degree_to_index = 
                                    static_cast<T>(SIN_COUNT) / 360;

    /* Fast sine table. */
    sf_inline std::array<T, SIN_COUNT> sintable = [] {
        std::array<T, SIN_COUNT> table;
        for (std::size_t i = 0; i < SIN_COUNT; ++i) {
             table[i] = 
             static_cast<T>(std::sin((i + 0.5f) / SIN_COUNT * math::TAU<T>));
        }

        table[0] = 0;
        table[static_cast<std::size_t>(90  * degree_to_index) & SIN_MASK] =  1;
        table[static_cast<std::size_t>(180 * degree_to_index) & SIN_MASK] =  0;
        table[static_cast<std::size_t>(270 * degree_to_index) & SIN_MASK] = -1;
        return table; 
    }();

public:

    constexpr sf_inline T 
    sin(const T& radians) {
        return sintable[static_cast<std::size_t>(radians * 
                                                 radian_to_index) & SIN_MASK];
    }

    constexpr sf_inline T 
    cos(const T& radians) {
        return sintable[static_cast<std::size_t>((radians + 
                                                  math::PI_DIV_2<T>) * 
                                                  radian_to_index) & SIN_MASK];
    }

};

template<typename T>
constexpr sf_inline T
sin(const T& x) {
    return math::fast_trig<T>().sin(x);
}

template<typename T>
constexpr sf_inline T
cos(const T& x) {
    return math::fast_trig<T>().cos(x);
}

/* Computes the tangent of the given value. */
template<typename T>
constexpr sf_inline T
tan(const T& a) {
    return math::sin(a) / math::cos(a);
}

/* Computes the arc tangent of the given value. */
template<typename T>
constexpr sf_inline T
atan(const T& a) {
    return std::atan(a);
}

/* Computes the arc sine of the given value. */
template<typename T>
constexpr sf_inline T
asin(const T& a) {
    return std::asin(a);
}

/* Computes the arc cosine of the given value. */
template<typename T>
constexpr sf_inline T
acos(const T& a) {
    return std::acos(a);
}

/* Computes the arc tangent of the given value. */
template<typename T>
constexpr sf_inline T
atan2(const T& a, const T& b) {
    return std::atan2(a, b);
}

/* Computes the hyperbolic sine of the given value. */
template<typename T>
constexpr sf_inline T
sinh(const T& a) {
    return std::sinh(a);
}

/* Computes the hyperbolic cosine of the given value. */
template<typename T>
constexpr sf_inline T
cosh(const T& a) {
    return std::cosh(a);
}

/* Computes the hyperbolic tangent of the given value. */
template<typename T>
constexpr sf_inline T
tanh(const T& a) {
    T exp_a     = math::exp(a);
    T exp_neg_a = math::exp(-a);
    return (exp_a - exp_neg_a) / 
           (exp_a + exp_neg_a);
}

/* Computes the hyperbolic arc sine of the given value. */
template<typename T>
constexpr sf_inline T
asinh(const T& a) {
    return math::log(a + math::sqrt(a * a + 1));
}

/* Computes the hyperbolic arc cosine of the given value. */
template<typename T>
constexpr sf_inline T
acosh(const T& a) {
    return math::log(a + math::sqrt(a * a - 1));
}

/* Computes the hyperbolic arc tangent of the given value. */
template<typename T>
constexpr sf_inline T
atanh(const T& a) {
    return math::log((1 + a) / (1 - a)) / 2;
}

/* Computes the hyperbolic secant of the given value. */
template<typename T>
constexpr sf_inline T
sech(const T& a) {
    return static_cast<T>(1) / math::cosh(a);
}

/* Computes the hyperbolic cosecant of the given value. */
template<typename T>
constexpr sf_inline T
csch(const T& a) {
    return static_cast<T>(1) / math::sinh(a);
}

/* Computes the hyperbolic cotangent of the given value. */
template<typename T>
constexpr sf_inline T
coth(const T& a) {
    return static_cast<T>(1) / math::tanh(a);
}

/*---------------------------*/
/*---------------------------*/
/* SIMD Elementary Functions */
/*---------------------------*/
/*---------------------------*/

/*----------------*/
/* Infinity & NaN */
/*----------------*/

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N> 
inf() {
    #if defined(__clang__)
        return scl::simd<T, N>(std::numeric_limits<T>::infinity());
    #else
        std::array<T, N> result;
        result.fill(std::numeric_limits<T>::infinity());
        return scl::simd<T, N>(result);
    #endif
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N> 
nan() {
    #if defined(__clang__)
        return scl::simd<T, N>(std::numeric_limits<T>::quiet_NaN());
    #else
        std::array<T, N> result;
        result.fill(std::numeric_limits<T>::quiet_NaN());
        return scl::simd<T, N>(result);
    #endif
}

/*-----------------*/
/* Basic Functions */
/*-----------------*/

/* Computes the absolute value of the given vector. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
abs(const scl::simd<T, N>& a) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::abs(a[i]);
    }
    return result;
}

/* Computes the minimum of two vectors using SIMD select. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
min(const scl::simd<T, N>& a, const scl::simd<T, N>& b) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::min(a[i], b[i]);
    }
    return result;
}

/* Computes the maximum of two vectors using SIMD select. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
max(const scl::simd<T, N>& a, const scl::simd<T, N>& b) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::max(a[i], b[i]);
    }
    return result;
}

/* Rounds the given vector to the nearest integer. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
round(const scl::simd<T, N>& a) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::round(a[i]);
    }
    return result;
}

/* Rounds the given vector to the nearest integer. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
floor(const scl::simd<T, N>& a) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::floor(a[i]);
    }
    return result;
}

/* Rounds the given vector to the nearest integer. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
ceil(const scl::simd<T, N>& a) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::ceil(a[i]);
    }
    return result;
}

/* Rounds the given vector to the nearest integer. */
template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
trunc(const scl::simd<T, N>& a) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = static_cast<T>(a[i]);
    }
    return result;
}

template<typename T, size_t N>
constexpr sf_inline scl::simd<T, N>
mod(const scl::simd<T, N>& a, const scl::simd<T, N>& b) {
    scl::simd<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = std::modf(a[i], b[i]);
    }
    return result;
}

/* Fused multiply-add operation. */
template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N> 
fma(const scl::simd<T, N>& a, const scl::simd<T, N>& b, const scl::simd<T, N>& c) {
    #if defined(__clang__)
        return scl::simd<T, N>(a.data * b.data + c.data);
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = std::fma(a.data[i], b.data[i], c.data[i]);
        }
        return scl::simd<T, N>(result);
    #endif
}

/* Multiply and inverse subtract. */
template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
fms(const scl::simd<T, N>& a, const scl::simd<T, N>& b, const scl::simd<T, N>& c) {
    #if defined(__clang__)
        return scl::simd<T, N>(a.data * b.data - c.data);
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = a.data[i] * b.data[i] - c.data[i];
        }
        return scl::simd<T, N>(result);
    #endif
}

/* Negated multiply-add */
template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N> 
fnma(const scl::simd<T, N>& a, const scl::simd<T, N>& b, const scl::simd<T, N>& c) {
    #if defined(__clang__)
        return scl::simd<T, N>(-(a.data * b.data) + c.data);
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = -(a.data[i] * b.data[i]) + c.data[i];
        }
        return scl::simd<T, N>(result);
    #endif
}

/* Negated multiply-subtract */
template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
fnms(const scl::simd<T, N>& a, const scl::simd<T, N>& b, const scl::simd<T, N>& c) {
    #if defined(__clang__)
        return scl::simd<T, N>(-(a.data * b.data) - c.data);
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = -(a.data[i] * b.data[i]) - c.data[i];
        }
        return scl::simd<T, N>(result);
    #endif
}


/*-----------------------*/
/* Polynomial Evaluation */
/*-----------------------*/

/* Polynomial evaluation using Horner's rule. */
template<typename T, std::size_t N, typename... Coeffs>
constexpr sf_inline scl::simd<T, N> 
polynomial(const scl::simd<T, N>& x, Coeffs... c) {
    std::array<T, sizeof...(c)> arr { static_cast<T>(c)... };
    scl::simd<T, N> result(arr.back());

    /* Horner's rule in descending order. */
    for (std::int32_t i  = static_cast<std::int32_t>(arr.size()) - 2; 
                      i >= 0; 
                    --i) {
        // result = result * x + arr[i];
        result = fma(result, x, scl::simd<T, N>(arr[i]));
    }
    return result;
}

/*-----------------------*/
/* Exponential Functions */
/*-----------------------*/

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
exp(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_exp_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_exp_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_exp_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_exp_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_exp_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_exp_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_exp_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_exp_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_exp_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_exp_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_exp_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_exp_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::exp(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
log(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_log_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_log_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_log_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_log_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_log_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_log_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_log_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_log_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_log_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_log_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_log_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_log_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::log(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
log2(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_log2_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_log2_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_log2_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_log2_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_log2_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_log2_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_log2_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_log2_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_log2_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_log2_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_log2_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_log2_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::log2(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
log10(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_log10_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_log10_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_log10_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_log10_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_log10_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_log10_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_log10_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_log10_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_log10_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_log10_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_log10_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_log10_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::log10(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
pow(const scl::simd<T, N>& x, const scl::simd<T, N>& y) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_pow_ps(x.data, y.data));
            #else
                simde__m128 input1 = std::bit_cast<simde__m128>(x.data);
                simde__m128 input2 = std::bit_cast<simde__m128>(y.data);
                simde__m128 result = simde_mm_pow_ps(input1, input2);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_pow_pd(x.data, y.data));
            #else
                simde__m128d input1 = std::bit_cast<simde__m128d>(x.data);
                simde__m128d input2 = std::bit_cast<simde__m128d>(y.data);
                simde__m128d result = simde_mm_pow_pd(input1, input2);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_pow_ps(x.data, y.data));
            #else
                simde__m256 input1 = std::bit_cast<simde__m256>(x.data);
                simde__m256 input2 = std::bit_cast<simde__m256>(y.data);
                simde__m256 result = simde_mm256_pow_ps(input1, input2);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_pow_pd(x.data, y.data));
            #else
                simde__m256d input1 = std::bit_cast<simde__m256d>(x.data);
                simde__m256d input2 = std::bit_cast<simde__m256d>(y.data);
                simde__m256d result = simde_mm256_pow_pd(input1, input2);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_pow_ps(x.data, y.data));
            #else
                simde__m512 input1 = std::bit_cast<simde__m512>(x.data);
                simde__m512 input2 = std::bit_cast<simde__m512>(y.data);
                simde__m512 result = simde_mm512_pow_ps(input1, input2);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_pow_pd(x.data, y.data));
            #else
                simde__m512d input1 = std::bit_cast<simde__m512d>(x.data);
                simde__m512d input2 = std::bit_cast<simde__m512d>(y.data);
                simde__m512d result = simde_mm512_pow_pd(input1, input2);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::pow(x[i], y[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
pow(const scl::simd<T, N>& x, const T& y) {
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::pow(x[i], y);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
sqrt(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_sqrt_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_sqrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_sqrt_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_sqrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_sqrt_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_sqrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_sqrt_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_sqrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_sqrt_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_sqrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_sqrt_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_sqrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::sqrt(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
cbrt(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_cbrt_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_cbrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_cbrt_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_cbrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_cbrt_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_cbrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_cbrt_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_cbrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_cbrt_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_cbrt_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_cbrt_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_cbrt_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::cbrt(x[i]);
    }
    return result;
}

/*-------------------------*/
/* Trigonometric Functions */
/*-------------------------*/

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
sin(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_sin_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_sin_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_sin_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_sin_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_sin_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_sin_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_sin_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_sin_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_sin_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_sin_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_sin_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_sin_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::sin(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
cos(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_cos_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_cos_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_cos_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_cos_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_cos_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_cos_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_cos_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_cos_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_cos_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_cos_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_cos_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_cos_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    // Fallback scalar implementation
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::cos(x[i]);
    }
    return result;
}

template<typename T, std::size_t N>
constexpr sf_inline scl::simd<T, N>
tan(const scl::simd<T, N>& x) {
    if constexpr (N <= 4) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm_tan_ps(x.data));
            #else
                simde__m128 input = std::bit_cast<simde__m128>(x.data);
                simde__m128 result = simde_mm_tan_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm_tan_pd(x.data));
            #else
                simde__m128d input = std::bit_cast<simde__m128d>(x.data);
                simde__m128d result = simde_mm_tan_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 8) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm256_tan_ps(x.data));
            #else
                simde__m256 input = std::bit_cast<simde__m256>(x.data);
                simde__m256 result = simde_mm256_tan_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm256_tan_pd(x.data));
            #else
                simde__m256d input = std::bit_cast<simde__m256d>(x.data);
                simde__m256d result = simde_mm256_tan_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    } else if constexpr (N <= 16) {
        if constexpr (std::is_same_v<T, f32>) {
            #if defined(__clang__)
                return scl::simd<f32, N>(simde_mm512_tan_ps(x.data));
            #else
                simde__m512 input = std::bit_cast<simde__m512>(x.data);
                simde__m512 result = simde_mm512_tan_ps(input);
                return scl::simd<f32, N>(std::bit_cast<std::array<f32, N>>(result));
            #endif
        } else if constexpr (std::is_same_v<T, f64>) {
            #if defined(__clang__)
                return scl::simd<f64, N>(simde_mm512_tan_pd(x.data));
            #else
                simde__m512d input = std::bit_cast<simde__m512d>(x.data);
                simde__m512d result = simde_mm512_tan_pd(input);
                return scl::simd<f64, N>(std::bit_cast<std::array<f64, N>>(result));
            #endif
        }
    }
    
    scl::simd<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = math::tan(x[i]);
    }
    return result;
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Generic Mathematical Utilities                                             */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/* Performs equality check using machine-epsilon. */
template<typename T>
constexpr sf_inline bool
is_equal(const T& a, const T& b) { 
    return math::abs(a - b) < std::numeric_limits<T>::epsilon();
}

/* Performs check for inequality using machine-epsilon. */
template<typename T>
constexpr sf_inline bool
is_not_equal(const T& a, const T& b) { 
    return math::abs(a - b) >= std::numeric_limits<T>::epsilon();
}

/* Performs check for less than using machine-epsilon. */
template<typename T>
constexpr sf_inline bool
is_less(const T& a, const T& b) { 
    return a < b && math::abs(a - b) >= std::numeric_limits<T>::epsilon();
}

/* Performs check for greater than using machine-epsilon. */
template<typename T>
constexpr sf_inline bool
is_greater(const T& a, const T& b) { 
    return a > b && math::abs(a - b) >= std::numeric_limits<T>::epsilon();
}

/* Radians to degrees. */
template<typename T>
constexpr sf_inline T
degrees(const T& x) { 
    return x * 180.0 / PI<T>; 
}

/* Degrees to radians. */
template<typename T>
constexpr sf_inline T
radians(const T& x) { 
    return x * PI<T> / 180.0; 
}

/* Saturates a value to the range [0, 1]. */
template<typename T>
constexpr sf_inline T
saturate(const T& x) { 
    return x < 0 ? 0 : (x > 1 ? 1 : x); 
}

/* Mutliplies a value by itself. */
template<typename T>
constexpr sf_inline T
square(const T& x) { 
    return x * x; 
}

/* Mutliplies a value by itself thrice. */
template<typename T>
constexpr sf_inline T
cubic(const T& x) { 
    return x * x * x; 
}

/* Mutliplies a value by itself frice. */
template<typename T>
constexpr sf_inline T
quartic(const T& x) { 
    return x * x * x * x; 
}

/* Mutliplies a value by itself five-rice. */
template<typename T>
constexpr sf_inline T
quintic(const T& x) { 
    return x * x * x * x * x; 
}

/* Clamps a value to the range [min, max]. */
template<typename T>
constexpr sf_inline T
clamp(const T& x, const T& min, const T& max) { 
    return x < min ? min : (x > max ? max : x); 
}

template<typename T, typename U>
constexpr sf_inline T
clamp(const T& x, const U& min, const U& max) { 
    return x < min ? min : (x > max ? max : x); 
}

/* Linearly interpolates between two values. */
template<typename T, typename U>
constexpr sf_inline T
lerp(const T& a, const T& b, const U& t) { 
    return a + (b - a) * t; 
}

/* Linearly interpolates between two values with smoothing. */
template<typename T>
constexpr sf_inline T
lerp_smooth(const T& a, const T& b, const T& t) { 
    return a + (b - a) * (3.0 - 2.0 * t) * t * t; 
}

/* Cubic interpolation between two values. */
template<typename T>
constexpr sf_inline T
cubic_interpolation(const T& a, const T& b, const T& c, 
                                const T& d, const T& t) { 
    T p = (d - c) - (a - b);
    T q = (a - b) - p;
    T r = c - a;
    T s = b;
    return p * t * t * t + q * t * t + r * t + s;
}

/* Cosine interpolation between two values. */
template<typename T>
constexpr sf_inline T
cosine_interpolation(const T& a, const T& b, const T& t) { 
    T ft = t * PI<T>;
    T f = (1.0 - math::cos(ft)) * 0.5;
    return a * (1.0 - f) + b * f;
}

/* Hermite interpolation between two values. */
template<typename T>
constexpr sf_inline T
hermite_interpolation(const T& a, const T& b, const T& c, 
                                  const T& d, const T& t) { 
    T t2 = t * t;
    T t3 = t2 * t;
    T m0 = (c - a) * 0.5;
    T m1 = (d - b) * 0.5;
    T a0 = 2.0 * t3 - 3.0 * t2 + 1.0;
    T a1 = t3 - 2.0 * t2 + t;
    T a2 = t3 - t2;
    T a3 = -2.0 * t3 + 3.0 * t2;
    return a0 * b + a1 * m0 + a2 * m1 + a3 * c;
}

/* Calculates the factorial of a given value. */
template<typename T>
constexpr sf_inline T
factorial(T x) { 
    return x == 0 ? 1 : x * math::factorial(x - 1); 
}

/* Calculates the binomial coefficient of a given value. */
template<typename T>
constexpr sf_inline T
binomial(T n, T k) { 
    return math::factorial(n) / 
          (math::factorial(k) * 
           math::factorial(n - k)); 
}

/* Calculates the greatest common divisor of two given values. */
template<typename T>
constexpr sf_inline T
gcd(T a, T b) { 
    return b == 0 ? a : math::gcd(b, a % b); 
}

/* Calculates the least common multiple of two given values. */
template<typename T>
constexpr sf_inline T
lcm(T a, T b) { 
    return (a * b) / math::gcd(a, b); 
}

/* Calculates the greatest common divisor of a given array of values. */
template<typename T>
constexpr sf_inline T
gcd_array(const T* a, std::size_t n) { 
    T result = a[0];
    for (std::size_t i = 1; i < n; ++i) {
        result = math::gcd(result, a[i]);
    }
    return result;
}

/* Calculates the least common multiple of a given array of values. */
template<typename T>
constexpr sf_inline T
lcm_array(const T* a, std::size_t n) { 
    T result = a[0];
    for (std::size_t i = 1; i < n; ++i) {
        result = math::lcm(result, a[i]);
    }
    return result;
}

/* Quadratic formula. */
template<typename T>
constexpr sf_inline T
quadratic(const T& a, const T& b, const T& c, const T& x) { 
    return ((b * b) - (4 * a * c)) / (2 * a * x);
}

/* Calculates the sum of a given array of n values. */
template<typename T>
constexpr sf_inline T
horizontal_sum(const T* a, std::size_t n) { 
    return n == 1 ? a[0] : a[0] + math::horizontal_sum(a + 1, n - 1);
}

/* Calculates the product of a given array of n values. */
template<typename T>
constexpr sf_inline T
horizontal_product(const T* a, std::size_t n) { 
    return n == 1 ? a[0] : a[0] * math::horizontal_product(a + 1, n - 1);
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Easing Functions                                                           */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

namespace easing {

/* Linear easing function. */
template<typename T>
constexpr sf_inline T
linear(const T& t) { 
    return t; 
}

template<typename T>
constexpr sf_inline T
quadratic_in(const T& t) {
    return math::square(t);
}

template<typename T>
constexpr sf_inline T
quadratic_out(const T& t) { 
    return t * (2.0 - t); 
}

template<typename T>
constexpr sf_inline T
quadratic_in_out(const T& t) {
    return (t < 0.5) ? (2.0 * math::square(t))
                     : (1.0 - math::square(-2.0 * t + 2.0) / 2.0);
}

template<typename T>
constexpr sf_inline T
cubic_in(const T& t) { 
    return math::cubic(t); 
}

template<typename T>
constexpr sf_inline T
cubic_out(const T& t) {
    return 1.0 - math::cubic(1.0 - t);
}

template<typename T>
constexpr sf_inline T
cubic_in_out(const T& t) {
    return (t < 0.5) ? (4.0 * math::cubic(t))
                     : (1.0 - math::cubic(-2.0 * t + 2.0) / 2.0);
}

template<typename T>
constexpr sf_inline T
quartic_in(const T& t) {
    return math::quartic(t);
}

template<typename T>
constexpr sf_inline T
quartic_out(const T& t) {
    return 1.0 - math::quartic(1.0 - t);
}

template<typename T>
constexpr sf_inline T
quartic_in_out(const T& t) {
    return (t < 0.5) ? (8.0 * math::quartic(t))
                     : (1.0 - math::quartic(-2.0 * t + 2.0) / 2.0);
}

template<typename T>
constexpr sf_inline T
quintic_in(const T& t) {
    return math::quintic(t);
}

template<typename T>
constexpr sf_inline T
quintic_out(const T& t) {
    return 1.0 - math::quintic(1.0 - t);
}

template<typename T>
constexpr sf_inline T
quintic_in_out(const T& t) {
    return (t < 0.5) ? (16.0 * math::quintic(t))
                     : (1.0  - math::quintic(-2.0 * t + 2.0) / 2.0);
}

template<typename T>
constexpr sf_inline T
sine_in(const T& t) {
    return 1.0 - math::cos((t * PI<T>) / 2.0);
}

template<typename T>
constexpr sf_inline T
sine_out(const T& t) {
    return math::sin((t * PI<T>) / 2.0);
}

template<typename T>
constexpr sf_inline T
sine_in_out(const T& t) {
    return (1.0 - math::cos(PI<T> * t)) / 2.0;
}

template<typename T>
constexpr sf_inline T
circular_in(const T& t) {
    return 1.0 - math::square(1.0 - math::square(t));
}

template<typename T>
constexpr sf_inline T
circular_out(const T& t) {
    return math::square(1.0 - math::square(t));
}

template<typename T>
constexpr sf_inline T
circular_in_out(const T& t) {
    return (t < 0.5) ? (1.0 - math::square(1.0 - 2.0 * t)) / 2.0
                     : (1.0 + math::square(2.0 * t - 1.0)) / 2.0;
}

template<typename T>
constexpr sf_inline T
exponential_in(const T& t) {
    return (t == 0.0) ? t : math::pow(2.0, 10.0 * (t - 1.0));
}

template<typename T>
constexpr sf_inline T
exponential_out(const T& t) {
    return (t == 1.0) ? t : 1.0 - math::pow(2.0, -10.0 * t);
}

template<typename T>
constexpr sf_inline T
exponential_in_out(const T& t) {
    return (t == 0.0 || 
            t == 1.0) ? 
            t : 
           (t <  0.5) ? 
        math::pow(2.0,  20.0 * t - 10.0)  / 2.0 : (2.0 - 
        math::pow(2.0, -20.0 * t + 10.0)) / 2.0;
}

template<typename T>
constexpr sf_inline T
elastic_in(const T& t) {
    return math::sin(13.0 * PI<T> / 2.0 * t) *
           math::pow(2.0, 10.0 * (t - 1.0));
}

template<typename T>
constexpr sf_inline T
elastic_out(const T& t) {
    return math::sin(-13.0 * PI<T> / 2.0 * (t + 1.0)) *
           math::pow(2.0, -10.0 * t) + 1.0;
}

template<typename T>
constexpr sf_inline T
elastic_in_out(const T& t) {
    const T pi_over_2       = PI<T> / 2.0;
    const T constant_a      = 13.0  * pi_over_2;
    const T constant_b      = 10.0;
    const T two_t           = 2.0   * t;
    const T two_t_minus_one = two_t - 1.0;
    if (t < 0.5) {
        return (math::sin(constant_a  * two_t)     * 
                math::pow(2.0, 
                          constant_b  * two_t_minus_one)) / 2.0;
    } else {
        return (math::sin(-constant_a * 
                          (two_t_minus_one + 1.0)) * 
                math::pow(2.0f, 
                          -constant_b * two_t_minus_one) + 2.0) / 2.0;
    }
}

template<typename T>
constexpr sf_inline T
back_in(const T& t) {
    return t * t * t - t * math::sin(t * PI<T>);
}

template<typename T>
constexpr sf_inline T
back_out(const T& t) {
    const T      one_minus_t = 1.0 - t;
    return 1.0 - one_minus_t * 
                 one_minus_t * 
                 one_minus_t -
                 one_minus_t * math::sin(one_minus_t * PI<T>);
}

template<typename T>
constexpr sf_inline T
back_in_out(const T& t) {
    const T two_t = 2.0 * t;
    const T back_constant1 = 2.70158;
    const T back_constant2 = 1.70158;
    if (t < 0.5) {
        return (math::pow(two_t, 2.0)    * 
                         (back_constant1 * two_t - back_constant2)) / 2.0;
    } else {
        const T two_t_minus_two = two_t  - 2.0;
        return (math::pow(two_t_minus_two, 2.0) * 
                         (back_constant1  * 
                          two_t_minus_two + 
                          back_constant2) + 2.0) / 2.0;
    }
}

template<typename T>
constexpr sf_inline T
bounce_out(T t) {
    const T bounce_constant  = 7.5625;
    const T first_threshold  = 1.0 / 2.75;
    const T second_threshold = 2.0 / 2.75;
    const T third_threshold  = 2.5 / 2.75;

    if (t < first_threshold) {
        return bounce_constant * t * t;
    } else if (t < second_threshold) {
        t -= 1.5 / 2.75;
        return bounce_constant * t * t + 0.75;
    } else if (t < third_threshold) {
        t -= 2.25 / 2.75;
        return bounce_constant * t * t + 0.9375;
    } else {
        t -= 2.625 / 2.75;
        return bounce_constant * t * t + 0.984375;
    }
}

template<typename T>
constexpr sf_inline T
bounce_in(const T& t) {
    return 1.0 - easing::bounce_out(1.0 - t);
}

template<typename T>
constexpr sf_inline T
bounce_in_out(const T& t) {
    return (t < 0.5) ? (1.0 - easing::bounce_out(1.0 - 2.0 * t)) / 2.0
                     : (1.0 + easing::bounce_out(2.0 * t - 1.0)) / 2.0;
}

} /* namespace easing */

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Random Number Generation                                                   */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* RNG - splitmix64 - Drop-in replacement for rand()                          */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/* Seed */
template<typename T>
sf_inline T 
rng_seed {
    []{
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<T> dis;
        return dis(gen);
    }()
};

/* Generates a random 64-bit integer. */
template<typename T>
constexpr sf_inline T
rand() {
    T z = (rng_seed<T>    += static_cast<T>(0x9E3779B97F4A7C15));
      z = (z ^ (z >> 30)) *  static_cast<T>(0xBF58476D1CE4E5B9);
      z = (z ^ (z >> 27)) *  static_cast<T>(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

/* Generates a random floating-point number within a specified range. */
template<typename T>
sf_inline T
rand_range(const T& min, const T& max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return static_cast<T>(dis(gen));
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 2D Vector                                                                  */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a 2D vector type and a number of functions for manipulating it.    */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/* 2D Vector */
template<typename T>
class vec2 {
public:

    union {
        struct { 
            T x, y; 
        };
        scl::simd<T, 2> data;
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor. Zero-initialize. */
    vec2() {
        this->x = 0;
        this->y = 0;
    }

    /* Initialize a vector with the given values. */
    vec2(const T& x, const T& y) {
        this->x = x;
        this->y = y;
    }

    /* Initialize a vector with the given value. */
    vec2(const T& value) {
        this->x = value;
        this->y = value;
    }

    /* Copy constructor */
    vec2(const vec2<T>& other) {
        this->x = other.x;
        this->y = other.y;
    }

    /* Construct a vector from an array. */
    vec2(const std::array<T, 2>& other) {
        this->x = other[0];
        this->y = other[1];
    }

    /* Construct a vector from a pointer. */
    vec2(const T* other) {
        this->x = other[0];
        this->y = other[1];
    }

    /*----------------------*/
    /*----------------------*/
    /* Assignment Operators */
    /*----------------------*/
    /*----------------------*/

    vec2&
    operator=(const vec2<T>& other) {
        this->x = other.x;
        this->y = other.y;
        return *this;
    }

    /*--------------------*/
    /*--------------------*/
    /* Indexing Operators */
    /*--------------------*/
    /*--------------------*/

    /* Access the vector using the [] operator. */
    T&
    operator[](std::size_t index) {
        return index == 0 ? this->x : this->y;
    }

    /* Access the vector using the [] operator. */
    const T&
    operator[](std::size_t index) const {
        return index == 0 ? this->x : this->y;
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/
   
    /* Returns the length of the vector. */
    constexpr T 
    length() const {
        return math::sqrt(this->x * this->x + this->y * this->y);
    }
    
    /* Normalize the vector. */
    constexpr vec2 
    normalize() const {
        T length = this->length();
        return vec2(this->x / length, this->y / length);
    }

    /* Returns the dot product of two vectors. */
    constexpr T
    dot_product(const vec2<T>& other) const {
        return this->x * other.x + this->y * other.y;
    }

    /* Rotate the vector by the given angle. */
    constexpr vec2
    rotate(const T& angle) const {
        T s = math::sin(angle);
        T c = math::cos(angle);
        return vec2(this->x * c - this->y * s, 
                    this->x * s + this->y * c);
    }

};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

/* Add two vectors. */
template<typename T>
constexpr sf_inline vec2<T>
operator+(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x + b.x, a.y + b.y);
}

/* Add a vector and a scalar. */
template<typename T>
constexpr sf_inline vec2<T>
operator+(const vec2<T>& a, const T& b) {
    return vec2<T>(a.x + b, a.y + b);
}

/* Add a scalar and a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator+(const T& a, const vec2<T>& b) {
    return vec2<T>(a + b.x, a + b.y);
}

/* Add and assign a vector to another vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator+=(vec2<T>& a, const vec2<T>& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

/* Add and assign a scalar to a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator+=(vec2<T>& a, const T& b) {
    a.x += b;
    a.y += b;
    return a;
}

/* Post-increment a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator++(vec2<T>& a, int) {
    vec2<T> tmp(a);
    a.x++;
    a.y++;
    return tmp;
}

/* Pre-increment a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator++(vec2<T>& a) {
    a.x++;
    a.y++;
    return a;
}

/* Subtract two vectors. */
template<typename T>
constexpr sf_inline vec2<T>
operator-(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x - b.x, a.y - b.y);
}

/* Subtract a vector and a scalar. */
template<typename T>
constexpr sf_inline vec2<T>
operator-(const vec2<T>& a, const T& b) {
    return vec2<T>(a.x - b, a.y - b);
}

/* Subtract a scalar and a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator-(const T& a, const vec2<T>& b) {
    return vec2<T>(a - b.x, a - b.y);
}

/* Subtract and assign a vector to another vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator-=(vec2<T>& a, const vec2<T>& b) {
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

/* Subtract and assign a scalar to a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator-=(vec2<T>& a, const T& b) {
    a.x -= b;
    a.y -= b;
    return a;
}

/* Post-decrement a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator--(vec2<T>& a, int) {
    vec2<T> tmp(a);
    a.x--;
    a.y--;
    return tmp;
}

/* Pre-decrement a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator--(vec2<T>& a) {
    a.x--;
    a.y--;
    return a;
}

/* Negate a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator-(const vec2<T>& a) {
    return vec2<T>(-a.x, -a.y);
}

/* Multiply two vectors. */
template<typename T>
constexpr sf_inline vec2<T>
operator*(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x * b.x, a.y * b.y);
}

/* Multiply a vector and a scalar. */
template<typename T>
constexpr sf_inline vec2<T>
operator*(const vec2<T>& a, const T& b) {
    return vec2<T>(a.x * b, a.y * b);
}

/* Multiply a scalar and a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator*(const T& a, const vec2<T>& b) {
    return vec2<T>(a * b.x, a * b.y);
}

/* Multiply and assign a vector to another vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator*=(vec2<T>& a, const vec2<T>& b) {
    a.x *= b.x;
    a.y *= b.y;
    return a;
}

/* Multiply and assign a scalar to a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator*=(vec2<T>& a, const T& b) {
    a.x *= b;
    a.y *= b;
    return a;
}

/* Divide two vectors. */
template<typename T>
constexpr sf_inline vec2<T>
operator/(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x / b.x, a.y / b.y);
}

/* Divide a vector and a scalar. */
template<typename T>
constexpr sf_inline vec2<T>
operator/(const vec2<T>& a, const T& b) {
    return vec2<T>(a.x / b, a.y / b);
}

/* Divide a scalar and a vector. */
template<typename T>
constexpr sf_inline vec2<T>
operator/(const T& a, const vec2<T>& b) {
    return vec2<T>(a / b.x, a / b.y);
}

/* Divide and assign a vector to another vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator/=(vec2<T>& a, const vec2<T>& b) {
    a.x /= b.x;
    a.y /= b.y;
    return a;
}

/* Divide and assign a scalar to a vector. */
template<typename T>
constexpr sf_inline vec2<T>&
operator/=(vec2<T>& a, const T& b) {
    a.x /= b;
    a.y /= b;
    return a;
}

/*----------------------*/
/*----------------------*/
/* Comparison Operators */ 
/*----------------------*/
/*----------------------*/

/* Compare two vectors for equality. */
template<typename T>
constexpr sf_inline bool
operator==(const vec2<T>& a, const vec2<T>& b) {
    return a.x == b.x && a.y == b.y;
}

/* Compare two vectors for inequality. */
template<typename T>
constexpr sf_inline bool
operator!=(const vec2<T>& a, const vec2<T>& b) {
    return a.x != b.x || a.y != b.y;
}

/* Less than operator. */
template<typename T>
constexpr sf_inline bool
operator<(const vec2<T>& a, const vec2<T>& b) {
    return a.x < b.x && a.y < b.y;
}

/* Greater than operator. */
template<typename T>
constexpr sf_inline bool
operator>(const vec2<T>& a, const vec2<T>& b) {
    return a.x > b.x && a.y > b.y;
}

/* Less than or equal to operator. */
template<typename T>
constexpr sf_inline bool
operator<=(const vec2<T>& a, const vec2<T>& b) {
    return a.x <= b.x && a.y <= b.y;
}

/* Greater than or equal to operator. */
template<typename T>
constexpr sf_inline bool
operator>=(const vec2<T>& a, const vec2<T>& b) {
    return a.x >= b.x && a.y >= b.y;
}

/* Logical negation operator. */
template<typename T>
constexpr sf_inline bool
operator!(const vec2<T>& a) {
    return !a.x && !a.y;
}

/*-------------------*/
/*-------------------*/
/* Bitwise Operators */
/*-------------------*/
/*-------------------*/

/* Bitwise AND operator. */
template<typename T>
constexpr sf_inline vec2<T>
operator&(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x & b.x, a.y & b.y);
}

/* Bitwise AND and assign operator. */
template<typename T>
constexpr sf_inline vec2<T>&
operator&=(vec2<T>& a, const vec2<T>& b) {
    a.x &= b.x;
    a.y &= b.y;
    return a;
}

/* Bitwise OR operator. */
template<typename T>
constexpr sf_inline vec2<T>
operator|(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x | b.x, a.y | b.y);
}

/* Bitwise OR and assign operator. */
template<typename T>
constexpr sf_inline vec2<T>&
operator|=(vec2<T>& a, const vec2<T>& b) {
    a.x |= b.x;
    a.y |= b.y;
    return a;
}

/* Bitwise XOR operator. */
template<typename T>
constexpr sf_inline vec2<T>
operator^(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T>(a.x ^ b.x, a.y ^ b.y);
}

/* Bitwise XOR and assign operator. */
template<typename T>
constexpr sf_inline vec2<T>&
operator^=(vec2<T>& a, const vec2<T>& b) {
    a.x ^= b.x;
    a.y ^= b.y;
    return a;
}

/* Bitwise NOT operator. */
template<typename T>
constexpr sf_inline vec2<T>
operator~(const vec2<T>& a) {
    return vec2<T>(~a.x, ~a.y);
}

/*------------------*/
/*------------------*/
/* Stream Operators */
/*------------------*/
/*------------------*/

/* Output stream operator. */
template<typename T>
std::ostream&
operator<<(std::ostream& os, const vec2<T>& v) {
    os << '{' << v.x << ", " << v.y << '}';
    return os;
}

/* Input stream operator. */
template<typename T>
std::istream&
operator>>(std::istream& is, vec2<T>& v) {
    is >> v.x >> v.y;
    return is;
}

/*------------------*/
/*------------------*/
/* Helper Functions */
/*------------------*/
/*------------------*/

/* Returns the length of a vector. */
template<typename T>
constexpr sf_inline T
length(const vec2<T>& v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

/* Returns the squared length of a vector. */
template<typename T>
constexpr sf_inline T
length_squared(const vec2<T>& v) {
    return v.x * v.x + v.y * v.y;
}

/* Returns the distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance(const vec2<T>& a, const vec2<T>& b) {
    return math::length(a - b);
}

/* Returns the squared distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance_squared(const vec2<T>& a, const vec2<T>& b) {
    return math::length_squared(a - b);
}

/* Returns the dot product of two vectors. */
template<typename T>
constexpr sf_inline T
dot_product(const vec2<T>& a, const vec2<T>& b) {
    return a.x * b.x + a.y * b.y;
}

/* Returns the cross product of two vectors. */
template<typename T>
constexpr sf_inline T
cross_product(const vec2<T>& a, const vec2<T>& b) {
    return a.x * b.y - a.y * b.x;
}

/* Returns the angle between two vectors. */
template<typename T>
constexpr sf_inline T
angle(const vec2<T>& a, const vec2<T>& b) {
    return std::acos(math::dot_product(a, b) / 
                    (math::length(a) * 
                     math::length(b)));
}

/* Linearly interpolates between two vectors. */
template<typename T>
constexpr sf_inline vec2<T>
lerp(const vec2<T>& a, const vec2<T>& b, const T& t) {
    return a + (b - a) * t;
}

/* Returns the normalized vector. */
template<typename T>
constexpr sf_inline vec2<T>
normalize(const vec2<T>& v) {
    return v / math::length(v);
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 3D Vector                                                                  */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a 3D vector type and a number of functions for manipulating it.    */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/* 3D Vector */
template<typename T>
class vec3 {
public:

    union {
        struct { 
            T x, y, z; 
        };
        std::array<T, 3> v;
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor. Zero-initialize. */
    vec3() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }

    /* Initialize a vector with the given values. */
    vec3(const T& x, const T& y, const T& z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    /* Initialize a vector with the given value. */
    vec3(const T& value) {
        this->x = value;
        this->y = value;
        this->z = value;
    }

    /* Copy constructor */
    vec3(const vec3& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
    }

    /* Construct a vector from an array. */
    vec3(const std::array<T, 3>& other) {
        this->x = other[0];
        this->y = other[1];
        this->z = other[2];
    }

    /*----------------------*/
    /*----------------------*/
    /* Assignment Operators */
    /*----------------------*/
    /*----------------------*/

    vec3&
    operator=(const vec3& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        return *this;
    }

    /*--------------------*/
    /*--------------------*/
    /* Indexing Operators */
    /*--------------------*/
    /*--------------------*/

    /* Access the vector by index. */
    T&
    operator[](std::size_t index) {
        return index == 0 ? this->x : 
               index == 1 ? this->y : 
                            this->z;
    }

    /* Access the vector by index. */
    const T&
    operator[](std::size_t index) const {
        return index == 0 ? this->x : 
               index == 1 ? this->y : 
                            this->z;
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    /* Horizontal addition of the vector. */
    constexpr T 
    horizontal_sum() const {
        return this->x + this->y + this->z;
    }

    /* Horizontal multiplication of the vector. */
    constexpr T 
    horizontal_product() const {
        return this->x * this->y * this->z;
    }

    constexpr T
    length() const {
        return math::sqrt(this->x * this->x + 
                          this->y * this->y + 
                          this->z * this->z);
    }
    
    constexpr vec3
    normalize() const {
        T length  = this->length();
        return vec3(this->x / length, 
                    this->y / length, 
                    this->z / length);
    }    

    constexpr vec3 
    rotate_x(const T& angle) {
        return vec3(this->x, 
                    this->y * math::cos(angle) - this->z * math::sin(angle), 
                    this->y * math::sin(angle) + this->z * math::cos(angle));
    }
    
    constexpr vec3 
    rotate_y(const T& angle) {
        return vec3(this->x * math::cos(angle) + this->z * math::sin(angle), 
                    this->y, 
                   -this->x * math::sin(angle) + this->z * math::cos(angle));
    }
    
    constexpr vec3 
    rotate_z(const T& angle) {
        return vec3(this->x * math::cos(angle) - this->y * math::sin(angle), 
                    this->x * math::sin(angle) + this->y * math::cos(angle), 
                    this->z);
    }
};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

/* Addition of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator+(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x + rhs.x, 
                   lhs.y + rhs.y, 
                   lhs.z + rhs.z);
}

/* Addition of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>
operator+(const vec3<T>& lhs, const T& rhs) {
    return vec3<T>(lhs.x + rhs, 
                   lhs.y + rhs, 
                   lhs.z + rhs);
}

/* Addition of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator+(const T& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs + rhs.x, 
                   lhs + rhs.y, 
                   lhs + rhs.z);
}

/* Post-increment of a vector. */
template<typename T>
constexpr sf_inline vec3<T>&
operator++(vec3<T>& lhs) {
    ++lhs.x;
    ++lhs.y;
    ++lhs.z;
    return lhs;
}

/* Pre-increment of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator++(vec3<T>& lhs, int) {
    vec3<T> tmp = lhs;
    ++lhs;
    return tmp;
}

/* Addition and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator+=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

/* Addition and assignment of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>&
operator+=(vec3<T>& lhs, const T& rhs) {
    lhs.x += rhs;
    lhs.y += rhs;
    lhs.z += rhs;
    return lhs;
}

/* Subtraction of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator-(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x - rhs.x, 
                   lhs.y - rhs.y, 
                   lhs.z - rhs.z);
}

/* Subtraction of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>
operator-(const vec3<T>& lhs, const T& rhs) {
    return vec3<T>(lhs.x - rhs, 
                   lhs.y - rhs, 
                   lhs.z - rhs);
}

/* Subtraction of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator-(const T& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs - rhs.x, 
                   lhs - rhs.y, 
                   lhs - rhs.z);
}

/* Post-decrement of a vector. */
template<typename T>
constexpr sf_inline vec3<T>&
operator--(vec3<T>& lhs) {
    --lhs.x;
    --lhs.y;
    --lhs.z;
    return lhs;
}

/* Pre-decrement of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator--(vec3<T>& lhs, int) {
    vec3<T> tmp = lhs;
    --lhs;
    return tmp;
}

/* Subtraction and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator-=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
}

/* Subtraction and assignment of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>&
operator-=(vec3<T>& lhs, const T& rhs) {
    lhs.x -= rhs;
    lhs.y -= rhs;
    lhs.z -= rhs;
    return lhs;
}

/* Negation of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator-(const vec3<T>& lhs) {
    return vec3<T>(-lhs.x, -lhs.y, -lhs.z);
}

/* Multiplication of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator*(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x * rhs.x, 
                   lhs.y * rhs.y, 
                   lhs.z * rhs.z);
}

/* Multiplication of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>
operator*(const vec3<T>& lhs, const T& rhs) {
    return vec3<T>(lhs.x * rhs, 
                   lhs.y * rhs, 
                   lhs.z * rhs);
}

/* Multiplication of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator*(const T& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs * rhs.x, 
                   lhs * rhs.y, 
                   lhs * rhs.z);
}

/* Multiplication and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator*=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    lhs.z *= rhs.z;
    return lhs;
}

/* Multiplication and assignment of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>&
operator*=(vec3<T>& lhs, const T& rhs) {
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;
    return lhs;
}

/* Division of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator/(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x / rhs.x, 
                   lhs.y / rhs.y, 
                   lhs.z / rhs.z);
}

/* Division of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>
operator/(const vec3<T>& lhs, const T& rhs) {
    return vec3<T>(lhs.x / rhs, 
                   lhs.y / rhs, 
                   lhs.z / rhs);
}

/* Division of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator/(const T& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs / rhs.x, 
                   lhs / rhs.y, 
                   lhs / rhs.z);
}

/* Division and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator/=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x /= rhs.x;
    lhs.y /= rhs.y;
    lhs.z /= rhs.z;
    return lhs;
}

/* Division and assignment of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec3<T>&
operator/=(vec3<T>& lhs, const T& rhs) {
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;
    return lhs;
}

/*----------------------*/
/*----------------------*/
/* Comparison Operators */
/*----------------------*/
/*----------------------*/

/* Equality of two vectors. */
template<typename T>
constexpr sf_inline bool
operator==(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x == rhs.x && 
           lhs.y == rhs.y && 
           lhs.z == rhs.z;
}

/* Inequality of two vectors. */
template<typename T>
constexpr sf_inline bool
operator!=(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x != rhs.x || 
           lhs.y != rhs.y || 
           lhs.z != rhs.z;
}

/* Less than of two vectors. */
template<typename T>
constexpr sf_inline bool
operator<(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x < rhs.x && 
           lhs.y < rhs.y && 
           lhs.z < rhs.z;
}

/* Greater than of two vectors. */
template<typename T>
constexpr sf_inline bool
operator>(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x > rhs.x && 
           lhs.y > rhs.y && 
           lhs.z > rhs.z;
}

/* Less than or equal to of two vectors. */
template<typename T>
constexpr sf_inline bool
operator<=(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x <= rhs.x && 
           lhs.y <= rhs.y && 
           lhs.z <= rhs.z;
}

/* Greater than or equal to of two vectors. */
template<typename T>
constexpr sf_inline bool
operator>=(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x >= rhs.x && 
           lhs.y >= rhs.y && 
           lhs.z >= rhs.z;
}

/* Logical and of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator&&(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x && rhs.x, 
                   lhs.y && rhs.y, 
                   lhs.z && rhs.z);
}

/* Logical or of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator||(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x || rhs.x, 
                   lhs.y || rhs.y, 
                   lhs.z || rhs.z);
}

/* Logical not of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator!(const vec3<T>& lhs) {
    return vec3<T>(!lhs.x, 
                   !lhs.y, 
                   !lhs.z);
}

/*-------------------*/
/*-------------------*/
/* Bitwise Operators */
/*-------------------*/
/*-------------------*/

/* Bitwise and of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator&(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x & rhs.x, 
                   lhs.y & rhs.y, 
                   lhs.z & rhs.z);
}

/* Bitwise and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator&=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x &= rhs.x;
    lhs.y &= rhs.y;
    lhs.z &= rhs.z;
    return lhs;
}

/* Bitwise or of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator|(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x | rhs.x, 
                   lhs.y | rhs.y, 
                   lhs.z | rhs.z);
}

/* Bitwise or assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator|=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x |= rhs.x;
    lhs.y |= rhs.y;
    lhs.z |= rhs.z;
    return lhs;
}

/* Bitwise xor of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
operator^(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x ^ rhs.x, 
                   lhs.y ^ rhs.y, 
                   lhs.z ^ rhs.z);
}

/* Bitwise xor assignment of two vectors. */
template<typename T>
constexpr sf_inline vec3<T>&
operator^=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x ^= rhs.x;
    lhs.y ^= rhs.y;
    lhs.z ^= rhs.z;
    return lhs;
}

/* Bitwise not of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator~(const vec3<T>& lhs) {
    return vec3<T>(~lhs.x, 
                   ~lhs.y, 
                   ~lhs.z);
}

/* Bitwise left shift of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator<<(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x << rhs.x, 
                   lhs.y << rhs.y, 
                   lhs.z << rhs.z);
}

/* Bitwise left shift assignment of a vector. */
template<typename T>
constexpr sf_inline vec3<T>&
operator<<=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x <<= rhs.x;
    lhs.y <<= rhs.y;
    lhs.z <<= rhs.z;
    return lhs;
}

/* Bitwise right shift of a vector. */
template<typename T>
constexpr sf_inline vec3<T>
operator>>(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.x >> rhs.x, 
                   lhs.y >> rhs.y, 
                   lhs.z >> rhs.z);
}

/* Bitwise right shift assignment of a vector. */
template<typename T>
constexpr sf_inline vec3<T>&
operator>>=(vec3<T>& lhs, const vec3<T>& rhs) {
    lhs.x >>= rhs.x;
    lhs.y >>= rhs.y;
    lhs.z >>= rhs.z;
    return lhs;
}

/*------------------*/
/*------------------*/
/* Stream Operators */
/*------------------*/
/*------------------*/

/* Output stream operator. */
template<typename T>
std::ostream&
operator<<(std::ostream& os, const vec3<T>& v) {
    os << "{" << v.x << ", " << v.y << ", " << v.z << "}";
    return os;
}

/* Input stream operator. */
template<typename T>
std::istream&
operator>>(std::istream& is, vec3<T>& v) {
    is >> v.x >> v.y >> v.z;
    return is;
}

/*------------------*/
/*------------------*/
/* Vector Functions */
/*------------------*/
/*------------------*/

/* Dot product of two vectors */
template<typename T>
constexpr sf_inline T
dot_product(const vec3<T>& lhs, const vec3<T>& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

/* Cross product of two vectors */
template<typename T>
constexpr sf_inline vec3<T>
cross_product(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(lhs.y * rhs.z - lhs.z * rhs.y, 
                   lhs.z * rhs.x - lhs.x * rhs.z, 
                   lhs.x * rhs.y - lhs.y * rhs.x);
}

/* Returns the length of a vector. */
template<typename T>
constexpr sf_inline T
length(const vec3<T>& v) {
    return math::sqrt(dot_product(v, v));
}

/* Returns the squared length of a vector. */
template<typename T>
constexpr sf_inline T
length_squared(const vec3<T>& v) {
    return dot_product(v, v);
}

/* Returns the distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance(const vec3<T>& lhs, const vec3<T>& rhs) {
    return length(lhs - rhs);
}

/* Returns the squared distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance_squared(const vec3<T>& lhs, const vec3<T>& rhs) {
    return length_squared(lhs - rhs);
}

/* Returns the normalized vector. */
template<typename T>
constexpr sf_inline vec3<T>
normalize(const vec3<T>& v) {
    return v / length(v);
}

/* Returns the angle between two vectors. */
template<typename T>
constexpr sf_inline T
angle(const vec3<T>& lhs, const vec3<T>& rhs) {
    return math::acos(dot_product(lhs, rhs) / (length(lhs) * length(rhs)));
}

/* Returns the reflection vector. */
template<typename T>
constexpr sf_inline vec3<T>
reflect(const vec3<T>& v, const vec3<T>& n) {
    return v - 2 * dot_product(v, n) * n;
}

/* Returns the refraction vector. */
template<typename T>
constexpr sf_inline vec3<T>
refract(const vec3<T>& v, const vec3<T>& n, const T& eta) {
    T k = 1 - eta * eta * (1 - dot_product(v, n) * dot_product(v, n));
    return k < 0 ? vec3<T>(0, 0, 0) : eta * v - 
        (eta * dot_product(v, n) + math::sqrt(k)) * n;
}

/* Returns the vector with the smallest components. */
template<typename T>
constexpr sf_inline vec3<T>
min(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(std::min(lhs.x, rhs.x), 
                   std::min(lhs.y, rhs.y), 
                   std::min(lhs.z, rhs.z));
}

/* Returns the vector with the largest components. */
template<typename T>
constexpr sf_inline vec3<T>
max(const vec3<T>& lhs, const vec3<T>& rhs) {
    return vec3<T>(std::max(lhs.x, rhs.x), 
                   std::max(lhs.y, rhs.y), 
                   std::max(lhs.z, rhs.z));
}

/* Returns the vector with the clamped components. */
template<typename T>
constexpr sf_inline vec3<T>
clamp(const vec3<T>& v, const vec3<T>& min, const vec3<T>& max) {
    return vec3<T>(math::clamp(v.x, min.x, max.x), 
                   math::clamp(v.y, min.y, max.y), 
                   math::clamp(v.z, min.z, max.z));
}

/* Linearly interpolates between two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
lerp(const vec3<T>& lhs, const vec3<T>& rhs, T t) {
    return (1 - t) * lhs + t * rhs;
}

/* Cosine interpolates between two vectors. */
template<typename T>
constexpr sf_inline vec3<T>
cosine_interpolate(const vec3<T>& lhs, const vec3<T>& rhs, const T& t) {
    T t2 = (1 - math::cos(t * PI<T>)) * 0.5;
    return (1 - t2) * lhs + t2 * rhs;
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 4D Vector                                                                  */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a 4D vector type and a number of functions for manipulating it.    */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

/* 4D Vector */
template<typename T>
class vec4 {
public:

    union {
        struct {
            T x;
            T y;
            T z;
            T w;
        };
        std::array<T, 4> v;
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor. Zero-initialize. */
    vec4() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
        this->w = 0;
    }

    /* Initialize a vector with the given values. */
    vec4(const T& x, const T& y, const T& z, const T& w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    /* Initialize a vector with the given value. */
    vec4(const T& value) {
        this->x = value;
        this->y = value;
        this->z = value;
        this->w = value;
    }

    /* Copy constructor */
    vec4(const vec4& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        this->w = other.w;
    }

    /* Construct a vector from an array. */
    vec4(const std::array<T, 4>& other) {
        this->x = other[0];
        this->y = other[1];
        this->z = other[2];
        this->w = other[3];
    }

    /* Convert from vec3 to vec4. */
    inline
    vec4(const vec3<T>& v, const T& w) {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
        this->w = 1.0f;
    }

    /*--------------------*/
    /*--------------------*/
    /* Indexing Operators */
    /*--------------------*/
    /*--------------------*/

    /* Access the vector by index. */
    T&
    operator[](std::size_t index) {
        return index == 0 ? this->x : 
               index == 1 ? this->y : 
               index == 2 ? this->z : 
                            this->w;
    }

    /* Access the vector by index. */
    const T&
    operator[](std::size_t index) const {
        return index == 0 ? this->x : 
               index == 1 ? this->y : 
               index == 2 ? this->z : 
                            this->w;
    }

    /*----------------------*/
    /*----------------------*/
    /* Assignment Operators */
    /*----------------------*/
    /*----------------------*/

    /* Assignment operator */
    vec4&
    operator=(const vec4& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        this->w = other.w;
        return *this;
    }

    /* Assign to a vec3. */
    vec4&
    operator=(const vec3<T>& v) {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
        this->w = 1.0f;
        return *this;
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    /* Horizontal addition of the vector. */
    constexpr T 
    horizontal_sum() const {
        return this->x + this->y + this->z + this->w;
    }

    /* Horizontal multiplication of the vector. */
    T 
    horizontal_product() const {
        return this->x * this->y * this->z * this->w;
    }
    
    /* Returns the length of the vector. */
    constexpr T 
    length() const {
        return math::sqrt(this->x * this->x + 
                          this->y * this->y + 
                          this->z * this->z + 
                          this->w * this->w);
    }
    
    /* Normalize the vector. */
    constexpr vec4<T>
    normalize() const {
        T len = this->length();
        return vec4<T>(this->x / len, this->y / len, 
                       this->z / len, this->w / len);
    }
    
    /* to vec3 */
    vec3<T>
    to_vec3() const {
        return vec3<T>(this->x, this->y, this->z);
    }
};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

/* Addition of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator+(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x + b.x, 
                   a.y + b.y, 
                   a.z + b.z, 
                   a.w + b.w);
}

/* Addition of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator+=(vec4<T>& a, const vec4<T>& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
    return a;
}

/* Addition of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>
operator+(const vec4<T>& a, const T& b) {
    return vec4<T>(a.x + b, 
                   a.y + b, 
                   a.z + b, 
                   a.w + b);
}

/* Addition of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator+(const T& a, const vec4<T>& b) {
    return vec4<T>(a + b.x, 
                   a + b.y, 
                   a + b.z, 
                   a + b.w);
}

/* Addition of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>&
operator+=(vec4<T>& a, const T& b) {
    a.x += b;
    a.y += b;
    a.z += b;
    a.w += b;
    return a;
}

/* Post-increment of a vector. */
template<typename T>
constexpr sf_inline vec4<T>&
operator++(vec4<T>& a) {
    a.x++;
    a.y++;
    a.z++;
    a.w++;
    return a;
}

/* Pre-increment of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator++(vec4<T>& a, int) {
    vec4<T> tmp(a);
    a.x++;
    a.y++;
    a.z++;
    a.w++;
    return tmp;
}

/* Subtraction of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator-(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x - b.x, 
                   a.y - b.y, 
                   a.z - b.z, 
                   a.w - b.w);
}

/* Subtraction of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator-=(vec4<T>& a, const vec4<T>& b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
    return a;
}

/* Subtraction of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>
operator-(const vec4<T>& a, const T& b) {
    return vec4<T>(a.x - b, 
                   a.y - b, 
                   a.z - b, 
                   a.w - b);
}

/* Subtraction of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>&
operator-=(vec4<T>& a, const T& b) {
    a.x -= b;
    a.y -= b;
    a.z -= b;
    a.w -= b;
    return a;
}

/* Subtraction of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator-(const T& a, const vec4<T>& b) {
    return vec4<T>(a - b.x, 
                   a - b.y, 
                   a - b.z, 
                   a - b.w);
}

/* Post-decrement of a vector. */
template<typename T>
constexpr sf_inline vec4<T>&
operator--(vec4<T>& a) {
    a.x--;
    a.y--;
    a.z--;
    a.w--;
    return a;
}

/* Pre-decrement of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator--(vec4<T>& a, int) {
    vec4<T> tmp(a);
    a.x--;
    a.y--;
    a.z--;
    a.w--;
    return tmp;
}

/* Negation of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator-(const vec4<T>& a) {
    return vec4<T>(-a.x, -a.y, -a.z, -a.w);
}

/* Multiplication of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator*(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x * b.x, 
                   a.y * b.y, 
                   a.z * b.z, 
                   a.w * b.w);
}

/* Multiplication of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator*=(vec4<T>& a, const vec4<T>& b) {
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
    a.w *= b.w;
    return a;
}

/* Multiplication of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>
operator*(const vec4<T>& a, const T& b) {
    return vec4<T>(a.x * b, 
                   a.y * b, 
                   a.z * b, 
                   a.w * b);
}

/* Multiplication of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator*(const T& a, const vec4<T>& b) {
    return vec4<T>(a * b.x, 
                   a * b.y, 
                   a * b.z, 
                   a * b.w);
}

/* Multiplication of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>&
operator*=(vec4<T>& a, const T& b) {
    a.x *= b;
    a.y *= b;
    a.z *= b;
    a.w *= b;
    return a;
}

/* Division of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator/(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x / b.x, 
                   a.y / b.y, 
                   a.z / b.z, 
                   a.w / b.w);
}
/* Division of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>
operator/(const vec4<T>& a, const T& b) {
    return vec4<T>(a.x / b, 
                   a.y / b, 
                   a.z / b, 
                   a.w / b);
}

/* Division of a vector and a scalar. */
template<typename T>
constexpr sf_inline vec4<T>&
operator/=(vec4<T>& a, const T& b) {
    a.x /= b;
    a.y /= b;
    a.z /= b;
    a.w /= b;
    return a;
}

/* Division of a scalar and a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator/(const T& a, const vec4<T>& b) {
    return vec4<T>(a / b.x, 
                   a / b.y, 
                   a / b.z, 
                   a / b.w);
}

/*----------------------*/
/*----------------------*/
/* Comparison Operators */
/*----------------------*/
/*----------------------*/

/* Equality comparison of two vectors. */
template<typename T>
constexpr sf_inline bool
operator==(const vec4<T>& a, const vec4<T>& b) {
    return (a.x == b.x) && 
           (a.y == b.y) && 
           (a.z == b.z) && 
           (a.w == b.w);
}

/* Inequality comparison of two vectors. */
template<typename T>
constexpr sf_inline bool
operator!=(const vec4<T>& a, const vec4<T>& b) {
    return (a.x != b.x) || 
           (a.y != b.y) || 
           (a.z != b.z) || 
           (a.w != b.w);
}

/* Less than of two vectors. */
template<typename T>
constexpr sf_inline bool
operator<(const vec4<T>& a, const vec4<T>& b) {
    return (a.x < b.x) && 
           (a.y < b.y) && 
           (a.z < b.z) && 
           (a.w < b.w);
}

/* Greater than of two vectors. */
template<typename T>
constexpr sf_inline bool
operator>(const vec4<T>& a, const vec4<T>& b) {
    return (a.x > b.x) && 
           (a.y > b.y) && 
           (a.z > b.z) && 
           (a.w > b.w);
}

/* Less than or equal to of two vectors. */
template<typename T>
constexpr sf_inline bool
operator<=(const vec4<T>& a, const vec4<T>& b) {
    return (a.x <= b.x) && 
           (a.y <= b.y) && 
           (a.z <= b.z) && 
           (a.w <= b.w);
}

/* Greater than or equal to of two vectors. */
template<typename T>
constexpr sf_inline bool
operator>=(const vec4<T>& a, const vec4<T>& b) {
    return (a.x >= b.x) && 
           (a.y >= b.y) && 
           (a.z >= b.z) && 
           (a.w >= b.w);
}

/* Logical and of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator&&(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x && b.x, 
                   a.y && b.y, 
                   a.z && b.z, 
                   a.w && b.w);
}

/* Logical or of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator||(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x || b.x, 
                   a.y || b.y, 
                   a.z || b.z, 
                   a.w || b.w);
}

/* Logical not of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator!(const vec4<T>& a) {
    return vec4<T>(!a.x, 
                   !a.y, 
                   !a.z, 
                   !a.w);
}

/*-------------------*/
/*-------------------*/
/* Bitwise Operators */
/*-------------------*/
/*-------------------*/

/* Bitwise and of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator&(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x & b.x, 
                   a.y & b.y, 
                   a.z & b.z, 
                   a.w & b.w);
}

/* Bitwise and assignment of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator&=(vec4<T>& a, const vec4<T>& b) {
    a.x &= b.x;
    a.y &= b.y;
    a.z &= b.z;
    a.w &= b.w;
    return a;
}

/* Bitwise or of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator|(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x | b.x, 
                   a.y | b.y, 
                   a.z | b.z, 
                   a.w | b.w);
}

/* Bitwise or assignment of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator|=(vec4<T>& a, const vec4<T>& b) {
    a.x |= b.x;
    a.y |= b.y;
    a.z |= b.z;
    a.w |= b.w;
    return a;
}

/* Bitwise xor of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
operator^(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x ^ b.x, 
                   a.y ^ b.y, 
                   a.z ^ b.z, 
                   a.w ^ b.w);
}

/* Bitwise xor assignment of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>&
operator^=(vec4<T>& a, const vec4<T>& b) {
    a.x ^= b.x;
    a.y ^= b.y;
    a.z ^= b.z;
    a.w ^= b.w;
    return a;
}

/* Bitwise not of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator~(const vec4<T>& a) {
    return vec4<T>(~a.x, 
                   ~a.y, 
                   ~a.z, 
                   ~a.w);
}

/* Bitwise left shift of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator<<(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x << b.x, 
                   a.y << b.y, 
                   a.z << b.z, 
                   a.w << b.w);
}

/* Bitwise left shift assignment of a vector. */
template<typename T>
constexpr sf_inline vec4<T>&
operator<<=(vec4<T>& a, const vec4<T>& b) {
    a.x <<= b.x;
    a.y <<= b.y;
    a.z <<= b.z;
    a.w <<= b.w;
    return a;
}

/* Bitwise right shift of a vector. */
template<typename T>
constexpr sf_inline vec4<T>
operator>>(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.x >> b.x, 
                   a.y >> b.y, 
                   a.z >> b.z, 
                   a.w >> b.w);
}

/* Bitwise right shift assignment of a vector. */
template<typename T>
constexpr sf_inline vec4<T>&
operator>>=(vec4<T>& a, const vec4<T>& b) {
    a.x >>= b.x;
    a.y >>= b.y;
    a.z >>= b.z;
    a.w >>= b.w;
    return a;
}

/*------------------*/
/*------------------*/
/* Stream Operators */
/*------------------*/
/*------------------*/

/* Output stream operator. */
template<typename T>
std::ostream&
operator<<(std::ostream& os, const vec4<T>& v) {
    os << "{" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "}";
    return os;
}

/* Input stream operator. */
template<typename T>
std::istream&
operator>>(std::istream& is, vec4<T>& v) {
    is >> v.x >> v.y >> v.z >> v.w;
    return is;
}

/*------------------*/
/*------------------*/
/* Vector Functions */
/*------------------*/
/*------------------*/

/* Returns the dot product of two vectors. */
template<typename T>
constexpr sf_inline T
dot_product(const vec4<T>& a, const vec4<T>& b) {
    return a.x * b.x + 
           a.y * b.y + 
           a.z * b.z + 
           a.w * b.w;
}

/* Returns the cross product of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
cross_product(const vec4<T>& a, const vec4<T>& b) {
    return vec4<T>(a.y * b.z - a.z * b.y, 
                   a.z * b.x - a.x * b.z, 
                   a.x * b.y - a.y * b.x, 
                   a.w * b.w - a.w * b.w);
}

/* Returns the length of a vector. */
template<typename T>
constexpr sf_inline T
length(const vec4<T>& v) {
    return math::sqrt(v.x * v.x + 
                      v.y * v.y + 
                      v.z * v.z + 
                      v.w * v.w);
}

/* Returns the squared length of a vector. */
template<typename T>
constexpr sf_inline T
length_squared(const vec4<T>& v) {
    return v.x * v.x + 
           v.y * v.y + 
           v.z * v.z + 
           v.w * v.w;
}

/* Returns the distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance(const vec4<T>& a, const vec4<T>& b) {
    return length(a - b);
}

/* Returns the squared distance between two vectors. */
template<typename T>
constexpr sf_inline T
distance_squared(const vec4<T>& a, const vec4<T>& b) {
    return length_squared(a - b);
}

/* Returns the normalized vector. */
template<typename T>
constexpr sf_inline vec4<T>
normalize(const vec4<T>& v) {
    return v / length(v);
}

/* Returns the reflection vector. */
template<typename T>
constexpr sf_inline vec4<T>
reflect(const vec4<T>& v, const vec4<T>& n) {
    return v - 2 * dot_product(v, n) * n;
}

/* Returns the refraction vector. */
template<typename T>
constexpr sf_inline vec4<T>
refract(const vec4<T>& v, const vec4<T>& n, T eta) {
    T k = 1 - eta * eta * (1 - dot_product(v, n) * dot_product(v, n));
    if (k < 0) {
        return vec4<T>(0, 0, 0, 0);
    }
    return eta * v - (eta * dot_product(v, n) + math::sqrt(k)) * n;
}

/* Returns the linear interpolation of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
lerp(const vec4<T>& a, const vec4<T>& b, const T& t) {
    return a + t * (b - a);
}

/* Cosine interpolation of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
cosine_interpolate(const vec4<T>& a, const vec4<T>& b, const T& t) {
    return math::lerp(a, b, (1 - math::cos(t * sf::math::PI<T>)) * 0.5);
}

/* Returns the smooth interpolation of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
smoothstep(const vec4<T>& a, const vec4<T>& b, T t) {
    return math::lerp(a, b, t * t * (3 - 2 * t));
}

/* Returns the smoother interpolation of two vectors. */
template<typename T>
constexpr sf_inline vec4<T>
smootherstep(const vec4<T>& a, const vec4<T>& b, T t) {
    return math::lerp(a, b, t * t * t * (t * (t * 6 - 15) + 10));
}

/* Returns the catmull-rom interpolation of four vectors. */
template<typename T>
constexpr sf_inline vec4<T>
catmull_rom(const vec4<T>& a, const vec4<T>& b, 
            const vec4<T>& c, const vec4<T>& d, const T& t) {
    return 0.5 * (2 * b + (c - a) * t + 
                  (2 * a - 5 * b + 4 * c - d) * t * t + 
                  (3 * b - a - 3 * c + d) * t * t * t);
}

/* Returns the hermite interpolation of four vectors. */
template<typename T>
constexpr sf_inline vec4<T>
hermite(const vec4<T>& a, const vec4<T>& b, 
        const vec4<T>& c, const vec4<T>& d, const T& t) {
    return (2 * t * t * t - 3 * t * t + 1) * a + 
           (t * t * t - 2 * t * t + t) * b + 
           (-2 * t * t * t + 3 * t * t) * c + 
           (t * t * t - t * t) * d;
}

/* Returns the cubic interpolation of four vectors. */
template<typename T>
constexpr sf_inline vec4<T>
cubic(const vec4<T>& a, const vec4<T>& b, 
      const vec4<T>& c, const vec4<T>& d, const T& t) {
    return b + 0.5 * t * (c - a + t * (2 * a - 5 * b + 4 * c - d + 
                                      t * (3 * (b - c) + d - a)));
}

/* Returns the quadratic interpolation of three vectors. */
template<typename T>
constexpr sf_inline vec4<T>
quadratic(const vec4<T>& a, const vec4<T>& b, 
          const vec4<T>& c, const T& t) {
    return a + t * (b - a + t * (2 * a - 5 * b + 4 * c - b));
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 2x2 Matrix                                                                 */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

template<typename T>
class mat2x2 {
public:

    union {
        struct {
            /* Matrix elements */
            T m00, m01;
            T m10, m11;
        };
        math::vec2<T> m[2];
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor is identity matrix. */
    constexpr
    mat2x2() {
        m00 = 1; m01 = 0;
        m10 = 0; m11 = 1;
    }

    constexpr
    mat2x2(const mat2x2& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m10 = m.m10; 
        this->m11 = m.m11;
    }

    constexpr
    mat2x2(const vec2<T>& m0, const vec2<T>& m1) {
        m00 = m0.x; m01 = m0.y;
        m10 = m1.x; m11 = m1.y;
    }

    constexpr
    mat2x2(const T& a00, const T& a01,
           const T& a10, const T& a11) {
        m00 = a00; m01 = a01;
        m10 = a10; m11 = a11;
    }

    /*---------------------------*/
    /*---------------------------*/
    /* Type Conversion Operators */
    /*---------------------------*/
    /*---------------------------*/

    template<typename U>
    constexpr inline
    operator mat2x2<U>() const {
        return mat2x2<U>(static_cast<U>(this->m00), 
                         static_cast<U>(this->m01), 
                         static_cast<U>(this->m10), 
                         static_cast<U>(this->m11));
    }

    /*------------------*/
    /*------------------*/
    /* Member Operators */
    /*------------------*/
    /*------------------*/

    constexpr mat2x2<T>&
    operator=(const mat2x2<T>& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m10 = m.m10; 
        this->m11 = m.m11;
        return *this;
    }

    constexpr vec2<T>&
    operator[](const std::size_t& i) {
        assert(i >= 0 && i < 2);
        return i ? m[1] : m[0];
    }

    constexpr vec2<T>
    operator[](const std::size_t& i) const {
        assert(i >= 0 && i < 2);
        return i ? m[1] : m[0];
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    constexpr mat2x2<T>
    transpose() const {
        return mat2x2<T>(m00, m10, m01, m11);
    }

    constexpr T
    determinant() const {
        return m00 * m11 - m01 * m10;
    }

    constexpr mat2x2<T>
    inverse() const {
        T det = determinant();
        return mat2x2<T>(m11 / det, -m01 / det, 
                        -m10 / det,  m00 / det);
    }

    constexpr mat2x2<T>
    adjugate() const {
        return mat2x2<T>(m11, -m01, 
                        -m10, m00);
    }

};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

template<typename T>
constexpr sf_inline mat2x2<T>
operator+(const mat2x2<T>& a, const mat2x2<T>& b) {
    return mat2x2<T>(a.m00 + b.m00, a.m01 + b.m01, 
                     a.m10 + b.m10, a.m11 + b.m11);
}

template<typename T>
constexpr sf_inline mat2x2<T>
operator-(const mat2x2<T>& a, const mat2x2<T>& b) {
    return mat2x2<T>(a.m00 - b.m00, a.m01 - b.m01, 
                     a.m10 - b.m10, a.m11 - b.m11);
}

template<typename T>
constexpr sf_inline mat2x2<T>
operator*(const mat2x2<T>& a, const mat2x2<T>& b) {
    return mat2x2<T>(a.m00 * b.m00 + a.m01 * b.m10, 
                     a.m00 * b.m01 + a.m01 * b.m11, 
                     a.m10 * b.m00 + a.m11 * b.m10, 
                     a.m10 * b.m01 + a.m11 * b.m11);
}

template<typename T>
constexpr sf_inline mat2x2<T>
operator*(const mat2x2<T>& m, const T& s) {
    return mat2x2<T>(m.m00 * s, m.m01 * s, 
                     m.m10 * s, m.m11 * s);
}

template<typename T>
constexpr sf_inline mat2x2<T>
operator*(const T& s, const mat2x2<T>& m) {
    return mat2x2<T>(m.m00 * s, m.m01 * s, 
                     m.m10 * s, m.m11 * s);
}

template<typename T>
constexpr sf_inline vec2<T>
operator*(const mat2x2<T>& m, const vec2<T>& v) {
    return vec2<T>(m.m00 * v.x + m.m01 * v.y, 
                   m.m10 * v.x + m.m11 * v.y);
}

template<typename T>
constexpr sf_inline mat2x2<T>
operator/(const mat2x2<T>& m, const T& s) {
    return mat2x2<T>(m.m00 / s, m.m01 / s, 
                     m.m10 / s, m.m11 / s);
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 3x3 Matrix                                                                 */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

template<typename T>
class mat3x3 {
public:

    union {
        struct {
            /* Matrix elements */
            T m00, m01, m02;
            T m10, m11, m12;
            T m20, m21, m22;
        };
        math::vec3<T> m[3];
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor is identity matrix. */
    constexpr
    mat3x3() {
        m00 = 1; m01 = 0; m02 = 0;
        m10 = 0; m11 = 1; m12 = 0;
        m20 = 0; m21 = 0; m22 = 1;
    }

    constexpr
    mat3x3(const mat3x3& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m02 = m.m02; 
        this->m10 = m.m10; 
        this->m11 = m.m11; 
        this->m12 = m.m12;
        this->m20 = m.m20;
        this->m21 = m.m21;
        this->m22 = m.m22;
    }

    constexpr
    mat3x3(const vec3<T>& m0, const vec3<T>& m1, const vec3<T>& m2) {
        m00 = m0.x; m01 = m0.y; m02 = m0.z;
        m10 = m1.x; m11 = m1.y; m12 = m1.z;
        m20 = m2.x; m21 = m2.y; m22 = m2.z;
    }

    constexpr
    mat3x3(const T& a00, const T& a01, const T& a02,
           const T& a10, const T& a11, const T& a12,
           const T& a20, const T& a21, const T& a22) {
        m00 = a00; m01 = a01; m02 = a02;
        m10 = a10; m11 = a11; m12 = a12;
        m20 = a20; m21 = a21; m22 = a22;
    }

    /*---------------------------*/
    /*---------------------------*/
    /* Type Conversion Operators */
    /*---------------------------*/
    /*---------------------------*/

    template<typename U>
    constexpr inline
    operator mat3x3<U>() const {
        return mat3x3<U>(static_cast<U>(this->m00), 
                         static_cast<U>(this->m01), 
                         static_cast<U>(this->m02),
                         static_cast<U>(this->m10),
                         static_cast<U>(this->m11),
                         static_cast<U>(this->m12),
                         static_cast<U>(this->m20),
                         static_cast<U>(this->m21),
                         static_cast<U>(this->m22));
    }

    /*------------------*/
    /*------------------*/
    /* Member Operators */
    /*------------------*/
    /*------------------*/

    constexpr mat3x3<T>&
    operator=(const mat3x3<T>& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m02 = m.m02; 
        this->m10 = m.m10; 
        this->m11 = m.m11; 
        this->m12 = m.m12;
        this->m20 = m.m20;
        this->m21 = m.m21;
        this->m22 = m.m22;
        return *this;
    }

    constexpr vec3<T>&
    operator[](const std::size_t& i) {
        assert(i >= 0 && i < 3);
        return i == 0 ? m[0] : (i == 1 ? m[1] : m[2]);
    }

    constexpr vec3<T>
    operator[](const std::size_t& i) const {
        assert(i >= 0 && i < 3);
        return i == 0 ? m[0] : (i == 1 ? m[1] : m[2]);
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    constexpr mat3x3<T>
    transpose() const {
        return mat3x3<T>(m00, m10, m20, 
                         m01, m11, m21, 
                         m02, m12, m22);
    }

    constexpr T
    determinant() const {
        return m00 * (m11 * m22 - m12 * m21) - 
               m01 * (m10 * m22 - m12 * m20) + 
               m02 * (m10 * m21 - m11 * m20);
    }

    constexpr mat3x3<T>
    inverse() const {
        T det = determinant();
        return mat3x3<T>((m11 * m22 - m12 * m21) / det, 
                        -(m01 * m22 - m02 * m21) / det, 
                         (m01 * m12 - m02 * m11) / det,
                        -(m10 * m22 - m12 * m20) / det, 
                         (m00 * m22 - m02 * m20) / det, 
                        -(m00 * m12 - m02 * m10) / det,
                         (m10 * m21 - m11 * m20) / det, 
                        -(m00 * m21 - m01 * m20) / det, 
                         (m00 * m11 - m01 * m10) / det);
    }

    constexpr mat3x3<T>
    adjugate() const {
        return mat3x3<T>( m11 * m22 - m12 * m21, 
                        -(m01 * m22 - m02 * m21), 
                          m01 * m12 - m02 * m11,
                        -(m10 * m22 - m12 * m20), 
                          m00 * m22 - m02 * m20, 
                        -(m00 * m12 - m02 * m10),
                          m10 * m21 - m11 * m20, 
                        -(m00 * m21 - m01 * m20), 
                          m00 * m11 - m01 * m10);
    }

};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

template<typename T>
constexpr sf_inline mat3x3<T>
operator+(const mat3x3<T>& a, const mat3x3<T>& b) {
    return mat3x3<T>(a.m00 + b.m00, a.m01 + b.m01, a.m02 + b.m02, 
                     a.m10 + b.m10, a.m11 + b.m11, a.m12 + b.m12, 
                     a.m20 + b.m20, a.m21 + b.m21, a.m22 + b.m22);
}

template<typename T>
constexpr sf_inline mat3x3<T>
operator-(const mat3x3<T>& a, const mat3x3<T>& b) {
    return mat3x3<T>(a.m00 - b.m00, a.m01 - b.m01, a.m02 - b.m02, 
                     a.m10 - b.m10, a.m11 - b.m11, a.m12 - b.m12, 
                     a.m20 - b.m20, a.m21 - b.m21, a.m22 - b.m22);
}

template<typename T>
constexpr sf_inline mat3x3<T>
operator*(const mat3x3<T>& a, const mat3x3<T>& b) {
    return mat3x3<T>(a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20, 
                     a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21, 
                     a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22, 
                     a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20, 
                     a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21, 
                     a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22, 
                     a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20, 
                     a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21, 
                     a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22);
}

template<typename T>
constexpr sf_inline mat3x3<T>
operator*(const mat3x3<T>& m, const T& s) {
    return mat3x3<T>(m.m00 * s, m.m01 * s, m.m02 * s, 
                     m.m10 * s, m.m11 * s, m.m12 * s, 
                     m.m20 * s, m.m21 * s, m.m22 * s);
}

template<typename T>
constexpr sf_inline mat3x3<T>
operator*(const T& s, const mat3x3<T>& m) {
    return mat3x3<T>(m.m00 * s, m.m01 * s, m.m02 * s, 
                     m.m10 * s, m.m11 * s, m.m12 * s, 
                     m.m20 * s, m.m21 * s, m.m22 * s);
}

template<typename T>
constexpr sf_inline vec3<T>
operator*(const mat3x3<T>& m, const vec3<T>& v) {
    return vec3<T>(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z, 
                   m.m10 * v.x + m.m11 * v.y + m.m12 * v.z, 
                   m.m20 * v.x + m.m21 * v.y + m.m22 * v.z);
}

template<typename T>
constexpr sf_inline mat3x3<T>
operator/(const mat3x3<T>& m, const T& s) {
    return mat3x3<T>(m.m00 / s, m.m01 / s, m.m02 / s, 
                     m.m10 / s, m.m11 / s, m.m12 / s, 
                     m.m20 / s, m.m21 / s, m.m22 / s);
}


/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* 4x4 Matrix                                                                 */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

template<typename T>
class mat4x4 {
public:

    union {
        struct {
            /* Matrix elements */
            T m00, m01, m02, m03;
            T m10, m11, m12, m13;
            T m20, m21, m22, m23;
            T m30, m31, m32, m33;
        };
        math::vec4<T> m[4];
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor is identity matrix. */
    constexpr
    mat4x4() {
        m00 = 1; m01 = 0; m02 = 0; m03 = 0;
        m10 = 0; m11 = 1; m12 = 0; m13 = 0;
        m20 = 0; m21 = 0; m22 = 1; m23 = 0;
        m30 = 0; m31 = 0; m32 = 0; m33 = 1;
    }

    constexpr
    mat4x4(const mat4x4& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m02 = m.m02; 
        this->m03 = m.m03; 
        this->m10 = m.m10; 
        this->m11 = m.m11; 
        this->m12 = m.m12; 
        this->m13 = m.m13;
        this->m20 = m.m20;
        this->m21 = m.m21;
        this->m22 = m.m22;
        this->m23 = m.m23;
        this->m30 = m.m30;
        this->m31 = m.m31;
        this->m32 = m.m32;
        this->m33 = m.m33;
    }

    constexpr
    mat4x4(const vec4<T>& m0, const vec4<T>& m1, 
           const vec4<T>& m2, const vec4<T>& m3) {
        m00 = m0.x; m01 = m0.y; m02 = m0.z; m03 = m0.w;
        m10 = m1.x; m11 = m1.y; m12 = m1.z; m13 = m1.w;
        m20 = m2.x; m21 = m2.y; m22 = m2.z; m23 = m2.w;
        m30 = m3.x; m31 = m3.y; m32 = m3.z; m33 = m3.w;
    }

    constexpr
    mat4x4(const T& a00, const T& a01, const T& a02, const T& a03,
           const T& a10, const T& a11, const T& a12, const T& a13,
           const T& a20, const T& a21, const T& a22, const T& a23,
           const T& a30, const T& a31, const T& a32, const T& a33) {
        m00 = a00; m01 = a01; m02 = a02; m03 = a03;
        m10 = a10; m11 = a11; m12 = a12; m13 = a13;
        m20 = a20; m21 = a21; m22 = a22; m23 = a23;
        m30 = a30; m31 = a31; m32 = a32; m33 = a33;
    }

    /*---------------------------*/
    /*---------------------------*/
    /* Type Conversion Operators */
    /*---------------------------*/
    /*---------------------------*/

    template<typename U>
    constexpr inline
    operator mat4x4<U>() const {
        return mat4x4<U>(static_cast<U>(this->m00), 
                         static_cast<U>(this->m01), 
                         static_cast<U>(this->m02), 
                         static_cast<U>(this->m03),
                         static_cast<U>(this->m10), 
                         static_cast<U>(this->m11), 
                         static_cast<U>(this->m12), 
                         static_cast<U>(this->m13),
                         static_cast<U>(this->m20), 
                         static_cast<U>(this->m21), 
                         static_cast<U>(this->m22), 
                         static_cast<U>(this->m23),
                         static_cast<U>(this->m30), 
                         static_cast<U>(this->m31), 
                         static_cast<U>(this->m32), 
                         static_cast<U>(this->m33));
    }

    /*------------------*/
    /*------------------*/
    /* Member Operators */
    /*------------------*/
    /*------------------*/

    constexpr mat4x4<T>&
    operator=(const mat4x4<T>& m) {
        this->m00 = m.m00; 
        this->m01 = m.m01; 
        this->m02 = m.m02; 
        this->m03 = m.m03; 
        this->m10 = m.m10; 
        this->m11 = m.m11; 
        this->m12 = m.m12; 
        this->m13 = m.m13;
        this->m20 = m.m20;
        this->m21 = m.m21;
        this->m22 = m.m22;
        this->m23 = m.m23;
        this->m30 = m.m30;
        this->m31 = m.m31;
        this->m32 = m.m32;
        this->m33 = m.m33;
        return *this;
    }

    /* Access the matrix by index. */
    constexpr T&
    operator[](const i32& index) {
        assert(index >= 0 && index < 16);
        return (&m00)[index];
    }

    /* Access the matrix by index. */
    constexpr const T&
    operator[](const i32& index) const {
        assert(index >= 0 && index < 16);
        return (&m00)[index];
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    /* Zero matrix. */
    void
    zero(void) {
        m00 = 0; m01 = 0; m02 = 0; m03 = 0;
        m10 = 0; m11 = 0; m12 = 0; m13 = 0;
        m20 = 0; m21 = 0; m22 = 0; m23 = 0;
        m30 = 0; m31 = 0; m32 = 0; m33 = 0;
    }

    /* Identity matrix. */
    void
    identity(void) {
        m00 = 1; m01 = 0; m02 = 0; m03 = 0;
        m10 = 0; m11 = 1; m12 = 0; m13 = 0;
        m20 = 0; m21 = 0; m22 = 1; m23 = 0;
        m30 = 0; m31 = 0; m32 = 0; m33 = 1;
    }

    /* Transpose the matrix. */
    constexpr mat4x4<T>
    transpose(void) const {
        return mat4x4<T>(m00, m10, m20, m30, 
                         m01, m11, m21, m31, 
                         m02, m12, m22, m32, 
                         m03, m13, m23, m33);
    }

    /* Trace of the matrix. */
    constexpr T
    trace(void) const {
        return m00 + m11 + m22 + m33;
    }

    /* Translate the matrix. */
    constexpr void
    translate(const T& tx, const T& ty, const T& tz) {
        m03 += tx;
        m13 += ty;
        m23 += tz;
    }

    /* Scale the matrix. */
    constexpr void
    scale(const T& sx, const T& sy, const T& sz) {
        m00 *= sx;
        m01 *= sx;
        m02 *= sx;
        m03 *= sx;
        m10 *= sy;
        m11 *= sy;
        m12 *= sy;
        m13 *= sy;
        m20 *= sz;
        m21 *= sz;
        m22 *= sz;
        m23 *= sz;
    }

    constexpr void
    scale(const vec3<T>& s) {
        m00 *= s.x;
        m01 *= s.x;
        m02 *= s.x;
        m03 *= s.x;
        m10 *= s.y;
        m11 *= s.y;
        m12 *= s.y;
        m13 *= s.y;
        m20 *= s.z;
        m21 *= s.z;
        m22 *= s.z;
        m23 *= s.z;
    }

    /* Determinant of the matrix. */
    constexpr T
    determinant(void) const {
        T a00 = this->m00;
        T a01 = this->m01;
        T a02 = this->m02;
        T a03 = this->m03;
        T a10 = this->m10;
        T a11 = this->m11;
        T a12 = this->m12;
        T a13 = this->m13;
        T a20 = this->m20;
        T a21 = this->m21;
        T a22 = this->m22;
        T a23 = this->m23;
        T a30 = this->m30;
        T a31 = this->m31;
        T a32 = this->m32;
        T a33 = this->m33;
        T b00 = a00 * a11 - a01 * a10;
        T b01 = a00 * a12 - a02 * a10;
        T b02 = a00 * a13 - a03 * a10;
        T b03 = a01 * a12 - a02 * a11;
        T b04 = a01 * a13 - a03 * a11;
        T b05 = a02 * a13 - a03 * a12;
        T b06 = a20 * a31 - a21 * a30;
        T b07 = a20 * a32 - a22 * a30;
        T b08 = a20 * a33 - a23 * a30;
        T b09 = a21 * a32 - a22 * a31;
        T b10 = a21 * a33 - a23 * a31;
        T b11 = a22 * a33 - a23 * a32;
        return b00 * b11 - b01 * b10 + b02 * b09 + 
               b03 * b08 - b04 * b07 + b05 * b06;
    }

    /* Invert the matrix. */
    mat4x4<T> 
    inverse(void) const {
        mat4x4 result;
        T a00 = this->m00;
        T a01 = this->m01;
        T a02 = this->m02;
        T a03 = this->m03;
        T a10 = this->m10;
        T a11 = this->m11;
        T a12 = this->m12;
        T a13 = this->m13;
        T a20 = this->m20;
        T a21 = this->m21;
        T a22 = this->m22;
        T a23 = this->m23;
        T a30 = this->m30;
        T a31 = this->m31;
        T a32 = this->m32;
        T a33 = this->m33;
        T b00 = a00 * a11 - a01 * a10;
        T b01 = a00 * a12 - a02 * a10;
        T b02 = a00 * a13 - a03 * a10;
        T b03 = a01 * a12 - a02 * a11;
        T b04 = a01 * a13 - a03 * a11;
        T b05 = a02 * a13 - a03 * a12;
        T b06 = a20 * a31 - a21 * a30;
        T b07 = a20 * a32 - a22 * a30;
        T b08 = a20 * a33 - a23 * a30;
        T b09 = a21 * a32 - a22 * a31;
        T b10 = a21 * a33 - a23 * a31;
        T b11 = a22 * a33 - a23 * a32;
        T inverse_determinant = 
            1.0 / (b00 * b11 - b01 * b10 + b02 * b09 + 
                   b03 * b08 - b04 * b07 + b05 * b06);
        result.m00 = ( a11 * b11 - a12 * b10 + a13 * b09) * inverse_determinant;
        result.m01 = (-a01 * b11 + a02 * b10 - a03 * b09) * inverse_determinant;
        result.m02 = ( a31 * b05 - a32 * b04 + a33 * b03) * inverse_determinant;
        result.m03 = (-a21 * b05 + a22 * b04 - a23 * b03) * inverse_determinant;
        result.m10 = (-a10 * b11 + a12 * b08 - a13 * b07) * inverse_determinant;
        result.m11 = ( a00 * b11 - a02 * b08 + a03 * b07) * inverse_determinant;
        result.m12 = (-a30 * b05 + a32 * b02 - a33 * b01) * inverse_determinant;
        result.m13 = ( a20 * b05 - a22 * b02 + a23 * b01) * inverse_determinant;
        result.m20 = ( a10 * b10 - a11 * b08 + a13 * b06) * inverse_determinant;
        result.m21 = (-a00 * b10 + a01 * b08 - a03 * b06) * inverse_determinant;
        result.m22 = ( a30 * b04 - a31 * b02 + a33 * b00) * inverse_determinant;
        result.m23 = (-a20 * b04 + a21 * b02 - a23 * b00) * inverse_determinant;
        result.m30 = (-a10 * b09 + a11 * b07 - a12 * b06) * inverse_determinant;
        result.m31 = ( a00 * b09 - a01 * b07 + a02 * b06) * inverse_determinant;
        result.m32 = (-a30 * b03 + a31 * b01 - a32 * b00) * inverse_determinant;
        result.m33 = ( a20 * b03 - a21 * b01 + a22 * b00) * inverse_determinant;
        return result;
    }
    
    /* Get x-rotation matrix. Angle must be provided in radians. */
    constexpr mat4x4<T>
    rotate_x(const T& angle) const {
        mat4x4<T> m;
        T c = math::cos(angle);
        T s = math::sin(angle);
        m.m11 =  c;
        m.m12 =  s;
        m.m21 = -s;
        m.m22 =  c;
        return m;
    }

    /* Get y-rotation matrix. Angle must be provided in radians. */
    constexpr mat4x4<T>
    rotate_y(const T& angle) const {
        mat4x4<T> m;
        T c = math::cos(angle);
        T s = math::sin(angle);
        m.m00 =  c;
        m.m02 = -s;
        m.m20 =  s;
        m.m22 =  c;
        return m;
    }

    /* Get z-rotation matrix. Angle must be provided in radians. */
    constexpr mat4x4<T>
    rotate_z(const T& angle) const {
        mat4x4<T> m;
        T c = math::cos(angle);
        T s = math::sin(angle);
        m.m00 =  c;
        m.m01 =  s;
        m.m10 = -s;
        m.m11 =  c;
        return m;
    }

};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

template<typename T>
constexpr sf_inline mat4x4<T>
operator+(const mat4x4<T>& a, const mat4x4<T>& b) {
    return mat4x4<T>(
        a.m00 + b.m00, a.m01 + b.m01, a.m02 + b.m02, a.m03 + b.m03,
        a.m10 + b.m10, a.m11 + b.m11, a.m12 + b.m12, a.m13 + b.m13,
        a.m20 + b.m20, a.m21 + b.m21, a.m22 + b.m22, a.m23 + b.m23,
        a.m30 + b.m30, a.m31 + b.m31, a.m32 + b.m32, a.m33 + b.m33
    );
}

template<typename T>
constexpr sf_inline mat4x4<T>
operator-(const mat4x4<T>& a, const mat4x4<T>& b) {
    return mat4x4<T>(
        a.m00 - b.m00, a.m01 - b.m01, a.m02 - b.m02, a.m03 - b.m03,
        a.m10 - b.m10, a.m11 - b.m11, a.m12 - b.m12, a.m13 - b.m13,
        a.m20 - b.m20, a.m21 - b.m21, a.m22 - b.m22, a.m23 - b.m23,
        a.m30 - b.m30, a.m31 - b.m31, a.m32 - b.m32, a.m33 - b.m33
    );
}

template<typename T>
constexpr sf_inline mat4x4<T>
operator*(const mat4x4<T>& a, const mat4x4<T>& b) {
    mat4x4<T> result;
    for (i32 i = 0; i < 4; ++i) {
        vec4<T> row = a.m[i];
        for (i32 j = 0; j < 4; ++j) {
             vec4<T> col    = vec4<T>(b.m[0][j], 
                                      b.m[1][j], 
                                      b.m[2][j], 
                                      b.m[3][j]);
             result.m[i][j] = dot_product(row, col);
        }
    }
    
    return result;
}

template<typename T>
constexpr sf_inline vec3<T>
operator*(const mat4x4<T>& m, const vec3<T>& v) {
    vec3<T> result;
    result.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
    result.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
    result.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
    return result;
}

template<typename T>
constexpr sf_inline vec4<T> 
operator*(const mat4x4<T>& m, const vec4<T>& v) {
    vec4<T> result;
    result.x = m.m00 * v.x + m.m01 * v.y + 
               m.m02 * v.z + m.m03 * v.w;
    result.y = m.m10 * v.x + m.m11 * v.y + 
               m.m12 * v.z + m.m13 * v.w;
    result.z = m.m20 * v.x + m.m21 * v.y + 
               m.m22 * v.z + m.m23 * v.w;
    result.w = m.m30 * v.x + m.m31 * v.y + 
               m.m32 * v.z + m.m33 * v.w;
    return result;
}

/*------------------*/
/*------------------*/
/* Stream Operators */
/*------------------*/
/*------------------*/

template<typename T>
sf_inline std::ostream&
operator<<(std::ostream& os, const mat4x4<T>& m) {
    os << "mat4x4(" << m.m00 << ", " << m.m01 << ", " 
                    << m.m02 << ", " << m.m03 << ", " 
                    << m.m10 << ", " << m.m11 << ", " 
                    << m.m12 << ", " << m.m13 << ", " 
                    << m.m20 << ", " << m.m21 << ", " 
                    << m.m22 << ", " << m.m23 << ", " 
                    << m.m30 << ", " << m.m31 << ", " 
                    << m.m32 << ", " << m.m33 << ")";
    return os;
}

/*------------------*/
/*------------------*/
/* Mat4x4 Functions */
/*------------------*/
/*------------------*/

template<typename T>
constexpr math::mat3x3<T> to_mat3x3(const math::mat4x4<T>& m) {
    return math::mat3x3<T>(
        m.m00, m.m01, m.m02,
        m.m10, m.m11, m.m12,
        m.m20, m.m21, m.m22
    );
}

template<typename T>
constexpr sf_inline mat4x4<T> 
zero(void) {
    mat4x4<T> m;
    m.m00 = 0;
    m.m01 = 0;
    m.m02 = 0;
    m.m03 = 0;
    m.m10 = 0;
    m.m11 = 0;
    m.m12 = 0;
    m.m13 = 0;
    m.m20 = 0;
    m.m21 = 0;
    m.m22 = 0;
    m.m23 = 0;
    m.m30 = 0;
    m.m31 = 0;
    m.m32 = 0;
    m.m33 = 0;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
identity(void) {
    mat4x4<T> m;
    m.m00 = 1;
    m.m11 = 1;
    m.m22 = 1;
    m.m33 = 1;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
scale(const T& sx, const T& sy, const T& sz) {
    mat4x4<T> m;
    m.m00 = sx;
    m.m11 = sy;
    m.m22 = sz;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
scale(const vec3<T>& v) {
    mat4x4<T> m;
    m.m00 = v.x;
    m.m11 = v.y;
    m.m22 = v.z;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
translate(const T& tx, const T& ty, const T& tz) {
    mat4x4<T> m;
    m.m03 = tx;
    m.m13 = ty;
    m.m23 = tz;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
translate(const vec3<T>& v) {
    mat4x4<T> m;
    m.m03 = v.x;
    m.m13 = v.y;
    m.m23 = v.z;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
transpose(const mat4x4<T>& m) {
    mat4x4<T> result;
    result.m00 = m.m00;
    result.m01 = m.m10;
    result.m02 = m.m20;
    result.m03 = m.m30;
    result.m10 = m.m01;
    result.m11 = m.m11;
    result.m12 = m.m21;
    result.m13 = m.m31;
    result.m20 = m.m02;
    result.m21 = m.m12;
    result.m22 = m.m22;
    result.m23 = m.m32;
    result.m30 = m.m03;
    result.m31 = m.m13;
    result.m32 = m.m23;
    result.m33 = m.m33;
    return result;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
rotate_x(const T& angle) {
    mat4x4<T> m;
    T c = math::cos(angle);
    T s = math::sin(angle);
    m.m11 =  c;
    m.m12 =  s;
    m.m21 = -s;
    m.m22 =  c;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
rotate_y(const T& angle) {
    mat4x4<T> m;
    T c = math::cos(angle);
    T s = math::sin(angle);
    m.m00 =  c;
    m.m02 = -s;
    m.m20 =  s;
    m.m22 =  c;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
rotate_z(const T& angle) {
    mat4x4<T> m;
    T c = math::cos(angle);
    T s = math::sin(angle);
    m.m00 =  c;
    m.m01 =  s;
    m.m10 = -s;
    m.m11 =  c;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
rotate_x(const mat4x4<T>& m, const T& angle) {
    mat4x4<T> rotation;
    T c = math::cos(angle);
    T s = math::sin(angle);
    rotation.m11 =  c;
    rotation.m12 =  s;
    rotation.m21 = -s;
    rotation.m22 =  c;
    return rotation * m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
rotate_y(const mat4x4<T>& m, const T& angle) {
    mat4x4<T> rotation;
    T c = math::cos(angle);
    T s = math::sin(angle);
    rotation.m00 =  c;
    rotation.m02 = -s;
    rotation.m20 =  s;
    rotation.m22 =  c;
    return rotation * m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
rotate_z(const mat4x4<T>& m, const T& angle) {
    mat4x4<T> rotation;
    T c = math::cos(angle);
    T s = math::sin(angle);
    rotation.m00 =  c;
    rotation.m01 =  s;
    rotation.m10 = -s;
    rotation.m11 =  c;
    return rotation * m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
rotate(const mat4x4<T>& m, const T& angle, const vec3<T>& axis) {
    T c = math::cos(angle);
    T s = math::sin(angle);
    T t = 1 - c;
    T x = axis.x;
    T y = axis.y;
    T z = axis.z;
    T tx = t * x;
    T ty = t * y;
    T tz = t * z;
    T sx = s * x;
    T sy = s * y;
    T sz = s * z;
    T txy = tx * y;
    T txz = tx * z;
    T tyz = ty * z;
    mat4x4<T> result;
    result.m00 = tx * x + c;
    result.m01 = txy + sz;
    result.m02 = txz - sy;
    result.m10 = txy - sz;
    result.m11 = ty * y + c;
    result.m12 = tyz + sx;
    result.m20 = txz + sy;
    result.m21 = tyz - sx;
    result.m22 = tz * z + c;
    result.m03 = result.m13 = result.m23 = 0;
    result.m30 = result.m31 = result.m32 = 0;
    result.m33 = 1;
    return result * m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
perspective_projection(const T& fov, 
                       const T& aspect, 
                       const T& z_near, 
                       const T& z_far) {
    mat4x4<T> m;
    m.zero();
    m.m00 = aspect * (1.0 / math::tan(fov / 2.0));
    m.m11 =           1.0 / math::tan(fov / 2.0);
    m.m22 =             z_far / (z_far - z_near);
    m.m23 = (-z_far * z_near) / (z_far - z_near);
    m.m32 = 1.0;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
orthographic_projection(const T& left, 
                        const T& right, 
                        const T& bottom, 
                        const T& top, 
                        const T& z_near, 
                        const T& z_far) {
    mat4x4<T> m;
    m.zero();
    m.m00 = 2.0 / (right - left);
    m.m11 = 2.0 / (top - bottom);
    m.m22 = 1.0 / (z_far - z_near);
    m.m23 = -z_near / (z_far - z_near);
    m.m30 = -(right + left) / (right - left);
    m.m31 = -(top + bottom) / (top - bottom);
    m.m32 = -z_near / (z_far - z_near);
    m.m33 = 1.0;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T> 
look_at(const vec3<T>& eye, const vec3<T>& target, const vec3<T>& up) {
    vec3<T> z = target - eye;
    vec3<T> x = cross_product(up, z.normalize());
    vec3<T> y = cross_product(z,  x.normalize());
    mat4x4<T> view_matrix;
    view_matrix.m00 = x.x;
    view_matrix.m01 = x.y;
    view_matrix.m02 = x.z;
    view_matrix.m03 = -dot_product(x, eye);
    view_matrix.m10 = y.x;
    view_matrix.m11 = y.y;
    view_matrix.m12 = y.z;
    view_matrix.m13 = -dot_product(y, eye);
    view_matrix.m20 = z.x;
    view_matrix.m21 = z.y;
    view_matrix.m22 = z.z;
    view_matrix.m23 = -dot_product(z, eye);
    view_matrix.m30 = 0.0;
    view_matrix.m31 = 0.0;
    view_matrix.m32 = 0.0;
    view_matrix.m33 = 1.0;
    return view_matrix;
}

template<typename T>
constexpr sf_inline mat4x4<T>
from_euler(const T& x, const T& y, const T& z) {
    mat4x4<T> m;
    T cx = math::cos(x);
    T sx = math::sin(x);
    T cy = math::cos(y);
    T sy = math::sin(y);
    T cz = math::cos(z);
    T sz = math::sin(z);
    m.m00 =  cy * cz;
    m.m01 =  cy * sz;
    m.m02 = -sy;
    m.m10 =  sx * sy * cz - cx * sz;
    m.m11 =  sx * sy * sz + cx * cz;
    m.m12 =  sx * cy;
    m.m20 =  cx * sy * cz + sx * sz;
    m.m21 =  cx * sy * sz - sx * cz;
    m.m22 =  cx * cy;
    return m;
}

template<typename T>
constexpr sf_inline mat4x4<T>
from_quaternion(const T& x, const T& y, const T& z, const T& w) {
    mat4x4<T> m;
    T xx = x * x;
    T xy = x * y;
    T xz = x * z;
    T xw = x * w;
    T yy = y * y;
    T yz = y * z;
    T yw = y * w;
    T zz = z * z;
    T zw = z * w;
    m.m00 = 1 - 2 * (yy + zz);
    m.m01 =     2 * (xy - zw);
    m.m02 =     2 * (xz + yw);
    m.m10 =     2 * (xy + zw);
    m.m11 = 1 - 2 * (xx + zz);
    m.m12 =     2 * (yz - xw);
    m.m20 =     2 * (xz - yw);
    m.m21 =     2 * (yz + xw);
    m.m22 = 1 - 2 * (xx + yy);
    return m;
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Quaternion                                                                 */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a quaternion type and a number of functions for manipulating it.   */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

template<typename T>
class quat {
public:

    union {
        struct {
            /* Quaternion components */
            T x, y, z, w;
        };
        /* Array of the quaternion's components. */
        std::array<T, 4> q = { 0, 0, 0, 1 };
        
    };

public:

    /*--------------*/
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /*--------------*/

    /* Default constructor is the identity quaternion. */
    constexpr
    quat() : 
        x(0), y(0), z(0), w(1) {}

    /* Constructor from components. */
    constexpr
    quat(const T& x, const T& y, const T& z, const T& w) : 
        x(x), y(y), z(z), w(w) {}

    /* Copy constructor. */
    constexpr
    quat(const quat& other) : 
        x(other.x), y(other.y), z(other.z), w(other.w) {}

    /* Initialize with matrix. */
    constexpr
    quat(const mat4x4<T>& mat) {
        *this = from_mat4x4(mat);
    }

    constexpr
    quat(T angle, const vec3<T>& axis) {
        T half_angle = angle / 2;
        T s = math::sin(half_angle);
        x = axis.x * s;
        y = axis.y * s;
        z = axis.z * s;
        w = math::cos(half_angle);
    }

    /*---------------------------*/
    /*---------------------------*/
    /* Type Conversion Operators */
    /*---------------------------*/
    /*---------------------------*/

    /* Convert to a 4x4 matrix. */
    constexpr operator mat4x4<T>() const {
        return to_mat4x4();
    }

    template<typename U>
    constexpr operator quat<U>() const {
        return quat<U>(static_cast<U>(x), 
                       static_cast<U>(y), 
                       static_cast<U>(z), 
                       static_cast<U>(w));
    }

    /*------------------*/
    /*------------------*/
    /* Member Operators */
    /*------------------*/
    /*------------------*/

    /* Assignment operator. */
    constexpr quat<T>&
    operator=(const quat<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        w = other.w;
        return *this;
    }

    /* Access the quaternion by index. */
    constexpr T&
    operator[](const int& index) {
        assert(index >= 0 && index < 4);
        return q[index];
    }

    /* Access the quaternion by index (const version). */
    constexpr const T&
    operator[](const int& index) const {
        assert(index >= 0 && index < 4);
        return q[index];
    }

    /*------------------*/
    /*------------------*/
    /* Member Functions */
    /*------------------*/
    /*------------------*/

    /* Identity quaternion. */
    constexpr quat<T>
    identity() {
        return quat<T>(0, 0, 0, 1);
    }

    /* Normalize the quaternion. */
    constexpr void
    normalize() {
        T norm = math::sqrt(x * x + y * y + z * z + w * w);
        x /= norm;
        y /= norm;
        z /= norm;
        w /= norm;
    }

    /* Conjugate of the quaternion. */
    constexpr quat<T>
    conjugate() const {
        return quat<T>(-x, -y, -z, w);
    }

    /* Inverse of the quaternion. */
    constexpr quat<T>
    inverse() const {
        quat<T> conj = conjugate();
        T norm = x * x + y * y + z * z + w * w;
        return quat<T>(conj.x / norm, 
                       conj.y / norm, 
                       conj.z / norm, 
                       conj.w / norm);
    }

    /* Magnitude of the quaternion. */
    constexpr T
    magnitude() const {
        return math::sqrt(x * x + y * y + z * z + w * w);
    }

    /* Convert the quaternion to a 4x4 matrix. */
    constexpr mat4x4<T>
    to_mat4x4() const {
        T xx = x * x;
        T yy = y * y;
        T zz = z * z;
        T ww = w * w;
        T xy = x * y;
        T xz = x * z;
        T xw = x * w;
        T yz = y * z;
        T yw = y * w;
        T zw = z * w;
        return mat4x4<T>(
            1 - 2 * (yy + zz), 2 * (xy - zw), 2 * (xz + yw), 0,
            2 * (xy + zw), 1 - 2 * (xx + zz), 2 * (yz - xw), 0,
            2 * (xz - yw), 2 * (yz + xw), 1 - 2 * (xx + yy), 0,
            0, 0, 0, 1
        );
    }

};

/*----------------------*/
/*----------------------*/
/* Arithmetic Operators */
/*----------------------*/
/*----------------------*/

template<typename T>
constexpr sf_inline quat<T>
operator+(const quat<T>& a, const quat<T>& b) {
    return quat<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

template<typename T>
constexpr sf_inline quat<T>&
operator+=(quat<T>& a, const quat<T>& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
    return a;
}

template<typename T>
constexpr sf_inline quat<T>
operator-(const quat<T>& a, const quat<T>& b) {
    return quat<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

template<typename T>
constexpr sf_inline quat<T>&
operator-=(quat<T>& a, const quat<T>& b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
    return a;
}

template<typename T>
constexpr sf_inline quat<T>
operator*(const quat<T>& a, const quat<T>& b) {
    return quat<T>(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}

template<typename T>
constexpr sf_inline math::vec3<T>
operator*(const math::vec3<T>& v, const quat<T>& q) {
    math::vec3<T> t = 2.0f * cross_product(math::vec3<T>(q.x, q.y, q.z), v);
    return     v + q.w * t + cross_product(math::vec3<T>(q.x, q.y, q.z), t);
}

template<typename T>
constexpr sf_inline quat<T>&
operator*=(quat<T>& a, const quat<T>& b) {
    a = a * b;
    return a;
}

template<typename T>
constexpr sf_inline quat<T>
operator*(const quat<T>& q, const T& scalar) {
    return quat<T>(q.x * scalar, q.y * scalar, q.z * scalar, q.w * scalar);
}

template<typename T>
constexpr sf_inline quat<T>
operator*(const T& scalar, const quat<T>& q) {
    return q * scalar;
}

template<typename T>
constexpr sf_inline quat<T>&
operator*=(quat<T>& q, const T& scalar) {
    q.x *= scalar;
    q.y *= scalar;
    q.z *= scalar;
    q.w *= scalar;
    return q;
}

template<typename T>
constexpr sf_inline quat<T>
operator/(const quat<T>& q, const T& scalar) {
    return quat<T>(q.x / scalar, q.y / scalar, q.z / scalar, q.w / scalar);
}

template<typename T>
constexpr sf_inline quat<T>&
operator/=(quat<T>& q, const T& scalar) {
    q.x /= scalar;
    q.y /= scalar;
    q.z /= scalar;
    q.w /= scalar;
    return q;
}

/*----------------------*/
/*----------------------*/
/* Comparison Operators */
/*----------------------*/
/*----------------------*/

template<typename T>
constexpr sf_inline bool
operator==(const quat<T>& a, const quat<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

template<typename T>
constexpr sf_inline bool
operator!=(const quat<T>& a, const quat<T>& b) {
    return !(a == b);
}

/*------------------*/
/*------------------*/
/* Stream Operators */
/*------------------*/
/*------------------*/

template<typename T>
std::ostream&
operator<<(std::ostream& os, const quat<T>& q) {
    os << "quat(" << q.x << ", " << q.y << ", " << q.z << ", " << q.w << ")";
    return os;
}

/*-------------------*/
/*-------------------*/
/* Utility Functions */
/*-------------------*/
/*-------------------*/

/* Create a quaternion from an axis and an angle. */
template<typename T>
constexpr sf_inline quat<T>
axis_angle(const vec3<T>& axis, const T& angle) {
    T half_angle = angle * 0.5f;
    T sin_half_angle = math::sin(half_angle);
    return quat<T>(axis.x * sin_half_angle,
                   axis.y * sin_half_angle,
                   axis.z * sin_half_angle,
                   math::cos(half_angle));
}

/* Dot product of two quaternions. */
template<typename T>
constexpr sf_inline T
dot_product(const quat<T>& a, const quat<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/* Normalize a quaternion. */
template<typename T>
constexpr sf_inline quat<T>
normalize(const quat<T>& q) {
    T norm = math::sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    return quat<T>(q.x / norm, q.y / norm, q.z / norm, q.w / norm);
}

/* Spherical linear interpolation between two quaternions. */
template<typename T>
constexpr sf_inline quat<T>
slerp(const quat<T>& a, const quat<T>& b, T t) {
    T cos_theta = dot(a, b);
    if (cos_theta < 0) {
        return slerp(-a, b, t);
    }

    if (cos_theta > 0.9995f) {
        return quat<T>(a.x + t * (b.x - a.x),
                       a.y + t * (b.y - a.y),
                       a.z + t * (b.z - a.z),
                       a.w + t * (b.w - a.w)).normalize();
    }

    T theta     = math::acos(cos_theta);
    T sin_theta = math::sin(theta);
    T w1        = math::sin((1 - t) * theta) / sin_theta;
    T w2        = math::sin(t * theta) / sin_theta;
    return quat<T>(w1 * a.x + w2 * b.x,
                   w1 * a.y + w2 * b.y,
                   w1 * a.z + w2 * b.z,
                   w1 * a.w + w2 * b.w);
}

/* Get a quaternion for a given rotation matrix. */
template<typename T>
constexpr sf_inline quat<T>
from_mat4x4(const mat4x4<T>& mat) {
    quat<T> result = { 0.0 };
    T four_w_squared_minus_1 = mat.m00 + mat.m11 + mat.m22;
    T four_x_squared_minus_1 = mat.m00 - mat.m11 - mat.m22;
    T four_y_squared_minus_1 = mat.m11 - mat.m00 - mat.m22;
    T four_z_squared_minus_1 = mat.m22 - mat.m00 - mat.m11;
    i32 largest_index = 0;
    T   largest_squared_minus_1 = four_w_squared_minus_1;

    if (four_x_squared_minus_1 > largest_squared_minus_1) {
        largest_squared_minus_1 = four_x_squared_minus_1;
        largest_index = 1;
    }

    if (four_y_squared_minus_1 > largest_squared_minus_1) {
        largest_squared_minus_1 = four_y_squared_minus_1;
        largest_index = 2;
    }

    if (four_z_squared_minus_1 > largest_squared_minus_1) {
        largest_squared_minus_1 = four_z_squared_minus_1;
        largest_index = 3;
    }

    T largest_value = math::sqrt(largest_squared_minus_1 + 1.0) * 0.5;
    T mult = 0.25 / largest_value;

    switch (largest_index) {
        case 0:
            result.w = largest_value;
            result.x = (mat.m21 - mat.m12) * mult;
            result.y = (mat.m02 - mat.m20) * mult;
            result.z = (mat.m10 - mat.m01) * mult;
            break;
        case 1:
            result.x = largest_value;
            result.w = (mat.m21 - mat.m12) * mult;
            result.y = (mat.m10 + mat.m01) * mult;
            result.z = (mat.m02 + mat.m20) * mult;
            break;
        case 2:
            result.y = largest_value;
            result.w = (mat.m02 - mat.m20) * mult;
            result.x = (mat.m10 + mat.m01) * mult;
            result.z = (mat.m21 + mat.m12) * mult;
            break;
        case 3:
            result.z = largest_value;
            result.w = (mat.m10 - mat.m01) * mult;
            result.x = (mat.m02 + mat.m20) * mult;
            result.y = (mat.m21 + mat.m12) * mult;
            break;
    }

    return result;
}

template<typename T>
constexpr math::vec3<T> quaternion_to_euler(const math::quat<T>& q) {
    math::vec3<T> euler;

    // Roll (x-axis rotation)
    T sinr_cosp = 2 * (q.w * q.x + q.y * q.z);
    T cosr_cosp = 1 - 2 * (q.x * q.x + q.y * q.y);
    euler.x = std::atan2(sinr_cosp, cosr_cosp);

    // Pitch (y-axis rotation)
    T sinp = 2 * (q.w * q.y - q.z * q.x);
    if (std::abs(sinp) >= 1)
        euler.y = std::copysign(math::PI<T> / 2.0, sinp); // use 90 degrees if out of range
    else
        euler.y = std::asin(sinp);

    // Yaw (z-axis rotation)
    T siny_cosp = 2 * (q.w * q.z + q.x * q.y);
    T cosy_cosp = 1 - 2 * (q.y * q.y + q.z * q.z);
    euler.z = std::atan2(siny_cosp, cosy_cosp);

    return euler;
}

template<typename T>
constexpr math::quat<T> euler_to_quaternion(const math::vec3<T>& euler) {
    T cy = math::cos(euler.z * 0.5);
    T sy = math::sin(euler.z * 0.5);
    T cp = math::cos(euler.y * 0.5);
    T sp = math::sin(euler.y * 0.5);
    T cr = math::cos(euler.x * 0.5);
    T sr = math::sin(euler.x * 0.5);

    math::quat<T> q;
    q.w = cr * cp * cy + sr * sp * sy;
    q.x = sr * cp * cy - cr * sp * sy;
    q.y = cr * sp * cy + sr * cp * sy;
    q.z = cr * cp * sy - sr * sp * cy;
    return q;
}

/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Simplex Noise                                                              */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Defines a number of functions for generating simplex noise. The simplex    */
/* noise implementation is based on the 2012 paper: Simplex Noise Demystified */
/* by Stefan Gustavson. The paper can be found in the reference material      */
/* folder.                                                                    */
/*                                                                            */
/*============================================================================*/
/*============================================================================*/

namespace noise {

/* Permutation table */
static u8 PERM[512] = { 
    208, 34,  231, 213, 32,  248, 233, 56,  
    245, 255, 247, 247, 40,  185, 248, 251, 
    245, 28,  124, 204, 204, 76,  36,  1,   
    107, 28,  234, 163, 202, 224, 245, 128, 
    167, 204, 9,   92,  217, 54,  239, 174, 
    173, 102, 193, 189, 190, 121, 100, 108,
    167, 44,  43,  77,  180, 204, 8,   81,  
    70,  223, 11,  38,  24,  254, 210, 210, 
    177, 32,  81,  195, 243, 125, 8,   169, 
    112, 32,  97,  53,  195, 13,  203, 9,   
    47,  104, 125, 117, 114, 124, 165, 203, 
    181, 235, 193, 206, 70,  180, 174, 0,   
    167, 181, 41,  164, 30,  116, 127, 198,
    245, 146, 87,  224, 149, 206, 57,  4,   
    192, 210, 65,  210, 129, 240, 178, 105, 
    228, 108, 245, 148, 140, 40,  35,  195, 
    161, 78,  24,  140, 71,  48,  140, 254, 
    38,  58,  65,  207, 215, 253, 65,  85,  
    208, 76,  62,  3,   237, 55,  89,  232, 
    90,  17,  212, 203, 149, 152, 140, 187,
    234, 177, 73,  174, 193, 100, 192, 143, 
    97,  53,  145, 135, 19,  103, 13,  90,  
    135, 151, 199, 91,  239, 247, 33,  39,  
    145, 101, 120, 99,  3,   186, 86,  99,  
    42,  30,  154, 120, 67,  87,  167, 135, 
    176, 183, 191, 253, 115, 184, 21,  233,
    58,  129, 233, 142, 39,  128, 211, 118, 
    41,  237, 203, 111, 79,  220, 135, 158, 
    137, 139, 255, 114, 20,  218, 113, 154, 
    58,  129, 233, 142, 39,  128, 211, 118, 
    27,  127, 246, 250, 1,   8,   198, 250, 
    209, 92,  222, 173, 21,  88,  102, 219,
    208, 34,  231, 213, 32,  248, 233, 56,  
    161, 78,  24,  140, 71,  48,  140, 254, 
    245, 28,  124, 204, 204, 76,  36,  1,   
    107, 28,  234, 163, 202, 224, 245, 128, 
    167, 204, 9,   92,  217, 54,  239, 174, 
    173, 102, 193, 189, 190, 121, 100, 108,
    167, 44,  43,  77,  180, 204, 8,   81,  
    70,  223, 11,  38,  24,  254, 210, 210, 
    245, 255, 247, 247, 40,  185, 248, 251, 
    177, 32,  81,  195, 243, 125, 8,   169, 
    50,  217, 64,  244, 157, 199, 121, 252, 
    112, 32,  97,  53,  195, 13,  203, 9,   
    47,  104, 125, 117, 114, 124, 165, 203, 
    181, 235, 193, 206, 70,  180, 174, 0,   
    245, 146, 87,  224, 149, 206, 57,  4,   
    192, 210, 65,  210, 129, 240, 178, 105, 
    228, 108, 245, 148, 140, 40,  35,  195, 
    38,  58,  65,  207, 215, 253, 65,  85,  
    208, 76,  62,  3,   237, 55,  89,  232, 
    50,  217, 64,  244, 157, 199, 121, 252, 
    167, 181, 41,  164, 30,  116, 127, 198,
    90,  17,  212, 203, 149, 152, 140, 187,
    234, 177, 73,  174, 193, 100, 192, 143, 
    97,  53,  145, 135, 19,  103, 13,  90,  
    145, 101, 120, 99,  3,   186, 86,  99,  
    41,  237, 203, 111, 79,  220, 135, 158, 
    42,  30,  154, 120, 67,  87,  167, 135, 
    176, 183, 191, 253, 115, 184, 21,  233,
    135, 151, 199, 91,  239, 247, 33,  39,  
    27,  127, 246, 250, 1,   8,   198, 250, 
    137, 139, 255, 114, 20,  218, 113, 154, 
    209, 92,  222, 173, 21,  88,  102, 219 
};

/* Constants */
template<typename T>
constexpr static T F2 = 0.366025403784438646763723170752936184;

template<typename T>
constexpr static T G2 = 0.211324865405187117745425609749021272;

template<typename T>
constexpr static T F3 = 0.333333333333333333333333333333333333;

template<typename T>
constexpr static T G3 = 0.166666666666666666666666666666666667;

/* 3D gradient vectors */
static const f32 GRADIENT_3D[16][3] = {
    {1, 1, 0}, {-1, 1, 0}, {1,-1, 0}, {-1,-1, 0},
    {1, 0, 1}, {-1, 0, 1}, {1, 0,-1}, {-1, 0,-1},
    {0, 1, 1}, { 0,-1, 1}, {0, 1,-1}, { 0,-1,-1},
    {1, 0,-1}, {-1, 0,-1}, {0,-1, 1}, { 0, 1, 1}
};

constexpr f32
simplex_2D(const f32& x, 
           const f32& y,
           const i32& octaves     = 1,
           const f32& persistence = 0.5f,
           const f32& lacunarity  = 2.0f) {
    if (octaves == 1) {
        f32 s = (x + y) * F2<f32>;
        f32 i = math::floor(x + s);
        f32 j = math::floor(y + s);
        f32 t = (i + j) * G2<f32>;
        f32 xy[3][2];
            xy[0][0] = x - (i - t);
            xy[0][1] = y - (j - t);
        i32 i1 = xy[0][0]  > xy[0][1] ? 1 : 0;
        i32 j1 = xy[0][0] <= xy[0][1] ? 1 : 0;
        xy[1][0] = xy[0][0] - i1 + G2<f32>;
        xy[1][1] = xy[0][1] - j1 + G2<f32>;
        xy[2][0] = xy[0][0] + G2<f32> * 2.0f - 1.0f;
        xy[2][1] = xy[0][1] + G2<f32> * 2.0f - 1.0f;

        i32 I = (i32)(i) & 255;
        i32 J = (i32)(j) & 255;
        i32 gi[3];

        gi[0] = PERM[I      + PERM[J]] % 12;
        gi[1] = PERM[I + i1 + PERM[J + j1]] % 12;
        gi[2] = PERM[I + 1  + PERM[J + 1]] % 12;

        f32 noise[3] = {0.0f, 0.0f, 0.0f};
        for (i32 c = 0; c <= 2; c++) {
             f32 f = 0.5 - xy[c][0] * xy[c][0] - 
                           xy[c][1] * xy[c][1];
            if (f > 0) {
                noise[c] = math::quartic(f) * 
                    (GRADIENT_3D[gi[c]][0] * xy[c][0] + 
                     GRADIENT_3D[gi[c]][1] * xy[c][1]);
            }
        }
        return (noise[0] + noise[1] + noise[2]) * 70.0f;
    }
    
    // Multiple octaves
    f32 total     = 0.0f;
    f32 frequency = 1.0f;
    f32 amplitude = 1.0f;
    f32 max_value = 0.0f;
    
    for (i32 i = 0; i < octaves; i++) {
        total += simplex_2D(x * frequency, y * frequency) * amplitude;
        max_value += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    return total / max_value;
}

constexpr f32
simplex_3D(const f32& x, 
           const f32& y,
           const f32& z,
           const i32& octaves     = 1,
           const f32& persistence = 0.5f,
           const f32& lacunarity  = 2.0f) {
    if (octaves == 1) {
        f32 s = (x + y + z) * F3<f32>;
        f32 i = math::floor(x + s);
        f32 j = math::floor(y + s);
        f32 k = math::floor(z + s);
        f32 t = (i + j + k) * G3<f32>;
        f32 xyz[4][3];
            xyz[0][0] = x - (i - t);
            xyz[0][1] = y - (j - t);
            xyz[0][2] = z - (k - t);
        i32 i1, j1, k1;
        i32 i2, j2, k2;
        if (xyz[0][0] >= xyz[0][1]) {
            if (xyz[0][1] >= xyz[0][2]) {
                i1 = 1; j1 = 0; k1 = 0;
                i2 = 1; j2 = 1; k2 = 0;
            } else if (xyz[0][0] >= xyz[0][2]) {
                i1 = 1; j1 = 0; k1 = 0;
                i2 = 1; j2 = 0; k2 = 1;
            } else {
                i1 = 0; j1 = 0; k1 = 1;
                i2 = 1; j2 = 0; k2 = 1;
            }
        } else {
            if (xyz[0][1] < xyz[0][2]) {
                i1 = 0; j1 = 0; k1 = 1;
                i2 = 0; j2 = 1; k2 = 1;
            } else if (xyz[0][0] < xyz[0][2]) {
                i1 = 0; j1 = 1; k1 = 0;
                i2 = 0; j2 = 1; k2 = 1;
            } else {
                i1 = 0; j1 = 1; k1 = 0;
                i2 = 1; j2 = 1; k2 = 0;
            }
        }
        xyz[1][0] = xyz[0][0] - i1 + G3<f32>;
        xyz[1][1] = xyz[0][1] - j1 + G3<f32>;
        xyz[1][2] = xyz[0][2] - k1 + G3<f32>;
        xyz[2][0] = xyz[0][0] - i2  + 2.0f * G3<f32>;
        xyz[2][1] = xyz[0][1] - j2  + 2.0f * G3<f32>;
        xyz[2][2] = xyz[0][2] - k2  + 2.0f * G3<f32>;
        xyz[3][0] = xyz[0][0] - 1.0f + 3.0f * G3<f32>;
        xyz[3][1] = xyz[0][1] - 1.0f + 3.0f * G3<f32>;
        xyz[3][2] = xyz[0][2] - 1.0f + 3.0f * G3<f32>;

        i32 I = (i32)(i) & 255;
        i32 J = (i32)(j) & 255;
        i32 K = (i32)(k) & 255;
        i32 gi[4];
        gi[0] = PERM[I      + PERM[J      + PERM[K]]]      % 12;
        gi[1] = PERM[I + i1 + PERM[J + j1 + PERM[K + k1]]] % 12;
        gi[2] = PERM[I + i2 + PERM[J + j2 + PERM[K + k2]]] % 12;
        gi[3] = PERM[I + 1  + PERM[J + 1  + PERM[K + 1]]]  % 12;
        
        f32 noise[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        for (i32 c = 0; c <= 3; c++) {
             f32 t = 0.6 - math::square(xyz[c][0]) - 
                           math::square(xyz[c][1]) -
                           math::square(xyz[c][2]); 
            if (t > 0) {
                t *= t;
                noise[c] = math::square(t) * (GRADIENT_3D[gi[c]][0] * xyz[c][0] +
                                              GRADIENT_3D[gi[c]][1] * xyz[c][1] +
                                              GRADIENT_3D[gi[c]][2] * xyz[c][2]);
            }
        }
        return (noise[0] + noise[1] + noise[2] + noise[3]) * 32.0f;
    }
    
    // Multiple octaves
    f32 total     = 0.0f;
    f32 frequency = 1.0f;
    f32 amplitude = 1.0f;
    f32 max_value = 0.0f;
    
    for (i32 i = 0; i < octaves; i++) {
        total += simplex_3D(x * frequency, y * frequency, z * frequency) * amplitude;
        max_value += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    return total / max_value;
}

} /* namespace noise */
} /* namespace math */
} /* namespace sf */
