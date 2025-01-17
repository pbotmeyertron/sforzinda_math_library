#pragma once

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* SIMD Class Library                                                         */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/

/* Standard Includes */
#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include <array>

/* Sforzinda Includes */
#include "sf_base.hpp"

namespace sf  {
namespace scl {

template <typename T, std::size_t N>
class simd {
public:

    static_assert(std::is_arithmetic<T>::value, 
                  "simd<T, N> requires arithmetic types");

    static_assert(N > 0, 
                  "simd<T, N> size N must be greater than zero");

public:

    using value_type        = T;

    using mask_element_type = 
    typename std::conditional<std::is_integral<T>::value,
             T, 
             std::int32_t>::type;

    #if defined(__clang__)
        using vector_type = T __attribute__((vector_size(sizeof(T) * N)));
        vector_type data;
    #else
        std::array<T, N> data;
    #endif

public:

    /*------------------------------------------*/
    /* Nested Mask Class for Logical Operations */
    /*------------------------------------------*/

    class mask_type {
    public:

        using element_type = mask_element_type;

        #if defined(__clang__)
            using vector_type = 
            element_type __attribute__((vector_size(sizeof(element_type) * N)));

            vector_type data;
        #else
            std::array<element_type, N> data;
        #endif

    public:

        /*--------------*/
        /* Constructors */
        /*--------------*/

        mask_type() = default;

        #if defined(__clang__)
            mask_type(const vector_type& v) : data(v) {}
        #else
            mask_type(const std::array<element_type, N>& arr) : data(arr) {}
        #endif

        /*----------------------*/
        /* Assignment Operators */
        /*----------------------*/

        #if defined(__clang__)
            constexpr mask_type& 
            operator=(const vector_type& v) {
                data = v;
                return *this;
            }
        #else
            constexpr mask_type& 
            operator=(const std::array<element_type, N>& arr) {
                data = arr;
                return *this;
            }
        #endif

        /*--------------------*/
        /* Indexing Operators */
        /*--------------------*/

        constexpr element_type 
        operator[](std::size_t idx) const {
            #if defined(__clang__)
                return data[idx];
            #else
                return data.at(idx);
            #endif
        }

        explicit 
        operator bool() const {
            for (std::size_t i = 0; i < N; ++i) {
                #if defined(__clang__)
                    if (data[i]) return true;
                #else
                    if (data.at(i)) return true;
                #endif
            }
            return false;
        }

        /*-------------------*/
        /* Logical Operators */
        /*-------------------*/

        mask_type 
        operator~() const {
            #if defined(__clang__)
                return mask_type{~data};
            #else
                std::array<element_type, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = ~data[i];
                }
                return mask_type{result};
            #endif
        }

        mask_type 
        operator&(const mask_type& rhs) const {
            #if defined(__clang__)
                return mask_type{data & rhs.data};
            #else
                std::array<element_type, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = data[i] & rhs.data[i];
                }
                return mask_type{result};
            #endif
        }

        mask_type 
        operator|(const mask_type& rhs) const {
            #if defined(__clang__)
                return mask_type{data | rhs.data};
            #else
                std::array<element_type, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = data[i] | rhs.data[i];
                }
                return mask_type{result};
            #endif
        }

        mask_type 
        operator^(const mask_type& rhs) const {
            #if defined(__clang__)
                return mask_type{data ^ rhs.data};
            #else
                std::array<element_type, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = data[i] ^ rhs.data[i];
                }
                return mask_type{result};
            #endif
        }

        mask_type& 
        operator&=(const mask_type& rhs) {
            #if defined(__clang__)
                data &= rhs.data;
            #else
                for (std::size_t i = 0; i < N; ++i) {
                     data[i] &= rhs.data[i];
                }
            #endif
            return *this;
        }

        mask_type& 
        operator|=(const mask_type& rhs) {
            #if defined(__clang__)
                data |= rhs.data;
            #else
                for (std::size_t i = 0; i < N; ++i) {
                     data[i] |= rhs.data[i];
                }
            #endif
            return *this;
        }

        mask_type& operator^=(const mask_type& rhs) {
            #if defined(__clang__)
                data ^= rhs.data;
            #else
                for (std::size_t i = 0; i < N; ++i) {
                     data[i] ^= rhs.data[i];
                }
            #endif
            return *this;
        }

        friend std::ostream& 
        operator<<(std::ostream& os, const mask_type& mask) {
            os << "{ ";
            for (std::size_t i = 0; i < N; ++i) {
                #if defined(__clang__)
                    os << mask.data[i];
                #else
                   os << mask.data.at(i);
                #endif
                   if (i < N - 1) os << ", ";
            }
            os << " }";
            return os;
        }

    }; // end of mask_type

public:

    /*--------------*/
    /* Constructors */
    /*--------------*/

    simd() = default;

    /* Initialize all elements with a scalar value. */
    explicit simd(T scalar) {
        #if defined(__clang__)
            data = vector_type{};
            data += scalar;
        #else
            data.fill(scalar);
        #endif
    }

    /* Construct with a vector type. */
    #if defined(__clang__)
        simd(const vector_type& vec) : data(vec) {}
    #else
        simd(const std::array<T, N>& arr) : data(arr) {}
    #endif

    /* Construct with an initializer list. */
    explicit 
    simd(std::initializer_list<T> init) {
        #if defined(__clang__)
            data = vector_type{};
        #endif
            std::size_t i = 0;
            for (auto it = init.begin(); it != init.end() && i < N; ++it, ++i) {
                 data[i] = *it;
            }
    }

    /* Construct with an array. */
    template <std::size_t M>
    explicit 
    simd(const std::array<T, M>& arr) {
        #if defined(__clang__)
            data = vector_type{};
        #endif
            std::size_t i = 0;
            for (auto it = arr.begin(); it != arr.end() && i < N; ++it, ++i) {
                 data[i] = *it;
            }
    }

    /*------------------------------*/
    /* Type Conversion Constructors */
    /*------------------------------*/

    explicit 
    operator bool() const {
        for (std::size_t i = 0; i < N; ++i) {
            #if defined(__clang__)
                if (data[i]) return true;
            #else
                if (data.at(i)) return true;
            #endif
        }
        return false;
    }

    explicit
    operator std::array<T, N>() const {
        #if defined(__clang__)
            return data;
        #else
            return data;
        #endif
    }

    /*----------------------*/
    /* Assignment Operators */
    /*----------------------*/

    simd& 
    operator=(const simd& rhs) {
        data = rhs.data;
        return *this;
    }

    simd& 
    operator=(T scalar) {
        #if defined(__clang__)
            data = vector_type{};
            data += scalar;
        #else
            data.fill(scalar);
        #endif
        return *this;
    }

    simd& 
    operator=(std::initializer_list<T> init) {
        #if defined(__clang__)
            data = vector_type{};
        #else
            data.fill(T{});
        #endif
            std::size_t i = 0;
            for (auto it = init.begin(); it != init.end() && i < N; ++it, ++i) {
                 data[i] = *it;
            }
            return *this;
    }

    template<std::size_t M>
    simd& 
    operator=(const std::array<T, M>& arr) {
        #if defined(__clang__)
            data = vector_type{};
        #else
            data.fill(T{});
        #endif
            std::size_t i = 0;
            for (auto it = arr.begin(); it != arr.end() && i < N; ++it, ++i) {
                 data[i] = *it;
            }
            return *this;
    }

    /*--------------------*/
    /* Indexing Operators */
    /*--------------------*/

    T 
    operator[](std::size_t idx) const {
        #if defined(__clang__)
            return data[idx];
        #else
            return data.at(idx);
        #endif
    }

    T& 
    operator[](std::size_t idx) {
        #if defined(__clang__)
            /* Union to access elements of vector type. */
            union {
                vector_type v;
                T           e[N];
            } u;

            static T temp;
                   u.v      = data;
                   temp     = u.e[idx];
                   u.e[idx] = temp;
                   data = u.v;
            return temp;
        #else
            return data.at(idx);
        #endif
    }

    constexpr void 
    load(const T* ptr) {
        #if defined(__clang__)
            std::memcpy(&data, ptr, sizeof(vector_type));
        #else
            std::memcpy(data.data(), ptr, sizeof(T) * N);
        #endif
    }

    constexpr void 
    store(T* ptr) const {
        #if defined(__clang__)
            std::memcpy(ptr, &data, sizeof(vector_type));
        #else
            std::memcpy(ptr, data.data(), sizeof(T) * N);
        #endif
    }

    constexpr void 
    store_reverse(T* ptr) const {
        for (std::size_t i = 0; i < N; ++i) {
            #if defined(__clang__)
                ptr[N - 1 - i] = data[i];
            #else
                ptr[N - 1 - i] = data.at(i);
            #endif
        }
    }

    /*------------------------------------------------------*/
    /* Arithmetic Unary and Arithmetic Assignment Operators */
    /*------------------------------------------------------*/

    /* Negate all elements. */
    simd 
    operator-() const {
        #if defined(__clang__)
            return simd(-data);
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = -data[i];
            }
            return simd(result);
        #endif
    }

    /* Ensure that the vector is positive. */
    simd 
    operator+() const { 
        return *this; 
    }

    simd& 
    operator+=(const T& rhs) {
        #if defined(__clang__)
            data += vector_type { rhs };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] += rhs;
            }
        #endif
        return *this;
    }

    simd& 
    operator-=(const T& rhs) {
        #if defined(__clang__)
            data -= vector_type { rhs };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] -= rhs;
            }
        #endif
        return *this;
    }

    simd& 
    operator*=(const T& rhs) {
        #if defined(__clang__)
            data *= vector_type { rhs };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] *= rhs;
            }
        #endif
        return *this;
    }

    simd& 
    operator/=(const T& rhs) {
        #if defined(__clang__)
            data /= vector_type { rhs };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] /= rhs;
            }
        #endif
        return *this;
    }

    simd& 
    operator+=(const simd& rhs) {
        #if defined(__clang__)
            data += rhs.data;
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] += rhs.data[i];
            }
        #endif
        return *this;
    }

    simd& 
    operator-=(const simd& rhs) {
        #if defined(__clang__)
            data -= rhs.data;
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] -= rhs.data[i];
            }
        #endif
        return *this;
    }

    simd& 
    operator*=(const simd& rhs) {
        #if defined(__clang__)
            data *= rhs.data;
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] *= rhs.data[i];
            }
        #endif
        return *this;
    }

    simd& 
    operator/=(const simd& rhs) {
        #if defined(__clang__)
            data /= rhs.data;
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 data[i] /= rhs.data[i];
            }
        #endif
        return *this;
    }

    simd& 
    operator++() {
        #if defined(__clang__)
            for (std::size_t i = 0; i < N; ++i) {
                data[i] += T(1);
            }
        #else
            for (std::size_t i = 0; i < N; ++i) {
                ++data[i];
            }
        #endif
        return *this;
    }

    simd 
    operator++(std::int32_t) {
        simd temp = *this;
        #if defined(__clang__)
            data += vector_type { T(1) };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 ++data[i];
            }
        #endif
        return temp;
    }

    simd& operator--() {
        #if defined(__clang__)
            for (std::size_t i = 0; i < N; ++i) {
                data[i] -= T(1);
            }
        #else
            for (std::size_t i = 0; i < N; ++i) {
                --data[i];
            }
        #endif
        return *this;
    }

    simd 
    operator--(std::int32_t) {
        simd temp = *this;
        #if defined(__clang__)
            data -= vector_type { T(1) };
        #else
            for (std::size_t i = 0; i < N; ++i) {
                 --data[i];
            }
        #endif
        return temp;
    }

    friend std::ostream& 
    operator<<(std::ostream& os, const simd& v) {
        os << "{ ";
        for (std::size_t i = 0; i < N; ++i) {
            #if defined(__clang__)
                os << v.data[i];
            #else
                os << v.data.at(i);
            #endif
            if (i < N - 1) os << ", ";
        }
        os << " }";
        return os;
    }

    /*----------------------*/
    /* Arithmetic Operators */
    /*----------------------*/

    friend simd 
    operator+(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data + rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] + rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator+(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data + vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] + rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator+(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } + rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs + rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator-(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data - rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] - rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator-(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data - vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] - rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator-(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } - rhs.data};
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs - rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator*(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data * rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] * rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator*(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data * vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] * rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator*(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } * rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs * rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator/(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data / rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] / rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator/(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data / vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] / rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator/(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } / rhs.data };
        #else
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = lhs / rhs.data[i];
                }
                return simd { result };
        #endif
    }

    /*----------------------*/
    /* Comparison Operators */
    /*----------------------*/

    friend mask_type 
    operator==(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data == rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] == rhs.data[i]           ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator==(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data == vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] == rhs                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator==(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } == rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs == rhs.data[i]                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator!=(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data != rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] != rhs.data[i]           ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator!=(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data != vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] != rhs                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator!=(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } != rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs != rhs.data[i]                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }

    friend mask_type 
    operator<(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            typename mask_type::vector_type result;
            for (std::size_t i = 0; i < N; ++i) {
                result[i] = lhs.data[i] < rhs.data[i]           ? 
                           ~typename mask_type::element_type(0) : 
                            typename mask_type::element_type(0);
            }
            return mask_type{result};
        #else
            std::array<typename mask_type::element_type, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                result[i] = lhs.data[i] < rhs.data[i]           ? 
                           ~typename mask_type::element_type(0) : 
                            typename mask_type::element_type(0);
            }
            return mask_type{result};
        #endif
    }
    
    friend mask_type 
    operator<(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data < vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] < rhs                    ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator<(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } < rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs < rhs.data[i]                    ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator>(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            typename mask_type::vector_type result;
            for (std::size_t i = 0; i < N; ++i) {
                result[i] = lhs.data[i] > rhs.data[i]           ? 
                           ~typename mask_type::element_type(0) : 
                            typename mask_type::element_type(0);
            }
            return mask_type{result};
        #else
            std::array<typename mask_type::element_type, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                result[i] = lhs.data[i] > rhs.data[i]           ? 
                           ~typename mask_type::element_type(0) : 
                            typename mask_type::element_type(0);
            }
            return mask_type{result};
        #endif
    }    

    friend mask_type 
    operator>(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data > vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] > rhs                    ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator>(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } > rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs > rhs.data[i]                    ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator<=(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data <= rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] <= rhs.data[i]           ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator<=(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data <= vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] <= rhs                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator<=(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } <= rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs <= rhs.data[i]                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator>=(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data >= rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] >= rhs.data[i]           ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator>=(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return mask_type { lhs.data >= vector_type { rhs } };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] >= rhs                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }
    
    friend mask_type 
    operator>=(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return mask_type { vector_type { lhs } >= rhs.data };
        #else
            std::array<typename mask_type::element_type, N> result{};
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs >= rhs.data[i]                   ? 
                             typename mask_type::element_type(~0) : 
                             typename mask_type::element_type(0);
            }
            return mask_type(result);
        #endif
    }    

    /*-------------------*/
    /* Bitwise Operators */
    /*-------------------*/

    friend simd
    operator<<(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data << rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] << rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator<<(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data << vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                result[i] = lhs.data[i] << rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator<<(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } << rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs << rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator>>(const simd& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { lhs.data >> rhs.data };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] >> rhs.data[i];
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator>>(const simd& lhs, const T& rhs) {
        #if defined(__clang__)
            return simd { lhs.data >> vector_type { rhs } };
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = lhs.data[i] >> rhs;
            }
            return simd { result };
        #endif
    }

    friend simd 
    operator>>(const T& lhs, const simd& rhs) {
        #if defined(__clang__)
            return simd { vector_type { lhs } >> rhs.data };
        #else
                std::array<T, N> result;
                for (std::size_t i = 0; i < N; ++i) {
                     result[i] = lhs >> rhs.data[i];
                }
                return simd { result };
        #endif
    }

    /*-------------------*/
    /* Utility Functions */
    /*-------------------*/

    constexpr T 
    horizontal_sum() const {
        T result = data[0];
        for (std::size_t i = 1; i < N; ++i) {
             result += data[i];
        }
        return result;
    }

    constexpr T 
    horizontal_product() const {
        T result = data[0];
        for (std::size_t i = 1; i < N; ++i) {
             result *= data[i];
        }
        return result;
    }

    /* Dot product of two vectors. */
    constexpr T 
    dot_product(const simd& lhs, const simd& rhs) {
        return (lhs * rhs).horizontal_sum();
    }

    constexpr simd 
    incremental_sequence() {
        #if defined(__clang__)
            simd result;
            for (std::size_t i = 0; i < N; ++i) {
                 result.data[i] = static_cast<T>(i);
            }
            return result;
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = static_cast<T>(i);
            }
            return simd { result };
        #endif
    }

    constexpr simd 
    incremental_sequence_reversed() {
        #if defined(__clang__)
            simd result;
            for (std::size_t i = 0; i < N; ++i) {
                 result.data[i] = static_cast<T>(N - 1 - i);
            }
            return result;
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = static_cast<T>(N - 1 - i);
            }
            return simd { result };
        #endif
    }

    /* Reverse the order of elements in the vector. */
    constexpr simd
    reverse() const {
        #if defined(__clang__)
            simd result;
            for (std::size_t i = 0; i < N; ++i) {
                 result.data[i] = data[N - 1 - i];
            }
            return result;
        #else
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i) {
                 result[i] = data[N - 1 - i];
            }
            return simd { result };
        #endif
    }

    /* Get low part of the vector. */
    template<std::size_t M>
    constexpr simd<T, M>
    get_low() const {
        static_assert(M <= N, "Low size must be less than or equal to simd size");
        simd<T, M> result;
        for (std::size_t i = 0; i < M; ++i) {
             result.data[i] = data[i];
        }
        return result;
    }

    /* Get high part of the vector. */
    template<std::size_t M>
    constexpr simd<T, M>
    get_high() const {
        static_assert(M <= N, "High size must be less than or equal to simd size");
        simd<T, M> result;
        for (std::size_t i = 0; i < M; ++i) {
             result.data[i] = data[N - M + i];
        }
        return result;
    }

};

/*-------------------------------------------------*/
/* Selection, Blending, Permutation, and Swizzling */
/*-------------------------------------------------*/

/* Cut off vector to n elements. The last elements are set to zero. */
template<typename T, std::size_t N>
constexpr sf_inline simd<T, N>
cutoff(const simd<T, N>& vector, std::size_t n) {
    #if defined(__clang__)
        simd<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result.data[i] = (i < n) ? vector.data[i] : T(0);
        }
        return result;
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result[i] = (i < n) ? vector.data[i] : T(0);
        }
        return simd<T, N> { result };
    #endif
}

/* Select between two simd vectors, element by element, based on mask */
template<typename T, std::size_t N>
constexpr sf_inline simd<T, N>
select(const typename simd<T, N>::mask_type& mask, 
       const          simd<T, N>&            a, 
       const          simd<T, N>&            b) {
    #if defined(__clang__)
        simd<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result.data[i] = mask.data[i] ? a.data[i] : b.data[i];
        }
        return result;
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result[i] = mask.data[i] ? a.data[i] : b.data[i];
        }
        return simd<T, N> { result };
    #endif
}

/* Blend two vectors according to immediate constant mask */
template<std::size_t... I, typename T, std::size_t N>
constexpr sf_inline simd<T, N>
blend(const simd<T, N>& a, const simd<T, N>& b) {
    #if defined(__clang__)
        simd<T, N> result;
        constexpr bool mask[N] = { ((I < N) ? true : false)... };
        for (std::size_t i = 0; i < N; ++i) {
             result.data[i] = mask[i] ? a.data[i] : b.data[i];
        }
        return result;
    #else
        std::array<T, N> result;
        constexpr bool mask[N] = { ((I < N) ? true : false)... };
        for (std::size_t i = 0; i < N; ++i) {
             result[i] = mask[i] ? a.data[i] : b.data[i];
        }
        return simd<T, N> { result };
    #endif
}

/* Permute elements in simd vector according to immediate indices */
template<std::size_t... I, typename T, std::size_t N>
constexpr sf_inline simd<T, sizeof...(I)>
permute(const simd<T, N>& vector) {
    static_assert(sizeof...(I) <= N, "Too many indices for simd vector size");
    
    #if defined(__clang__)
        using result_type = simd<T, sizeof...(I)>;
        typename result_type::vector_type result;
        constexpr std::size_t indices[] = { I... };
        for (std::size_t i = 0; i < sizeof...(I); ++i) {
             result[i] = vector.data[indices[i]];
        }
        return result_type { result };
    #else
        std::array<T, sizeof...(I)> result;
        constexpr std::size_t indices[] = { I... };
        for (std::size_t i = 0; i < sizeof...(I); ++i) {
             result[i] = vector.data[indices[i]];
        }
        return simd<T, sizeof...(I)> { result };
    #endif
}

/* Shuffle elements from two vectors according to immediate indices */
template<std::size_t... I, typename T, std::size_t N>
constexpr sf_inline simd<T, sizeof...(I)>
shuffle(const simd<T, N>& a, const simd<T, N>& b) {
    static_assert(sizeof...(I) <= 2*N, 
                  "Too many indices for simd vector size");
    
    #if defined(__clang__)
        using result_type = simd<T, sizeof...(I)>;
        typename result_type::vector_type result;
        constexpr std::size_t indices[] = { I... };
        for (std::size_t i = 0; i < sizeof...(I); ++i) {
            if (indices[i] < N) {
                result[i] = a.data[indices[i]];
            } else {
                result[i] = b.data[indices[i] - N];
            }
        }
        return result_type{result};
    #else
        std::array<T, sizeof...(I)> result;
        constexpr std::size_t indices[] = { I... };
        for (std::size_t i = 0; i < sizeof...(I); ++i) {
            if (indices[i] < N) {
                result[i] = a.data[indices[i]];
            } else {
                result[i] = b.data[indices[i] - N];
            }
        }
        return simd<T, sizeof...(I)> { result };
    #endif
}

/* Swap two simd vectors */
template<typename T, std::size_t N>
constexpr sf_inline void
swap(simd<T, N>& a, simd<T, N>& b) {
    simd<T, N> temp = a;
               a    = b;
               b    = temp;
}

/* Sign combine two simd vectors */
template<typename T, std::size_t N>
constexpr sf_inline simd<T, N>
sign_combine(const simd<T, N>& a, const simd<T, N>& b) {
    #if defined(__clang__)
        return simd<T, N> { a.data ^ 
                           (b.data & (T(1) << (sizeof(T) * 8 - 1))) };
    #else
        std::array<T, N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result[i] = a.data[i] ^ 
                        (b.data[i] & (T(1) << (sizeof(T) * 8 - 1)));
        }
        return simd<T, N> { result };
    #endif
}

/* Split simd vector into two simd vectors. */
template<typename T, std::size_t N>
constexpr sf_inline std::pair<simd<T, N/2>, simd<T, N/2>>
split(const simd<T, N>& vector) {
    static_assert(N % 2 == 0, "Split size must be even");
    #if defined(__clang__)
        simd<T, N/2> a;
        simd<T, N/2> b;
        for (std::size_t i = 0; i < N/2; ++i) {
             a.data[i] = vector.data[i];
             b.data[i] = vector.data[N/2 + i];
        }
        return std::make_pair(a, b);
    #else
        std::array<T, N/2> a;
        std::array<T, N/2> b;
        for (std::size_t i = 0; i < N/2; ++i) {
             a[i] = vector.data[i];
             b[i] = vector.data[N/2 + i];
        }
        return std::make_pair(simd<T, N/2> { a }, simd<T, N/2> { b });
    #endif
}

/* Merge two simd vectors into one simd vector. */
template<typename T, std::size_t N>
constexpr sf_inline simd<T, 2*N>
merge(const simd<T, N>& a, const simd<T, N>& b) {
    #if defined(__clang__)
        simd<T, 2*N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result.data[i] = a.data[i];
             result.data[N + i] = b.data[i];
        }
        return result;
    #else
        std::array<T, 2*N> result;
        for (std::size_t i = 0; i < N; ++i) {
             result[i] = a.data[i];
             result[N + i] = b.data[i];
        }
        return simd<T, 2*N> { result };
    #endif
}

/*-----------------------------*/
/* Logical and State Functions */
/*-----------------------------*/

/* Returns true if all bits are 1. */
template<typename T, std::size_t N>
constexpr sf_inline bool
horizontal_and(const typename simd<T, N>::mask_type& mask) {
    bool result = true;
    for (std::size_t i = 0; i < N; ++i) {
         result &= mask.data[i];
    }
    return result;
}

/* Returns true if any bit is 1. */
template<typename T, std::size_t N>
constexpr sf_inline bool
horizontal_or(const typename simd<T, N>::mask_type& mask) {
    bool result = false;
    for (std::size_t i = 0; i < N; ++i) {
         result |= mask.data[i];
    }
    return result;
}

/* Returns true if all bits are 0. */
template<typename T, std::size_t N>
constexpr sf_inline bool
horizontal_not(const typename simd<T, N>::mask_type& mask) {
    bool result = true;
    for (std::size_t i = 0; i < N; ++i) {
         result &= !mask.data[i];
    }
    return result;
}

/* Convert boolean vector to integer bitfield. */
template<typename T, std::size_t N>
constexpr sf_inline std::size_t
to_bitfield(const typename simd<T, N>::mask_type& mask) {
    std::size_t result = 0;
    for (std::size_t i = 0; i < N; ++i) {
         result |= (mask.data[i] ? 1 : 0) << i;
    }
    return result;
}

/* Convert integer bitfield to boolean vector. */
template<typename T, std::size_t N>
constexpr sf_inline typename simd<T, N>::mask_type
to_mask(std::size_t bitfield) {
    typename simd<T, N>::mask_type result;
    for (std::size_t i = 0; i < N; ++i) {
         result.data[i] = (bitfield & (1 << i)) != 0;
    }
    return result;
}

} /* namespace scl */
} /* namespace sf */
