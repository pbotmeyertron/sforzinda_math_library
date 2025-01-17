#pragma once

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Sforzinda Game Engine Base Configuration                                   */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/

/* Standard Includes */
#include <cstdint>

/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Compiler Helpers - Compiler-Specific Attributes and Macros Used Globally   */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* The engine currently supports the following compilers:                     */
/*     - Clang                                                                */
/*     - GCC                                                                  */
/*     - MSVC                                                                 */
/*                                                                            */
/*============================================================================*/

/* Forces alignment of a variable or struct member to n bytes */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_align(n)         __attribute__((aligned(n)))
#elif defined(_MSC_VER)
   #define sf_align(n)         __declspec(align(n))
#else
   #define sf_align(n)
#endif

/* Minimizes structure padding to reduce memory usage */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_packed           __attribute__((packed))
#elif defined(_MSC_VER)
   #define sf_packed           __pragma(pack(1))
#else
   #define sf_packed
#endif

/* Standard inline function declaration */
#define sf_inline               static inline

/* Forces function to be inlined, overriding compiler's decision */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_always_inline    static __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
   #define sf_always_inline    static __forceinline
#else
   #define sf_always_inline    static inline
#endif

/* Prevents function from being inlined */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_no_inline        __attribute__((noinline))
#elif defined(_MSC_VER)
   #define sf_no_inline        __declspec(noinline)
#else
   #define sf_no_inline
#endif

/* Indicates function never returns to caller */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_noreturn         __attribute__((noreturn))
#elif defined(_MSC_VER)
   #define sf_noreturn         __declspec(noreturn)
#else
   #define sf_noreturn
#endif

/* Controls symbol visibility in shared libraries */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_public           __attribute__((visibility("default")))
   #define sf_private          __attribute__((visibility("hidden")))
#elif defined(_MSC_VER)
   #define sf_public           __declspec(dllexport)
   #define sf_private
#else
   #define sf_public
   #define sf_private
#endif

/* Hints to compiler about likely/unlikely code paths for optimization */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_expect_yes(x)    __builtin_expect(!!(x), 1)
   #define sf_expect_no(x)     __builtin_expect((x), 0)
   #define sf_likely(x)        __builtin_expect(!!(x), 1)
   #define sf_unlikely(x)      __builtin_expect((x), 0)
#else
   #define sf_expect_yes(x)    (x)
   #define sf_expect_no(x)     (x)
   #define sf_likely(x)        (x)
   #define sf_unlikely(x)      (x)
#endif

/* Indicates code path is unreachable, allowing compiler optimizations */
#if defined(__clang__) || defined(__GNUC__)
   #define sf_unreachable()    __builtin_unreachable()
#elif defined(_MSC_VER)
   #define sf_unreachable()    __assume(0)
#else
   #define sf_unreachable()
#endif

/* Suppresses unused parameter warnings */
#define sf_unused_parameter(x)  (void)(x)

/* Safely computes offset of a member within a struct */
#if defined(__clang__) || defined(__GNUC__) || defined(_MSC_VER)
   #define sf_offset_of(s,m)   __builtin_offsetof(s,m)
#else
   #define sf_offset_of(s,m)   ((size_t)&(((s*)0)->m))
#endif

/*============================================================================*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/* Fundamental Types                                                          */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Includes the fundamental types used by the Sforzinda library. All          */
/* primitive types are typedef'd to their respective standard C99 types.      */
/*                                                                            */
/*============================================================================*/

namespace sf {

/* Signed 8-Bit Integer */
using i8  = std::int8_t;
/* Signed 16-Bit Integer */
using i16 = std::int16_t;
/* Signed 32-Bit Integer */
using i32 = std::int32_t;
/* Signed 64-Bit Integer */
using i64 = std::int64_t;

/* Unsigned 8-Bit Integer */
using u8  = std::uint8_t;
/* Unsigned 16-Bit Integer */
using u16 = std::uint16_t;
/* Unsigned 32-Bit Integer */
using u32 = std::uint32_t;
/* Unsigned 64-Bit Integer */
using u64 = std::uint64_t;

/* 32-Bit Floating-Point Number */
using f32  = float;
/* 64-Bit Floating-Point Number */
using f64  = double;

} /* namespace sf */
