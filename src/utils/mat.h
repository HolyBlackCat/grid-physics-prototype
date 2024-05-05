// mat.h
// Vector and matrix math
// Version 3.25
// Generated, don't touch.

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <iterator>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>


#ifndef IMP_MATH_IS_CONSTANT
#  ifndef _MSC_VER
#    define IMP_MATH_IS_CONSTANT(...) __builtin_constant_p(__VA_ARGS__)
#  else
#    define IMP_MATH_IS_CONSTANT(...) false
#  endif
#endif

#ifndef IMP_MATH_UNREACHABLE
#  ifndef _MSC_VER
#    define IMP_MATH_UNREACHABLE(...) __builtin_unreachable()
#  else
#    define IMP_MATH_UNREACHABLE(...) __assume(false)
#  endif
#endif

#ifndef IMP_MATH_SMALL_FUNC
#  ifndef _MSC_VER
#    define IMP_MATH_SMALL_FUNC   __attribute__((__always_inline__, __artificial__)) inline // Need explicit inline, otherwise `artificial` complains, even on implicitly inline functions.
#    define IMP_MATH_SMALL_LAMBDA __attribute__((__always_inline__, __artificial__))
#  else
#    define IMP_MATH_SMALL_FUNC   [[msvc::forceinline]]
#    define IMP_MATH_SMALL_LAMBDA [[msvc::forceinline]] // There is also `__forceinline`, but it doesn't work on lambdas.
#  endif
#endif

// Vectors and matrices

namespace Math
{
    inline namespace Utility // Scalar concepts
    {
        template <typename T> concept cvref_unqualified = std::is_same_v<T, std::remove_cvref_t<T>>;

        // Whether a type is a scalar.
        template <typename T> struct helper_is_scalar : std::is_arithmetic<T> {}; // Not `std::is_scalar`, because that includes pointers.
        template <typename T> concept scalar = cvref_unqualified<T> && helper_is_scalar<T>::value;
        template <typename T> concept scalar_maybe_const = scalar<std::remove_const_t<T>>;
    }

    inline namespace Vector // Declarations
    {
        template <int D, scalar T> struct vec;
        template <int D, scalar T> struct rect;
        template <int W, int H, scalar T> struct mat;
    }

    inline namespace Alias // Short type aliases
    {
        template <scalar T> using vec2 = vec<2,T>; template <scalar T> using vec3 = vec<3,T>; template <scalar T> using vec4 = vec<4,T>;
        template <scalar T> using rect2 = rect<2,T>; template <scalar T> using rect3 = rect<3,T>; template <scalar T> using rect4 = rect<4,T>;
        template <scalar T> using mat2x2 = mat<2,2,T>; template <scalar T> using mat3x2 = mat<3,2,T>; template <scalar T> using mat4x2 = mat<4,2,T>;
        template <scalar T> using mat2x3 = mat<2,3,T>; template <scalar T> using mat3x3 = mat<3,3,T>; template <scalar T> using mat4x3 = mat<4,3,T>;
        template <scalar T> using mat2x4 = mat<2,4,T>; template <scalar T> using mat3x4 = mat<3,4,T>; template <scalar T> using mat4x4 = mat<4,4,T>;
        template <scalar T> using mat2 = mat2x2<T>; template <scalar T> using mat3 = mat3x3<T>; template <scalar T> using mat4 = mat4x4<T>;

        template <int D> using bvec = vec<D,bool>;
        template <int D> using brect = rect<D,bool>;
        template <int W, int H> using bmat = mat<W,H,bool>;
        using bvec2 = vec<2,bool>; using bvec3 = vec<3,bool>; using bvec4 = vec<4,bool>;
        using brect2 = rect<2,bool>; using brect3 = rect<3,bool>; using brect4 = rect<4,bool>;
        using bmat2x2 = mat<2,2,bool>; using bmat3x2 = mat<3,2,bool>; using bmat4x2 = mat<4,2,bool>;
        using bmat2x3 = mat<2,3,bool>; using bmat3x3 = mat<3,3,bool>; using bmat4x3 = mat<4,3,bool>;
        using bmat2x4 = mat<2,4,bool>; using bmat3x4 = mat<3,4,bool>; using bmat4x4 = mat<4,4,bool>;
        using bmat2 = bmat2x2; using bmat3 = bmat3x3; using bmat4 = bmat4x4;

        template <int D> using cvec = vec<D,char>;
        template <int D> using crect = rect<D,char>;
        template <int W, int H> using cmat = mat<W,H,char>;
        using cvec2 = vec<2,char>; using cvec3 = vec<3,char>; using cvec4 = vec<4,char>;
        using crect2 = rect<2,char>; using crect3 = rect<3,char>; using crect4 = rect<4,char>;
        using cmat2x2 = mat<2,2,char>; using cmat3x2 = mat<3,2,char>; using cmat4x2 = mat<4,2,char>;
        using cmat2x3 = mat<2,3,char>; using cmat3x3 = mat<3,3,char>; using cmat4x3 = mat<4,3,char>;
        using cmat2x4 = mat<2,4,char>; using cmat3x4 = mat<3,4,char>; using cmat4x4 = mat<4,4,char>;
        using cmat2 = cmat2x2; using cmat3 = cmat3x3; using cmat4 = cmat4x4;

        template <int D> using ucvec = vec<D,unsigned char>;
        template <int D> using ucrect = rect<D,unsigned char>;
        template <int W, int H> using ucmat = mat<W,H,unsigned char>;
        using ucvec2 = vec<2,unsigned char>; using ucvec3 = vec<3,unsigned char>; using ucvec4 = vec<4,unsigned char>;
        using ucrect2 = rect<2,unsigned char>; using ucrect3 = rect<3,unsigned char>; using ucrect4 = rect<4,unsigned char>;
        using ucmat2x2 = mat<2,2,unsigned char>; using ucmat3x2 = mat<3,2,unsigned char>; using ucmat4x2 = mat<4,2,unsigned char>;
        using ucmat2x3 = mat<2,3,unsigned char>; using ucmat3x3 = mat<3,3,unsigned char>; using ucmat4x3 = mat<4,3,unsigned char>;
        using ucmat2x4 = mat<2,4,unsigned char>; using ucmat3x4 = mat<3,4,unsigned char>; using ucmat4x4 = mat<4,4,unsigned char>;
        using ucmat2 = ucmat2x2; using ucmat3 = ucmat3x3; using ucmat4 = ucmat4x4;

        template <int D> using scvec = vec<D,signed char>;
        template <int D> using screct = rect<D,signed char>;
        template <int W, int H> using scmat = mat<W,H,signed char>;
        using scvec2 = vec<2,signed char>; using scvec3 = vec<3,signed char>; using scvec4 = vec<4,signed char>;
        using screct2 = rect<2,signed char>; using screct3 = rect<3,signed char>; using screct4 = rect<4,signed char>;
        using scmat2x2 = mat<2,2,signed char>; using scmat3x2 = mat<3,2,signed char>; using scmat4x2 = mat<4,2,signed char>;
        using scmat2x3 = mat<2,3,signed char>; using scmat3x3 = mat<3,3,signed char>; using scmat4x3 = mat<4,3,signed char>;
        using scmat2x4 = mat<2,4,signed char>; using scmat3x4 = mat<3,4,signed char>; using scmat4x4 = mat<4,4,signed char>;
        using scmat2 = scmat2x2; using scmat3 = scmat3x3; using scmat4 = scmat4x4;

        template <int D> using svec = vec<D,short>;
        template <int D> using srect = rect<D,short>;
        template <int W, int H> using smat = mat<W,H,short>;
        using svec2 = vec<2,short>; using svec3 = vec<3,short>; using svec4 = vec<4,short>;
        using srect2 = rect<2,short>; using srect3 = rect<3,short>; using srect4 = rect<4,short>;
        using smat2x2 = mat<2,2,short>; using smat3x2 = mat<3,2,short>; using smat4x2 = mat<4,2,short>;
        using smat2x3 = mat<2,3,short>; using smat3x3 = mat<3,3,short>; using smat4x3 = mat<4,3,short>;
        using smat2x4 = mat<2,4,short>; using smat3x4 = mat<3,4,short>; using smat4x4 = mat<4,4,short>;
        using smat2 = smat2x2; using smat3 = smat3x3; using smat4 = smat4x4;

        template <int D> using usvec = vec<D,unsigned short>;
        template <int D> using usrect = rect<D,unsigned short>;
        template <int W, int H> using usmat = mat<W,H,unsigned short>;
        using usvec2 = vec<2,unsigned short>; using usvec3 = vec<3,unsigned short>; using usvec4 = vec<4,unsigned short>;
        using usrect2 = rect<2,unsigned short>; using usrect3 = rect<3,unsigned short>; using usrect4 = rect<4,unsigned short>;
        using usmat2x2 = mat<2,2,unsigned short>; using usmat3x2 = mat<3,2,unsigned short>; using usmat4x2 = mat<4,2,unsigned short>;
        using usmat2x3 = mat<2,3,unsigned short>; using usmat3x3 = mat<3,3,unsigned short>; using usmat4x3 = mat<4,3,unsigned short>;
        using usmat2x4 = mat<2,4,unsigned short>; using usmat3x4 = mat<3,4,unsigned short>; using usmat4x4 = mat<4,4,unsigned short>;
        using usmat2 = usmat2x2; using usmat3 = usmat3x3; using usmat4 = usmat4x4;

        template <int D> using ivec = vec<D,int>;
        template <int D> using irect = rect<D,int>;
        template <int W, int H> using imat = mat<W,H,int>;
        using ivec2 = vec<2,int>; using ivec3 = vec<3,int>; using ivec4 = vec<4,int>;
        using irect2 = rect<2,int>; using irect3 = rect<3,int>; using irect4 = rect<4,int>;
        using imat2x2 = mat<2,2,int>; using imat3x2 = mat<3,2,int>; using imat4x2 = mat<4,2,int>;
        using imat2x3 = mat<2,3,int>; using imat3x3 = mat<3,3,int>; using imat4x3 = mat<4,3,int>;
        using imat2x4 = mat<2,4,int>; using imat3x4 = mat<3,4,int>; using imat4x4 = mat<4,4,int>;
        using imat2 = imat2x2; using imat3 = imat3x3; using imat4 = imat4x4;

        template <int D> using uvec = vec<D,unsigned int>;
        template <int D> using urect = rect<D,unsigned int>;
        template <int W, int H> using umat = mat<W,H,unsigned int>;
        using uvec2 = vec<2,unsigned int>; using uvec3 = vec<3,unsigned int>; using uvec4 = vec<4,unsigned int>;
        using urect2 = rect<2,unsigned int>; using urect3 = rect<3,unsigned int>; using urect4 = rect<4,unsigned int>;
        using umat2x2 = mat<2,2,unsigned int>; using umat3x2 = mat<3,2,unsigned int>; using umat4x2 = mat<4,2,unsigned int>;
        using umat2x3 = mat<2,3,unsigned int>; using umat3x3 = mat<3,3,unsigned int>; using umat4x3 = mat<4,3,unsigned int>;
        using umat2x4 = mat<2,4,unsigned int>; using umat3x4 = mat<3,4,unsigned int>; using umat4x4 = mat<4,4,unsigned int>;
        using umat2 = umat2x2; using umat3 = umat3x3; using umat4 = umat4x4;

        template <int D> using lvec = vec<D,long>;
        template <int D> using lrect = rect<D,long>;
        template <int W, int H> using lmat = mat<W,H,long>;
        using lvec2 = vec<2,long>; using lvec3 = vec<3,long>; using lvec4 = vec<4,long>;
        using lrect2 = rect<2,long>; using lrect3 = rect<3,long>; using lrect4 = rect<4,long>;
        using lmat2x2 = mat<2,2,long>; using lmat3x2 = mat<3,2,long>; using lmat4x2 = mat<4,2,long>;
        using lmat2x3 = mat<2,3,long>; using lmat3x3 = mat<3,3,long>; using lmat4x3 = mat<4,3,long>;
        using lmat2x4 = mat<2,4,long>; using lmat3x4 = mat<3,4,long>; using lmat4x4 = mat<4,4,long>;
        using lmat2 = lmat2x2; using lmat3 = lmat3x3; using lmat4 = lmat4x4;

        template <int D> using ulvec = vec<D,unsigned long>;
        template <int D> using ulrect = rect<D,unsigned long>;
        template <int W, int H> using ulmat = mat<W,H,unsigned long>;
        using ulvec2 = vec<2,unsigned long>; using ulvec3 = vec<3,unsigned long>; using ulvec4 = vec<4,unsigned long>;
        using ulrect2 = rect<2,unsigned long>; using ulrect3 = rect<3,unsigned long>; using ulrect4 = rect<4,unsigned long>;
        using ulmat2x2 = mat<2,2,unsigned long>; using ulmat3x2 = mat<3,2,unsigned long>; using ulmat4x2 = mat<4,2,unsigned long>;
        using ulmat2x3 = mat<2,3,unsigned long>; using ulmat3x3 = mat<3,3,unsigned long>; using ulmat4x3 = mat<4,3,unsigned long>;
        using ulmat2x4 = mat<2,4,unsigned long>; using ulmat3x4 = mat<3,4,unsigned long>; using ulmat4x4 = mat<4,4,unsigned long>;
        using ulmat2 = ulmat2x2; using ulmat3 = ulmat3x3; using ulmat4 = ulmat4x4;

        template <int D> using llvec = vec<D,long long>;
        template <int D> using llrect = rect<D,long long>;
        template <int W, int H> using llmat = mat<W,H,long long>;
        using llvec2 = vec<2,long long>; using llvec3 = vec<3,long long>; using llvec4 = vec<4,long long>;
        using llrect2 = rect<2,long long>; using llrect3 = rect<3,long long>; using llrect4 = rect<4,long long>;
        using llmat2x2 = mat<2,2,long long>; using llmat3x2 = mat<3,2,long long>; using llmat4x2 = mat<4,2,long long>;
        using llmat2x3 = mat<2,3,long long>; using llmat3x3 = mat<3,3,long long>; using llmat4x3 = mat<4,3,long long>;
        using llmat2x4 = mat<2,4,long long>; using llmat3x4 = mat<3,4,long long>; using llmat4x4 = mat<4,4,long long>;
        using llmat2 = llmat2x2; using llmat3 = llmat3x3; using llmat4 = llmat4x4;

        template <int D> using ullvec = vec<D,unsigned long long>;
        template <int D> using ullrect = rect<D,unsigned long long>;
        template <int W, int H> using ullmat = mat<W,H,unsigned long long>;
        using ullvec2 = vec<2,unsigned long long>; using ullvec3 = vec<3,unsigned long long>; using ullvec4 = vec<4,unsigned long long>;
        using ullrect2 = rect<2,unsigned long long>; using ullrect3 = rect<3,unsigned long long>; using ullrect4 = rect<4,unsigned long long>;
        using ullmat2x2 = mat<2,2,unsigned long long>; using ullmat3x2 = mat<3,2,unsigned long long>; using ullmat4x2 = mat<4,2,unsigned long long>;
        using ullmat2x3 = mat<2,3,unsigned long long>; using ullmat3x3 = mat<3,3,unsigned long long>; using ullmat4x3 = mat<4,3,unsigned long long>;
        using ullmat2x4 = mat<2,4,unsigned long long>; using ullmat3x4 = mat<3,4,unsigned long long>; using ullmat4x4 = mat<4,4,unsigned long long>;
        using ullmat2 = ullmat2x2; using ullmat3 = ullmat3x3; using ullmat4 = ullmat4x4;

        template <int D> using fvec = vec<D,float>;
        template <int D> using frect = rect<D,float>;
        template <int W, int H> using fmat = mat<W,H,float>;
        using fvec2 = vec<2,float>; using fvec3 = vec<3,float>; using fvec4 = vec<4,float>;
        using frect2 = rect<2,float>; using frect3 = rect<3,float>; using frect4 = rect<4,float>;
        using fmat2x2 = mat<2,2,float>; using fmat3x2 = mat<3,2,float>; using fmat4x2 = mat<4,2,float>;
        using fmat2x3 = mat<2,3,float>; using fmat3x3 = mat<3,3,float>; using fmat4x3 = mat<4,3,float>;
        using fmat2x4 = mat<2,4,float>; using fmat3x4 = mat<3,4,float>; using fmat4x4 = mat<4,4,float>;
        using fmat2 = fmat2x2; using fmat3 = fmat3x3; using fmat4 = fmat4x4;

        template <int D> using dvec = vec<D,double>;
        template <int D> using drect = rect<D,double>;
        template <int W, int H> using dmat = mat<W,H,double>;
        using dvec2 = vec<2,double>; using dvec3 = vec<3,double>; using dvec4 = vec<4,double>;
        using drect2 = rect<2,double>; using drect3 = rect<3,double>; using drect4 = rect<4,double>;
        using dmat2x2 = mat<2,2,double>; using dmat3x2 = mat<3,2,double>; using dmat4x2 = mat<4,2,double>;
        using dmat2x3 = mat<2,3,double>; using dmat3x3 = mat<3,3,double>; using dmat4x3 = mat<4,3,double>;
        using dmat2x4 = mat<2,4,double>; using dmat3x4 = mat<3,4,double>; using dmat4x4 = mat<4,4,double>;
        using dmat2 = dmat2x2; using dmat3 = dmat3x3; using dmat4 = dmat4x4;

        template <int D> using ldvec = vec<D,long double>;
        template <int D> using ldrect = rect<D,long double>;
        template <int W, int H> using ldmat = mat<W,H,long double>;
        using ldvec2 = vec<2,long double>; using ldvec3 = vec<3,long double>; using ldvec4 = vec<4,long double>;
        using ldrect2 = rect<2,long double>; using ldrect3 = rect<3,long double>; using ldrect4 = rect<4,long double>;
        using ldmat2x2 = mat<2,2,long double>; using ldmat3x2 = mat<3,2,long double>; using ldmat4x2 = mat<4,2,long double>;
        using ldmat2x3 = mat<2,3,long double>; using ldmat3x3 = mat<3,3,long double>; using ldmat4x3 = mat<4,3,long double>;
        using ldmat2x4 = mat<2,4,long double>; using ldmat3x4 = mat<3,4,long double>; using ldmat4x4 = mat<4,4,long double>;
        using ldmat2 = ldmat2x2; using ldmat3 = ldmat3x3; using ldmat4 = ldmat4x4;

        template <int D> using i8vec = vec<D,std::int8_t>;
        template <int D> using i8rect = rect<D,std::int8_t>;
        template <int W, int H> using i8mat = mat<W,H,std::int8_t>;
        using i8vec2 = vec<2,std::int8_t>; using i8vec3 = vec<3,std::int8_t>; using i8vec4 = vec<4,std::int8_t>;
        using i8rect2 = rect<2,std::int8_t>; using i8rect3 = rect<3,std::int8_t>; using i8rect4 = rect<4,std::int8_t>;
        using i8mat2x2 = mat<2,2,std::int8_t>; using i8mat3x2 = mat<3,2,std::int8_t>; using i8mat4x2 = mat<4,2,std::int8_t>;
        using i8mat2x3 = mat<2,3,std::int8_t>; using i8mat3x3 = mat<3,3,std::int8_t>; using i8mat4x3 = mat<4,3,std::int8_t>;
        using i8mat2x4 = mat<2,4,std::int8_t>; using i8mat3x4 = mat<3,4,std::int8_t>; using i8mat4x4 = mat<4,4,std::int8_t>;
        using i8mat2 = i8mat2x2; using i8mat3 = i8mat3x3; using i8mat4 = i8mat4x4;

        template <int D> using u8vec = vec<D,std::uint8_t>;
        template <int D> using u8rect = rect<D,std::uint8_t>;
        template <int W, int H> using u8mat = mat<W,H,std::uint8_t>;
        using u8vec2 = vec<2,std::uint8_t>; using u8vec3 = vec<3,std::uint8_t>; using u8vec4 = vec<4,std::uint8_t>;
        using u8rect2 = rect<2,std::uint8_t>; using u8rect3 = rect<3,std::uint8_t>; using u8rect4 = rect<4,std::uint8_t>;
        using u8mat2x2 = mat<2,2,std::uint8_t>; using u8mat3x2 = mat<3,2,std::uint8_t>; using u8mat4x2 = mat<4,2,std::uint8_t>;
        using u8mat2x3 = mat<2,3,std::uint8_t>; using u8mat3x3 = mat<3,3,std::uint8_t>; using u8mat4x3 = mat<4,3,std::uint8_t>;
        using u8mat2x4 = mat<2,4,std::uint8_t>; using u8mat3x4 = mat<3,4,std::uint8_t>; using u8mat4x4 = mat<4,4,std::uint8_t>;
        using u8mat2 = u8mat2x2; using u8mat3 = u8mat3x3; using u8mat4 = u8mat4x4;

        template <int D> using i16vec = vec<D,std::int16_t>;
        template <int D> using i16rect = rect<D,std::int16_t>;
        template <int W, int H> using i16mat = mat<W,H,std::int16_t>;
        using i16vec2 = vec<2,std::int16_t>; using i16vec3 = vec<3,std::int16_t>; using i16vec4 = vec<4,std::int16_t>;
        using i16rect2 = rect<2,std::int16_t>; using i16rect3 = rect<3,std::int16_t>; using i16rect4 = rect<4,std::int16_t>;
        using i16mat2x2 = mat<2,2,std::int16_t>; using i16mat3x2 = mat<3,2,std::int16_t>; using i16mat4x2 = mat<4,2,std::int16_t>;
        using i16mat2x3 = mat<2,3,std::int16_t>; using i16mat3x3 = mat<3,3,std::int16_t>; using i16mat4x3 = mat<4,3,std::int16_t>;
        using i16mat2x4 = mat<2,4,std::int16_t>; using i16mat3x4 = mat<3,4,std::int16_t>; using i16mat4x4 = mat<4,4,std::int16_t>;
        using i16mat2 = i16mat2x2; using i16mat3 = i16mat3x3; using i16mat4 = i16mat4x4;

        template <int D> using u16vec = vec<D,std::uint16_t>;
        template <int D> using u16rect = rect<D,std::uint16_t>;
        template <int W, int H> using u16mat = mat<W,H,std::uint16_t>;
        using u16vec2 = vec<2,std::uint16_t>; using u16vec3 = vec<3,std::uint16_t>; using u16vec4 = vec<4,std::uint16_t>;
        using u16rect2 = rect<2,std::uint16_t>; using u16rect3 = rect<3,std::uint16_t>; using u16rect4 = rect<4,std::uint16_t>;
        using u16mat2x2 = mat<2,2,std::uint16_t>; using u16mat3x2 = mat<3,2,std::uint16_t>; using u16mat4x2 = mat<4,2,std::uint16_t>;
        using u16mat2x3 = mat<2,3,std::uint16_t>; using u16mat3x3 = mat<3,3,std::uint16_t>; using u16mat4x3 = mat<4,3,std::uint16_t>;
        using u16mat2x4 = mat<2,4,std::uint16_t>; using u16mat3x4 = mat<3,4,std::uint16_t>; using u16mat4x4 = mat<4,4,std::uint16_t>;
        using u16mat2 = u16mat2x2; using u16mat3 = u16mat3x3; using u16mat4 = u16mat4x4;

        template <int D> using i32vec = vec<D,std::int32_t>;
        template <int D> using i32rect = rect<D,std::int32_t>;
        template <int W, int H> using i32mat = mat<W,H,std::int32_t>;
        using i32vec2 = vec<2,std::int32_t>; using i32vec3 = vec<3,std::int32_t>; using i32vec4 = vec<4,std::int32_t>;
        using i32rect2 = rect<2,std::int32_t>; using i32rect3 = rect<3,std::int32_t>; using i32rect4 = rect<4,std::int32_t>;
        using i32mat2x2 = mat<2,2,std::int32_t>; using i32mat3x2 = mat<3,2,std::int32_t>; using i32mat4x2 = mat<4,2,std::int32_t>;
        using i32mat2x3 = mat<2,3,std::int32_t>; using i32mat3x3 = mat<3,3,std::int32_t>; using i32mat4x3 = mat<4,3,std::int32_t>;
        using i32mat2x4 = mat<2,4,std::int32_t>; using i32mat3x4 = mat<3,4,std::int32_t>; using i32mat4x4 = mat<4,4,std::int32_t>;
        using i32mat2 = i32mat2x2; using i32mat3 = i32mat3x3; using i32mat4 = i32mat4x4;

        template <int D> using u32vec = vec<D,std::uint32_t>;
        template <int D> using u32rect = rect<D,std::uint32_t>;
        template <int W, int H> using u32mat = mat<W,H,std::uint32_t>;
        using u32vec2 = vec<2,std::uint32_t>; using u32vec3 = vec<3,std::uint32_t>; using u32vec4 = vec<4,std::uint32_t>;
        using u32rect2 = rect<2,std::uint32_t>; using u32rect3 = rect<3,std::uint32_t>; using u32rect4 = rect<4,std::uint32_t>;
        using u32mat2x2 = mat<2,2,std::uint32_t>; using u32mat3x2 = mat<3,2,std::uint32_t>; using u32mat4x2 = mat<4,2,std::uint32_t>;
        using u32mat2x3 = mat<2,3,std::uint32_t>; using u32mat3x3 = mat<3,3,std::uint32_t>; using u32mat4x3 = mat<4,3,std::uint32_t>;
        using u32mat2x4 = mat<2,4,std::uint32_t>; using u32mat3x4 = mat<3,4,std::uint32_t>; using u32mat4x4 = mat<4,4,std::uint32_t>;
        using u32mat2 = u32mat2x2; using u32mat3 = u32mat3x3; using u32mat4 = u32mat4x4;

        template <int D> using i64vec = vec<D,std::int64_t>;
        template <int D> using i64rect = rect<D,std::int64_t>;
        template <int W, int H> using i64mat = mat<W,H,std::int64_t>;
        using i64vec2 = vec<2,std::int64_t>; using i64vec3 = vec<3,std::int64_t>; using i64vec4 = vec<4,std::int64_t>;
        using i64rect2 = rect<2,std::int64_t>; using i64rect3 = rect<3,std::int64_t>; using i64rect4 = rect<4,std::int64_t>;
        using i64mat2x2 = mat<2,2,std::int64_t>; using i64mat3x2 = mat<3,2,std::int64_t>; using i64mat4x2 = mat<4,2,std::int64_t>;
        using i64mat2x3 = mat<2,3,std::int64_t>; using i64mat3x3 = mat<3,3,std::int64_t>; using i64mat4x3 = mat<4,3,std::int64_t>;
        using i64mat2x4 = mat<2,4,std::int64_t>; using i64mat3x4 = mat<3,4,std::int64_t>; using i64mat4x4 = mat<4,4,std::int64_t>;
        using i64mat2 = i64mat2x2; using i64mat3 = i64mat3x3; using i64mat4 = i64mat4x4;

        template <int D> using u64vec = vec<D,std::uint64_t>;
        template <int D> using u64rect = rect<D,std::uint64_t>;
        template <int W, int H> using u64mat = mat<W,H,std::uint64_t>;
        using u64vec2 = vec<2,std::uint64_t>; using u64vec3 = vec<3,std::uint64_t>; using u64vec4 = vec<4,std::uint64_t>;
        using u64rect2 = rect<2,std::uint64_t>; using u64rect3 = rect<3,std::uint64_t>; using u64rect4 = rect<4,std::uint64_t>;
        using u64mat2x2 = mat<2,2,std::uint64_t>; using u64mat3x2 = mat<3,2,std::uint64_t>; using u64mat4x2 = mat<4,2,std::uint64_t>;
        using u64mat2x3 = mat<2,3,std::uint64_t>; using u64mat3x3 = mat<3,3,std::uint64_t>; using u64mat4x3 = mat<4,3,std::uint64_t>;
        using u64mat2x4 = mat<2,4,std::uint64_t>; using u64mat3x4 = mat<3,4,std::uint64_t>; using u64mat4x4 = mat<4,4,std::uint64_t>;
        using u64mat2 = u64mat2x2; using u64mat3 = u64mat3x3; using u64mat4 = u64mat4x4;

        template <int D> using xvec = vec<D,std::ptrdiff_t>;
        template <int D> using xrect = rect<D,std::ptrdiff_t>;
        template <int W, int H> using xmat = mat<W,H,std::ptrdiff_t>;
        using xvec2 = vec<2,std::ptrdiff_t>; using xvec3 = vec<3,std::ptrdiff_t>; using xvec4 = vec<4,std::ptrdiff_t>;
        using xrect2 = rect<2,std::ptrdiff_t>; using xrect3 = rect<3,std::ptrdiff_t>; using xrect4 = rect<4,std::ptrdiff_t>;
        using xmat2x2 = mat<2,2,std::ptrdiff_t>; using xmat3x2 = mat<3,2,std::ptrdiff_t>; using xmat4x2 = mat<4,2,std::ptrdiff_t>;
        using xmat2x3 = mat<2,3,std::ptrdiff_t>; using xmat3x3 = mat<3,3,std::ptrdiff_t>; using xmat4x3 = mat<4,3,std::ptrdiff_t>;
        using xmat2x4 = mat<2,4,std::ptrdiff_t>; using xmat3x4 = mat<3,4,std::ptrdiff_t>; using xmat4x4 = mat<4,4,std::ptrdiff_t>;
        using xmat2 = xmat2x2; using xmat3 = xmat3x3; using xmat4 = xmat4x4;

        template <int D> using zvec = vec<D,std::size_t>;
        template <int D> using zrect = rect<D,std::size_t>;
        template <int W, int H> using zmat = mat<W,H,std::size_t>;
        using zvec2 = vec<2,std::size_t>; using zvec3 = vec<3,std::size_t>; using zvec4 = vec<4,std::size_t>;
        using zrect2 = rect<2,std::size_t>; using zrect3 = rect<3,std::size_t>; using zrect4 = rect<4,std::size_t>;
        using zmat2x2 = mat<2,2,std::size_t>; using zmat3x2 = mat<3,2,std::size_t>; using zmat4x2 = mat<4,2,std::size_t>;
        using zmat2x3 = mat<2,3,std::size_t>; using zmat3x3 = mat<3,3,std::size_t>; using zmat4x3 = mat<4,3,std::size_t>;
        using zmat2x4 = mat<2,4,std::size_t>; using zmat3x4 = mat<3,4,std::size_t>; using zmat4x4 = mat<4,4,std::size_t>;
        using zmat2 = zmat2x2; using zmat3 = zmat3x3; using zmat4 = zmat4x4;
    }

    namespace Custom // Customization points
    {
        // Specializing this adds corresponding constructors and conversion operators to vectors and matrices.
        template <typename From, typename To>
        struct Convert
        {
            // To operator()(const From &) const {...}
            // static constexpr bool is_explicit = false;
        };

        template <typename From, typename To>
        concept convertible = requires(const Convert<From, To> conv, const From from)
        {
            { conv(from) } -> std::same_as<To>;
        };
    }

    inline namespace Utility // Helper templates
    {
        // Some of the concept definitions here are redundant.
        // In some cases this sanitizes user specializations. In some cases it should help with subsumption.

        // Check if `T` is a vector type.
        template <typename T> struct helper_is_vector : std::false_type {};
        template <int D, typename T> struct helper_is_vector<vec<D,T>> : std::true_type {};
        template <typename T> concept vector = cvref_unqualified<T>/*redundant*/ && helper_is_vector<T>::value;
        template <typename T> concept vector_maybe_const = vector<std::remove_const_t<T>>;

        template <typename T> concept vector_or_scalar = scalar<T> || vector<T>;
        template <typename T> concept vector_or_scalar_maybe_const = scalar_maybe_const<T> || vector_maybe_const<T>;

        // Checks if any of `P...` are vector types.
        template <typename ...P> inline constexpr bool any_vectors_v = (vector<P> || ...);

        // Check if `T` is a matrix type.
        template <typename T> struct helper_is_matrix : std::false_type {};
        template <int W, int H, typename T> struct helper_is_matrix<mat<W,H,T>> : std::true_type {};
        template <typename T> concept matrix = cvref_unqualified<T>/*redundant*/ && helper_is_matrix<T>::value;
        template <typename T> concept square_matrix = matrix<T> && T::width == T::height;
        template <typename T> concept matrix_maybe_const = matrix<std::remove_const_t<T>>;
        template <typename T> concept square_matrix_maybe_const = square_matrix<std::remove_const_t<T>>;

        // For vectors returns their element type, for scalars returns them unchanged.
        template <typename T> struct helper_vec_base {using type = T;};
        template <int D, typename T> struct helper_vec_base<      vec<D,T>> {using type =       T;};
        template <int D, typename T> struct helper_vec_base<const vec<D,T>> {using type = const T;};
        template <vector_or_scalar_maybe_const T> using vec_base_t = typename helper_vec_base<T>::type;
        // This version accepts any type, and returns unknown types unchanged.
        template <typename T> using vec_base_weak_t = typename helper_vec_base<T>::type;

        // Whether `T` is a vector with the base type `U`.
        template <typename T, typename U> concept vector_with_base = vector<T> && std::same_as<U, vec_base_t<T>>;
        // Whether `T` is a vector or scalar with the base type `U`.
        template <typename T, typename U> concept vector_or_scalar_with_base = vector_or_scalar<T> && std::same_as<U, vec_base_t<T>>;

        // For vectors returns the number of elements, for scalars returns 1.
        template <typename T> struct helper_vec_size : std::integral_constant<int, 1> {};
        template <int D, typename T> struct helper_vec_size<      vec<D,T>> : std::integral_constant<int, D> {};
        template <int D, typename T> struct helper_vec_size<const vec<D,T>> : std::integral_constant<int, D> {};
        template <vector_or_scalar_maybe_const T> inline constexpr int vec_size_v = helper_vec_size<T>::value;
        template <typename T> inline constexpr int vec_size_weak_v = helper_vec_size<T>::value;

        // If `D == 1` or `T == void`, returns `T`. Otherwise returns `vec<D,T>`.
        template <int D, typename T> struct helper_ver_or_scalar {using type = vec<D,T>;};
        template <int D, typename T> requires(D == 1 || std::is_void_v<T>) struct helper_ver_or_scalar<D, T> {using type = T;};
        template <int D, typename T> using vec_or_scalar_t = typename helper_ver_or_scalar<D,T>::type;

        // If the set {D...} is either {N} or {1,N}, returns `N`.
        // If the set {D...} is empty, returns `1`.
        // Otherwise returns 0.
        template <int ...D> inline constexpr int common_vec_size_or_zero_v = []{
            int ret = 1;
            bool ok = ((D == 1 ? true : ret == 1 || ret == D ? (void(ret = D), true) : false) && ...);
            return ok * ret;
        }();

        template <int ...D> concept have_common_vec_size = common_vec_size_or_zero_v<D...> != 0;

        // If the set {D...} is either {N} or {1,N}, returns `N`.
        // If the set {D...} is empty, returns `1`.
        // Otherwise causes a soft error.
        template <int ...D> requires have_common_vec_size<D...>
        inline constexpr int common_vec_size_v = common_vec_size_or_zero_v<D...>;

        // If `A` is a vector, changes its element type to `B`. If `A` is scalar, returns `B`.
        // In any case, preserves constness of `A`.
        template <typename A, typename B> struct helper_change_vec_base {using type = B;};
        template <typename A, typename B> struct helper_change_vec_base<const A,B> {using type = const typename helper_change_vec_base<A, B>::type;};
        template <int D, typename A, typename B> struct helper_change_vec_base<vec<D,A>,B> {using type = vec<D,B>;};
        template <vector_or_scalar_maybe_const A, scalar B> using change_vec_base_t = typename helper_change_vec_base<A,B>::type;
        // This version accepts any types, and treats them as scalars.
        template <typename A, typename B> using change_vec_base_weak_t = typename helper_change_vec_base<A,B>::type;

        // Whether `T` is a floating-point type, or a vector of such.
        template <typename T> struct helper_is_floating_point_scalar : std::is_floating_point<T> {};
        template <typename T> concept floating_point_scalar = scalar<T> && helper_is_floating_point_scalar<T>::value;
        template <typename T> concept floating_point_vector = vector<T>/*reject const types*/ && floating_point_scalar<vec_base_t<T>>;
        template <typename T> concept floating_point_vector_or_scalar = floating_point_scalar<T> || floating_point_vector<T>;

        // Whether `T` is an integral type, or a vector of such.
        template <typename T> struct helper_is_integral_scalar : std::is_integral<T> {};
        template <typename T> concept integral_scalar = scalar<T> && helper_is_integral_scalar<T>::value;
        template <typename T> concept integral_vector = vector<T> && integral_scalar<vec_base_t<T>>;
        template <typename T> concept integral_vector_or_scalar = integral_scalar<T> || integral_vector<T>;

        // Whether `T` is a signed/unsigned integral type, or a vector of such.
        template <typename T> struct helper_is_unsigned_integral_scalar : std::is_unsigned<T> {};
        template <typename T> concept   signed_integral_scalar = integral_scalar<T> && !helper_is_unsigned_integral_scalar<T>::value;
        template <typename T> concept unsigned_integral_scalar = integral_scalar<T> &&  helper_is_unsigned_integral_scalar<T>::value;
        template <typename T> concept   signed_integral_vector = integral_vector<T> &&   signed_integral_scalar<vec_base_t<T>>;
        template <typename T> concept unsigned_integral_vector = integral_vector<T> && unsigned_integral_scalar<vec_base_t<T>>;
        template <typename T> concept   signed_integral_vector_or_scalar =   signed_integral_scalar<T> ||   signed_integral_vector<T>;
        template <typename T> concept unsigned_integral_vector_or_scalar = unsigned_integral_scalar<T> || unsigned_integral_vector<T>;

        template <typename T> concept signed_maybe_floating_point_scalar = signed_integral_scalar<T> || floating_point_scalar<T>;
        template <typename T> concept signed_maybe_floating_point_vector = signed_integral_vector<T> || floating_point_vector<T>;
        template <typename T> concept signed_maybe_floating_point_vector_or_scalar = signed_integral_vector_or_scalar<T> || floating_point_vector_or_scalar<T>;

        // Returns a reasonable 'floating-point counterpart' for a type.
        // Currently if the type is not floating-point, returns `float`. Otherwise returns the same type.
        // If `T` is a vector, it's base type is changed according to the same rules.
        template <vector_or_scalar T> using floating_point_t = std::conditional_t<floating_point_vector_or_scalar<T>, T, change_vec_base_t<T, float>>;

        // 3-way compares two scalar or vector types to determine which one is 'larger'.
        // Considers the types equivalent only if they are the same.
        template <cvref_unqualified A, cvref_unqualified B> inline constexpr std::partial_ordering compare_types_v = []{
            if constexpr (std::is_same_v<A, B>)
                return std::partial_ordering::equivalent;
            else if constexpr (!vector_or_scalar<A> || !vector_or_scalar<B>)
                return std::partial_ordering::unordered;
            else if constexpr (vec_size_v<A> != vec_size_v<B>)
                return std::partial_ordering::unordered;
            else if constexpr (floating_point_vector_or_scalar<A> < floating_point_vector_or_scalar<B>)
                return std::partial_ordering::less;
            else if constexpr (floating_point_vector_or_scalar<A> > floating_point_vector_or_scalar<B>)
                return std::partial_ordering::greater;
            else if constexpr (signed_integral_vector_or_scalar<A> != signed_integral_vector_or_scalar<B>)
                return std::partial_ordering::unordered;
            else if constexpr (sizeof(vec_base_t<A>) < sizeof(vec_base_t<B>))
                return std::partial_ordering::less;
            else if constexpr (sizeof(vec_base_t<A>) > sizeof(vec_base_t<B>))
                return std::partial_ordering::greater;
            else
                return std::partial_ordering::unordered;
        }();

        // Internal, see below for the public interface.
        // Given a list of scalar and vector types, determines the "larger' type among them according to `compare_types_v`.
        // Returns `void` on failure.
        // If vector types are present, all of them must have the same size, and the resulting type will also be a vector.
        template <typename ...P> struct helper_larger {};
        template <typename T> struct helper_larger<T> {using type = T;};
        template <typename A, typename B, typename C, typename ...P> requires requires{typename helper_larger<B,C,P...>::type;} struct helper_larger<A,B,C,P...> {using type = typename helper_larger<A, typename helper_larger<B,C,P...>::type>::type;};
        template <typename A, typename B> requires(compare_types_v<A,B> == std::partial_ordering::equivalent) struct helper_larger<A,B> {using type = A;};
        template <typename A, typename B> requires(compare_types_v<A,B> == std::partial_ordering::less      ) struct helper_larger<A,B> {using type = B;};
        template <typename A, typename B> requires(compare_types_v<A,B> == std::partial_ordering::greater   ) struct helper_larger<A,B> {using type = A;};
        // Causes a soft error if there's no larger type.
        template <cvref_unqualified ...P> using larger_t = vec_or_scalar_t<common_vec_size_v<vec_size_weak_v<P>...>, typename helper_larger<std::remove_cv_t<vec_base_weak_t<P>>...>::type>;

        // Checks if it's possible to determine the 'larger' type among `P`.
        template <typename ...P> concept have_larger_type = requires{typename larger_t<P...>;};

        // Whether the conversion of `A` to `B` is not narrowing. Doesn't fail when there's no conversion, should return false in that case.
        template <typename A, typename B> concept safely_convertible_to = std::is_same_v<larger_t<A, B>, B>;

        struct uninit {}; // A constructor tag to leave a vector/matrix uninitialized.

        // Wrappers for different kinds of comparisons.
        template <vector_or_scalar T> struct compare_any {const T &value; [[nodiscard]] explicit constexpr compare_any(const T &value) : value(value) {}};
        template <vector_or_scalar T> struct compare_all {const T &value; [[nodiscard]] explicit constexpr compare_all(const T &value) : value(value) {}};
        template <vector_or_scalar T> struct compare_none {const T &value; [[nodiscard]] explicit constexpr compare_none(const T &value) : value(value) {}};
        template <vector_or_scalar T> struct compare_not_all {const T &value; [[nodiscard]] explicit constexpr compare_not_all(const T &value) : value(value) {}};
        template <vector_or_scalar T> struct compare_elemwise {const T &value; [[nodiscard]] explicit constexpr compare_elemwise(const T &value) : value(value) {}};
        // Tags for different kinds of comparisons.
        struct compare_any_tag {template <vector_or_scalar T> [[nodiscard]] constexpr compare_any<T> operator()(const T &value) const {return compare_any(value);}};
        struct compare_all_tag {template <vector_or_scalar T> [[nodiscard]] constexpr compare_all<T> operator()(const T &value) const {return compare_all(value);}};
        struct compare_none_tag {template <vector_or_scalar T> [[nodiscard]] constexpr compare_none<T> operator()(const T &value) const {return compare_none(value);}};
        struct compare_not_all_tag {template <vector_or_scalar T> [[nodiscard]] constexpr compare_not_all<T> operator()(const T &value) const {return compare_not_all(value);}};
        struct compare_elemwise_tag {template <vector_or_scalar T> [[nodiscard]] constexpr compare_elemwise<T> operator()(const T &value) const {return compare_elemwise(value);}};
    }

    inline namespace Utility // Helpers for operators
    {
        // Returns i-th vector element. For other types ignores the index.
        template <typename T>
        [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr decltype(auto) vec_elem(int i, T &&vec)
        {
            if constexpr (std::is_lvalue_reference_v<T>)
            {
                if constexpr (!vector<std::remove_cvref_t<T>>)
                    return vec;
                else
                    return vec[i];
            }
            else
            {
                if constexpr (!vector<std::remove_cvref_t<T>>)
                    return std::move(vec);
                else
                    return std::move(vec[i]);
            }
        }

        // Helper for applying a function to one or several scalars or vectors.
        // Mixing scalars and vectors is allowed, but vectors must have the same size.
        // If at least one vector is passed, the result is also a vector.
        // If `D != 1`, forces the result to be the vector of this size, or causes a hard error if not possible.
        template <int D = 1, typename F, typename ...P, typename = std::enable_if_t<(vector_or_scalar_maybe_const<std::remove_reference_t<P>> && ...)>> // Trying to put this condition into `requires` crashes Clang 14.
        IMP_MATH_SMALL_FUNC constexpr auto apply_elementwise(F &&func, P &&... params) -> vec_or_scalar_t<common_vec_size_v<D, vec_size_v<std::remove_reference_t<P>>...>, decltype(std::declval<F>()(vec_elem(0, std::declval<P>())...))>
        {
            constexpr int size = common_vec_size_v<D, vec_size_v<std::remove_reference_t<P>>...>;
            using R = vec_or_scalar_t<size, decltype(std::declval<F>()(vec_elem(0, std::declval<P>())...))>;

            if constexpr (std::is_void_v<R>)
            {
                for (int i = 0; i < size; i++)
                    func(vec_elem(i, params)...); // No forwarding to prevent moving.
                return void();
            }
            else
            {
                R ret{};
                for (int i = 0; i < size; i++)
                    vec_elem(i, ret) = func(vec_elem(i, params)...); // No forwarding to prevent moving.
                return ret;
            }
        }

        template <vector_or_scalar T> [[nodiscard]] constexpr bool any_nonzero_elements(const T &value)
        {
            if constexpr (vector<T>)
                return value.any();
            else
                return bool(value);
        }
        template <vector_or_scalar T> [[nodiscard]] constexpr bool all_nonzero_elements(const T &value)
        {
            if constexpr (vector<T>)
                return value.all();
            else
                return bool(value);
        }
        template <vector_or_scalar T> [[nodiscard]] constexpr bool none_nonzero_elements(const T &value)
        {
            if constexpr (vector<T>)
                return value.none();
            else
                return !bool(value);
        }
        template <vector_or_scalar T> [[nodiscard]] constexpr bool not_all_nonzero_elements(const T &value)
        {
            if constexpr (vector<T>)
                return value.not_all();
            else
                return !bool(value);
        }
    }

    inline namespace Vector // Operators
    {
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator+(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() + std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a + b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator-(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() - std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a - b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator*(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() * std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a * b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator/(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() / std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a / b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator%(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() % std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a % b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator^(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() ^ std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a ^ b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator&(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() & std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a & b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator|(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() | std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a | b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator<<(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() << std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a << b;}, a, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator>>(const A &a, const B &b) -> vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, decltype(std::declval<vec_base_t<A>>() >> std::declval<vec_base_t<B>>())> {return apply_elementwise([](vec_base_t<A> a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {return a >> b;}, a, b);}
        template <vector V> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator~(const V &v) -> change_vec_base_t<V, decltype(~v.x)> {return apply_elementwise([](vec_base_t<V> v) IMP_MATH_SMALL_LAMBDA {return ~v;}, v);}
        template <vector V> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator+(const V &v) -> change_vec_base_t<V, decltype(+v.x)> {return apply_elementwise([](vec_base_t<V> v) IMP_MATH_SMALL_LAMBDA {return +v;}, v);}
        template <vector V> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator-(const V &v) -> change_vec_base_t<V, decltype(-v.x)> {return apply_elementwise([](vec_base_t<V> v) IMP_MATH_SMALL_LAMBDA {return -v;}, v);}
        template <vector_or_scalar_with_base<bool> V> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr auto operator!(const V &v) -> change_vec_base_t<V, decltype(!v.x)> {return apply_elementwise([](vec_base_t<V> v) IMP_MATH_SMALL_LAMBDA {return !v;}, v);}
        template <vector V> IMP_MATH_SMALL_FUNC constexpr V &operator++(V &v) {apply_elementwise([](vec_base_t<V> &v) IMP_MATH_SMALL_LAMBDA {++v;}, v); return v;}
        template <vector V> IMP_MATH_SMALL_FUNC constexpr V operator++(V &v, int) {V ret = v; apply_elementwise([](vec_base_t<V> &v) IMP_MATH_SMALL_LAMBDA {++v;}, v); return ret;}
        template <vector V> IMP_MATH_SMALL_FUNC constexpr V &operator--(V &v) {apply_elementwise([](vec_base_t<V> &v) IMP_MATH_SMALL_LAMBDA {--v;}, v); return v;}
        template <vector V> IMP_MATH_SMALL_FUNC constexpr V operator--(V &v, int) {V ret = v; apply_elementwise([](vec_base_t<V> &v) IMP_MATH_SMALL_LAMBDA {--v;}, v); return ret;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator+=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() += std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a += b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator-=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() -= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a -= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator*=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() *= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a *= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator/=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() /= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a /= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator%=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() %= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a %= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator^=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() ^= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a ^= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator&=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() &= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a &= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator|=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() |= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a |= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator<<=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() <<= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a <<= b;}, a, b); return a;}
        template <vector A, safely_convertible_to<A> B> IMP_MATH_SMALL_FUNC constexpr auto operator>>=(A &a, const B &b) -> decltype(std::enable_if_t<vector<A> && vector_or_scalar<B>>(), void(std::declval<vec_base_t<A> &>() >>= std::declval<vec_base_t<B>>()), std::declval<A &>()) {apply_elementwise([](vec_base_t<A> &a, vec_base_t<B> b) IMP_MATH_SMALL_LAMBDA {a >>= b;}, a, b); return a;}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) < b; else return a < compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::less{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::less{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::less{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::less{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::less{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::less{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::less{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::less{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::less{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::less{}, a, b.value);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) > b; else return a > compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::greater{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::greater{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::greater{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::greater{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::greater{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::greater{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::greater{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::greater{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::greater{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::greater{}, a, b.value);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<=(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) <= b; else return a <= compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::less_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::less_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::less_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::less_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::less_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::less_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::less_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator<=(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::less_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<=(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::less_equal{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator<=(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::less_equal{}, a, b.value);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>=(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) >= b; else return a >= compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::greater_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::greater_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::greater_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::greater_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::greater_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::greater_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::greater_equal{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator>=(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::greater_equal{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>=(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::greater_equal{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator>=(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::greater_equal{}, a, b.value);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(const A &a, const B &b) {if constexpr (vector<A>) return compare_all(a) == b; else return a == compare_all(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator==(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator==(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::equal_to{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator==(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::equal_to{}, a, b.value);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(const A &a, const B &b) {if constexpr (vector<A>) return compare_any(a) != b; else return a != compare_any(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::not_equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::not_equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::not_equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::not_equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::not_equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::not_equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::not_equal_to{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator!=(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::not_equal_to{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator!=(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::not_equal_to{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator!=(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::not_equal_to{}, a, b.value);}
        template <vector_or_scalar_with_base<bool> A, vector_or_scalar_with_base<bool> B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator&&(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) && b; else return a && compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::logical_and{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::logical_and{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::logical_and{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::logical_and{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::logical_and{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::logical_and{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::logical_and{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator&&(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::logical_and{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator&&(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::logical_and{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator&&(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::logical_and{}, a, b.value);}
        template <vector_or_scalar_with_base<bool> A, vector_or_scalar_with_base<bool> B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator||(const A &a, const B &b) {if constexpr (vector<A>) return compare_elemwise(a) || b; else return a || compare_elemwise(b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(compare_any<A> &&a, const B &b) {return any_nonzero_elements(apply_elementwise(std::logical_or{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(const A &a, compare_any<B> &&b) {return any_nonzero_elements(apply_elementwise(std::logical_or{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(compare_all<A> &&a, const B &b) {return all_nonzero_elements(apply_elementwise(std::logical_or{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(const A &a, compare_all<B> &&b) {return all_nonzero_elements(apply_elementwise(std::logical_or{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(compare_none<A> &&a, const B &b) {return none_nonzero_elements(apply_elementwise(std::logical_or{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(const A &a, compare_none<B> &&b) {return none_nonzero_elements(apply_elementwise(std::logical_or{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(compare_not_all<A> &&a, const B &b) {return not_all_nonzero_elements(apply_elementwise(std::logical_or{}, a.value, b));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool operator||(const A &a, compare_not_all<B> &&b) {return not_all_nonzero_elements(apply_elementwise(std::logical_or{}, a, b.value));}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator||(compare_elemwise<A> &&a, const B &b) {return apply_elementwise(std::logical_or{}, a.value, b);}
        template <vector_or_scalar A, vector_or_scalar B> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<common_vec_size_v<vec_size_v<A>, vec_size_v<B>>, bool> operator||(const A &a, compare_elemwise<B> &&b) {return apply_elementwise(std::logical_or{}, a, b.value);}

        //{ input/output
        template <typename A, typename B, int D, typename T> std::basic_ostream<A,B> &operator<<(std::basic_ostream<A,B> &s, const vec<D,T> &v)
        {
            s.width(0);
            s << '[';
            for (int i = 0; i < D; i++)
            {
                if (i != 0)
                    s << ',';
                s << v[i];
            }
            s << ']';
            return s;
        }
        template <typename A, typename B, int W, int H, typename T> std::basic_ostream<A,B> &operator<<(std::basic_ostream<A,B> &s, const mat<W,H,T> &v)
        {
            s.width(0);
            s << '[';
            for (int y = 0; y < H; y++)
            {
                if (y != 0)
                    s << ';';
                for (int x = 0; x < W; x++)
                {
                    if (x != 0)
                        s << ',';
                    s << v[x][y];
                }
            }
            s << ']';
            return s;
        }
        template <typename A, typename B, int D, typename T> std::basic_istream<A,B> &operator>>(std::basic_istream<A,B> &s, vec<D,T> &v)
        {
            s.width(0);
            for (int i = 0; i < D; i++)
                s >> v[i];
            return s;
        }
        template <typename A, typename B, int W, int H, typename T> std::basic_istream<A,B> &operator>>(std::basic_istream<A,B> &s, mat<W,H,T> &v)
        {
            s.width(0);
            for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++)
                s >> v[x][y];
            return s;
        }
        //} input/output

        //{ matrix multiplication
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const mat2x2<A> &a, const vec2<B> &b) {return {a.x.x*b.x + a.y.x*b.y, a.x.y*b.x + a.y.y*b.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const mat3x2<A> &a, const vec3<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const mat4x2<A> &a, const vec4<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z + a.w.x*b.w, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z + a.w.y*b.w};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const mat2x3<A> &a, const vec2<B> &b) {return {a.x.x*b.x + a.y.x*b.y, a.x.y*b.x + a.y.y*b.y, a.x.z*b.x + a.y.z*b.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const mat3x3<A> &a, const vec3<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z, a.x.z*b.x + a.y.z*b.y + a.z.z*b.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const mat4x3<A> &a, const vec4<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z + a.w.x*b.w, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z + a.w.y*b.w, a.x.z*b.x + a.y.z*b.y + a.z.z*b.z + a.w.z*b.w};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const mat2x4<A> &a, const vec2<B> &b) {return {a.x.x*b.x + a.y.x*b.y, a.x.y*b.x + a.y.y*b.y, a.x.z*b.x + a.y.z*b.y, a.x.w*b.x + a.y.w*b.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const mat3x4<A> &a, const vec3<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z, a.x.z*b.x + a.y.z*b.y + a.z.z*b.z, a.x.w*b.x + a.y.w*b.y + a.z.w*b.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const mat4x4<A> &a, const vec4<B> &b) {return {a.x.x*b.x + a.y.x*b.y + a.z.x*b.z + a.w.x*b.w, a.x.y*b.x + a.y.y*b.y + a.z.y*b.z + a.w.y*b.w, a.x.z*b.x + a.y.z*b.y + a.z.z*b.z + a.w.z*b.w, a.x.w*b.x + a.y.w*b.y + a.z.w*b.z + a.w.w*b.w};}
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const vec2<A> &a, const mat2x2<B> &b) {return {a.x*b.x.x + a.y*b.x.y, a.x*b.y.x + a.y*b.y.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const vec3<A> &a, const mat2x3<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec2<larger_t<A,B>> operator*(const vec4<A> &a, const mat2x4<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z + a.w*b.x.w, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z + a.w*b.y.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x2<larger_t<A,B>> operator*(const mat2x2<A> &a, const mat2x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x2<larger_t<A,B>> operator*(const mat3x2<A> &a, const mat2x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x2<larger_t<A,B>> operator*(const mat4x2<A> &a, const mat2x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x3<larger_t<A,B>> operator*(const mat2x3<A> &a, const mat2x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x3<larger_t<A,B>> operator*(const mat3x3<A> &a, const mat2x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x3<larger_t<A,B>> operator*(const mat4x3<A> &a, const mat2x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x4<larger_t<A,B>> operator*(const mat2x4<A> &a, const mat2x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y, a.x.w*b.x.x + a.y.w*b.x.y, a.x.w*b.y.x + a.y.w*b.y.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x4<larger_t<A,B>> operator*(const mat3x4<A> &a, const mat2x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat2x4<larger_t<A,B>> operator*(const mat4x4<A> &a, const mat2x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z + a.w.w*b.x.w, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z + a.w.w*b.y.w};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const vec2<A> &a, const mat3x2<B> &b) {return {a.x*b.x.x + a.y*b.x.y, a.x*b.y.x + a.y*b.y.y, a.x*b.z.x + a.y*b.z.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const vec3<A> &a, const mat3x3<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z, a.x*b.z.x + a.y*b.z.y + a.z*b.z.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec3<larger_t<A,B>> operator*(const vec4<A> &a, const mat3x4<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z + a.w*b.x.w, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z + a.w*b.y.w, a.x*b.z.x + a.y*b.z.y + a.z*b.z.z + a.w*b.z.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x2<larger_t<A,B>> operator*(const mat2x2<A> &a, const mat3x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x2<larger_t<A,B>> operator*(const mat3x2<A> &a, const mat3x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x2<larger_t<A,B>> operator*(const mat4x2<A> &a, const mat3x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x3<larger_t<A,B>> operator*(const mat2x3<A> &a, const mat3x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y, a.x.z*b.z.x + a.y.z*b.z.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x3<larger_t<A,B>> operator*(const mat3x3<A> &a, const mat3x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x3<larger_t<A,B>> operator*(const mat4x3<A> &a, const mat3x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z + a.w.z*b.z.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x4<larger_t<A,B>> operator*(const mat2x4<A> &a, const mat3x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y, a.x.z*b.z.x + a.y.z*b.z.y, a.x.w*b.x.x + a.y.w*b.x.y, a.x.w*b.y.x + a.y.w*b.y.y, a.x.w*b.z.x + a.y.w*b.z.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x4<larger_t<A,B>> operator*(const mat3x4<A> &a, const mat3x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z, a.x.w*b.z.x + a.y.w*b.z.y + a.z.w*b.z.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat3x4<larger_t<A,B>> operator*(const mat4x4<A> &a, const mat3x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z + a.w.z*b.z.w, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z + a.w.w*b.x.w, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z + a.w.w*b.y.w, a.x.w*b.z.x + a.y.w*b.z.y + a.z.w*b.z.z + a.w.w*b.z.w};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const vec2<A> &a, const mat4x2<B> &b) {return {a.x*b.x.x + a.y*b.x.y, a.x*b.y.x + a.y*b.y.y, a.x*b.z.x + a.y*b.z.y, a.x*b.w.x + a.y*b.w.y};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const vec3<A> &a, const mat4x3<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z, a.x*b.z.x + a.y*b.z.y + a.z*b.z.z, a.x*b.w.x + a.y*b.w.y + a.z*b.w.z};}
        template <typename A, typename B> [[nodiscard]] constexpr vec4<larger_t<A,B>> operator*(const vec4<A> &a, const mat4x4<B> &b) {return {a.x*b.x.x + a.y*b.x.y + a.z*b.x.z + a.w*b.x.w, a.x*b.y.x + a.y*b.y.y + a.z*b.y.z + a.w*b.y.w, a.x*b.z.x + a.y*b.z.y + a.z*b.z.z + a.w*b.z.w, a.x*b.w.x + a.y*b.w.y + a.z*b.w.z + a.w*b.w.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x2<larger_t<A,B>> operator*(const mat2x2<A> &a, const mat4x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.x*b.w.x + a.y.x*b.w.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y, a.x.y*b.w.x + a.y.y*b.w.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x2<larger_t<A,B>> operator*(const mat3x2<A> &a, const mat4x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x2<larger_t<A,B>> operator*(const mat4x2<A> &a, const mat4x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z + a.w.x*b.w.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z + a.w.y*b.w.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x3<larger_t<A,B>> operator*(const mat2x3<A> &a, const mat4x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.x*b.w.x + a.y.x*b.w.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y, a.x.y*b.w.x + a.y.y*b.w.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y, a.x.z*b.z.x + a.y.z*b.z.y, a.x.z*b.w.x + a.y.z*b.w.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x3<larger_t<A,B>> operator*(const mat3x3<A> &a, const mat4x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z, a.x.z*b.w.x + a.y.z*b.w.y + a.z.z*b.w.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x3<larger_t<A,B>> operator*(const mat4x3<A> &a, const mat4x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z + a.w.x*b.w.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z + a.w.y*b.w.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z + a.w.z*b.z.w, a.x.z*b.w.x + a.y.z*b.w.y + a.z.z*b.w.z + a.w.z*b.w.w};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x4<larger_t<A,B>> operator*(const mat2x4<A> &a, const mat4x2<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y, a.x.x*b.y.x + a.y.x*b.y.y, a.x.x*b.z.x + a.y.x*b.z.y, a.x.x*b.w.x + a.y.x*b.w.y, a.x.y*b.x.x + a.y.y*b.x.y, a.x.y*b.y.x + a.y.y*b.y.y, a.x.y*b.z.x + a.y.y*b.z.y, a.x.y*b.w.x + a.y.y*b.w.y, a.x.z*b.x.x + a.y.z*b.x.y, a.x.z*b.y.x + a.y.z*b.y.y, a.x.z*b.z.x + a.y.z*b.z.y, a.x.z*b.w.x + a.y.z*b.w.y, a.x.w*b.x.x + a.y.w*b.x.y, a.x.w*b.y.x + a.y.w*b.y.y, a.x.w*b.z.x + a.y.w*b.z.y, a.x.w*b.w.x + a.y.w*b.w.y};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x4<larger_t<A,B>> operator*(const mat3x4<A> &a, const mat4x3<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z, a.x.z*b.w.x + a.y.z*b.w.y + a.z.z*b.w.z, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z, a.x.w*b.z.x + a.y.w*b.z.y + a.z.w*b.z.z, a.x.w*b.w.x + a.y.w*b.w.y + a.z.w*b.w.z};}
        template <typename A, typename B> [[nodiscard]] constexpr mat4x4<larger_t<A,B>> operator*(const mat4x4<A> &a, const mat4x4<B> &b) {return {a.x.x*b.x.x + a.y.x*b.x.y + a.z.x*b.x.z + a.w.x*b.x.w, a.x.x*b.y.x + a.y.x*b.y.y + a.z.x*b.y.z + a.w.x*b.y.w, a.x.x*b.z.x + a.y.x*b.z.y + a.z.x*b.z.z + a.w.x*b.z.w, a.x.x*b.w.x + a.y.x*b.w.y + a.z.x*b.w.z + a.w.x*b.w.w, a.x.y*b.x.x + a.y.y*b.x.y + a.z.y*b.x.z + a.w.y*b.x.w, a.x.y*b.y.x + a.y.y*b.y.y + a.z.y*b.y.z + a.w.y*b.y.w, a.x.y*b.z.x + a.y.y*b.z.y + a.z.y*b.z.z + a.w.y*b.z.w, a.x.y*b.w.x + a.y.y*b.w.y + a.z.y*b.w.z + a.w.y*b.w.w, a.x.z*b.x.x + a.y.z*b.x.y + a.z.z*b.x.z + a.w.z*b.x.w, a.x.z*b.y.x + a.y.z*b.y.y + a.z.z*b.y.z + a.w.z*b.y.w, a.x.z*b.z.x + a.y.z*b.z.y + a.z.z*b.z.z + a.w.z*b.z.w, a.x.z*b.w.x + a.y.z*b.w.y + a.z.z*b.w.z + a.w.z*b.w.w, a.x.w*b.x.x + a.y.w*b.x.y + a.z.w*b.x.z + a.w.w*b.x.w, a.x.w*b.y.x + a.y.w*b.y.y + a.z.w*b.y.z + a.w.w*b.y.w, a.x.w*b.z.x + a.y.w*b.z.y + a.z.w*b.z.z + a.w.w*b.z.w, a.x.w*b.w.x + a.y.w*b.w.y + a.z.w*b.w.z + a.w.w*b.w.w};}

        template <typename A, typename B, int D> constexpr vec<D,A> &operator*=(vec<D,A> &a, const mat<D,D,B> &b) {a = a * b; return a;}
        template <typename A, typename B, int W, int H> constexpr mat<W,H,A> &operator*=(mat<W,H,A> &a, const mat<W,W,B> &b) {a = a * b; return a;}
        //} matrix multiplication
    }

    inline namespace Utility // Low-level helper functions
    {
        //{ Custom operators
        struct op_type_dot {};
        struct op_type_cross {};

        template <typename A> struct op_expr_type_dot
        {
            A &&a;
            template <typename B> [[nodiscard]] constexpr decltype(auto) operator/(B &&b) {return std::forward<A>(a).dot(std::forward<B>(b));}
            template <typename B> constexpr decltype(auto) operator/=(B &&b) {a = std::forward<A>(a).dot(std::forward<B>(b)); return std::forward<A>(a);}
        };
        template <typename A> struct op_expr_type_cross
        {
            A &&a;
            template <typename B> [[nodiscard]] constexpr decltype(auto) operator/(B &&b) {return std::forward<A>(a).cross(std::forward<B>(b));}
            template <typename B> constexpr decltype(auto) operator/=(B &&b) {a = std::forward<A>(a).cross(std::forward<B>(b)); return std::forward<A>(a);}
        };

        template <typename T> inline constexpr op_expr_type_dot<T> operator/(T &&param, op_type_dot) {return {std::forward<T>(param)};}
        template <typename T> inline constexpr op_expr_type_cross<T> operator/(T &&param, op_type_cross) {return {std::forward<T>(param)};}
        //} Custom operators

        //{ Ranges
        template <integral_vector_or_scalar T> class vector_range_t
        {
            T vec_begin{};
            T vec_end{};

          public:
            class iterator
            {
                friend class vector_range_t<T>;

                T vec_begin{};
                T vec_end{};
                T vec_cur{};
                bool finished = true;

                iterator(T vec_begin, T vec_end) : vec_begin(vec_begin), vec_end(vec_end), vec_cur(vec_begin), finished(compare_any(vec_begin) >= vec_end) {}

              public:
                using difference_type   = std::ptrdiff_t;
                using value_type        = T;
                using pointer           = const T *;
                using reference         = const T &;
                using iterator_category = std::forward_iterator_tag;

                iterator() {}

                iterator &operator++()
                {
                    for (int i = 0; i < vec_size_v<T>; i++)
                    {
                        auto &elem = vec_elem(i, vec_cur);
                        elem++;
                        if (elem < vec_elem(i, vec_end))
                            break;
                        elem = vec_elem(i, vec_begin);
                        if (i == vec_size_v<T> - 1)
                            finished = true;
                    }

                    return *this;
                }
                iterator operator++(int)
                {
                    iterator ret = *this;
                    ++(*this);
                    return ret;
                }

                reference operator*() const
                {
                    return vec_cur;
                }
                pointer operator->() const
                {
                    return &vec_cur;
                }

                bool operator==(const iterator &other) const
                {
                    if (finished != other.finished)
                        return false;
                    if (finished && other.finished)
                        return true;
                    return vec_cur == other.vec_cur;
                }
            };

            vector_range_t() {}
            vector_range_t(T vec_begin, T vec_end) : vec_begin(vec_begin), vec_end(vec_end) {}

            iterator begin() const
            {
                return iterator(vec_begin, vec_end);
            }

            iterator end() const
            {
                return {};
            }

            [[nodiscard]] friend vector_range_t operator+(const vector_range_t &range, std::same_as<T> auto offset)
            {
                return vector_range_t(range.vec_begin + offset, range.vec_end + offset);
            }
            [[nodiscard]] friend vector_range_t operator+(std::same_as<T> auto offset, const vector_range_t &range)
            {
                return range + offset;
            }
        };

        template <integral_vector_or_scalar T> class vector_range_halfbound
        {
            T vec_begin{};

          public:
            vector_range_halfbound(T vec_begin) : vec_begin(vec_begin) {}

            [[nodiscard]] friend vector_range_t<T> operator<(const vector_range_halfbound &range, std::same_as<T> auto point)
            {
                return vector_range_t<T>(range.vec_begin, point);
            }
            [[nodiscard]] friend vector_range_t<T> operator<=(const vector_range_halfbound &range, std::same_as<T> auto point)
            {
                return range < point+1;
            }
        };

        struct vector_range_factory
        {
            template <integral_vector_or_scalar T> vector_range_t<T> operator()(T size) const
            {
                return vector_range_t<T>(T(0), size);
            }

            template <int D, typename T> vector_range_t<vec<D,T>> operator()(rect<D,T> r) const
            {
                return vector_range_t<vec<D,T>>(r.a, r.size());
            }

            template <integral_vector_or_scalar T> friend vector_range_halfbound<T> operator<=(T point, vector_range_factory)
            {
                return {point};
            }
            template <integral_vector_or_scalar T> friend vector_range_halfbound<T> operator<(T point, vector_range_factory)
            {
                return point+1 <= vector_range_factory{};
            }
        };
        //} Ranges
    }

    inline namespace Common // Common functions
    {
        // Named operators.
        inline constexpr op_type_dot dot;
        inline constexpr op_type_cross cross;

        // Comparison tags.
        inline constexpr compare_any_tag any;
        inline constexpr compare_all_tag all;
        inline constexpr compare_none_tag none;
        inline constexpr compare_not_all_tag not_all;
        inline constexpr compare_elemwise_tag elemwise;

        // Helper class for writing nested loops.
        // Example usage:
        //   for (auto v : vec_a <= vector_range <= vec_b) // `<` are also allowed, in one or both positions.
        //   for (auto v : vector_range(vec_a)) // Equivalent to `vec..(0) <= vector_range < vec_a`.
        inline constexpr vector_range_factory vector_range;

        // The value of pi.
        template <scalar T> [[nodiscard]] constexpr T pi() {return T(3.14159265358979323846l);}
        constexpr float       f_pi  = pi<float>();
        constexpr double      d_pi  = pi<double>();
        constexpr long double ld_pi = pi<long double>();

        // Conversions between degrees and radians.
        template <vector_or_scalar T> [[nodiscard]] constexpr auto to_rad(T in)
        {
            using fp_t = floating_point_t<T>;
            return in * pi<fp_t>() / fp_t(180);
        }
        template <vector_or_scalar T> [[nodiscard]] constexpr auto to_deg(T in)
        {
            using fp_t = floating_point_t<T>;
            return in * fp_t(180) / pi<fp_t>();
        }

        // Returns the sign of the argument as `int` or `ivecN`.
        template <vector_or_scalar T> [[nodiscard]] constexpr change_vec_base_t<T,int> sign(T val)
        {
            // Works on scalars and vectors.
            return (val > 0) - (val < 0);
        }
        // Returns the sign of `a - b`. Unlike `sign(a - b)`, not affected by overflow.
        // Refuses to work if one of the arguments is a signed integer, and the other is unsigned.
        template <vector_or_scalar A, vector_or_scalar B> requires have_larger_type<A, B>
        [[nodiscard]] constexpr auto diffsign(A a, B b) -> vec_or_scalar_t<common_vec_size_v<vec_size_v<A>,vec_size_v<B>>,int>
        {
            // Works on scalars and vectors.
            return (a > b) - (a < b);
        }

        // `clamp[_var][_min|_max|_abs] (value, min, max)`.
        // Clamps scalars or vectors.
        // `_var` functions modify the first parameter instead of returning the result.
        // `_min` functions don't have a `max` parameter, and vice versa.
        // `_abs` functions don't have a `min` parameter, they use `-max` as `min`.
        // If both `min` and `max` are omitted, 0 and 1 are assumed.
        // If bounds contradict each other, only the `max` bound is used.

        template <vector_or_scalar A, safely_convertible_to<A> B>
        constexpr void clamp_var_min(A &var, B min)
        {
            if constexpr (!any_vectors_v<A,B>)
            {
                if (!(var >= min)) // The condition is written like this to catch NaNs, they always compare to false.
                    var = min;
            }
            else
            {
                apply_elementwise(clamp_var_min<vec_base_t<A>, vec_base_t<B>>, var, min);
            }
        }

        template <vector_or_scalar A, safely_convertible_to<A> B>
        constexpr void clamp_var_max(A &var, B max)
        {
            if constexpr (!any_vectors_v<A,B>)
            {
                if (!(var <= max)) // The condition is written like this to catch NaNs, they always compare to false.
                    var = max;
            }
            else
            {
                apply_elementwise(clamp_var_max<vec_base_t<A>, vec_base_t<B>>, var, max);
            }
        }

        template <vector_or_scalar A, safely_convertible_to<A> B, safely_convertible_to<A> C>
        constexpr void clamp_var(A &var, B min, C max)
        {
            clamp_var_min(var, min);
            clamp_var_max(var, max);
        }

        template <vector_or_scalar A, safely_convertible_to<A> B> requires signed_maybe_floating_point_vector_or_scalar<B>
        constexpr void clamp_var_abs(A &var, B abs_max)
        {
            clamp_var(var, -abs_max, abs_max);
        }

        template <vector_or_scalar A, safely_convertible_to<A> B>
        [[nodiscard]] constexpr A clamp_min(A val, B min)
        {
            clamp_var_min(val, min);
            return val;
        }

        template <vector_or_scalar A, safely_convertible_to<A> B>
        [[nodiscard]] constexpr A clamp_max(A val, B max)
        {
            clamp_var_max(val, max);
            return val;
        }

        template <vector_or_scalar A, safely_convertible_to<A> B, safely_convertible_to<A> C>
        [[nodiscard]] constexpr A clamp(A val, B min, C max)
        {
            clamp_var(val, min, max);
            return val;
        }

        template <vector_or_scalar A, safely_convertible_to<A> B> requires signed_maybe_floating_point_vector_or_scalar<B>
        [[nodiscard]] constexpr A clamp_abs(A val, B abs_max)
        {
            clamp_var_abs(val, abs_max);
            return val;
        }

        template <vector_or_scalar A> [[nodiscard]] constexpr A clamp(A val) {return clamp(val, 0, 1);}
        template <vector_or_scalar A> [[nodiscard]] constexpr A clamp_min(A val) {return clamp_min(val, 0);}
        template <vector_or_scalar A> [[nodiscard]] constexpr A clamp_max(A val) {return clamp_max(val, 1);}
        template <vector_or_scalar A> [[nodiscard]] constexpr A clamp_abs(A val) {return clamp_abs(val, 1);}
        template <vector_or_scalar A> constexpr void clamp_var(A &var) {clamp_var(var, 0, 1);}
        template <vector_or_scalar A> constexpr void clamp_var_min(A &var) {clamp_var_min(var, 0);}
        template <vector_or_scalar A> constexpr void clamp_var_max(A &var) {clamp_var_max(var, 1);}
        template <vector_or_scalar A> constexpr void clamp_var_abs(A &var) {clamp_var_abs(var, 1);}

        // Rounds a floating-point scalar or vector.
        // Returns an integral type (`int` by default).
        template <signed_integral_scalar I = int, floating_point_vector_or_scalar F>
        [[nodiscard]] change_vec_base_t<F,I> iround(F x)
        {
            if constexpr (!any_vectors_v<F>)
            {
                // This seems to be faster than `std::lround()`.
                return I(std::round(x));
            }
            else
            {
                return apply_elementwise(iround<I, vec_base_t<F>>, x);
            }
        }

        // Various useful functions.
        // Some of them are imported from `std` and extended to operate on vectors. Some are custom.

        using std::abs;
        template <vector T>
        [[nodiscard]] T abs(T x)
        {
            return apply_elementwise([](auto val){return std::abs(val);}, x);
        }

        using std::round;
        template <floating_point_vector T>
        [[nodiscard]] T round(T x)
        {
            return apply_elementwise([](auto val){return std::round(val);}, x);
        }

        using std::floor;
        template <floating_point_vector T>
        [[nodiscard]] T floor(T x)
        {
            return apply_elementwise([](auto val){return std::floor(val);}, x);
        }

        using std::ceil;
        template <floating_point_vector T>
        [[nodiscard]] T ceil(T x)
        {
            return apply_elementwise([](auto val){return std::ceil(val);}, x);
        }

        using std::trunc;
        template <floating_point_vector T>
        [[nodiscard]] T trunc(T x)
        {
            return apply_elementwise([](auto val){return std::trunc(val);}, x);
        }

        template <floating_point_vector T>
        [[nodiscard]] T round_maxabs(T x) // Round away from zero.
        {
            return apply_elementwise([](auto val){return val < 0 ? std::floor(val) : std::ceil(val);}, x);
        }

        template <floating_point_vector T>
        [[nodiscard]] T frac(T x)
        {
            if constexpr (!any_vectors_v<T>)
                return std::modf(x, 0);
            else
                return apply_elementwise(frac<vec_base_t<T>>, x);
        }

        using std::nextafter;
        template <floating_point_vector_or_scalar A, floating_point_vector_or_scalar B>
        requires any_vectors_v<A, B> && std::is_same_v<vec_base_t<A>, vec_base_t<B>> && have_larger_type<A, B>
        [[nodiscard]] A nextafter(A a, B b)
        {
            return apply_elementwise([](auto a, auto b){return std::nextafter(a, b);}, a, b);
        }

        // Integer division, slightly changed to behave nicely for negative values of the left operand:
        //           i : -4  -3  -2  -1  0  1  2  3  4
        // div_ex(i,2) : -2  -2  -1  -1  0  0  1  1  2
        template <integral_vector_or_scalar A, integral_vector_or_scalar B>
        [[nodiscard]] constexpr A div_ex(A a, B b)
        {
            if constexpr (!any_vectors_v<A,B>)
            {
                if (a >= 0)
                    return a / b;
                else
                    return (a + 1) / b - sign(b);
            }
            else
            {
                return apply_elementwise(div_ex<vec_base_t<A>, vec_base_t<B>>, a, b);
            }
        }

        // True integral modulo that remains periodic for negative values of the left operand.
        template <integral_vector_or_scalar A, integral_vector_or_scalar B>
        [[nodiscard]] constexpr A mod_ex(A a, B b)
        {
            if constexpr (!any_vectors_v<A,B>)
            {
                if (a >= 0)
                    return a % b;
                else
                    return abs(b) - 1 + (a + 1) % b;
            }
            else
            {
                return apply_elementwise(mod_ex<vec_base_t<A>, vec_base_t<B>>, a, b);
            }
        }

        // Divide `a / b`, rounding away from zero.
        // Supports both integers and floating-point numbers, including vectors.
        template <signed_maybe_floating_point_vector_or_scalar A, signed_maybe_floating_point_vector_or_scalar B>
        [[nodiscard]] constexpr larger_t<A, B> div_maxabs(A a, B b)
        {
            if constexpr (!any_vectors_v<A, B>)
            {
                if constexpr (integral_scalar<A> && integral_scalar<B>)
                {
                    return (a + (abs(b) - 1) * sign(a)) / b;
                }
                else
                {
                    using T = larger_t<A, B>;
                    T ret = T(a) / T(b);
                    return round_maxabs(ret);
                }
            }
            else
            {
                return apply_elementwise(div_maxabs<vec_base_t<A>, vec_base_t<B>>, a, b);
            }
        }

        // A simple implementation of `pow` for non-negative integral powers.
        template <vector_or_scalar A, integral_scalar B>
        [[nodiscard]] constexpr A ipow(A a, B b)
        {
            A ret = 1;
            while (b > 0)
            {
                if (b & 1)
                ret *= a;
                a *= a;
                b >>= 1;
            }
            return ret;
        }

        using std::pow;
        template <vector_or_scalar A, vector_or_scalar B>
        requires any_vectors_v<A, B>
        [[nodiscard]] auto pow(A a, B b)
        {
            return apply_elementwise([](auto val_a, auto val_b){return std::pow(val_a, val_b);}, a, b);
        }

        // Computes the smooth step function. Doesn't clamp `x`.
        template <floating_point_vector_or_scalar T>
        [[nodiscard]] constexpr T smoothstep(T x)
        {
            // No special handling required for `T` being a vector.
            return (3 - 2*x) * x*x;
        }

        // Performs linear interpolation. Returns `a * (1-factor) + b * factor`.
        template <floating_point_scalar F, vector_or_scalar A, vector_or_scalar B>
        requires have_larger_type<A, B>
        [[nodiscard]] constexpr auto mix(F factor, A a, B b)
        {
            // No special handling required for the parameters being vectors.
            using type = larger_t<A, B>;
            return type(a) * (1-factor) + type(b) * factor;
        }

        // Returns a `min` or `max` value of the parameters.
        template <typename ...P> [[nodiscard]] constexpr larger_t<P...> min(P ... params)
        {
            if constexpr (!any_vectors_v<P...>)
                return std::min({larger_t<P...>(params)...});
            else
                return apply_elementwise(min<vec_base_t<P>...>, params...);
        }
        template <typename ...P> [[nodiscard]] constexpr larger_t<P...> max(P ... params)
        {
            if constexpr (!any_vectors_v<P...>)
                return std::max({larger_t<P...>(params)...});
            else
                return apply_elementwise(max<vec_base_t<P>...>, params...);
        }

        // Returns `[min(a,b), max(a,b)]`. Like `std::minmax`, but returns by value and can handle vectors.
        template <typename A, typename B> [[nodiscard]] constexpr std::pair<larger_t<A, B>, larger_t<A, B>> sort_two(A a, B b)
        {
            using T = larger_t<A, B>;
            std::pair<T, T> ret;
            for (int i = 0; i < vec_size_weak_v<T>; i++)
            {
                auto a_elem = vec_elem(i, a);
                auto b_elem = vec_elem(i, b);
                if (b_elem < a_elem)
                    vec_elem(i, ret.first) = b_elem, vec_elem(i, ret.second) = a_elem;
                else
                    vec_elem(i, ret.first) = a_elem, vec_elem(i, ret.second) = b_elem;
            }
            return ret;
        }
        // Sorts `{a,b}` in place. Sorts vectors element-wise.
        template <typename T> constexpr void sort_two_var(T &a, T &b)
        {
            if constexpr (!any_vectors_v<T>)
            {
                if (b < a)
                    std::swap(a, b);
            }
            else
            {
                apply_elementwise(sort_two_var<vec_base_t<T>>, a, b);
            }
        }
    }

    inline namespace Misc // Misc functions
    {
        // A functor that performs linear mapping on scalars or vectors.
        template <floating_point_vector_or_scalar T>
        struct linear_mapping
        {
            T scale = T(1), offset = T(0);

            constexpr linear_mapping() {}

            constexpr linear_mapping(T src_a, T src_b, T dst_a, T dst_b)
            {
                T factor = 1 / (src_a - src_b);
                scale = (dst_a - dst_b) * factor;
                offset = (dst_b * src_a - dst_a * src_b) * factor;
            }

            constexpr T operator()(T x) const
            {
                return x * scale + offset;
            }

            using matrix_t = mat<vec_size_v<T>+1, vec_size_v<T>+1, vec_base_t<T>>;
            constexpr matrix_t matrix() const
            {
                matrix_t ret{};
                for (int i = 0; i < vec_size_v<T>; i++)
                {
                    ret[i][i] = scale[i];
                    ret[vec_size_v<T>][i] = offset[i];
                }
                return ret;
            }
        };

        // Like `nextafter()`, but works with integers as well.
        template <vector_or_scalar A, vector_or_scalar B>
        [[nodiscard]] larger_t<A, B> next_value_towards(A value, B target)
        {
            using type = larger_t<A, B>;
            if constexpr (floating_point_vector_or_scalar<type>)
                return nextafter(type(value), type(target));
            else
                return type(value) + diffsign(type(target), type(value)); // The plain `sign()` could overflow here.
        }
        // Returns the next or previous representable value.
        // Refuses to increment the largest representable value, and returns it unchanged.
        // If asked to increment infinity in either direction, returns the closest representable value.
        // If given NaN, returns NaN.
        template <bool Prev, vector_or_scalar T>
        [[nodiscard]] T next_or_prev_value(T value)
        {
            return next_value_towards(value, Prev ? std::numeric_limits<vec_base_t<T>>::lowest() : std::numeric_limits<vec_base_t<T>>::max());
        }
        template <vector_or_scalar T> [[nodiscard]] T next_value(T value) {return next_or_prev_value<false>(value);}
        template <vector_or_scalar T> [[nodiscard]] T prev_value(T value) {return next_or_prev_value<true >(value);}

        // Shrinks a vector as little as possible to give it specific proportions.
        // Always returns a floating-point type.
        template <vector A, vector B> requires have_larger_type<A, B>
        [[nodiscard]] constexpr auto shrink_to_proportions(A value, B proportions)
        {
            using type = larger_t<floating_point_t<A>,floating_point_t<B>>;
            return (type(value) / type(proportions)).min() * type(proportions);
        }
        // Expands a vector as little as possible to give it specific proportions.
        // Always returns a floating-point type.
        template <vector A, vector B> requires have_larger_type<A, B>
        [[nodiscard]] constexpr auto expand_to_proportions(A value, B proportions)
        {
            using type = larger_t<floating_point_t<A>,floating_point_t<B>>;
            return (type(value) / type(proportions)).max() * type(proportions);
        }

        // Given two lines, each defined by a point `p?` an a direction `d?`, returns the intersection point offset.
        // If multiplied by `d1` and added to `p1`, that gives the intersection point.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr T point_dir_intersection_factor(vec2<T> p1, vec2<T> d1, vec2<T> p2, vec2<T> d2)
        {
            // This was solved programmatically, not sure how exactly it works.
            return ((p2 - p1) /cross/ d2) / (d1 /cross/ d2);
        }
        // Same, but returns the factors for both lines.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr std::array<T, 2> point_dir_intersection_factor_two_way(vec2<T> p1, vec2<T> d1, vec2<T> p2, vec2<T> d2)
        {
            // This was solved programmatically, not sure how exactly it works.
            vec2<T> p = p2 - p1;
            T d = d1 /cross/ d2;
            return {(p /cross/ d2) / d, (p /cross/ d1) / d};
        }
        // Given two lines, each defined by a point `p?` an a direction `d?`, returns the intersection point.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr vec2<T> point_dir_intersection(vec2<T> p1, vec2<T> d1, vec2<T> p2, vec2<T> d2)
        {
            return p1 + d1 * (point_dir_intersection_factor)(p1, d1, p2, d2);
        }
        // Given two lines, each defined by two points, returns the intersection point.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr vec2<T> line_intersection(vec2<T> a1, vec2<T> b1, vec2<T> a2, vec2<T> b2)
        {
            return (point_dir_intersection)(a1, b1 - a1, a2, b2 - a2);
        }

        // Finds an intersection point of a line and a plane.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr vec3<T> line_plane_intersection(vec3<T> line_point, vec3<T> line_dir, vec3<T> plane_point, vec3<T> plane_normal)
        {
            return (plane_point - line_point).dot(plane_normal) / line_dir.dot(plane_normal) * line_dir + line_point;
        }

        // Projects a point onto a line. `dir` is assumed to be normalized.
        template <vector T>
        [[nodiscard]] constexpr T project_onto_line_norm(T point, T dir)
        {
            return dir * point.dot(dir);
        }
        // Projects a point onto a line.
        template <floating_point_vector T>
        [[nodiscard]] constexpr T project_onto_line(T point, T dir)
        {
            return project_onto_line_norm(point, dir.norm());
        }

        // Projects a point onto a plane. `plane_normal` is assumed to be normalized.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr vec3<T> project_onto_plane_norm(vec3<T> point, vec3<T> plane_normal)
        {
            return point - project_onto_line_norm(point, plane_normal);
        }
        // Projects a point onto a plane.
        template <floating_point_scalar T>
        [[nodiscard]] constexpr vec3<T> project_onto_plane(vec3<T> point, vec3<T> plane_normal)
        {
            return project_onto_plane_norm(point, plane_normal.norm());
        }

        // Compares the angles of `a` and `b` without doing any trigonometry. Works with integers too.
        // The assumed angles are in range [0;2pi), with +X having angle 0.
        // Zero vectors are considered to be greater than everything else.
        template <scalar T>
        [[nodiscard]] constexpr bool less_positively_rotated(vec2<T> a, vec2<T> b)
        {
            // This check makes (0,0) worse than any other vector,
            // and doesn't seem to affect the result if zero vectors are not involved.
            if (int d = (a == vec2<T>()) - (b == vec2<T>()))
                return d < 0;

            if (int d = (a.y < 0) - (b.y < 0))
                return d < 0;
            if (int d = (a.y == 0 && a.x < 0) - (b.y == 0 && b.x < 0))
                return d < 0;

            return a.x * b.y > b.x * a.y;
        }

        // Same, but angle 0 is mapped to `dir` instead of +X.
        template <scalar T>
        [[nodiscard]] constexpr bool less_positively_rotated(vec2<T> dir, vec2<T> a, vec2<T> b)
        {
            mat2<T> mat(dir, dir.rot90());
            return less_positively_rotated(a * mat, b * mat);
        }

        // Rounds `value` to type `I`, with compensation: `comp` is added to it before rounding, then updated to the difference between rounded and unrounded value.
        // This makes the average return value converge to `value`.
        template <integral_scalar I = int, floating_point_vector_or_scalar F>
        [[nodiscard]] constexpr change_vec_base_t<F,I> round_with_compensation(F value, F &comp)
        {
            // Works on scalars and vectors.
            change_vec_base_t<F,I> ret = iround<I>(value += comp);
            comp = value - ret;
            return ret;
        }

        // Produces points to fill a cuboid (line, rect, cube, and so on), either entirely or only the borders.
        // `a` and `b` are the corners, inclusive. `step` is the step, the sign is ignored.
        // `pred` lets you select what parts of the cuboid to output. It's is either `nullptr` (output everything)
        // or `bool pred(unsigned int mask)`, where the mask receives all combinations of N bits, where N is `vec_size_v<T>`.
        // If `pred` returns true, the corresponding region is emitted using repeated calls to `func`, which is `bool func(T &&point)`.
        // If `func` returns true, the function stops immediately and also returns true. Otherwise returns false when done.
        // The number of `1`s in the mask (`std::popcount(mask)`) describes the dimensions of the region: 0 = points, 1 = lines, 2 = rects, and so on.
        // If the i-th bit is set, the region extends in i-th dimension. Each mask corresponds to a set of parallel lines/planes/etc,
        // and the zero mask corresponds to the corners of the cuboid.
        template <signed_maybe_floating_point_vector_or_scalar T, typename F1 = std::nullptr_t, typename F2>
        bool for_each_cuboid_point(T a, T b, T step, F1 &&pred, F2 &&func)
        {
            // Fix the sign of the `step`.
            for (int i = 0; i < vec_size_v<T>; i++)
            {
                vec_elem(i, step) *= sign(vec_elem(i, b) - vec_elem(i, a)) * sign(vec_elem(i, step));
                // We don't want zero step.
                if (vec_elem(i, step) == 0) vec_elem(i, step) = 1;
            }

            using int_vec = change_vec_base_t<T, int>;
            int_vec count = abs(div_maxabs(b - a, step)) - 1;

            if constexpr (std::is_null_pointer_v<std::remove_cvref_t<F1>>)
            {
                // A simple algorithm to fill the whole cuboid.
                for (int_vec pos : vector_range(count + 2))
                {
                    T value;
                    for (int i = 0; i < vec_size_v<T>; i++)
                        vec_elem(i, value) = vec_elem(i, pos) == vec_elem(i, count) + 1 ? vec_elem(i, b) : vec_elem(i, a) + vec_elem(i, step) * vec_elem(i, pos);
                    if (func(std::move(value)))
                        return true;
                }
            }
            else
            {
                // A more advanced algorithm to control separate regions.
                for (unsigned int i = 0; i < 1u << vec_size_v<T>; i++)
                {
                    // Stop early if we don't want this region.
                    // The casts stop `pred` from doing weird things.
                    if (!bool(pred((unsigned int)i)))
                        continue;

                    // Get the number of points in the region, in each dimension.
                    bool bad_region = false;
                    int_vec region_size;
                    for (int j = 0; j < vec_size_v<T>; j++)
                    {
                        if (i & 1u << j)
                        {
                            if ((vec_elem(j, region_size) = vec_elem(j, count)) <= 0)
                            {
                                bad_region = true;
                                break;
                            }
                        }
                        else
                        {
                            vec_elem(j, region_size) = vec_elem(j, a) == vec_elem(j, b) ? 1 : 2;
                        }
                    }
                    if (bad_region)
                        continue; // A degenerate region.

                    // Output points.
                    for (int_vec pos : vector_range(region_size))
                    {
                        T value;
                        for (int j = 0; j < vec_size_v<T>; j++)
                        {
                            if (!(i & 1u << j))
                                vec_elem(j, value) = vec_elem(j, vec_elem(j, pos) ? b : a);
                            else
                                vec_elem(j, value) = vec_elem(j, a) + (vec_elem(j, pos) + 1) * vec_elem(j, step);
                        }
                        if (func(std::move(value)))
                            return true;
                    }
                }
            }

            return false;
        }

        // Produces points to fill a cuboid (line, rect, cube, and so on), either entirely or only the borders. Writes the points of type `T` to `*iter++`.
        // `a` and `b` are the corners, inclusive. `step` is the step, the sign is ignored.
        // `D` is the dimensions of the output. `D == -1` and `D == vec_size_v<T>` mean "fill the whole cuboid".
        // `D == 0` only outputs the corner points, `D == 1` outputs lines, `D == 2` outputs planes, and so on.
        template <int D = -1, signed_maybe_floating_point_vector_or_scalar T, typename I>
        requires(D >= -1 && D <= vec_size_v<T>)
        void make_cuboid(T a, T b, T step, I iter)
        {
            if constexpr (D == -1 || D == vec_size_v<T>)
            {
                for_each_cuboid_point(a, b, step, nullptr, [&](T &&point)
                {
                    *iter++ = std::move(point);
                    return false;
                });
            }
            else
            {
                for_each_cuboid_point(a, b, step, [](unsigned int mask)
                {
                    return std::popcount(mask) <= D;
                },
                [&](T &&point)
                {
                    *iter++ = std::move(point);
                    return false;
                });
            }
        }

        // Same, but writes the output to a container.
        template <typename C, int D = -1, signed_maybe_floating_point_vector_or_scalar T>
        [[nodiscard]] C make_cuboid(T a, T b, T step)
        {
            C ret;
            make_cuboid(a, b, step, std::back_inserter(ret));
            return ret;
        }
    }

    inline namespace Vector // Definitions
    {
        //{ Vectors
        template <typename T> struct vec<2,T> // vec2
        {
            using type = T;
            using rect_type = rect2<T>;
            static constexpr int size = 2;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            type x, y;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            IMP_MATH_SMALL_FUNC constexpr vec() : x{}, y{} {}
            IMP_MATH_SMALL_FUNC constexpr vec(uninit) {}
            IMP_MATH_SMALL_FUNC constexpr vec(type x, type y) : x(x), y(y) {}
            IMP_MATH_SMALL_FUNC explicit constexpr vec(type obj) : x(obj), y(obj) {}
            template <scalar U> IMP_MATH_SMALL_FUNC explicit(!safely_convertible_to<U,type>) constexpr vec(vec2<U> obj) : x(obj.x), y(obj.y) {}
            template <typename U> requires Custom::convertible<U, vec> explicit(Custom::Convert<U, vec>::is_explicit) constexpr vec(const U &obj) {*this = Custom::Convert<U, vec>{}(obj);}
            template <typename U> requires Custom::convertible<vec, U> explicit(Custom::Convert<vec, U>::is_explicit) operator U() const {return Custom::Convert<vec, U>{}(*this);}
            template <scalar U> [[nodiscard]] constexpr vec2<U> to() const {return vec2<U>(U(x), U(y));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      type *)((      char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const type *)((const char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x;}
            [[nodiscard]] explicit constexpr operator bool() const requires(!std::is_same_v<type, bool>) {return any();} // Use the explicit methods below for vectors of bool.
            [[nodiscard]] constexpr bool any() const {return x || y;}
            [[nodiscard]] constexpr bool all() const {return x && y;}
            [[nodiscard]] constexpr bool none() const {return !any();}
            [[nodiscard]] constexpr bool not_all() const {return !all();}
            [[nodiscard]] constexpr auto sum() const {return x + y;}
            [[nodiscard]] constexpr auto diff() const {return x - y;}
            [[nodiscard]] constexpr auto prod() const {return x * y;}
            [[nodiscard]] constexpr auto ratio() const {return floating_point_t<type>(x) / floating_point_t<type>(y);}
            [[nodiscard]] constexpr type min() const {return std::min({x,y});}
            [[nodiscard]] constexpr type max() const {return std::max({x,y});}
            [[nodiscard]] constexpr vec abs() const {return vec(std::abs(x), std::abs(y));}
            template <typename C> [[nodiscard]] constexpr auto index(C &&container) const -> vec<common_vec_size_v<2,vec_size_v<std::decay_t<decltype(container[x])>>>,vec_base_t<std::decay_t<decltype(container[x])>>> {return {Math::vec_elem(0,container[x]), Math::vec_elem(1,container[y])};}
            [[nodiscard]] constexpr vec3<type> to_vec3(type nz) const {return {x, y, nz};}
            [[nodiscard]] constexpr vec4<type> to_vec4(type nz, type nw) const {return {x, y, nz, nw};}
            [[nodiscard]] constexpr vec3<type> to_vec3() const {return {x, y, 0};}
            [[nodiscard]] constexpr vec4<type> to_vec4() const {return {x, y, 0, 0};}
            [[nodiscard]] constexpr auto len_sq() const {return x*x + y*y;}
            [[nodiscard]] constexpr floating_point_t<type> len() const {return std::sqrt(floating_point_t<type>(len_sq()));}
            [[nodiscard]] constexpr floating_point_t<vec> norm() const {if (auto l = len()) return *this / l; else return vec(0);}
            [[nodiscard]] constexpr auto approx_len() const {return floating_point_t<type>(len_sq() + 1) / 2;} // Accurate only around `len()==1`. Always greater or equal than the true length.
            [[nodiscard]] constexpr auto approx_inv_len() const {return 2 / floating_point_t<type>(len_sq() + 1);}
            [[nodiscard]] constexpr auto approx_norm() const {return *this * approx_inv_len();} // Guaranteed to converge to `len()==1` eventually, when starting from any finite `len_sq()`.
            [[nodiscard]] static constexpr vec axis(int a, type len = 1) {vec ret{}; ret[mod_ex(a,2)] = len; return ret;}
            [[nodiscard]] constexpr vec only_component(int a) {vec ret{}; a = mod_ex(a,2); ret[a] = (*this)[a]; return ret;}
            [[nodiscard]] static constexpr vec dir(type angle, type len = 1) requires is_floating_point {return vec(std::cos(angle) * len, std::sin(angle) * len);}
            template <scalar U = floating_point_t<type>> [[nodiscard]] constexpr U angle() const {return std::atan2(U(y), U(x));}
            [[nodiscard]] constexpr vec rot90(int steps = 1) const {switch (steps & 3) {default: return *this; case 1: return {-y,x}; case 2: return -*this; case 3: return {y,-x};}}
            [[nodiscard]] static constexpr vec dir4(int index, type len = 1) {return vec(len,0).rot90(index);}
            [[nodiscard]] static constexpr vec dir4_diag(int index, type len = 1) {return vec(len,len).rot90(index);}
            [[nodiscard]] static constexpr vec dir8(int index, type len = 1) {vec array[8]{vec(len,0),vec(len,len),vec(0,len),vec(-len,len),vec(-len,0),vec(-len,-len),vec(0,-len),vec(len,-len)}; return array[index & 7];}
            [[nodiscard]] constexpr int angle4_round() const {type s = sum(); type d = diff(); return d<0&&s>=0?1:x<0&&d<=0?2:y<0&&s<=0?3:0;} // Non-cardinal directions round to the closest one, diagnoals round backwards, (0,0) returns zero.
            [[nodiscard]] constexpr int angle4_floor() const {return y>0&&x<=0?1:x<0?2:y<0?3:0;}
            [[nodiscard]] constexpr int angle8_sign() const {return y>0?(x>0?1:x==0?2:3):y<0?(x<0?5:x==0?6:7):(x<0?4:0);} // Non-cardinal directions count as diagonals, (0,0) returns zero.
            [[nodiscard]] constexpr int angle8_floor() const {type s = sum(); type d = diff(); return y<0&&d>=0?(x<0?5:s<0?6:7):x<=0&&d<0?(y<=0?4:s<=0?3:2):y>0&&d<=0?1:0;}
            template <typename U> [[nodiscard]] constexpr auto dot(const vec2<U> &o) const {return x * o.x + y * o.y;}
            template <typename U> [[nodiscard]] constexpr auto cross(const vec2<U> &o) const {return x * o.y - y * o.x;}
            [[nodiscard]] constexpr auto tie() & {return std::tie(x,y);}
            [[nodiscard]] constexpr auto tie() const & {return std::tie(x,y);}
            template <int I> [[nodiscard]] constexpr type &get() & {return std::get<I>(tie());}
            template <int I> [[nodiscard]] constexpr const type &get() const & {return std::get<I>(tie());}
            [[nodiscard]] friend constexpr std::array<T,2> format_as(const vec &v) {return {v.x,v.y};}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_any<vec> operator()(compare_any_tag) const {return compare_any(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_all<vec> operator()(compare_all_tag) const {return compare_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_none<vec> operator()(compare_none_tag) const {return compare_none(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_not_all<vec> operator()(compare_not_all_tag) const {return compare_not_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_elemwise<vec> operator()(compare_elemwise_tag) const {return compare_elemwise(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr mat2<type> to_rotation_matrix() const {return mat2<type>(*this, rot90());}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<type> tiny_rect() const {return rect_to(next_value(*this));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> rect_to(vec2<U> b) const {rect2<larger_t<type,U>> ret; ret.a = *this; ret.b = b; return ret;}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> rect_size(vec2<U> b) const {return rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> rect_size(U b) const {return rect_size(vec2<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> centered_rect_size(vec2<U> b) const {return (*this - b/2).rect_size(b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> centered_rect_size(U b) const {return centered_rect_size(vec2<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> centered_rect_halfsize(vec2<U> b) const {return (*this - b).rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect2<larger_t<type,U>> centered_rect_halfsize(U b) const {return centered_rect_halfsize(vec2<U>(b));}
        };

        template <typename T> struct vec<3,T> // vec3
        {
            using type = T;
            using rect_type = rect3<T>;
            static constexpr int size = 3;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            type x, y, z;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            IMP_MATH_SMALL_FUNC constexpr vec() : x{}, y{}, z{} {}
            IMP_MATH_SMALL_FUNC constexpr vec(uninit) {}
            IMP_MATH_SMALL_FUNC constexpr vec(type x, type y, type z) : x(x), y(y), z(z) {}
            IMP_MATH_SMALL_FUNC explicit constexpr vec(type obj) : x(obj), y(obj), z(obj) {}
            template <scalar U> IMP_MATH_SMALL_FUNC explicit(!safely_convertible_to<U,type>) constexpr vec(vec3<U> obj) : x(obj.x), y(obj.y), z(obj.z) {}
            template <typename U> requires Custom::convertible<U, vec> explicit(Custom::Convert<U, vec>::is_explicit) constexpr vec(const U &obj) {*this = Custom::Convert<U, vec>{}(obj);}
            template <typename U> requires Custom::convertible<vec, U> explicit(Custom::Convert<vec, U>::is_explicit) operator U() const {return Custom::Convert<vec, U>{}(*this);}
            template <scalar U> [[nodiscard]] constexpr vec3<U> to() const {return vec3<U>(U(x), U(y), U(z));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      type *)((      char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const type *)((const char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x;}
            [[nodiscard]] explicit constexpr operator bool() const requires(!std::is_same_v<type, bool>) {return any();} // Use the explicit methods below for vectors of bool.
            [[nodiscard]] constexpr bool any() const {return x || y || z;}
            [[nodiscard]] constexpr bool all() const {return x && y && z;}
            [[nodiscard]] constexpr bool none() const {return !any();}
            [[nodiscard]] constexpr bool not_all() const {return !all();}
            [[nodiscard]] constexpr auto sum() const {return x + y + z;}
            [[nodiscard]] constexpr auto prod() const {return x * y * z;}
            [[nodiscard]] constexpr type min() const {return std::min({x,y,z});}
            [[nodiscard]] constexpr type max() const {return std::max({x,y,z});}
            [[nodiscard]] constexpr vec abs() const {return vec(std::abs(x), std::abs(y), std::abs(z));}
            template <typename C> [[nodiscard]] constexpr auto index(C &&container) const -> vec<common_vec_size_v<3,vec_size_v<std::decay_t<decltype(container[x])>>>,vec_base_t<std::decay_t<decltype(container[x])>>> {return {Math::vec_elem(0,container[x]), Math::vec_elem(1,container[y]), Math::vec_elem(2,container[z])};}
            [[nodiscard]] constexpr vec2<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr vec4<type> to_vec4(type nw) const {return {x, y, z, nw};}
            [[nodiscard]] constexpr vec4<type> to_vec4() const {return {x, y, z, 0};}
            [[nodiscard]] constexpr auto len_sq() const {return x*x + y*y + z*z;}
            [[nodiscard]] constexpr floating_point_t<type> len() const {return std::sqrt(floating_point_t<type>(len_sq()));}
            [[nodiscard]] constexpr floating_point_t<vec> norm() const {if (auto l = len()) return *this / l; else return vec(0);}
            [[nodiscard]] constexpr auto approx_len() const {return floating_point_t<type>(len_sq() + 1) / 2;} // Accurate only around `len()==1`. Always greater or equal than the true length.
            [[nodiscard]] constexpr auto approx_inv_len() const {return 2 / floating_point_t<type>(len_sq() + 1);}
            [[nodiscard]] constexpr auto approx_norm() const {return *this * approx_inv_len();} // Guaranteed to converge to `len()==1` eventually, when starting from any finite `len_sq()`.
            [[nodiscard]] static constexpr vec axis(int a, type len = 1) {vec ret{}; ret[mod_ex(a,3)] = len; return ret;}
            [[nodiscard]] constexpr vec only_component(int a) {vec ret{}; a = mod_ex(a,3); ret[a] = (*this)[a]; return ret;}
            template <typename U> [[nodiscard]] constexpr auto dot(const vec3<U> &o) const {return x * o.x + y * o.y + z * o.z;}
            template <typename U> [[nodiscard]] constexpr auto cross(const vec3<U> &o) const -> vec3<decltype(x * o.x - x * o.x)> {return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};}
            [[nodiscard]] constexpr auto tie() & {return std::tie(x,y,z);}
            [[nodiscard]] constexpr auto tie() const & {return std::tie(x,y,z);}
            template <int I> [[nodiscard]] constexpr type &get() & {return std::get<I>(tie());}
            template <int I> [[nodiscard]] constexpr const type &get() const & {return std::get<I>(tie());}
            [[nodiscard]] friend constexpr std::array<T,3> format_as(const vec &v) {return {v.x,v.y,v.z};}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_any<vec> operator()(compare_any_tag) const {return compare_any(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_all<vec> operator()(compare_all_tag) const {return compare_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_none<vec> operator()(compare_none_tag) const {return compare_none(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_not_all<vec> operator()(compare_not_all_tag) const {return compare_not_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_elemwise<vec> operator()(compare_elemwise_tag) const {return compare_elemwise(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<type> tiny_rect() const {return rect_to(next_value(*this));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> rect_to(vec3<U> b) const {rect3<larger_t<type,U>> ret; ret.a = *this; ret.b = b; return ret;}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> rect_size(vec3<U> b) const {return rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> rect_size(U b) const {return rect_size(vec3<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> centered_rect_size(vec3<U> b) const {return (*this - b/2).rect_size(b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> centered_rect_size(U b) const {return centered_rect_size(vec3<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> centered_rect_halfsize(vec3<U> b) const {return (*this - b).rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect3<larger_t<type,U>> centered_rect_halfsize(U b) const {return centered_rect_halfsize(vec3<U>(b));}
        };

        template <typename T> struct vec<4,T> // vec4
        {
            using type = T;
            using rect_type = rect4<T>;
            static constexpr int size = 4;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            type x, y, z, w;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &a() {return w;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &a() const {return w;}
            IMP_MATH_SMALL_FUNC constexpr vec() : x{}, y{}, z{}, w{} {}
            IMP_MATH_SMALL_FUNC constexpr vec(uninit) {}
            IMP_MATH_SMALL_FUNC constexpr vec(type x, type y, type z, type w) : x(x), y(y), z(z), w(w) {}
            IMP_MATH_SMALL_FUNC explicit constexpr vec(type obj) : x(obj), y(obj), z(obj), w(obj) {}
            template <scalar U> IMP_MATH_SMALL_FUNC explicit(!safely_convertible_to<U,type>) constexpr vec(vec4<U> obj) : x(obj.x), y(obj.y), z(obj.z), w(obj.w) {}
            template <typename U> requires Custom::convertible<U, vec> explicit(Custom::Convert<U, vec>::is_explicit) constexpr vec(const U &obj) {*this = Custom::Convert<U, vec>{}(obj);}
            template <typename U> requires Custom::convertible<vec, U> explicit(Custom::Convert<vec, U>::is_explicit) operator U() const {return Custom::Convert<vec, U>{}(*this);}
            template <scalar U> [[nodiscard]] constexpr vec4<U> to() const {return vec4<U>(U(x), U(y), U(z), U(w));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      type *)((      char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const type *)((const char *)this + sizeof(type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x;}
            [[nodiscard]] explicit constexpr operator bool() const requires(!std::is_same_v<type, bool>) {return any();} // Use the explicit methods below for vectors of bool.
            [[nodiscard]] constexpr bool any() const {return x || y || z || w;}
            [[nodiscard]] constexpr bool all() const {return x && y && z && w;}
            [[nodiscard]] constexpr bool none() const {return !any();}
            [[nodiscard]] constexpr bool not_all() const {return !all();}
            [[nodiscard]] constexpr auto sum() const {return x + y + z + w;}
            [[nodiscard]] constexpr auto prod() const {return x * y * z * w;}
            [[nodiscard]] constexpr type min() const {return std::min({x,y,z,w});}
            [[nodiscard]] constexpr type max() const {return std::max({x,y,z,w});}
            [[nodiscard]] constexpr vec abs() const {return vec(std::abs(x), std::abs(y), std::abs(z), std::abs(w));}
            template <typename C> [[nodiscard]] constexpr auto index(C &&container) const -> vec<common_vec_size_v<4,vec_size_v<std::decay_t<decltype(container[x])>>>,vec_base_t<std::decay_t<decltype(container[x])>>> {return {Math::vec_elem(0,container[x]), Math::vec_elem(1,container[y]), Math::vec_elem(2,container[z]), Math::vec_elem(3,container[w])};}
            [[nodiscard]] constexpr vec2<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr vec3<type> to_vec3() const {return {x, y, z};}
            [[nodiscard]] constexpr auto len_sq() const {return x*x + y*y + z*z + w*w;}
            [[nodiscard]] constexpr floating_point_t<type> len() const {return std::sqrt(floating_point_t<type>(len_sq()));}
            [[nodiscard]] constexpr floating_point_t<vec> norm() const {if (auto l = len()) return *this / l; else return vec(0);}
            [[nodiscard]] constexpr auto approx_len() const {return floating_point_t<type>(len_sq() + 1) / 2;} // Accurate only around `len()==1`. Always greater or equal than the true length.
            [[nodiscard]] constexpr auto approx_inv_len() const {return 2 / floating_point_t<type>(len_sq() + 1);}
            [[nodiscard]] constexpr auto approx_norm() const {return *this * approx_inv_len();} // Guaranteed to converge to `len()==1` eventually, when starting from any finite `len_sq()`.
            [[nodiscard]] static constexpr vec axis(int a, type len = 1) {vec ret{}; ret[mod_ex(a,4)] = len; return ret;}
            [[nodiscard]] constexpr vec only_component(int a) {vec ret{}; a = mod_ex(a,4); ret[a] = (*this)[a]; return ret;}
            template <typename U> [[nodiscard]] constexpr auto dot(const vec4<U> &o) const {return x * o.x + y * o.y + z * o.z + w * o.w;}
            [[nodiscard]] constexpr auto tie() & {return std::tie(x,y,z,w);}
            [[nodiscard]] constexpr auto tie() const & {return std::tie(x,y,z,w);}
            template <int I> [[nodiscard]] constexpr type &get() & {return std::get<I>(tie());}
            template <int I> [[nodiscard]] constexpr const type &get() const & {return std::get<I>(tie());}
            [[nodiscard]] friend constexpr std::array<T,4> format_as(const vec &v) {return {v.x,v.y,v.z,v.w};}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_any<vec> operator()(compare_any_tag) const {return compare_any(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_all<vec> operator()(compare_all_tag) const {return compare_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_none<vec> operator()(compare_none_tag) const {return compare_none(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_not_all<vec> operator()(compare_not_all_tag) const {return compare_not_all(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr compare_elemwise<vec> operator()(compare_elemwise_tag) const {return compare_elemwise(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<type> tiny_rect() const {return rect_to(next_value(*this));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> rect_to(vec4<U> b) const {rect4<larger_t<type,U>> ret; ret.a = *this; ret.b = b; return ret;}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> rect_size(vec4<U> b) const {return rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> rect_size(U b) const {return rect_size(vec4<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> centered_rect_size(vec4<U> b) const {return (*this - b/2).rect_size(b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> centered_rect_size(U b) const {return centered_rect_size(vec4<U>(b));}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> centered_rect_halfsize(vec4<U> b) const {return (*this - b).rect_to(*this + b);}
            template <scalar U = type> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect4<larger_t<type,U>> centered_rect_halfsize(U b) const {return centered_rect_halfsize(vec4<U>(b));}
        };

        template <typename ...P> requires(sizeof...(P) >= 2 && sizeof...(P) <= 4) vec(P...) -> vec<sizeof...(P), larger_t<P...>>;
        //} Vectors

        //{ Matrices
        template <typename T> struct mat<2,2,T> // mat2x2
        {
            using type = T;
            using member_type = vec2<T>;
            static constexpr int width = 2, height = 2;
            static constexpr int size = 2;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            constexpr mat() : mat(1,0,0,1) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y) : x(x), y(y) {}
            constexpr mat(type xx, type yx, type xy, type yy) : x(xx,xy), y(yx,yy) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat2x2<U> &obj) : x(obj.x), y(obj.y) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y;}
            template <scalar U> [[nodiscard]] constexpr mat2x2<U> to() const {return mat2x2<U>(U(x.x), U(y.x), U(x.y), U(y.y));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,2> format_as(const mat &m) {return {m.x,m.y};}
            [[nodiscard]] constexpr mat3x2<type> to_vec3(const member_type &nz) const {return {x, y, nz};}
            [[nodiscard]] constexpr mat4x2<type> to_vec4(const member_type &nz, const member_type &nw) const {return {x, y, nz, nw};}
            [[nodiscard]] constexpr mat3x2<type> to_vec3() const {return to_vec3({});}
            [[nodiscard]] constexpr mat4x2<type> to_vec4() const {return to_vec4({}, {});}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,0,x.y,y.y,0};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,0,0,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,0,x.y,y.y,0,0,0,1};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,0,0,x.y,y.y,0,0,0,0,1,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,0,0,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,0,x.y,y.y,0,0,0,1,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,0,0,x.y,y.y,0,0,0,0,1,0,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat2x2<T> transpose() const {return {x.x,x.y,y.x,y.y};}
            [[nodiscard]] constexpr mat inverse() requires is_floating_point
            {
                mat ret{};

                ret.x.x =  y.y;
                ret.y.x = -y.x;

                type d = x.x * ret.x.x + x.y * ret.y.x;
                if (d == 0) return {};
                d = 1 / d;
                ret.x.x *= d;
                ret.y.x *= d;

                ret.x.y = -x.y * d;
                ret.y.y =  x.x * d;

                return ret;
            }
            [[nodiscard]] static constexpr mat scale(vec2<type> v)
            {
                return { v.x , 0   ,
                         0   , v.y };
            }
            [[nodiscard]] static constexpr mat rotate(type angle) requires is_floating_point
            {
                type c = std::cos(angle);
                type s = std::sin(angle);
                return { c, -s ,
                         s, c  };
            }
        };

        template <typename T> struct mat<2,3,T> // mat2x3
        {
            using type = T;
            using member_type = vec3<T>;
            static constexpr int width = 2, height = 3;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            constexpr mat() : mat(1,0,0,1,0,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y) : x(x), y(y) {}
            constexpr mat(type xx, type yx, type xy, type yy, type xz, type yz) : x(xx,xy,xz), y(yx,yy,yz) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat2x3<U> &obj) : x(obj.x), y(obj.y) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y;}
            template <scalar U> [[nodiscard]] constexpr mat2x3<U> to() const {return mat2x3<U>(U(x.x), U(y.x), U(x.y), U(y.y), U(x.z), U(y.z));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,2> format_as(const mat &m) {return {m.x,m.y};}
            [[nodiscard]] constexpr mat3x3<type> to_vec3(const member_type &nz) const {return {x, y, nz};}
            [[nodiscard]] constexpr mat4x3<type> to_vec4(const member_type &nz, const member_type &nw) const {return {x, y, nz, nw};}
            [[nodiscard]] constexpr mat3x3<type> to_vec3() const {return to_vec3({});}
            [[nodiscard]] constexpr mat4x3<type> to_vec4() const {return to_vec4({}, {});}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,0,x.y,y.y,0};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,0,0,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,0,x.y,y.y,0,x.z,y.z,1};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,0,0,x.y,y.y,0,0,x.z,y.z,1,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,x.z,y.z,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,0,x.y,y.y,0,x.z,y.z,1,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,0,0,x.y,y.y,0,0,x.z,y.z,1,0,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat3x2<T> transpose() const {return {x.x,x.y,x.z,y.x,y.y,y.z};}
        };

        template <typename T> struct mat<2,4,T> // mat2x4
        {
            using type = T;
            using member_type = vec4<T>;
            static constexpr int width = 2, height = 4;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            constexpr mat() : mat(1,0,0,1,0,0,0,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y) : x(x), y(y) {}
            constexpr mat(type xx, type yx, type xy, type yy, type xz, type yz, type xw, type yw) : x(xx,xy,xz,xw), y(yx,yy,yz,yw) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat2x4<U> &obj) : x(obj.x), y(obj.y) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y;}
            template <scalar U> [[nodiscard]] constexpr mat2x4<U> to() const {return mat2x4<U>(U(x.x), U(y.x), U(x.y), U(y.y), U(x.z), U(y.z), U(x.w), U(y.w));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,2> format_as(const mat &m) {return {m.x,m.y};}
            [[nodiscard]] constexpr mat3x4<type> to_vec3(const member_type &nz) const {return {x, y, nz};}
            [[nodiscard]] constexpr mat4x4<type> to_vec4(const member_type &nz, const member_type &nw) const {return {x, y, nz, nw};}
            [[nodiscard]] constexpr mat3x4<type> to_vec3() const {return to_vec3({});}
            [[nodiscard]] constexpr mat4x4<type> to_vec4() const {return to_vec4({}, {});}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,0,x.y,y.y,0};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,0,0,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,x.z,y.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,0,x.y,y.y,0,x.z,y.z,1};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,0,0,x.y,y.y,0,0,x.z,y.z,1,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,0,x.y,y.y,0,x.z,y.z,1,x.w,y.w,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,0,0,x.y,y.y,0,0,x.z,y.z,1,0,x.w,y.w,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat4x2<T> transpose() const {return {x.x,x.y,x.z,x.w,y.x,y.y,y.z,y.w};}
        };

        template <typename T> struct mat<3,2,T> // mat3x2
        {
            using type = T;
            using member_type = vec2<T>;
            static constexpr int width = 3, height = 2;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            constexpr mat() : mat(1,0,0,0,1,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z) : x(x), y(y), z(z) {}
            constexpr mat(type xx, type yx, type zx, type xy, type yy, type zy) : x(xx,xy), y(yx,yy), z(zx,zy) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat3x2<U> &obj) : x(obj.x), y(obj.y), z(obj.z) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z;}
            template <scalar U> [[nodiscard]] constexpr mat3x2<U> to() const {return mat3x2<U>(U(x.x), U(y.x), U(z.x), U(x.y), U(y.y), U(z.y));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,3> format_as(const mat &m) {return {m.x,m.y,m.z};}
            [[nodiscard]] constexpr mat2x2<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat4x2<type> to_vec4(const member_type &nw) const {return {x, y, z, nw};}
            [[nodiscard]] constexpr mat4x2<type> to_vec4() const {return to_vec4({});}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,z.x,x.y,y.y,z.y,0,0,1};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,0,0,1,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,0,0,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,z.x,x.y,y.y,z.y,0,0,1,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,0,0,1,0,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat2x3<T> transpose() const {return {x.x,x.y,y.x,y.y,z.x,z.y};}
        };

        template <typename T> struct mat<3,3,T> // mat3x3
        {
            using type = T;
            using member_type = vec3<T>;
            static constexpr int width = 3, height = 3;
            static constexpr int size = 3;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            constexpr mat() : mat(1,0,0,0,1,0,0,0,1) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z) : x(x), y(y), z(z) {}
            constexpr mat(type xx, type yx, type zx, type xy, type yy, type zy, type xz, type yz, type zz) : x(xx,xy,xz), y(yx,yy,yz), z(zx,zy,zz) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat3x3<U> &obj) : x(obj.x), y(obj.y), z(obj.z) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z;}
            template <scalar U> [[nodiscard]] constexpr mat3x3<U> to() const {return mat3x3<U>(U(x.x), U(y.x), U(z.x), U(x.y), U(y.y), U(z.y), U(x.z), U(y.z), U(z.z));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,3> format_as(const mat &m) {return {m.x,m.y,m.z};}
            [[nodiscard]] constexpr mat2x3<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat4x3<type> to_vec4(const member_type &nw) const {return {x, y, z, nw};}
            [[nodiscard]] constexpr mat4x3<type> to_vec4() const {return to_vec4({});}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,z.x,x.y,y.y,z.y};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,x.z,y.z};}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,x.z,y.z,z.z,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,x.z,y.z,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,x.z,y.z,z.z,0,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat3x3<T> transpose() const {return {x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z};}
            [[nodiscard]] constexpr mat inverse() const requires is_floating_point
            {
                mat ret{};

                ret.x.x =  y.y * z.z - z.y * y.z;
                ret.y.x = -y.x * z.z + z.x * y.z;
                ret.z.x =  y.x * z.y - z.x * y.y;

                type d = x.x * ret.x.x + x.y * ret.y.x + x.z * ret.z.x;
                if (d == 0) return {};
                d = 1 / d;
                ret.x.x *= d;
                ret.y.x *= d;
                ret.z.x *= d;

                ret.x.y = (-x.y * z.z + z.y * x.z) * d;
                ret.y.y = ( x.x * z.z - z.x * x.z) * d;
                ret.z.y = (-x.x * z.y + z.x * x.y) * d;
                ret.x.z = ( x.y * y.z - y.y * x.z) * d;
                ret.y.z = (-x.x * y.z + y.x * x.z) * d;
                ret.z.z = ( x.x * y.y - y.x * x.y) * d;

                return ret;
            }
            [[nodiscard]] static constexpr mat scale(vec2<type> v) {return mat2<T>::scale(v).to_mat3();}
            [[nodiscard]] static constexpr mat scale(vec3<type> v)
            {
                return { v.x , 0   , 0   ,
                         0   , v.y , 0   ,
                         0   , 0   , v.z };
            }
            [[nodiscard]] static constexpr mat ortho(vec2<type> min, vec2<type> max) requires is_floating_point
            {
                return { 2 / (max.x - min.x) , 0                   , (min.x + max.x) / (min.x - max.x) ,
                         0                   , 2 / (max.y - min.y) , (min.y + max.y) / (min.y - max.y) ,
                         0                   , 0                   , 1                                 };
            }
            [[nodiscard]] static constexpr mat translate(vec2<type> v)
            {
                return { 1, 0, v.x ,
                         0, 1, v.y ,
                         0, 0, 1   };
            }
            [[nodiscard]] static constexpr mat rotate(type angle) {return mat2<T>::rotate(angle).to_mat3();}
            [[nodiscard]] static constexpr mat rotate_with_normalized_axis(vec3<type> axis, type angle)
            {
                type c = std::cos(angle);
                type s = std::sin(angle);
                return { axis.x * axis.x * (1 - c) + c          , axis.x * axis.y * (1 - c) - axis.z * s , axis.x * axis.z * (1 - c) + axis.y * s,
                         axis.y * axis.x * (1 - c) + axis.z * s , axis.y * axis.y * (1 - c) + c          , axis.y * axis.z * (1 - c) - axis.x * s,
                         axis.x * axis.z * (1 - c) - axis.y * s , axis.y * axis.z * (1 - c) + axis.x * s , axis.z * axis.z * (1 - c) + c         };
            }
            [[nodiscard]] static constexpr mat rotate(vec3<type> axis, type angle) requires is_floating_point
            {
                return rotate_with_normalized_axis(axis.norm(), angle);
            }
        };

        template <typename T> struct mat<3,4,T> // mat3x4
        {
            using type = T;
            using member_type = vec4<T>;
            static constexpr int width = 3, height = 4;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            constexpr mat() : mat(1,0,0,0,1,0,0,0,1,0,0,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z) : x(x), y(y), z(z) {}
            constexpr mat(type xx, type yx, type zx, type xy, type yy, type zy, type xz, type yz, type zz, type xw, type yw, type zw) : x(xx,xy,xz,xw), y(yx,yy,yz,yw), z(zx,zy,zz,zw) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat3x4<U> &obj) : x(obj.x), y(obj.y), z(obj.z) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z;}
            template <scalar U> [[nodiscard]] constexpr mat3x4<U> to() const {return mat3x4<U>(U(x.x), U(y.x), U(z.x), U(x.y), U(y.y), U(z.y), U(x.z), U(y.z), U(z.z), U(x.w), U(y.w), U(z.w));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,3> format_as(const mat &m) {return {m.x,m.y,m.z};}
            [[nodiscard]] constexpr mat2x4<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat4x4<type> to_vec4(const member_type &nw) const {return {x, y, z, nw};}
            [[nodiscard]] constexpr mat4x4<type> to_vec4() const {return to_vec4({});}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,z.x,x.y,y.y,z.y};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,x.z,y.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,x.z,y.z,z.z,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,x.z,y.z,x.w,y.w};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,z.x,0,x.y,y.y,z.y,0,x.z,y.z,z.z,0,x.w,y.w,z.w,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat4x3<T> transpose() const {return {x.x,x.y,x.z,x.w,y.x,y.y,y.z,y.w,z.x,z.y,z.z,z.w};}
        };

        template <typename T> struct mat<4,2,T> // mat4x2
        {
            using type = T;
            using member_type = vec2<T>;
            static constexpr int width = 4, height = 2;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z, w;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &a() {return w;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &a() const {return w;}
            constexpr mat() : mat(1,0,0,0,0,1,0,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}), w(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z, const member_type &w) : x(x), y(y), z(z), w(w) {}
            constexpr mat(type xx, type yx, type zx, type wx, type xy, type yy, type zy, type wy) : x(xx,xy), y(yx,yy), z(zx,zy), w(wx,wy) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat4x2<U> &obj) : x(obj.x), y(obj.y), z(obj.z), w(obj.w) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z && a.w==b.w;}
            template <scalar U> [[nodiscard]] constexpr mat4x2<U> to() const {return mat4x2<U>(U(x.x), U(y.x), U(z.x), U(w.x), U(x.y), U(y.y), U(z.y), U(w.y));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,4> format_as(const mat &m) {return {m.x,m.y,m.z,m.w};}
            [[nodiscard]] constexpr mat2x2<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat3x2<type> to_vec3() const {return {x, y, z};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,z.x,x.y,y.y,z.y};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,0,0};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,z.x,x.y,y.y,z.y,0,0,1};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y,0,0,1,0};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,0,0,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,z.x,x.y,y.y,z.y,0,0,1,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y,0,0,1,0,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat2x4<T> transpose() const {return {x.x,x.y,y.x,y.y,z.x,z.y,w.x,w.y};}
        };

        template <typename T> struct mat<4,3,T> // mat4x3
        {
            using type = T;
            using member_type = vec3<T>;
            static constexpr int width = 4, height = 3;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z, w;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &a() {return w;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &a() const {return w;}
            constexpr mat() : mat(1,0,0,0,0,1,0,0,0,0,1,0) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}), w(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z, const member_type &w) : x(x), y(y), z(z), w(w) {}
            constexpr mat(type xx, type yx, type zx, type wx, type xy, type yy, type zy, type wy, type xz, type yz, type zz, type wz) : x(xx,xy,xz), y(yx,yy,yz), z(zx,zy,zz), w(wx,wy,wz) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat4x3<U> &obj) : x(obj.x), y(obj.y), z(obj.z), w(obj.w) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z && a.w==b.w;}
            template <scalar U> [[nodiscard]] constexpr mat4x3<U> to() const {return mat4x3<U>(U(x.x), U(y.x), U(z.x), U(w.x), U(x.y), U(y.y), U(z.y), U(w.y), U(x.z), U(y.z), U(z.z), U(w.z));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,4> format_as(const mat &m) {return {m.x,m.y,m.z,m.w};}
            [[nodiscard]] constexpr mat2x3<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat3x3<type> to_vec3() const {return {x, y, z};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,z.x,x.y,y.y,z.y};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,x.z,y.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,x.z,y.z,0,0};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z,0,0,0};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4x4() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y,x.z,y.z,z.z,w.z,0,0,0,1};}
            [[nodiscard]] constexpr mat4x4<type> to_mat4() const {return to_mat4x4();}
            [[nodiscard]] constexpr mat3x4<T> transpose() const {return {x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z,w.x,w.y,w.z};}
        };

        template <typename T> struct mat<4,4,T> // mat4x4
        {
            using type = T;
            using member_type = vec4<T>;
            static constexpr int width = 4, height = 4;
            static constexpr int size = 4;
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            member_type x, y, z, w;
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &r() {return x;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &r() const {return x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &g() {return y;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &g() const {return y;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &b() {return z;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &b() const {return z;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr type &a() {return w;} [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const type &a() const {return w;}
            constexpr mat() : mat(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1) {}
            constexpr mat(uninit) : x(uninit{}), y(uninit{}), z(uninit{}), w(uninit{}) {}
            constexpr mat(const member_type &x, const member_type &y, const member_type &z, const member_type &w) : x(x), y(y), z(z), w(w) {}
            constexpr mat(type xx, type yx, type zx, type wx, type xy, type yy, type zy, type wy, type xz, type yz, type zz, type wz, type xw, type yw, type zw, type ww) : x(xx,xy,xz,xw), y(yx,yy,yz,yw), z(zx,zy,zz,zw), w(wx,wy,wz,ww) {}
            template <scalar U> explicit(!safely_convertible_to<U,T>) constexpr mat(const mat4x4<U> &obj) : x(obj.x), y(obj.y), z(obj.z), w(obj.w) {}
            template <typename U> requires Custom::convertible<U, mat> explicit(Custom::Convert<U, mat>::is_explicit) constexpr mat(const U &obj) {*this = Custom::Convert<U, mat>{}(obj);}
            template <typename U> requires Custom::convertible<mat, U> explicit(Custom::Convert<mat, U>::is_explicit) operator U() const {return Custom::Convert<mat, U>{}(*this);}
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(const mat &a, const mat &b) {return a.x==b.x && a.y==b.y && a.z==b.z && a.w==b.w;}
            template <scalar U> [[nodiscard]] constexpr mat4x4<U> to() const {return mat4x4<U>(U(x.x), U(y.x), U(z.x), U(w.x), U(x.y), U(y.y), U(z.y), U(w.y), U(x.z), U(y.z), U(z.z), U(w.z), U(x.w), U(y.w), U(z.w), U(w.w));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr       member_type &operator[](int i)       {if (!IMP_MATH_IS_CONSTANT(i)) return *(      member_type *)((      char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr const member_type &operator[](int i) const {if (!IMP_MATH_IS_CONSTANT(i)) return *(const member_type *)((const char *)this + sizeof(member_type)*i); else if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else if (i == 3) return w; IMP_MATH_UNREACHABLE();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC type *as_array() {return &x.x;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC const type *as_array() const {return &x.x;}
            [[nodiscard]] friend constexpr std::array<T,4> format_as(const mat &m) {return {m.x,m.y,m.z,m.w};}
            [[nodiscard]] constexpr mat2x4<type> to_vec2() const {return {x, y};}
            [[nodiscard]] constexpr mat3x4<type> to_vec3() const {return {x, y, z};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2x2() const {return {x.x,y.x,x.y,y.y};}
            [[nodiscard]] constexpr mat2x2<type> to_mat2() const {return to_mat2x2();}
            [[nodiscard]] constexpr mat3x2<type> to_mat3x2() const {return {x.x,y.x,z.x,x.y,y.y,z.y};}
            [[nodiscard]] constexpr mat4x2<type> to_mat4x2() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y};}
            [[nodiscard]] constexpr mat2x3<type> to_mat2x3() const {return {x.x,y.x,x.y,y.y,x.z,y.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3x3() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z};}
            [[nodiscard]] constexpr mat3x3<type> to_mat3() const {return to_mat3x3();}
            [[nodiscard]] constexpr mat4x3<type> to_mat4x3() const {return {x.x,y.x,z.x,w.x,x.y,y.y,z.y,w.y,x.z,y.z,z.z,w.z};}
            [[nodiscard]] constexpr mat2x4<type> to_mat2x4() const {return {x.x,y.x,x.y,y.y,x.z,y.z,x.w,y.w};}
            [[nodiscard]] constexpr mat3x4<type> to_mat3x4() const {return {x.x,y.x,z.x,x.y,y.y,z.y,x.z,y.z,z.z,x.w,y.w,z.w};}
            [[nodiscard]] constexpr mat4x4<T> transpose() const {return {x.x,x.y,x.z,x.w,y.x,y.y,y.z,y.w,z.x,z.y,z.z,z.w,w.x,w.y,w.z,w.w};}
            [[nodiscard]] constexpr mat inverse() const requires is_floating_point
            {
                mat ret;

                ret.x.x =  y.y * z.z * w.w - y.y * z.w * w.z - z.y * y.z * w.w + z.y * y.w * w.z + w.y * y.z * z.w - w.y * y.w * z.z;
                ret.y.x = -y.x * z.z * w.w + y.x * z.w * w.z + z.x * y.z * w.w - z.x * y.w * w.z - w.x * y.z * z.w + w.x * y.w * z.z;
                ret.z.x =  y.x * z.y * w.w - y.x * z.w * w.y - z.x * y.y * w.w + z.x * y.w * w.y + w.x * y.y * z.w - w.x * y.w * z.y;
                ret.w.x = -y.x * z.y * w.z + y.x * z.z * w.y + z.x * y.y * w.z - z.x * y.z * w.y - w.x * y.y * z.z + w.x * y.z * z.y;

                type d = x.x * ret.x.x + x.y * ret.y.x + x.z * ret.z.x + x.w * ret.w.x;
                if (d == 0) return {};
                d = 1 / d;
                ret.x.x *= d;
                ret.y.x *= d;
                ret.z.x *= d;
                ret.w.x *= d;

                ret.x.y = (-x.y * z.z * w.w + x.y * z.w * w.z + z.y * x.z * w.w - z.y * x.w * w.z - w.y * x.z * z.w + w.y * x.w * z.z) * d;
                ret.y.y = ( x.x * z.z * w.w - x.x * z.w * w.z - z.x * x.z * w.w + z.x * x.w * w.z + w.x * x.z * z.w - w.x * x.w * z.z) * d;
                ret.z.y = (-x.x * z.y * w.w + x.x * z.w * w.y + z.x * x.y * w.w - z.x * x.w * w.y - w.x * x.y * z.w + w.x * x.w * z.y) * d;
                ret.w.y = ( x.x * z.y * w.z - x.x * z.z * w.y - z.x * x.y * w.z + z.x * x.z * w.y + w.x * x.y * z.z - w.x * x.z * z.y) * d;
                ret.x.z = ( x.y * y.z * w.w - x.y * y.w * w.z - y.y * x.z * w.w + y.y * x.w * w.z + w.y * x.z * y.w - w.y * x.w * y.z) * d;
                ret.y.z = (-x.x * y.z * w.w + x.x * y.w * w.z + y.x * x.z * w.w - y.x * x.w * w.z - w.x * x.z * y.w + w.x * x.w * y.z) * d;
                ret.z.z = ( x.x * y.y * w.w - x.x * y.w * w.y - y.x * x.y * w.w + y.x * x.w * w.y + w.x * x.y * y.w - w.x * x.w * y.y) * d;
                ret.w.z = (-x.x * y.y * w.z + x.x * y.z * w.y + y.x * x.y * w.z - y.x * x.z * w.y - w.x * x.y * y.z + w.x * x.z * y.y) * d;
                ret.x.w = (-x.y * y.z * z.w + x.y * y.w * z.z + y.y * x.z * z.w - y.y * x.w * z.z - z.y * x.z * y.w + z.y * x.w * y.z) * d;
                ret.y.w = ( x.x * y.z * z.w - x.x * y.w * z.z - y.x * x.z * z.w + y.x * x.w * z.z + z.x * x.z * y.w - z.x * x.w * y.z) * d;
                ret.z.w = (-x.x * y.y * z.w + x.x * y.w * z.y + y.x * x.y * z.w - y.x * x.w * z.y - z.x * x.y * y.w + z.x * x.w * y.y) * d;
                ret.w.w = ( x.x * y.y * z.z - x.x * y.z * z.y - y.x * x.y * z.z + y.x * x.z * z.y + z.x * x.y * y.z - z.x * x.z * y.y) * d;

                return ret;
            }
            [[nodiscard]] static constexpr mat scale(vec3<type> v) {return mat3<T>::scale(v).to_mat4();}
            [[nodiscard]] static constexpr mat ortho(vec2<type> min, vec2<type> max, type near, type far) requires is_floating_point
            {
                return { 2 / (max.x - min.x) , 0                   , 0                , (min.x + max.x) / (min.x - max.x) ,
                         0                   , 2 / (max.y - min.y) , 0                , (min.y + max.y) / (min.y - max.y) ,
                         0                   , 0                   , 2 / (near - far) , (near + far) / (near - far)       ,
                         0                   , 0                   , 0                , 1                                 };
            }
            [[nodiscard]] static constexpr mat look_at(vec3<type> src, vec3<type> dst, vec3<type> local_up) requires is_floating_point
            {
                vec3<type> v3 = (src-dst).norm();
                vec3<type> v1 = local_up.cross(v3).norm();
                vec3<type> v2 = v3.cross(v1);
                return { v1.x , v1.y , v1.z , -src.x*v1.x-src.y*v1.y-src.z*v1.z ,
                         v2.x , v2.y , v2.z , -src.x*v2.x-src.y*v2.y-src.z*v2.z ,
                         v3.x , v3.y , v3.z , -src.x*v3.x-src.y*v3.y-src.z*v3.z ,
                         0    , 0    , 0    , 1                                 };
            }
            [[nodiscard]] static constexpr mat translate(vec3<type> v)
            {
                return { 1 , 0 , 0 , v.x ,
                         0 , 1 , 0 , v.y ,
                         0 , 0 , 1 , v.z ,
                         0 , 0 , 0 , 1   };
            }
            [[nodiscard]] static constexpr mat rotate_with_normalized_axis(vec3<type> axis, type angle) {return mat3<T>::rotate_with_normalized_axis(axis, angle).to_mat4();}
            [[nodiscard]] static constexpr mat rotate(vec3<type> axis, type angle) {return mat3<T>::rotate(axis, angle).to_mat4();}
            [[nodiscard]] static constexpr mat perspective(type wh_aspect, type y_fov, type near, type far) requires is_floating_point
            {
                y_fov = type(1) / std::tan(y_fov / 2);
                return { y_fov / wh_aspect , 0     , 0                           , 0                             ,
                         0                 , y_fov , 0                           , 0                             ,
                         0                 , 0     , (near + far) / (near - far) , 2 * near * far / (near - far) ,
                         0                 , 0     , -1                          , 0                             };
            }
        };

        template <scalar ...P> requires(sizeof...(P) == 4) mat(P...) -> mat<2, 2, larger_t<P...>>;
        template <scalar ...P> requires(sizeof...(P) == 9) mat(P...) -> mat<3, 3, larger_t<P...>>;
        template <scalar ...P> requires(sizeof...(P) == 16) mat(P...) -> mat<4, 4, larger_t<P...>>;
        template <typename ...P> requires(sizeof...(P) >= 2 && sizeof...(P) <= 4 && ((vec_size_v<P> == 2) && ...)) mat(P...) -> mat<sizeof...(P), 2, larger_t<typename P::type...>>;
        template <typename ...P> requires(sizeof...(P) >= 2 && sizeof...(P) <= 4 && ((vec_size_v<P> == 3) && ...)) mat(P...) -> mat<sizeof...(P), 3, larger_t<typename P::type...>>;
        template <typename ...P> requires(sizeof...(P) >= 2 && sizeof...(P) <= 4 && ((vec_size_v<P> == 4) && ...)) mat(P...) -> mat<sizeof...(P), 4, larger_t<typename P::type...>>;
        //} Matrices

        //{ Rects
        template <int D, scalar T> struct rect
        {
            using type = T;
            using vec_type = vec<D,T>;
            static constexpr int dim = D; // `size` is already used as a function name.
            static constexpr bool is_floating_point = floating_point_scalar<type>;
            vec_type a, b; // `a` is inclusive, `b` is exclusive.
            IMP_MATH_SMALL_FUNC constexpr rect() {} // No fancy constructors, use helpers in `vec`.
            IMP_MATH_SMALL_FUNC constexpr rect(uninit) : a(uninit{}), b(uninit{}) {}
            template <scalar U> IMP_MATH_SMALL_FUNC explicit(!safely_convertible_to<U,T>) constexpr rect(rect<D,U> r) : a(r.a), b(r.b) {}
            template <scalar U> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,U> to() const {return vec<D,U>(a).rect_to(vec<D,U>(b));}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec_type size() const {return b - a;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec_type center() const {return a + size() / 2;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool has_length() const {return (b > a).any();}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool has_area() const {return (b > a).sum() >= 2;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool has_volume() const requires(D >= 3) {return (b > a).sum() >= 3;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool has_4d_volume() const requires(D >= 4) {return (b > a).sum() >= 4;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect fix() const {rect ret = *this; sort_two_var(ret.a, ret.b); return ret;} // Swap components of `a` and `b` to order them correctly.
            // Offsetting.
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset_a(vec<D,U> x) const {rect<D,larger_t<T,U>> ret = *this; ret.a += x; return ret;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset_a(U        x) const {rect<D,larger_t<T,U>> ret = *this; ret.a += x; return ret;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset_b(vec<D,U> x) const {rect<D,larger_t<T,U>> ret = *this; ret.b += x; return ret;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset_b(U        x) const {rect<D,larger_t<T,U>> ret = *this; ret.b += x; return ret;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset  (vec<D,U> x) const {return offset_a(x).offset_b(x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> offset  (U        x) const {return offset_a(x).offset_b(x);}
            // Operators. Those apply to both corners. Don't forget to `.fix()` the rect if you negate it.
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect operator+() const {return *this;}
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect operator-() const {return (-a).rect_to(-b);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+(rect r, vec<D,U> x) {return (r.a + x).rect_to(r.b + x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+(rect r, U        x) {return (r.a + x).rect_to(r.b + x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+(vec<D,U> x, rect r) {return (x + r.a).rect_to(x + r.b);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+(U        x, rect r) {return (x + r.a).rect_to(x + r.b);}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+=(rect &r, vec<D,U> x) {r = r + x; return r;}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator+=(rect &r, U        x) {r = r + x; return r;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-(rect r, vec<D,U> x) {return (r.a - x).rect_to(r.b - x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-(rect r, U        x) {return (r.a - x).rect_to(r.b - x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-(vec<D,U> x, rect r) {return (x - r.a).rect_to(x - r.b);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-(U        x, rect r) {return (x - r.a).rect_to(x - r.b);}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-=(rect &r, vec<D,U> x) {r = r - x; return r;}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator-=(rect &r, U        x) {r = r - x; return r;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*(rect r, vec<D,U> x) {return (r.a * x).rect_to(r.b * x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*(rect r, U        x) {return (r.a * x).rect_to(r.b * x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*(vec<D,U> x, rect r) {return (x * r.a).rect_to(x * r.b);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*(U        x, rect r) {return (x * r.a).rect_to(x * r.b);}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*=(rect &r, vec<D,U> x) {r = r * x; return r;}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator*=(rect &r, U        x) {r = r * x; return r;}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/(rect r, vec<D,U> x) {return (r.a / x).rect_to(r.b / x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/(rect r, U        x) {return (r.a / x).rect_to(r.b / x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/(vec<D,U> x, rect r) {return (x / r.a).rect_to(x / r.b);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/(U        x, rect r) {return (x / r.a).rect_to(x / r.b);}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/=(rect &r, vec<D,U> x) {r = r / x; return r;}
            template <safely_convertible_to<T> U = T> IMP_MATH_SMALL_FUNC friend constexpr rect<D,larger_t<T,U>> operator/=(rect &r, U        x) {r = r / x; return r;}
            // Expanding and shrinking.
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> expand(vec<D,U> x) const {return offset_a(-x).offset_b(x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> expand(U        x) const {return offset_a(-x).offset_b(x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> shrink(vec<D,U> x) const {return offset_a(x).offset_b(-x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> shrink(U        x) const {return offset_a(x).offset_b(-x);}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> expand_dir(vec<D,U> x) const {return offset_a(min(x,larger_t<T,U>{})).offset_b(max(x,larger_t<T,U>{}));}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> expand_dir(U        x) const {return expand_dir(vec<D,U>(x));}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> shrink_dir(vec<D,U> x) const {return offset_a(max(x,larger_t<T,U>{})).offset_b(min(x,larger_t<T,U>{}));}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> shrink_dir(U        x) const {return shrink_dir(vec<D,U>(x));}
            // Checking collisions.
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool contains(vec<D,U> p) const {return (p >= a).all() && (p </*sic*/ b).all();}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool contains(rect<D,U> r) const {return (r.a >= a).all() && (r.b <= b).all();}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr bool touches(rect r) const {return (r.a < b).all() && (r.b > a).all();}
            // Modifying the rect.
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> combine(vec<D,U>  p) const {return combine(p.tiny_rect());}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> combine(rect<D,U> r) const {return min(a, r.a).rect_to(max(b, r.b));}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr rect<D,larger_t<T,U>> intersect(rect<D,U> r) const {return max(a, r.a).rect_to(min(b, r.b));}
            template <scalar U = T> [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec<D,larger_t<T,U>> clamp(vec<D,U> p) const {return min(max(p, a), prev_value(b)); }
            // Constructing the contour.
            [[nodiscard]] IMP_MATH_SMALL_FUNC constexpr vec_type corner(int i) const requires(dim==2) {return vec_type((i+1)&2?b.x:a.x, i&2?b.y:a.y);}
            [[nodiscard]] constexpr std::array<vec_type, 4> to_contour() const requires(dim==2) {std::array<vec_type, 4> ret; for (int i=0;i<4;i++) ret[i]=corner(i); return ret;}
            // Comparisons.
            [[nodiscard]] IMP_MATH_SMALL_FUNC friend constexpr bool operator==(rect x, rect y) {return x.a == y.a && x.b == y.b;}
        };

        // input/output
        template <typename A, typename B, int D, typename T> constexpr std::basic_ostream<A,B> &operator<<(std::basic_ostream<A, B> &s, const rect<D, T> &r)
        {
            return s << r.a << ".." << r.b;
        }
        template <typename A, typename B, int D, typename T> constexpr std::basic_istream<A,B> &operator>>(std::basic_istream<A, B> &s, rect<D, T> &r)
        {
            return s >> r.a >> r.b;
        }
        //} Rects
    }

    namespace Export
    {
        using Vector::vec; // Vector and matrix definitions. We use this instead of `using namespace Vector` to avoid bringing...
        using Vector::mat; // ...the overloaded operators into the global namespace, mostly for better error messages and build speed.
        using namespace Alias; // Convenient type aliases.
        using namespace Common; // Common functions.

        // Common types.
        using std::int8_t;
        using std::uint8_t;
        using std::int16_t;
        using std::uint16_t;
        using std::int32_t;
        using std::uint32_t;
        using std::int64_t;
        using std::uint64_t;
        using std::size_t;
        using std::ptrdiff_t;
        using std::intptr_t;
        using std::uintptr_t;

        // Common standard functions.
        using std::sqrt;
        using std::cos;
        using std::sin;
        using std::tan;
        using std::acos;
        using std::asin;
        using std::atan;
        using std::atan2;
    }
}

namespace std
{
    template <int D, typename T> struct less         <Math::vec<D,T>> {constexpr bool operator()(const Math::vec<D,T> &a, const Math::vec<D,T> &b) const {return a.tie() <  b.tie();}};
    template <int D, typename T> struct greater      <Math::vec<D,T>> {constexpr bool operator()(const Math::vec<D,T> &a, const Math::vec<D,T> &b) const {return a.tie() >  b.tie();}};
    template <int D, typename T> struct less_equal   <Math::vec<D,T>> {constexpr bool operator()(const Math::vec<D,T> &a, const Math::vec<D,T> &b) const {return a.tie() <= b.tie();}};
    template <int D, typename T> struct greater_equal<Math::vec<D,T>> {constexpr bool operator()(const Math::vec<D,T> &a, const Math::vec<D,T> &b) const {return a.tie() >= b.tie();}};

    template <int D, typename T> struct hash<Math::vec<D,T>>
    {
        std::size_t operator()(const Math::vec<D,T> &v) const
        {
            std::size_t ret = std::hash<decltype(v.x)>{}(v.x);
            for (int i = 1; i < D; i++)
                ret ^= std::hash<decltype(v.x)>{}(v[i]) + std::size_t(0x9E3779B97F4A7C16) + (ret << 6) + (ret >> 2); // Boost uses something similar.
            return ret;
        }
    };
}

// Quaternions

namespace Math
{
    inline namespace Quat // Quaternions.
    {
        template <floating_point_scalar T> struct quat
        {
            using type = T;
            using vec3_t = vec3<T>;
            using vec4_t = vec4<T>;
            using mat3_t = mat3<T>;
            type x = 0, y = 0, z = 0, w = 1; // This represents zero rotation.

            constexpr quat() {}
            constexpr quat(type x, type y, type z, type w) : x(x), y(y), z(z), w(w) {}
            explicit constexpr quat(const vec4_t &vec) : x(vec.x), y(vec.y), z(vec.z), w(vec.w) {}

            // Normalizes the axis. If it's already normalized, use `with_normalized_axis()` instead.
            constexpr quat(vec3_t axis, type angle) {*this = with_normalized_axis(axis.norm(), angle);}
            [[nodiscard]] static constexpr quat with_normalized_axis(vec3_t axis, type angle) {angle *= type(0.5); return quat((axis * std::sin(angle)).to_vec4(std::cos(angle)));}

            [[nodiscard]] constexpr vec4_t as_vec() const {return {x, y, z, w};}
            [[nodiscard]] constexpr vec3_t xyz() const {return {x, y, z};}
            [[nodiscard]] type *as_array() {return &x;}
            [[nodiscard]] const type *as_array() const {return &x;}

            [[nodiscard]] constexpr quat norm() const {return quat(as_vec().norm());}
            [[nodiscard]] constexpr quat approx_norm() const {return quat(as_vec().approx_norm());}

            [[nodiscard]] constexpr vec3_t axis_denorm() const { return xyz(); }
            [[nodiscard]] constexpr vec3_t axis_norm() const { return xyz().norm(); }
            [[nodiscard]] constexpr float angle() const { return 2 * std::atan2(xyz().len(), w); }

            // Negates the rotation. Not strictly an inversion in the mathematical sense, since the length stays unchanged (while it's supposed to become `1 / old_length`).
            [[nodiscard]] constexpr quat inverse() const {return quat(xyz().to_vec4(-w));}
            // Negates the three imaginary parts of the quaternion, `xyz`. Effectively inverts the rotation, but works slower than `inverse()`. Useful only for low-level quaternion things.
            [[nodiscard]] constexpr quat conjugate() const {return quat((-xyz()).to_vec4(w));}

            // Uses iterative normalization to keep denormalization from accumulating due to lack of precision.
            template <typename U> [[nodiscard]] constexpr quat<larger_t<T,U>> operator*(const quat<U> &other) const {return mult_without_norm(other).approx_norm();}
            constexpr quat &operator*=(const quat &other) {return *this = *this * other;}

            // Simple quaternion multiplication, without any normalization.
            template <typename U> [[nodiscard]] constexpr quat<larger_t<T,U>> mult_without_norm(const quat<U> &other) const
            {
                return quat<larger_t<T,U>>(vec4<larger_t<T,U>>(
                    x * other.w + w * other.x - z * other.y + y * other.z,
                    y * other.w + z * other.x + w * other.y - x * other.z,
                    z * other.w - y * other.x + x * other.y + w * other.z,
                    w * other.w - x * other.x - y * other.y - z * other.z
                ));
            }

            // Transforms a vector by this quaternion. Only makes sense if the quaternion is normalized.
            template <typename U> [[nodiscard]] constexpr vec3<larger_t<T,U>> operator*(const vec3<U> &other) const
            {
                // This is called the "Euler-Rodrigues formula".
                // We could also use `*this * other * this->conjugate()`, but that looks less optimized.
                vec3<larger_t<T,U>> tmp = xyz().cross(other);
                return other + 2 * w * tmp + 2 * xyz().cross(tmp);
            }

            // Transforms a vector by this quaternion, inversed. Mimics a similar matrix operation.
            template <typename U> [[nodiscard]] friend constexpr vec3<larger_t<T,U>> operator*(const vec3<U> &v, const quat &q)
            {
                return q.inverse() * v;
            }

            // Returns a rotation matrix for this quaternion. Only makes sense if the quaternion is normalized.
            [[nodiscard]] constexpr mat3_t matrix() const
            {
                return mat3_t(
                    1 - (2*y*y + 2*z*z), 2*x*y - 2*z*w, 2*x*z + 2*y*w,
                    2*x*y + 2*z*w, 1 - (2*x*x + 2*z*z), 2*y*z - 2*x*w,
                    2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - (2*x*x + 2*y*y)
                );
            }

            // Returns a rotation matrix for this quaternion. Works even if the quaternion is not normalized.
            [[nodiscard]] constexpr mat3_t matrix_from_denorm() const
            {
                type f = 1 / as_vec().len_sq();
                mat3_t m = matrix();
                return mat3_t(m.x * f, m.y * f, m.z * f);
            }
        };

        using fquat = quat<float>;
        using dquat = quat<double>;
        using ldquat = quat<long double>;

        template <typename A, typename B, typename T> constexpr std::basic_ostream<A,B> &operator<<(std::basic_ostream<A,B> &s, const quat<T> &q)
        {
            s.width(0);
            if (q.axis_denorm() == vec3<T>(0))
                s << "[angle=0";
            else
                s << "[axis=" << q.axis_denorm()/q.axis_denorm().max() << " angle=" << to_deg(q.angle()) << "(deg)";
            return s << " len=" << q.as_vec().len() << ']';
        }

        template <typename A, typename B, typename T> constexpr std::basic_istream<A,B> &operator>>(std::basic_istream<A,B> &s, quat<T> &q)
        {
            vec4<T> vec;
            s >> vec;
            q = quat(vec);
            return s;
        }
    }

    inline namespace Utility
    {
        // Check if `T` is a quaternion type (possibly const).
        template <typename T> struct is_quat_impl : std::false_type {};
        template <typename T> struct is_quat_impl<      quat<T>> : std::true_type {};
        template <typename T> struct is_quat_impl<const quat<T>> : std::true_type {};
        template <typename T> inline constexpr bool is_quat_v = is_quat_impl<T>::value;
    }

    namespace Export
    {
        using namespace Quat;
    }
}

namespace std
{
    template <typename T> struct less         <Math::quat<T>> {constexpr bool operator()(const Math::quat<T> &a, const Math::quat<T> &b) const {return a.as_vec().tie() <  b.as_vec().tie();}};
    template <typename T> struct greater      <Math::quat<T>> {constexpr bool operator()(const Math::quat<T> &a, const Math::quat<T> &b) const {return a.as_vec().tie() >  b.as_vec().tie();}};
    template <typename T> struct less_equal   <Math::quat<T>> {constexpr bool operator()(const Math::quat<T> &a, const Math::quat<T> &b) const {return a.as_vec().tie() <= b.as_vec().tie();}};
    template <typename T> struct greater_equal<Math::quat<T>> {constexpr bool operator()(const Math::quat<T> &a, const Math::quat<T> &b) const {return a.as_vec().tie() >= b.as_vec().tie();}};

    template <typename T> struct hash<Math::quat<T>>
    {
        std::size_t operator()(const Math::quat<T> &q) const
        {
            return std::hash<Math::vec4<T>>{}(q.as_vec());
        }
    };
}

using namespace Math::Export;
