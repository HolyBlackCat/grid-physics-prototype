#pragma once

#include <algorithm>
#include <cstddef>

#include "meta/common.h"

// This file offers compile-time strings that can be used as template parameters.
// Example 1:
//     template <Meta::ConstString Name> void foo() {std::cout << Name.str << '\n';}
//     foo<"123">();
// Example 2:
//     template <Meta::ConstString Name> void foo(Meta::ConstStringTag<Name>) {std::cout << Name.str << '\n';}
//     foo("123"_const);

namespace Meta
{
    namespace impl
    {
        // Does nothing, but causes an error if called from a `consteval` function.
        inline void ExpectedNullTerminatedArray() {}
    }

    // A string that can be used as a template parameter.
    template <std::size_t N>
    struct ConstString
    {
        char str[N]{};

        static constexpr std::size_t size = N - 1;

        consteval ConstString() {}
        consteval ConstString(const char (&new_str)[N])
        {
            if (new_str[N-1] != '\0')
                impl::ExpectedNullTerminatedArray();
            std::copy_n(new_str, size, str);
        }

        [[nodiscard]] constexpr std::string_view view() const &
        {
            return {str, str + size};
        }
        [[nodiscard]] constexpr std::string_view view() const && = delete;
    };

    template <std::size_t A, std::size_t B>
    [[nodiscard]] constexpr ConstString<A + B - 1> operator+(const ConstString<A> &a, const ConstString<B> &b)
    {
        ConstString<A + B - 1> ret;
        std::copy_n(a.str, a.size, ret.str);
        std::copy_n(b.str, b.size, ret.str + a.size);
        return ret;
    }

    template <std::size_t A, std::size_t B>
    [[nodiscard]] constexpr ConstString<A + B - 1> operator+(const ConstString<A> &a, const char (&b)[B])
    {
        return a + ConstString<B>(b);
    }

    template <std::size_t A, std::size_t B>
    [[nodiscard]] constexpr ConstString<A + B - 1> operator+(const char (&a)[A], const ConstString<B> &b)
    {
        return ConstString<A>(a) + b;
    }


    // A tag structure returned by `operator""_const` below.
    template <Meta::ConstString S>
    struct ConstStringTag
    {
        static constexpr Meta::ConstString value = S;
    };

    // Returns a string encoded into a template parameter of a tag structure `ConstStringTag`.
    template <Meta::ConstString S>
    [[nodiscard]] constexpr ConstStringTag<S> operator""_const()
    {
        return {};
    }
}

using Meta::operator""_const;
