#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

#include "stream/input.h"

namespace Refl::Utils
{
    // Skips whitespace and c++-style comments (`//` and `/* */`).
    // Returns true if skipped at least one character.
    // Throws if there is an unterminated `/*` comment.
    inline bool SkipWhitespaceAndComments(Stream::Input &input)
    {
        bool first = true;

        while (1)
        {
            if (!input.MoreData())
                break;

            // Skip any whitespace.
            bool have_whitespace = input.Discard<Stream::any>(Stream::Char::IsWhitespace{});

            bool have_comment = false;

            // Check for a slash.
            if (input.MoreData() && input.PeekChar() == '/')
            {
                input.SkipOne();
                char ch = input.MoreData() ? input.PeekChar() : '\0';
                if (ch != '*' && ch != '/')
                {
                    // This is not a comment, unget the slash.
                    input.Seek(-1, Stream::relative);
                }
                else
                {
                    // This is a comment, skip it.
                    have_comment = 1;
                    input.SkipOne();
                    if (ch == '/') // One-line comment.
                    {
                        while (input.MoreData())
                        {
                            char ch = input.PeekChar();
                            if (ch == '\r' || ch == '\n')
                            {
                                input.SkipOne();
                                break;
                            }
                            input.SkipOne();
                        }
                    }
                    else // ch == '*', multi-line comment.
                    {
                        auto pos = input.Position();

                        char prev = 0;
                        while (1)
                        {
                            if (!input.MoreData())
                            {
                                input.Seek(pos-2, Stream::absolute); // Move cursor to the beginning of the comment, for better error reporting.
                                throw std::runtime_error(input.GetExceptionPrefix() + "Unterminated comment.");
                            }
                            char ch = input.PeekChar();
                            if (ch == '/' && prev == '*')
                            {
                                input.SkipOne();
                                break;
                            }
                            prev = ch;
                            input.SkipOne();
                        }
                    }
                }
            }

            if (!have_whitespace && !have_comment)
                break;

            first = false;
        }

        return !first;
    }


    // Constexpr alternative to `std::strcmp`.
    constexpr int const_strcmp(const char *a, const char *b)
    {
        while (*a == *b)
        {
            if (*a == '\0')
                return 0;
            a++;
            b++;
        }
        return (unsigned char)*a - (unsigned char)*b;
    }

    struct NameIndexPair
    {
        const char *name = nullptr;
        std::size_t index = 0;

        constexpr bool operator==(const NameIndexPair &other) const
        {
            return *this <=> other == 0;
        }

        constexpr std::weak_ordering operator<=>(const NameIndexPair &other) const
        {
            return *this <=> other.name;
        }

        constexpr std::weak_ordering operator<=>(const char *other_name) const
        {
            return const_strcmp(name, other_name) <=> 0;
        }
    };

    // An universal function to look up strings in immutable lists.
    // `F` is a pointer to a constexpr function that returns an array of names: `std::array<const char *, N> (*)(auto index)`.
    // `name` is a name that we're looking for. If it's not found, -1 is returned.
    // Avoid using lambdas as `F`. If you do that in a header, you will most likely get an ODR violation.
    template <auto F> std::size_t GetStringIndex(const char *name)
    {
        static constexpr auto array = []
        {
            auto name_array = F();
            std::array<NameIndexPair, name_array.size()> array{};
            for (std::size_t i = 0; i < array.size(); i++)
            {
                array[i].name = name_array[i];
                array[i].index = i;
            }
            std::sort(array.begin(), array.end());
            return array;
        }();
        static_assert(std::adjacent_find(array.begin(), array.end()) == array.end(), "Duplicate string in a static list.");
        auto it = std::lower_bound(array.begin(), array.end(), name);
        if (it == array.end() || *it != NameIndexPair{name, 0})
            return -1;
        return it->index;
    }
}
