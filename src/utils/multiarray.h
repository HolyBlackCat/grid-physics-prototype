#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <utility>

#include "meta/lists.h"
#include "program/errors.h"
#include "strings/format.h"
#include "utils/mat.h"


template <int D, typename T, std::signed_integral Index = std::ptrdiff_t>
class MultiArray
{
  public:
    static constexpr int dimensions = D;
    static_assert(dimensions >= 2, "Arrays with less than 2 dimensions are not supported.");
    static_assert(dimensions <= 4, "Arrays with more than 4 dimensions are not supported.");

    using type = T;
    using index_t = Index;
    using index_vec_t = vec<D, index_t>;

    struct ReflHelper; // Our reflection metadata uses this to access private fields.

  private:
    index_vec_t size_vec{};
    std::vector<type> storage;

  public:
    constexpr MultiArray() {}

    MultiArray(index_vec_t size_vec) : size_vec(size_vec), storage(size_vec.prod())
    {
        ASSERT(size_vec.min() >= 0, "Invalid multiarray size.");
        if (size_vec(any) <= 0)
            size_vec = {};
    }
    MultiArray(index_vec_t size_vec, const T &init) : size_vec(size_vec), storage(size_vec.prod(), init)
    {
        ASSERT(size_vec.min() >= 0, "Invalid multiarray size.");
        if (size_vec(any) <= 0)
            size_vec = {};
    }
    template <typename A, A ...I>
    MultiArray(Meta::value_list<I...>, const std::array<type, index_vec_t(I...).prod()> &data) : size_vec(I...), storage(data.begin(), data.end())
    {
        static_assert(std::is_integral_v<A>, "Indices must be integral.");
        static_assert(((I >= 0) && ...), "Invalid multiarray size.");
        if (size_vec(any) <= 0)
            size_vec = {};
    }

    [[nodiscard]] index_vec_t size() const
    {
        return size_vec;
    }

    [[nodiscard]] typename index_vec_t::rect_type bounds() const
    {
        return index_vec_t{}.rect_size(size_vec);
    }

    [[nodiscard]] bool pos_in_range(index_vec_t pos) const
    {
        return bounds().contains(pos);
    }

    [[nodiscard]] type &unsafe_at(index_vec_t pos)
    {
        ASSERT(pos_in_range(pos), STR("Multiarray indices out of range. Indices are ", (pos), " but the array size is ", (size_vec), "."));

        index_t index = 0;
        index_t factor = 1;

        for (int i = 0; i < dimensions; i++)
        {
            index += factor * pos[i];
            factor *= size_vec[i];
        }

        return storage[index];
    }
    [[nodiscard]] type &safe_throwing_at(index_vec_t pos)
    {
        if (!pos_in_range(pos))
            throw std::runtime_error(FMT("Multiarray index {} is out of range. The array size is {}.", pos, size_vec));
        return unsafe_at(pos);
    }
    [[nodiscard]] type &safe_nonthrowing_at(index_vec_t pos)
    {
        if (!pos_in_range(pos))
            Program::HardError(FMT("Multiarray index {} is out of range. The array size is {}.", pos, size_vec));
        return unsafe_at(pos);
    }
    [[nodiscard]] type &clamped_at(index_vec_t pos)
    {
        clamp_var(pos, 0, size_vec-1);
        return unsafe_at(pos);
    }
    [[nodiscard]] type try_get(index_vec_t pos)
    {
        if (!pos_in_range(pos))
            return {};
        return unsafe_at(pos);
    }
    void try_set(index_vec_t pos, const type &obj)
    {
        if (!pos_in_range(pos))
            return;
        unsafe_at(pos) = obj;
    }
    void try_set(index_vec_t pos, type &&obj)
    {
        if (!pos_in_range(pos))
            return;
        unsafe_at(pos) = std::move(obj);
    }

    [[nodiscard]] const type &unsafe_at(index_vec_t pos) const
    {
        return const_cast<MultiArray *>(this)->unsafe_at(pos);
    }
    [[nodiscard]] const type &safe_throwing_at(index_vec_t pos) const
    {
        return const_cast<MultiArray *>(this)->safe_throwing_at(pos);
    }
    [[nodiscard]] const type &safe_nonthrowing_at(index_vec_t pos) const
    {
        return const_cast<MultiArray *>(this)->safe_nonthrowing_at(pos);
    }
    [[nodiscard]] const type &clamped_at(index_vec_t pos) const
    {
        return const_cast<MultiArray *>(this)->clamped_at(pos);
    }
    [[nodiscard]] type try_get(index_vec_t pos) const
    {
        return const_cast<MultiArray *>(this)->try_get(pos);
    }
    void try_set(index_vec_t pos, const type &obj) const
    {
        return const_cast<MultiArray *>(this)->try_set(pos, obj);
    }
    void try_set(index_vec_t pos, type &&obj) const
    {
        return const_cast<MultiArray *>(this)->try_set(pos, std::move(obj));
    }

    [[nodiscard]] index_t element_count() const
    {
        return storage.size();
    }
    [[nodiscard]] type *elements()
    {
        return storage.data();
    }
    [[nodiscard]] const type *elements() const
    {
        return storage.data();
    }

    // Resizes the array and/or offsets it by the specified amount.
    // Any out-of-range elements are destroyed.
    void resize(index_vec_t new_size, index_vec_t offset = {})
    {
        if (new_size == size_vec && offset == 0)
            return;
        *this = resize_copy(new_size, offset);
    }

    // Resizes a copy of the array and/or offsets it by the specified amount.
    // Any out-of-range elements are destroyed.
    [[nodiscard]] MultiArray resize_copy(index_vec_t new_size, index_vec_t offset = {}) const
    {
        if (new_size(any) == 0)
            return {}; // Target is empty, stop early.

        if (new_size == size_vec && offset == 0)
            return *this; // No changes are needed.

        MultiArray ret(new_size);
        new_size = ret.size(); // This sanitizes the size.

        if (size_vec(any) == 0)
            return ret; // Source is empty, return zeroed array.

        index_vec_t source_start = clamp_min(-offset, 0);
        index_vec_t source_end = clamp_max(new_size - offset, size_vec);

        for (index_vec_t pos : source_start <= vector_range < source_end)
            ret.safe_nonthrowing_at(pos + offset) = std::move(safe_nonthrowing_at(pos));

        return ret;
    }
};

template <typename T, typename Index = std::ptrdiff_t> using Array2D = MultiArray<2, T, Index>;
template <typename T, typename Index = std::ptrdiff_t> using Array3D = MultiArray<3, T, Index>;
template <typename T, typename Index = std::ptrdiff_t> using Array4D = MultiArray<4, T, Index>;
