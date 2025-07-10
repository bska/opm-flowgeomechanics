/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_GEOMECH_DUNE_UTILITIES_HPP_INCLUDED
#define OPM_GEOMECH_DUNE_UTILITIES_HPP_INCLUDED

#include <dune/grid/common/datahandleif.hh>

#include <algorithm>
#include <cstddef>
#include <set>
#include <vector>

namespace Dune
{
template <class GridView, class Vector>
class EntityVectorVectorDataHandle
    : public Dune::CommDataHandleIF<EntityVectorVectorDataHandle<GridView, Vector>,
                                    typename Vector::value_type>
{
public:
    /// \brief the data type we send
    using DataType = typename Vector::value_type;

    /// \brief Constructor
    /// \param data The vector of data vectors
    /// \param gridView The gridview the data is attached to.
    EntityVectorVectorDataHandle(Vector& data, const GridView& gridView, int codim)
        : data_(data)
        , gridView_(gridView)
        , codim_(codim)
    {
    }

    bool contains(int /* dim */, int codim) const
    {
        return codim == codim_;
    }

#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
    bool fixedsize(int /* dim */, int /* codim */) const
    {
        return true;
    }
#else

    bool fixedSize(int /* dim */, int /* codim */) const
    {
        return true;
    }
#endif

    template <class EntityType>
    std::size_t size(const EntityType /* entity */) const
    {
        return 1;
    }

    template <class BufferType, class EntityType>
    void gather(BufferType& buffer, const EntityType& e) const
    {

        buffer.write(data_[gridView_.indexSet().index(e)]);
    }

    template <class BufferType, class EntityType>
    void scatter(BufferType& buffer, const EntityType& e, [[maybe_unused]] std::size_t n)
    {
        assert(n == 1);
        buffer.read(data_[gridView_.indexSet().index(e)]);
    }

private:
    Vector& data_;
    const GridView& gridView_;
    int codim_;
};

template <class GridView, class Vector>
class MaxEntityVectorVectorDataHandle
    : public Dune::CommDataHandleIF<MaxEntityVectorVectorDataHandle<GridView, Vector>,
                                    typename Vector::value_type>
{
public:
    /// \brief the data type we send
    using DataType = typename Vector::value_type;

    /// \brief Constructor
    /// \param data The vector of data vectors
    /// \param gridView The gridview the data is attached to.
    MaxEntityVectorVectorDataHandle(Vector& data, const GridView& gridView, int codim)
        : data_(data)
        , gridView_(gridView)
        , codim_(codim)
    {
        all_data_.resize(gridView.size(codim));
    }

    bool contains(int /* dim */, int codim) const
    {
        return codim == codim_;
    }

#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
    bool fixedsize(int /* dim */, int /* codim */) const
    {
        return true;
    }
#else
    bool fixedSize(int /* dim */, int /* codim */) const
    {
        return true;
    }
#endif

    template <class EntityType>
    std::size_t size(const EntityType /* entity */) const
    {
        return 1;
    }

    template <class BufferType, class EntityType>
    void gather(BufferType& buffer, const EntityType& e) const
    {

        buffer.write(data_[gridView_.indexSet().index(e)]);
    }

    template <class BufferType, class EntityType>
    void scatter(BufferType& buffer, const EntityType& e, [[maybe_unused]] std::size_t n)
    {
        assert(n == 1);

        DataType value;
        value = 5;
        buffer.read(value);

        const auto ix = this->gridView_.indexSet().index(e);

        data_[ix] = std::max(data_[ix], value);
        all_data_[ix].push_back(value);
    }

    const std::vector<std::vector<int>>& all_data() const
    {
        return all_data_;
    }

private:
    Vector& data_;
    std::vector<std::vector<int>> all_data_;
    const GridView& gridView_;
    int codim_;
};

template <typename T>
struct MaxGatherScatter
{
    using V = typename CommPolicy<T>::IndexedType;

    static V gather(const T& a, std::size_t i)
    {
        return a[i];
    }

    static void scatter(T& a, V v, std::size_t i)
    {
        a[i] = std::max(a[i], v);
    }
};

template <>
struct CommPolicy<std::vector<std::vector<int>>>
{
    using Type = std::vector<std::vector<int>>;
    using IndexedType = int;
    using IndexedTypeFlag = VariableSize;

    static const void* getAddress(const Type& v, const int index)
    {
        return &(v[index][0]);
    }

    static std::size_t getSize(const Type& v, const int index)
    {
        return v[index].size();
    }
};

struct VariableVectorAdderGatherScatter
{
    using Container = std::vector<std::vector<int>>;
    using IndexedType = typename CommPolicy<Container>::IndexedType;

    static IndexedType gather(const Container& cont, const std::size_t i, const std::size_t j)
    {
        return cont[i][j];
    }

    static void
    scatter(Container& cont, const IndexedType& data, const std::size_t i, const std::size_t j)
    {
        std::set<int> tmp(cont[i].begin(), cont[i].end());

        const auto tmp_size = tmp.size();
        tmp.insert(data);

        if ((tmp_size != tmp.size()) && (tmp.size() > cont[i].size())) {
            std::cout << " size error" << std::endl;
            std::cout << "tmp" << std::endl;

            for (const auto& t : tmp) {
                std::cout << t << " ";
            }

            std::cout << std::endl;
            std::cout << "cont" << std::endl;
            for (const auto& t : cont[i]) {
                std::cout << t << " ";
            }

            assert(false);
        }

        std::copy(tmp.begin(), tmp.end(), cont[i].begin());
    }
};

/** \brief gather/scatter callback for communcation */
template <typename T>
struct FreeCopyGatherScatter
{
    using V = typename CommPolicy<T>::IndexedType;

    static V gather(const T& a, std::size_t i)
    {
        return a[i];
    }

    static void scatter(T& a, V v, std::size_t i)
    {
        a[i] = v;
    }
};

template <typename T>
struct FreeAddGatherScatter
{
    using V = typename CommPolicy<T>::IndexedType;

    static V gather(const T& a, const std::size_t i)
    {
        return a[i];
    }

    static void scatter(T& a, V v, const std::size_t i)
    {
        a[i] += v;
    }
};

} // namespace Dune

#endif // OPM_GEOMECH_DUNE_UTILITIES_HPP_INCLUDED
