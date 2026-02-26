//
// Simple bidirectional map to replace boost::bimap
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_BIMAP_HH
#define OPAL_BIMAP_HH

#include <map>
#include <stdexcept>

/// \brief Simple bidirectional map with lookup in both directions.
/// \tparam Left Key type for the left map.
/// \tparam Right Key type for the right map.
template <typename Left, typename Right>
class BiMap {
public:
    using left_map  = std::map<Left, Right>;
    using right_map = std::map<Right, Left>;

private:
    // Declare maps first so they're initialized before left/right references
    left_map left_map_;
    right_map right_map_;

public:
    /// \brief Pair relation used for initialization.
    struct relation {
        Left left;
        Right right;
        /// \brief Construct a relation from left/right values.
        /// \param l Input: left value.
        /// \param r Input: right value.
        relation(const Left& l, const Right& r) : left(l), right(r) {}
    };

    // Access to underlying maps for compatibility (using references)
    /// \brief Left-map view with find/at helpers.
    struct left_view {
        left_map& map_;
        /// \brief Construct a left view from a map reference.
        /// \param m Input: left map reference.
        left_view(left_map& m) : map_(m) {}
        /// \brief Find iterator for \p key in the left map.
        /// \param key Input: left key.
        /// \return Output: iterator to entry or \c end().
        typename left_map::iterator find(const Left& key) const { return map_.find(key); }
        /// \brief End iterator for the left map.
        /// \return Output: end iterator.
        typename left_map::iterator end() const { return map_.end(); }
        /// \brief Get mapped value by key (throws if missing).
        /// \param key Input: left key.
        /// \return Output: mapped right value.
        const Right& at(const Left& key) const { return map_.at(key); }
    };

    /// \brief Right-map view with find/at helpers.
    struct right_view {
        right_map& map_;
        /// \brief Construct a right view from a map reference.
        /// \param m Input: right map reference.
        right_view(right_map& m) : map_(m) {}
        /// \brief Find iterator for \p key in the right map.
        /// \param key Input: right key.
        /// \return Output: iterator to entry or \c end().
        typename right_map::iterator find(const Right& key) const { return map_.find(key); }
        /// \brief End iterator for the right map.
        /// \return Output: end iterator.
        typename right_map::iterator end() const { return map_.end(); }
        /// \brief Get mapped value by key (throws if missing).
        /// \param key Input: right key.
        /// \return Output: mapped left value.
        const Left& at(const Right& key) const { return map_.at(key); }
    };

    /// \brief Construct an empty bimap.
    BiMap() : left_map_(), right_map_(), left(left_map_), right(right_map_) {}

    /// \brief Left view accessor.
    left_view left;
    /// \brief Right view accessor.
    right_view right;

    /// \brief Insert or overwrite a left/right association.
    /// \param left Input: left key.
    /// \param right Input: right key.
    void insert(const Left& left, const Right& right) {
        left_map_[left]   = right;
        right_map_[right] = left;
    }

    /// \brief Find an entry by left key.
    /// \param key Input: left key.
    /// \return Output: iterator to entry or \c end().
    typename left_map::iterator left_find(const Left& key) { return left_map_.find(key); }

    /// \brief Find an entry by right key.
    /// \param key Input: right key.
    /// \return Output: iterator to entry or \c end().
    typename right_map::iterator right_find(const Right& key) { return right_map_.find(key); }

    /// \brief End iterator for left map.
    /// \return Output: end iterator.
    typename left_map::iterator left_end() { return left_map_.end(); }
    /// \brief End iterator for right map.
    /// \return Output: end iterator.
    typename right_map::iterator right_end() { return right_map_.end(); }

    /// \brief Get mapped right value by left key (throws if missing).
    /// \param key Input: left key.
    /// \return Output: mapped right value.
    const Right& left_at(const Left& key) const {
        auto it = left_map_.find(key);
        if (it == left_map_.end()) {
            throw std::out_of_range("Key not found in left map");
        }
        return it->second;
    }

    /// \brief Get mapped left value by right key (throws if missing).
    /// \param key Input: right key.
    /// \return Output: mapped left value.
    const Left& right_at(const Right& key) const {
        auto it = right_map_.find(key);
        if (it == right_map_.end()) {
            throw std::out_of_range("Key not found in right map");
        }
        return it->second;
    }
};

/// \brief Helper function to create a BiMap from an initializer list.
/// \param relations Input: list of left/right relations.
/// \return Output: initialized bimap instance.
template <typename Left, typename Right>
BiMap<Left, Right> make_bimap(
    std::initializer_list<typename BiMap<Left, Right>::relation> relations) {
    BiMap<Left, Right> bimap;
    for (const auto& rel : relations) {
        bimap.insert(rel.left, rel.right);
    }
    return bimap;
}

#endif
