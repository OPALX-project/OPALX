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

template<typename Left, typename Right>
class BiMap {
public:
    using left_map = std::map<Left, Right>;
    using right_map = std::map<Right, Left>;
    
private:
    // Declare maps first so they're initialized before left/right references
    left_map left_map_;
    right_map right_map_;
    
public:
    struct relation {
        Left left;
        Right right;
        relation(const Left& l, const Right& r) : left(l), right(r) {}
    };
    
    // Access to underlying maps for compatibility (using references)
    struct left_view {
        left_map& map_;
        left_view(left_map& m) : map_(m) {}
        typename left_map::iterator find(const Left& key) const { return map_.find(key); }
        typename left_map::iterator end() const { return map_.end(); }
        const Right& at(const Left& key) const { return map_.at(key); }
    };
    
    struct right_view {
        right_map& map_;
        right_view(right_map& m) : map_(m) {}
        typename right_map::iterator find(const Right& key) const { return map_.find(key); }
        typename right_map::iterator end() const { return map_.end(); }
        const Left& at(const Right& key) const { return map_.at(key); }
    };
    
    BiMap() : left_map_(), right_map_(), left(left_map_), right(right_map_) {}
    
    left_view left;
    right_view right;
    
    void insert(const Left& left, const Right& right) {
        left_map_[left] = right;
        right_map_[right] = left;
    }
    
    typename left_map::iterator left_find(const Left& key) {
        return left_map_.find(key);
    }
    
    typename right_map::iterator right_find(const Right& key) {
        return right_map_.find(key);
    }
    
    typename left_map::iterator left_end() { return left_map_.end(); }
    typename right_map::iterator right_end() { return right_map_.end(); }
    
    const Right& left_at(const Left& key) const {
        auto it = left_map_.find(key);
        if (it == left_map_.end()) {
            throw std::out_of_range("Key not found in left map");
        }
        return it->second;
    }
    
    const Left& right_at(const Right& key) const {
        auto it = right_map_.find(key);
        if (it == right_map_.end()) {
            throw std::out_of_range("Key not found in right map");
        }
        return it->second;
    }
};

// Helper function to create a BiMap from initializer list
template<typename Left, typename Right>
BiMap<Left, Right> make_bimap(std::initializer_list<typename BiMap<Left, Right>::relation> relations) {
    BiMap<Left, Right> bimap;
    for (const auto& rel : relations) {
        bimap.insert(rel.left, rel.right);
    }
    return bimap;
}

#endif

