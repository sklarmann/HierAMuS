// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <unordered_map>
#include <iterator>
#include <memory>

namespace HierAMuS::Geometry {

template<class T>
class GeometryListSingle {

public:
  GeometryListSingle() : m_number_of_objects(0), m_max_object_id(0) {};
  virtual ~GeometryListSingle(){};

  void add_element(indexType id) {
    m_object_map[id].set_id(id);
    if (id > m_max_object_id)
      m_max_object_id = id;
    ++m_number_of_objects;
  };

  void add_element(indexType id, T &&element) {
    if (m_object_map.find(id) == m_object_map.end()) {
      m_object_map.emplace(id, element);
    }
  }

  auto get_element_reference(indexType id) -> T & { 
    return m_object_map.at(id);
  };

  auto get_element_pointer(indexType id) -> T * {
    return &this->get_element_reference(id);
  };

  auto lastElement() -> indexType { return m_max_object_id; };
  auto getNumberOfElements() -> indexType { return m_number_of_objects; };

  // Iterator functions
  auto begin() -> auto {
    return m_object_map.begin();
  };
  auto end() -> auto {
    return m_object_map.end();
  };

  void clear() {
    m_max_object_id = 0;
    m_number_of_objects = 0;
    m_object_map.clear();
  };

private:
  indexType m_max_object_id;
  indexType m_number_of_objects;
  std::unordered_map<indexType, T> m_object_map;
};

} // namespace HierAMuS::Geometry