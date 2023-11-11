// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <unordered_map>
#include <iterator>
#include <memory>

namespace HierAMuS::Geometry {

template<class T>
class GeometryList {

public:
  GeometryList() : m_number_of_objects(0), m_max_object_id(0) {};
  virtual ~GeometryList(){};

  template<class D>
  void add_element(indexType id) {
    if (m_object_map.find(id) == m_object_map.end()) {
      m_object_map[id] = std::make_shared<D>();
      m_object_map[id]->set_id(id);
      if (id > m_max_object_id)
        m_max_object_id = id;
      ++m_number_of_objects;
    }
  };
  template<class D>
  void add_element(indexType id, D &&element) {
    if (m_object_map.find(id) == m_object_map.end()) {
      element.setId(id);
      m_object_map[id].make_shared(element);
      if (id > m_max_object_id)
        m_max_object_id = id;
      ++m_number_of_objects;
    }
  }

  auto get_element(indexType id) -> std::shared_ptr<T> { return m_object_map.at(id); };
  auto get_element_reference(indexType id) -> T & { return (*m_object_map.at(id)); };
  auto get_element_pointer(indexType id) -> T * {
    return &(*m_object_map.at(id));
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

  auto has_element(indexType id) -> bool {
    if (m_object_map.find(id) != m_object_map.end()) {
      return true;
    }
    return false;
  }

private:
  indexType m_max_object_id;
  indexType m_number_of_objects;
  std::map<indexType, std::shared_ptr<T>> m_object_map;
};

} // namespace HierAMuS::Geometry