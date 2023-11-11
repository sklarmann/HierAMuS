#pragma once

#include <algorithm>

template <typename C, typename T> bool contains(C &&c, T e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};

template<size_t N> 
constexpr size_t minArraySize() {
  if constexpr (N < 0) {
    return 0;
  } else {
    return N;
  }
}

template <size_t N, typename T, T value>
constexpr std::array<T, N> initArray() {
  constexpr size_t Nt = minArraySize<N>();
  std::array<T, Nt> a;
  for (auto i = 0; i < Nt; ++i) {
    a[i] = value;
  }
  return a;
};

template <class U, size_t N, typename T, T value>
constexpr U initArray2() {
  using temp = typename U::value_type;
  if constexpr (std::is_same<U, std::vector<temp>>::value) {
    return U();
  } else {
    constexpr size_t Nt = minArraySize<N>();
    std::array<T, Nt> a;
    for (auto i = 0; i < Nt; ++i) {
      a[i] = value;
    }
    return a;
  }
};
