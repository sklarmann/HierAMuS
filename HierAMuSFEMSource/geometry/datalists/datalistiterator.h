// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

template<class T, typename indexType, typename data>
class datalistiterator {
public:
  datalistiterator(T *obj, indexType ind) : main(obj), currentIndex(ind) {};

  auto operator!=(const datalistiterator &second) -> bool {return this->currentIndex!=second.currentIndex;};
  auto operator==(const datalistiterator &second) -> bool {return this->currentIndex==second.currentIndex;};

  auto operator++() -> datalistiterator & {++this->currentIndex; return(*this);};

  auto operator*() -> data{return this->main->getItemLocalNumber(this->currentIndex);};

  auto operator->() -> data{return this->main->getItemLocalNumber(this->currentIndex);};

private:
  T *main;
  indexType currentIndex;

};

