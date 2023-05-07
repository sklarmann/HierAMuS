#pragma once

#include <utility>

namespace HierAMuS {
template <class master, typename indx, typename data>
class rangedVectorIterator {
public:
  rangedVectorIterator(master *ml, indx pl) : main(ml), pos(pl){};

  bool operator==(const rangedVectorIterator &second) const {
    return (second.pos == this->pos);
  };
  bool operator!=(const rangedVectorIterator &second) const {
    return (second.pos != this->pos);
  };

  rangedVectorIterator &operator++() {
    this->main->nextElement(this->pos);
    return (*this);
  };

  data operator->(){return (*this->main)[this->pos];};

  std::pair<indx &, data &> operator*() {

    return std::pair<indx &, data &>(this->pos, (*this->main)[this->pos]);
  };

private:
  master *main;
  indx pos;
};
} // namespace HierAMuS