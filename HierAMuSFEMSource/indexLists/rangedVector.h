
#pragma once

#include <iostream>
#include <type_traits>
#include <vector>

#include <indexLists/rangedVectorIterator.h>

namespace HierAMuS {
template <typename indx, class returnval> class rangedVectorBase {
  typedef rangedVectorIterator<rangedVectorBase<indx, returnval>, indx,
                               returnval>
      iterator;
  typedef rangedVectorIterator<rangedVectorBase<indx, returnval>, indx,
                               const returnval>
      const_iterator;

public:
  returnval &operator[](const indx &idx) {
    return data[idx - this->startindex];
  };
  indx getLastElementsIndex() {
    indx retidx = this->startindex + static_cast<indx>(this->data.size()) - 1;
    // this->data.size() > 0 ? retidx = this->startindex + this->data.size() - 1
    //                       : retidx = 0;
    return retidx;
  };

  indx getNumberOfElements() { return this->numberOfElements; };

  bool inList(const indx &idx);

  void nextElement(indx &elem) {
    indx pos = elem - this->startindex;
    ++pos;
    ++elem;
    if (pos < this->data.size()) {
      bool search = !this->dataInit[pos];
      while (search) {
        ++pos;
        ++elem;
        search = !this->dataInit[pos];
      }
    }
  };

  iterator begin() { return iterator(this, this->startindex); };
  iterator end() {
    return iterator(this, this->startindex + static_cast<indx>(this->data.size()));
  };

protected:
  bool inListRange(const indx &idx);

  /* data */
  indx startindex = 0, numberOfElements = 0;
  std::vector<returnval> data;
  std::vector<bool> dataInit;
};

template <typename indx, class returnval>
bool rangedVectorBase<indx, returnval>::inListRange(const indx &idx) {
  if (this->data.size() == 0) {
    return false;
  } else if (idx < this->startindex ||
             idx >= this->startindex + this->data.size()) {
    return false;
  }
  return true;
}
template <typename indx, class returnval>
bool rangedVectorBase<indx, returnval>::inList(const indx &idx) {
  if (this->data.size() == 0) {
    return false;
  } else if (idx < this->startindex ||
             idx >= this->startindex + this->data.size()) {
    return false;
  } else {
    indx pos = idx - this->startindex;
    return this->dataInit[pos];
  }
}

template <typename indx, class returnval>
class rangedVector : public rangedVectorBase<indx, returnval> {

public:
  rangedVector(/* args */);
  ~rangedVector();

  void insertAtIndex(const indx &idx, const returnval &toAdd);

protected:
};

template <typename indx, class returnval>
rangedVector<indx, returnval>::rangedVector(/* args */) = default;

template <typename indx, class returnval>
rangedVector<indx, returnval>::~rangedVector() = default;

template <typename indx, class returnval>
void rangedVector<indx, returnval>::insertAtIndex(const indx &idx,
                                                  const returnval &toAdd) {
  if (this->inListRange(idx)) {
    indx pos = idx - this->startindex;
    this->data[pos] = toAdd;
    if(!this->dataInit[pos]) ++this->numberOfElements;
    this->dataInit[pos] = true;

  } else if (this->data.size() == 0) {
    this->startindex = idx;
    this->data.push_back(toAdd);
    this->dataInit.push_back(true);
    ++this->numberOfElements;
  } else {
    if (idx < this->startindex) {
      indx numToAdd = this->startindex - idx;
      this->data.insert(this->data.begin(), numToAdd, toAdd);
      this->dataInit.insert(this->dataInit.begin(), numToAdd, false);
      this->startindex = idx;
      this->dataInit[0] = true;
    } else {
      indx numToAdd =
          idx - this->startindex - static_cast<indx>(this->data.size()) + 1;
      this->data.insert(this->data.end(), numToAdd, toAdd);
      this->dataInit.insert(this->dataInit.end(), numToAdd, false);
      this->dataInit[idx - this->startindex] = true;
    }
    ++this->numberOfElements;
  }
}

/*
 *  Pointer Implementation
 */

template <typename indx, class returnval>
class rangedVector<indx, returnval *>
    : public rangedVectorBase<indx, returnval *> {
public:
  ~rangedVector() {
    for (auto it = this->data.begin(); it != this->data.end(); ++it) {
      delete *it;
    }
  };

  void insertAtIndex(const indx &idx, returnval *toAdd);
};

template <typename indx, class returnval>
void rangedVector<indx, returnval *>::insertAtIndex(const indx &idx,
                                                    returnval *toAdd) {
  if (this->inListRange(idx)) {
    indx pos = idx - this->startindex;
    if (this->dataInit[pos]){
        delete this->data[pos];
    } else {
        ++this->numberOfElements;
    }
      
    this->data[pos] = toAdd;
    this->dataInit[pos] = true;
  } else if (this->data.size() == 0) {
    this->startindex = idx;
    this->data.push_back(toAdd);
    this->dataInit.push_back(true);
    ++this->numberOfElements;
  } else {
    if (idx < this->startindex) {
      indx numToAdd = this->startindex - idx;
      this->data.insert(this->data.begin(), numToAdd, NULL);
      this->dataInit.insert(this->dataInit.begin(), numToAdd, false);
      this->startindex = idx;
      this->dataInit[0] = true;
      this->data[0] = toAdd;
    } else {
      indx numToAdd =
          idx - this->startindex - static_cast<indx>(this->data.size()) + 1;
      this->data.insert(this->data.end(), numToAdd, NULL);
      this->dataInit.insert(this->dataInit.end(), numToAdd, false);
      this->dataInit[idx - this->startindex] = true;
      this->data[idx - this->startindex] = toAdd;
    }
    ++this->numberOfElements;
  }
}

} // namespace HierAMuS