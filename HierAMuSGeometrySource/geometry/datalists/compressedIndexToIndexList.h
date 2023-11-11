// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <vector>
#include "datatypes.h"
namespace HierAMuS {



class compressedIndexToIndexList {
public:
    compressedIndexToIndexList() : m_initialized(false){};
    ~compressedIndexToIndexList(){};


    /** @brief Initializes the compressedIndexToIndexList with the given start and end index.
     *
     *  Sets the start index and the end index of the compressedIndexToIndexList. After this, the number of entries at each index needs to be counted.
     *
     *  @param [in] startIndex The start index of the compressedIndexToIndexList.
     *  @param [in] endIndex The end index of the compressedIndexToIndexList.
     */
    void initialize(indexType startIndex, indexType endIndex);

    /** @brief Increases the number of entries at the given index by one.
     *
     *  Counts the number of entries at the given index. The index must be between the start and end index of the compressedIndexToIndexList.
     *
     *  @param [in] index The index at which the number of entries is increased by one.
     */
    void countAtIndex(indexType index);

    /** @brief Finalizes the setup after all entries are counted.
     *
     *  After all entries are counted, this function needs to be called to finalize the setup of the compressedIndexToIndexList.
     *
     */
    void finalizeSetup();

    /** @brief Adds a value to the given index.
     *
     *
     *
     *  @param [in] index The index at which the value is added.
     *  @param [in] value The value to add the index to.
     */
    void add(indexType index, indexType value);

    /** @brief Returns a vector of values associated with the given index.
     *
     *
     *  @param [in] index The index at which the entries are returned.
     *  @return std::vector<indexType> with values associated to the index.
     */
    auto get(indexType index) -> std::vector<indexType>;



private:
    bool m_initialized;
    indexType m_startIndex;
    std::vector<indexType> m_index;
    std::vector<indexType> m_data;
};

} // namespace HierAMuS