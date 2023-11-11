// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "compressedIndexToIndexList.h"


#include <sstream>
#include <vector>


namespace HierAMuS {

void compressedIndexToIndexList::initialize(indexType startIndex, indexType endIndex){
    if(m_initialized){
        std::stringstream msg;
        msg << "Trying to initialize a compressedIndexToIndexList twice!";
        throw std::runtime_error(msg.str());
    }
    m_startIndex = startIndex;
    m_index.resize(endIndex-startIndex+2);
    std::fill(m_index.begin(), m_index.end(), 0);
    m_initialized = true;
}

void compressedIndexToIndexList::countAtIndex(indexType index){

    m_index[index-m_startIndex] += 1;
}


void compressedIndexToIndexList::finalizeSetup(){
    indexType pos = 0;
    indexType sum = 0;
    indexType tmp;
    for(auto i=0;i<m_index.size();++i){
        tmp = m_index[i];
        m_index[i] = pos;
        pos += tmp;
    }

    sum = pos + tmp;
    m_data.resize(sum);
    std::fill(m_data.begin(), m_data.end(), -1);
}

void compressedIndexToIndexList::add(indexType index, indexType value){
    bool search = true;
    indexType pos = m_index[index-m_startIndex];
    while(search){
        if(m_data[pos] == -1){
            m_data[pos] = value;
            search = false;
        }
        pos++;
    }
    if(pos > m_index[index+1]){
        throw std::runtime_error("Trying to add more entries than counted!");
    }
}

auto compressedIndexToIndexList::get(indexType index) -> std::vector<indexType> {
    std::vector<indexType> ret(&m_data[m_index[index-m_startIndex]], &m_data[m_index[index-m_startIndex+1]]);
    return ret;
}


}