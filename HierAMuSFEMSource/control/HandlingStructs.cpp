// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HandlingStructs.h"

namespace HierAMuS {
auto trim(const std::string &str, const std::string &whitespace)
    -> std::string 
   {
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content

  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}
InfoData::InfoData() {

}
} // namespace HierAMuS
