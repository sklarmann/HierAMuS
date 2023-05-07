// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#ifndef FILEHANDLING_H_
#define FILEHANDLING_H_

#ifndef __LINUX_

#include <map>
#include <string>
#include <fstream>
#include <control/HandlingStructs.h>

namespace HierAMuS{
	void set_files(std::map<FileHandling, std::string> &fileinformation);
}
#endif

#endif /* FILEHANDLING_H_ */
