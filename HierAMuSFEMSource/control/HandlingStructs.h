// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


#include "control/OutputHandler.h"

#include <map>
#include <string>


namespace HierAMuS {

	/**
	* @brief Contains informational datas.
	* @param ifile Input file for reading.
	* @param ofile Output file for reading
	* @param fileNames Map which contains data about the files.
	*/
	enum class FileHandling {infile, outfile, directory};


	auto trim(const std::string& str,
            const std::string& whitespace = " \t") -> std::string;

	/**
	 * @brief InfoData class which contains information of the current input file and the output object.
	 * Contains a map which maps the strings of infile, outfile, directory to the corresponding enums.
	 */
	struct InfoData {
	  InfoData();
		std::map<FileHandling, std::string> fileNames;
		OutputHandler Log;
		bool runProg;
		std::string envName = "CPPFEM";

	};


} /*End of namespace*/

