set(indexType "long long" CACHE STRING "Set the index data type")
set_property(CACHE indexType PROPERTY STRINGS int64_t int32_t "long long" )

set(Precision double CACHE STRING "Precision setting with which the calculation will be performed")
set_property(CACHE Precision PROPERTY STRINGS  double single)

set(prepstring "#pragma once\n")
string(APPEND prepstring "#include <vector>\n")
string(APPEND prepstring "#include <string>\n")
string(APPEND prepstring "\n")
string(APPEND prepstring "typedef std::vector<std::string> stringVector;\n")

if("${Precision}" STREQUAL "single")
	string(APPEND prepstring "typedef float prec;\n")
endif()

if("${Precision}" STREQUAL "double")
	string(APPEND prepstring "typedef double prec;\n")
endif()


if("${indexType}" STREQUAL "int64_t")
	string(APPEND prepstring "typedef int64_t indexType;\n")
endif()

if("${indexType}" STREQUAL "int32_t")
	string(APPEND prepstring "typedef int32_t indexType;\n")
endif()

if("${indexType}" STREQUAL "long long")
	string(APPEND prepstring "typedef long long indexType;\n")
endif()


set(BOOST_REQUIRED OFF)

add_definitions("-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=${indexType}")

configure_file(${CMAKE_SOURCE_DIR}/datatypes.h.in ${CMAKE_SOURCE_DIR}/datatypes.h @ONLY)
