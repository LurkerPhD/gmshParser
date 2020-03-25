/*=============================================================================
*
*   Copyright (C) 2020 All rights reserved.
*
*   Filename:		gmsh2json.cpp
*
*   Author: Wang Zhecheng - wangzhecheng@yeah.net
*
*   Date: 2020-02-14 04:49
*
*   Last Editors: Wang Zhecheng - wangzhecheng@yeah.net
*
*   Last modified:	2020-03-21 01:14
*
*   Description:
*
=============================================================================*/
#include "gmshParser.hpp"
#include <filesystem>

void parse_directory(const std::filesystem::path& path) {
  gmshParser<gmsh2json> parser;
  if(path.empty())
    return;
  for(auto& p : std::filesystem::recursive_directory_iterator(path))
    if(p.path().extension() == ".msh")
      parser.parseGmshFile(p);

  return;
}

int main(int argc, char const* argv[]) {
  std::string filename;

  if(argc > 1)
    filename = argv[1];
  else
    // filename = "/Users/lurker.phd/Documents/geoxfem/gmshParser";
    // filename = "/Users/lurker.phd/Documents/geoxfem/demo";
    //  filename = "/Users/lurker.phd/Documents/geoxfem/test/test_crack";
     filename = "/Users/lurker.phd/Documents/geoxfem/test/test_sif";

  parse_directory(filename.c_str());

  return 0;
}
