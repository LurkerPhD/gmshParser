#include "tinyxml2.h"
#include <gmsh.h>
#include <memory>
#include <string>
#include <vector>

using namespace tinyxml2;

std::pair<bool, std::string> parseName(const std::string &_src, const std::string &_patternStr, const std::size_t &_index = 1);

class gmshParser
{
private:
    class impl;
    std::unique_ptr<impl> m_pImpl;

public:
    /// default constructor
    gmshParser();

    /// default destructor
    ~gmshParser();

    /// get class name
    std::string getClassName() const { return "gmshParser"; }

    ///----------------------------------------------------------------------------
    ///----------------------------------------------------------------------------
    ///----------------------------------------------------------------------------

    void readFromGmshFile(const std::string &_filename);
    ///----------------------------------------------------------------------------
};
