#include "tinyxml2.h"
#include <gmsh.h>
#include <memory>
#include <string>
#include <vector>

using namespace tinyxml2;

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
