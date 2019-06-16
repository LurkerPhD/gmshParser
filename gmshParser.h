#include "tinyxml2.h"
#include <gmsh.h>
#include <string>
#include <vector>

using namespace tinyxml2;

class gmshParser
{
private:
    ///
    XMLDocument *doc;

    std::string filename;

    ///
    int dimension;

public:
    /// default constructor
    gmshParser() {}

    /// default destructor
    ~gmshParser() {}

    /// get class name
    std::string getClassName() const { return "gmshParser"; }

    ///----------------------------------------------------------------------------
    ///----------------------------------------------------------------------------
    ///----------------------------------------------------------------------------

    void readFromGmshFile(const std::string &_filename);
    ///----------------------------------------------------------------------------

private:
    void buildMainFrame();

    XMLElement *addChildElement(const std::string &_name, XMLNode *root);

    XMLComment *addChildComment(const std::string &_content, XMLNode *root);

    XMLText *addChildText(const std::string &_content, XMLNode *root);
    ///----------------------------------------------------------------------------

    void addDescribe();

    void addBreakPoint();

    void addOutput();

    void addCalculationConfiguration();

    void addTopology();

    void addMaterials();

    void addInitCondition();

    void addPhases();

    void addNodes();

    void addElements();

    void addMaterial(std::vector<std::string> &tokens, const int &id);

    void addBoundaries();

    void initialPhase(XMLElement *root, const int &id);

    void insertDirichlet(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags);

    void insertNeumann(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags);
};
