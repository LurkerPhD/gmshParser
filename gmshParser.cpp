#include <algorithm>
#include <gmsh.h>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>

#include "tinyxml2.h"

using namespace tinyxml2;

int XMLCreater(XMLDocument *doc);

void insertElemMetaData(const char *typeID, XMLElement *elemMetaData);

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id);

void insertBC(XMLElement *root, std::vector<std::string> &name, const int &id, const std::vector<int> &nodeTags);

void insertMat(XMLElement *root, std::vector<std::string> &tokens, const int &id);

const char *elemType(const int &dim, const int &nSize);

const std::string elemNameRule(const int &dim, const int &nEdge, const char *interpolation);

XMLElement *QueryElementByAttribute(XMLElement *root, const std::string &Attri_Name, const std::string &value);

int main(int argc, char *argv[])
{
  /* code */
  char *filename;
  if (argc > 1)
    filename = argv[1];
  else
    filename = "demo/test1.msh";
  XMLDocument *doc = new XMLDocument;
  XMLCreater(doc);
  std::cout << "xml file created.\n";

  gmsh::initialize();
  if (filename == nullptr)
    filename = (char *)"demo/test1.msh";
  gmsh::open(filename);
  std::cout << "gmsh file imported.\n";
  int dimension = gmsh::model::getDimension();
  gmsh::model::mesh::renumberNodes();
  gmsh::model::mesh::renumberElements();
  // get all nodes
  {
    std::vector<int> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);
    // std::cout << "number of nodes: " << nodeTags.size() << "\n";
    XMLElement *nodes = doc->FirstChildElement("topology")
                            ->FirstChildElement("nodes");
    for (int inode : nodeTags)
    {
      int nodeNum = nodes->IntAttribute("number") + 1;
      nodes->SetAttribute("number", nodeNum);
      XMLElement *newNode = doc->NewElement("node");
      nodes->InsertEndChild(newNode);
      newNode->SetAttribute("id", nodeNum);
      // newNode->SetAttribute("type", ("femNode" + std::to_string(dimension) + "D").c_str());
      std::vector<int> temp = {(inode - 1) * 3 + 0, (inode - 1) * 3 + 1, (inode - 1) * 3 + 2};
      std::vector<double> temp2 = {coord[(inode - 1) * 3 + 0], coord[(inode - 1) * 3 + 1], coord[(inode - 1) * 3 + 2]};
      newNode->SetAttribute("x", coord[(inode - 1) * 3 + 0]);
      if (dimension > 1)
        newNode->SetAttribute("y", coord[(inode - 1) * 3 + 1]);
      if (dimension > 2)
        newNode->SetAttribute("z", coord[(inode - 1) * 3 + 2]);
    }
  }
  std::cout << "node information parsed.\n";
  //get all elements
  {
    XMLElement *elements = doc->FirstChildElement("topology")
                               ->FirstChildElement("elements");
    // std::cout << "step 0.\n";

    XMLElement *materials = doc->FirstChildElement("topology")
                                ->FirstChildElement("materials");
    // std::cout << "step 1.\n";
    XMLElement *elemMetaData = doc->FirstChildElement("topology")
                                   ->FirstChildElement("elemMetaData");
    std::vector<int> elementTypes;
    std::vector<std::vector<int>> elementTags;
    std::vector<std::vector<int>> nodeTags;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, dimension);
    for (auto x : dimTags)
    {
      std::vector<int> PhyTags;
      std::string PhyName;
      int mat_ID;
      gmsh::model::getPhysicalGroupsForEntity(x.first, x.second, PhyTags);
      switch (PhyTags.size())
      {
      case 0:
      {
        std::cout << "Entity" << x.second << "has no material!";
        std::cin.get();
        break;
      }
      case 1:
      {
        elementTypes.clear();
        elementTags.clear();
        nodeTags.clear();
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, x.first, x.second);
        gmsh::model::getPhysicalName(dimension, PhyTags[0], PhyName);
        // parsing physical names
        std::vector<std::string> tokens;
        std::stringstream ss;
        ss.str(PhyName);
        while (ss.good())
        {
          std::string substr;
          std::getline(ss, substr, '_');
          tokens.push_back(substr);
        }
        int mat_pos = tokens[0].find("mat");
        if (mat_pos != std::string::npos)
        {
          std::string mat_ID_str = tokens[0].substr(mat_pos + 3, 1);
          mat_ID = std::stoi(mat_ID_str);
          insertMat(materials, tokens, mat_ID);
        }
        int elemNum = 0;
        for (int i = 0; i < elementTypes.size(); i++)
        {
          std::string elemName;
          int order;
          int numNodes;
          std::vector<double> parametricCoord;
          gmsh::model::mesh::getElementProperties(elementTypes[i], elemName, dimension, order, numNodes, parametricCoord);
          for (int j = 0; j < elementTags[i].size(); j++)
          {
            elemNum++;
            elements->SetAttribute("number", elemNum);
            XMLElement *newElem = doc->NewElement("elem");
            elements->InsertEndChild(newElem);
            newElem->SetAttribute("id", elemNum);
            // std::string type = "C" + std::to_string(dimension) + "D" + std::to_string(numNodes) + "PE";
            std::string type;
            type = elemNameRule(dimension, numNodes, "Lin");
            newElem->SetAttribute("type", type.c_str());
            insertElemMetaData(type.c_str(), elemMetaData);
            newElem->SetAttribute("mat", mat_ID);
            for (int k = 0; k < numNodes; k++)
              newElem->SetAttribute(("v" + std::to_string(k + 1)).c_str(),
                                    nodeTags[i][j * numNodes + k]);
          }
        }
        break;
      }
      default:
      {
        std::cout << "Entity" << x.second << "contains more than one materials!";
        std::cin.get();
      }
      }
    }
  }
  std::cout << "element information parsed.\n";

  //get all boundary condition information
  {
    XMLElement *phases = doc
                             ->FirstChildElement("phases");
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);
    int Ux_id = 0;
    int Uy_id = 0;
    int Uz_id = 0;
    int pore_id = 0;
    int Fx_id = 0;
    int Fy_id = 0;
    int Fz_id = 0;
    int flow_id = 0;
    for (auto x : dimTags)
    {
      std::string PhyName;
      std::vector<int> nodeTags;
      std::vector<double> coords;
      gmsh::model::getPhysicalName(x.first, x.second, PhyName);
      // parsing physical names
      std::vector<std::string> tokens;
      std::stringstream ss;
      ss.str(PhyName);
      while (ss.good())
      {
        std::string substr;
        std::getline(ss, substr, '_');
        tokens.push_back(substr);
      }
      int P_pos = tokens[0].find("P");
      if (P_pos != std::string::npos)
      {
        // when phase detected
        std::string phase_ID = PhyName.substr(P_pos + 1, 1);
        XMLElement *phase = phases->QueryElementByAttribute("phase", "id", phase_ID.c_str());
        if (phase == nullptr)
        {
          phase = doc->NewElement("phase");
          phases->InsertEndChild(phase);
          int nPhase = phases->IntAttribute("number");
          phases->SetAttribute("number", ++nPhase);
          initialPhase(doc, phase, std::stoi(phase_ID));
          phase->SetAttribute("id", phase_ID.c_str());
        }
        // Dirichlet detected
        if (tokens[1].find("Ux") != std::string::npos)
        {
          // Ux boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Ux_id, nodeTags);
          continue;
        }
        if (tokens[1].find("Uy") != std::string::npos)
        {
          // Uy boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Uy_id, nodeTags);
          continue;
        }
        if (tokens[1].find("Uz") != std::string::npos)
        {
          // Uz boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Uz_id, nodeTags);
          continue;
        }
        if (tokens[1].find("pore") != std::string::npos)
        {
          // pore boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++pore_id, nodeTags);
          continue;
        }
        // Neumann detected
        if (tokens[1].find("Fx") != std::string::npos)
        {
          // Fx boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fx_id, nodeTags);
          continue;
        }
        if (tokens[1].find("Fy") != std::string::npos)
        {
          // Fy boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fy_id, nodeTags);
          continue;
        }
        if (tokens[1].find("Fz") != std::string::npos)
        {
          // Fz boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fz_id, nodeTags);
          continue;
        }
        if (tokens[1].find("flow") != std::string::npos)
        {
          // flow boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++flow_id, nodeTags);
          continue;
        }
      }
    }
  }
  std::cout << "boundary condition information parsed.\n";
  doc->SaveFile(strcat(filename, ".xml"));
  std::cout << "xml file has been written.\n";
  delete doc;
  return 0;
}

void insertElemMetaData(const char *MetaDataID, XMLElement *elemMetaData)
{
  if (!elemMetaData->FirstChildElement(MetaDataID))
  {
    XMLElement *newElemMetaData = elemMetaData->GetDocument()->NewElement(MetaDataID);
    elemMetaData->InsertEndChild(newElemMetaData);
    newElemMetaData->SetAttribute("nGauss", 4);
    int n = (elemMetaData->IntAttribute("number"));
    elemMetaData->SetAttribute("number", ++n);
  }
  return;
}

int XMLCreater(XMLDocument *doc)
{
  const char *declaration =
      "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
  doc->Parse(declaration); //会覆盖xml所有内容
  // root node Describe
  XMLElement *describe = doc->NewElement("describe");
  doc->InsertEndChild(describe);
  {
    XMLComment *describeComment = doc->NewComment("This is an example of input with xml format");
    describe->InsertEndChild(describeComment);
    // XMLText *describeText =
    //     doc->NewText("This is a example of input with xml formate");
    // describe->InsertEndChild(describeText);
    describe->SetAttribute("method", "fem");
    describe->SetAttribute("problem", "static");
    describe->SetAttribute("couple", "solid");
    describe->SetAttribute("equation", "implicit");
    describe->SetAttribute("solver", "default");
    describe->SetAttribute("dim", 2);
  }
  // root node start
  XMLElement *breakpoint = doc->NewElement("breakpoint");
  doc->InsertEndChild(breakpoint);
  breakpoint->SetAttribute("file", "");

  // root node Problem
  // output node in doc
  XMLElement *output = doc->NewElement("output");
  doc->InsertEndChild(output);
  output->SetAttribute("path", "");
  output->SetAttribute("plot", "mesh+gauss+bc");
  output->SetAttribute("dispAmplifier", 200);

  // calConfig node in doc
  XMLElement *calConfig = doc->NewElement("calConfig");
  doc->InsertEndChild(calConfig);
  {
    // iteration node in calConfig
    XMLElement *iteration = doc->NewElement("iteration");
    calConfig->InsertEndChild(iteration);
    iteration->SetAttribute("maxLimit", 100);
    // loadApply node in calConfig
    XMLElement *loadApply = doc->NewElement("loadApply");
    calConfig->InsertEndChild(loadApply);
    loadApply->SetAttribute("incremental", true);
    // largeDeform node in calConfig
    XMLElement *largeDeform = doc->NewElement("largeDeform");
    calConfig->InsertEndChild(largeDeform);
    largeDeform->SetAttribute("on", false);
    // kayMat node in calConfig
    XMLElement *kayMat = doc->NewElement("kayMat");
    calConfig->InsertEndChild(kayMat);
    kayMat->SetAttribute("constant", false); 
    // gravity node in calConfig
    XMLElement *gravity = doc->NewElement("gravity");
    calConfig->InsertEndChild(gravity);
    gravity->SetAttribute("on", false);
    // xfem node in calConfig
    XMLElement *xfem = doc->NewElement("xfem");
    calConfig->InsertEndChild(xfem);
    {
      xfem->SetAttribute("active", true);
      // crackSegLength node in xfem
      XMLElement *crackSegLength = doc->NewElement("crackSegment");
      xfem->InsertEndChild(crackSegLength);
      crackSegLength->SetAttribute("dLength", 0.03);
      // SIF node in xfem
      XMLElement *SIF = doc->NewElement("SIF");
      xfem->InsertEndChild(SIF);
      SIF->SetAttribute("calMethod", "Jint");
      {
        // radiusFactor node in SIF
        XMLElement *radiusFactor = doc->NewElement("integralRadius");
        SIF->InsertEndChild(radiusFactor);
        radiusFactor->SetAttribute("factor", 9);
        // qShape node in SIF
        XMLElement *qShape = doc->NewElement("qShape");
        SIF->InsertEndChild(qShape);
        qShape->SetAttribute("exponent", 1);
      }
      // initCrack node in xfem
      XMLElement *initCrack = doc->NewElement("initCrack");
      xfem->InsertEndChild(initCrack);
      initCrack->SetAttribute("number", 0);
    } // end xfem
  }   // end calConfig

  // mesh node in doc
  XMLElement *mesh = doc->NewElement("topology");
  doc->InsertEndChild(mesh);
  {
    // elemMetaData node in doc
    XMLElement *elemMetaData = doc->NewElement("elemMetaData");
    mesh->InsertEndChild(elemMetaData);
    elemMetaData->SetAttribute("number", 0);
    // materials node in doc
    XMLElement *materials = doc->NewElement("materials");
    mesh->InsertEndChild(materials);
    materials->SetAttribute("number", 0);
    // nodes node in mesh
    XMLElement *nodes = doc->NewElement("nodes");
    mesh->InsertEndChild(nodes);
    nodes->SetAttribute("number", 0);
    // elements node in mesh
    XMLElement *elements = doc->NewElement("elements");
    mesh->InsertEndChild(elements);
    elements->SetAttribute("number", 0);
  }

  // phases node in doc
  XMLElement *phases = doc->NewElement("phases");
  doc->InsertEndChild(phases);
  phases->SetAttribute("number", 0);

  return 0;
}

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id)
{
  phase->SetAttribute("id", id);
  phase->SetAttribute("totalTime", "1000");
  phase->SetAttribute("deltaTime", "10");
  phase->SetAttribute("nElem", "9");
  // phase->SetAttribute("BCnumber", "0");
  // {
  //   // Dirichlet node in phase
  //   XMLElement *Dirichlet = doc->NewElement("Dirichlet");
  //   phase->InsertEndChild(Dirichlet);
  //   Dirichlet->SetAttribute("number", 0);
  //   XMLElement *Neumann = doc->NewElement("Neumann");
  //   phase->InsertEndChild(Neumann);
  //   Neumann->SetAttribute("number", 0);
  // }
}

void insertBC(XMLElement *root, std::vector<std::string> &tokens, const int &id, const std::vector<int> &nodeTags)
{
  XMLDocument *doc = root->GetDocument();
  XMLElement *newBC = doc->NewElement((tokens[1] + "BC").c_str());
  root->InsertEndChild(newBC);
  int BCsize = root->IntAttribute("BCsize");
  root->SetAttribute("BCsize", ++BCsize);
  newBC->SetAttribute("id", id);
  newBC->SetAttribute("type", "linear");
  for (int node : nodeTags)
  {
    XMLElement *newNode = doc->NewElement(tokens[1].c_str());
    newBC->InsertEndChild(newNode);
    newNode->SetAttribute("node", node);
    newNode->SetAttribute("start", (tokens[0] + "_" + tokens[1] + tokens[2] + "_start").c_str());
    newNode->SetAttribute("end", (tokens[0] + "_" + tokens[1] + tokens[2] + "_end").c_str());
  }
}

void insertMat(XMLElement *root, std::vector<std::string> &tokens, const int &id)
{
  if (root->QueryElementByAttribute("mat", "id", std::to_string(id).c_str()))
    return;
  XMLElement *mat = root->FindOrCreatChildElement("mat", "id", std::to_string(id).c_str());
  int nMat = std::stoi(root->Attribute("number"));
  root->SetAttribute("number", ++nMat);
  std::string type = tokens[1];
  XMLDocument *doc = root->GetDocument();
  XMLElement *rho = doc->NewElement("rho");
  mat->InsertEndChild(rho);
  rho->SetAttribute("value", "1");
  XMLElement *fricAngle = doc->NewElement("fricAngle");
  mat->InsertEndChild(fricAngle);
  fricAngle->SetAttribute("value", "2");
  XMLElement *coh = doc->NewElement("coh");
  mat->InsertEndChild(coh);
  coh->SetAttribute("value", "3");
  XMLElement *kx = doc->NewElement("kx");
  mat->InsertEndChild(kx);
  kx->SetAttribute("value", "4");
  XMLElement *ky = doc->NewElement("ky");
  mat->InsertEndChild(ky);
  ky->SetAttribute("value", "5");
  XMLElement *ft = doc->NewElement("ft");
  mat->InsertEndChild(ft);
  ft->SetAttribute("value", "6");
  XMLElement *precon = doc->NewElement("precon");
  mat->InsertEndChild(precon);
  precon->SetAttribute("value", "7");
  XMLElement *Kc = doc->NewElement("Kc");
  mat->InsertEndChild(Kc);
  Kc->SetAttribute("value", "8");
  XMLElement *tModule = doc->NewElement("tModule");
  mat->InsertEndChild(tModule);
  tModule->SetAttribute("value", "9");
  XMLElement *tPoisson = doc->NewElement("tPoisson");
  mat->InsertEndChild(tPoisson);
  tPoisson->SetAttribute("value", "10");
  XMLElement *eModule = doc->NewElement("eModule");
  mat->InsertEndChild(eModule);
  eModule->SetAttribute("value", "11");
  XMLElement *ePoisson = doc->NewElement("ePoisson");
  mat->InsertEndChild(ePoisson);
  ePoisson->SetAttribute("value", "12");
  if (type == "elastic")
  {
    XMLElement *eModule = doc->NewElement("eModule");
    mat->InsertEndChild(eModule);
    eModule->SetAttribute("value", "101");
    XMLElement *ePoisson = doc->NewElement("ePoisson");
    mat->InsertEndChild(ePoisson);
    ePoisson->SetAttribute("value", "102");
  }
  else if (type == "DuncanChang")
  {
    XMLElement *Rf = doc->NewElement("Rf");
    mat->InsertEndChild(Rf);
    Rf->SetAttribute("value", "201");
    XMLElement *k = doc->NewElement("k");
    mat->InsertEndChild(k);
    k->SetAttribute("value", "202");
    XMLElement *n = doc->NewElement("n");
    mat->InsertEndChild(n);
    n->SetAttribute("value", "203");
    XMLElement *G = doc->NewElement("G");
    mat->InsertEndChild(G);
    G->SetAttribute("value", "204");
    XMLElement *F = doc->NewElement("F");
    mat->InsertEndChild(F);
    F->SetAttribute("value", "205");
    XMLElement *d = doc->NewElement("d");
    mat->InsertEndChild(d);
    d->SetAttribute("value", "206");
  }
  else if (type == "CamClay")
  {
    ;
  }
  else if (type == "MohrColumb")
  {
    ;
  }
  else if (type == "DruckPrager")
  {
    ;
  }
  mat->SetAttribute("type", type.c_str());
  return;
}

const char *elemType(const int &dim, const int &nSize)
{
  switch (dim)
  {
  case 1:
    return "truss";
    break;
  case 2:
    switch (nSize)
    {
    case 3:
      return "triangle";
      break;
    case 4:
      return "quadrilateral";
      break;
    default:
      std::cout << "Invalid nSize in 2D element!" << std::endl;
      std::cin.get();
      break;
    }
    break;
  case 3:
    switch (nSize)
    {
    case 4:
      return "tetrahedron";
      break;
    case 6:
      return "prism";
      break;
    case 8:
      return "hexahedron";
      break;
    default:
      break;
    }
    break;
  default:
    std::cout << "Invalid dimensions:" << dim << std::endl;
    std::cin.get();
    break;
  }
  return "";
}

const std::string elemNameRule(const int &dim, const int &nEdge, const char *interpolation)
{
  std::string dimention;
  switch (dim)
  {
  case 1:
    dimention = "T";
    break;
  case 2:
    dimention = "PE";
    break;
  case 3:
    dimention = "C";
    break;
  default:
    std::cout << "Error dimension!\n";
    break;
  }
  return dimention + std::to_string(nEdge) + "B" + interpolation;
}