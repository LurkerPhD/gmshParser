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

void insertElemType(const char *typeID, XMLElement *elemType);

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id);

void insertBC(XMLElement *root, std::vector<std::string> &name, const int &id, const std::vector<int> &nodeTags);

void insertMat(XMLElement *root, std::vector<std::string> &tokens, const int &id);

const char *elemType(const int &dim, const int &nSize);

XMLElement *
QueryElementByAttribute(XMLElement *root, const std::string &Attri_Name, const std::string &value);

int main(int argc, char *argv[])
{
  /* code */
  char *filename;
  filename = argv[1];
  XMLDocument *doc = new XMLDocument;
  XMLCreater(doc);
  std::cout << "Step 0\n";

  gmsh::initialize();
  if (filename == nullptr)
    filename = (char *)"t3.msh";
  gmsh::open(filename);
  std::cout << "Step 1\n";
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
      int nodeNum = nodes->IntAttribute("size") + 1;
      nodes->SetAttribute("size", nodeNum);
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
  //get all elements
  {
    XMLElement *elements = doc->FirstChildElement("topology")
                               ->FirstChildElement("elements");
    XMLElement *materials = doc->FirstChildElement("topology")
                                ->FirstChildElement("materials");
    XMLElement *elemType = doc->FirstChildElement("topology")
                               ->FirstChildElement("elemType");
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
            elements->SetAttribute("size", elemNum);
            XMLElement *newElem = doc->NewElement("elem");
            elements->InsertEndChild(newElem);
            newElem->SetAttribute("id", elemNum);
            std::string type = "C" + std::to_string(dimension) + "D" + std::to_string(numNodes) + "PE";
            newElem->SetAttribute("type", type.c_str());
            insertElemType(type.c_str(), elemType);
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
    elements->SetAttribute("order", 1);
  }
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
          int nPhase = phases->IntAttribute("size");
          phases->SetAttribute("size", ++nPhase);
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
          break;
        }
        if (tokens[1].find("Uy") != std::string::npos)
        {
          // Uy boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Uy_id, nodeTags);
          break;
        }
        if (tokens[1].find("Uz") != std::string::npos)
        {
          // Uz boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Uz_id, nodeTags);
          break;
        }
        if (tokens[1].find("pore") != std::string::npos)
        {
          // pore boundary condition
          //XMLElement *entry = phase->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++pore_id, nodeTags);
          break;
        }
        // Neumann detected
        if (tokens[1].find("Fx") != std::string::npos)
        {
          // Fx boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fx_id, nodeTags);
          break;
        }
        if (tokens[1].find("Fy") != std::string::npos)
        {
          // Fy boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fy_id, nodeTags);
          break;
        }
        if (tokens[1].find("Fz") != std::string::npos)
        {
          // Fz boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++Fz_id, nodeTags);
          break;
        }
        if (tokens[1].find("flow") != std::string::npos)
        {
          // flow boundary condition
          //XMLElement *entry = phase->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(phase, tokens, ++flow_id, nodeTags);
          break;
        }
      }
    }
  }
  std::cout << "done\n";
  doc->SaveFile(strcat(filename, ".xml"));
  delete doc;
  return 0;
}

void insertElemType(const char *typeID, XMLElement *elemType)
{
  if (!elemType->FirstChildElement(typeID))
  {
    XMLElement *newElemType = elemType->GetDocument()->NewElement(typeID);
    elemType->InsertEndChild(newElemType);
    newElemType->SetAttribute("nGauss", 0);
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
    describe->SetAttribute("model", "femSolid");
    describe->SetAttribute("dim", 2);
  }
  // root node start
  XMLElement *start = doc->NewElement("start");
  doc->InsertEndChild(start);
  start->SetAttribute("file", "");

  // root node Problem
  // output node in doc
  XMLElement *output = doc->NewElement("output");
  doc->InsertEndChild(output);
  output->SetAttribute("path", "");
  output->SetAttribute("dispAmplifier", 200);

  // modelConfig node in doc
  XMLElement *modelConfig = doc->NewElement("modelConfig");
  doc->InsertEndChild(modelConfig);
  {
    // iterLimit node in modelConfig
    XMLElement *iterLimit = doc->NewElement("iterLimit");
    modelConfig->InsertEndChild(iterLimit);
    iterLimit->SetAttribute("value", 100);
    // gravity node in modelConfig
    XMLElement *gravity = doc->NewElement("gravity");
    modelConfig->InsertEndChild(gravity);
    gravity->SetAttribute("value", false);
    // stress_state node in modelConfig
    XMLElement *stress_state = doc->NewElement("stress_state");
    modelConfig->InsertEndChild(stress_state);
    stress_state->SetAttribute("value", "plain_strain");
    // xfem node in modelConfig
    XMLElement *xfem = doc->NewElement("xfem");
    modelConfig->InsertEndChild(xfem);
    {
      xfem->SetAttribute("active", true);
      // crackSegLength node in xfem
      XMLElement *crackSegLength = doc->NewElement("crackSegLength");
      xfem->InsertEndChild(crackSegLength);
      crackSegLength->SetAttribute("value", 0.03);
      // SIFcal node in xfem
      XMLElement *SIFcal = doc->NewElement("SIFcal");
      xfem->InsertEndChild(SIFcal);
      SIFcal->SetAttribute("value", "Jint");
      // SIFfactor node in xfem
      XMLElement *SIFfactor = doc->NewElement("SIFfactor");
      xfem->InsertEndChild(SIFfactor);
      SIFfactor->SetAttribute("value", 9);
      // qShape node in xfem
      XMLElement *qShape = doc->NewElement("qShape");
      xfem->InsertEndChild(qShape);
      qShape->SetAttribute("value", 1);
      // initialSegment node in xfem
      XMLElement *initialSegment = doc->NewElement("initialSegment");
      xfem->InsertEndChild(initialSegment);
      initialSegment->SetAttribute("size", 0);
    } // end xfem
  }   // end modelConfig

  // phases node in doc
  XMLElement *phases = doc->NewElement("phases");
  doc->InsertEndChild(phases);
  phases->SetAttribute("size", 0);

  // mesh node in doc
  XMLElement *mesh = doc->NewElement("topology");
  doc->InsertEndChild(mesh);
  {
    // elemType node in doc
    XMLElement *elemType = doc->NewElement("elemType");
    mesh->InsertEndChild(elemType);
    elemType->SetAttribute("size", 0);
    // materials node in doc
    XMLElement *materials = doc->NewElement("materials");
    mesh->InsertEndChild(materials);
    materials->SetAttribute("size", 0);
    // nodes node in mesh
    XMLElement *nodes = doc->NewElement("nodes");
    mesh->InsertEndChild(nodes);
    nodes->SetAttribute("size", 0);
    // elements node in mesh
    XMLElement *elements = doc->NewElement("elements");
    mesh->InsertEndChild(elements);
    elements->SetAttribute("size", 0);
  }
  return 0;
}

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id)
{
  phase->SetAttribute("id", id);
  phase->SetAttribute("totalTime", "1000");
  phase->SetAttribute("deltaTime", "10");
  phase->SetAttribute("nElem", "9");
  phase->SetAttribute("BCsize", "0");
  // {
  //   // Dirichlet node in phase
  //   XMLElement *Dirichlet = doc->NewElement("Dirichlet");
  //   phase->InsertEndChild(Dirichlet);
  //   Dirichlet->SetAttribute("size", 0);
  //   XMLElement *Neumann = doc->NewElement("Neumann");
  //   phase->InsertEndChild(Neumann);
  //   Neumann->SetAttribute("size", 0);
  // }
}

void insertBC(XMLElement *root, std::vector<std::string> &tokens, const int &id, const std::vector<int> &nodeTags)
{
  XMLDocument *doc = root->GetDocument();
  XMLElement *newBC = doc->NewElement(tokens[1].c_str());
  root->InsertEndChild(newBC);
  int BCsize = root->IntAttribute("BCsize");
  root->SetAttribute("BCsize", ++BCsize);
  newBC->SetAttribute("id", id);
  newBC->SetAttribute("type", "progressive");
  for (int node : nodeTags)
  {
    XMLElement *newNode = doc->NewElement("BCnode");
    newBC->InsertEndChild(newNode);
    newNode->SetAttribute("id", node);
    newNode->SetAttribute("start", (tokens[0] + "_" + tokens[1] + tokens[2] + "_start").c_str());
    newNode->SetAttribute("end", (tokens[0] + "_" + tokens[1] + tokens[2] + "_end").c_str());
  }
}

void insertMat(XMLElement *root, std::vector<std::string> &tokens, const int &id)
{
  if (root->QueryElementByAttribute("mat", "id", std::to_string(id).c_str()))
    return;
  XMLElement *mat = root->FindOrCreatChildElement("mat", "id", std::to_string(id).c_str());
  int nMat = std::stoi(root->Attribute("size"));
  root->SetAttribute("size", ++nMat);
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
  XMLElement *reModule = doc->NewElement("reModule");
  mat->InsertEndChild(reModule);
  reModule->SetAttribute("value", "11");
  XMLElement *rePoisson = doc->NewElement("rePoisson");
  mat->InsertEndChild(rePoisson);
  rePoisson->SetAttribute("value", "12");
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