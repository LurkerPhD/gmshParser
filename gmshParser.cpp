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

void initialPhase(XMLDocument *doc, XMLElement *problem, const int &id);

void insertBC(XMLElement *root, std::vector<std::string> &name, const int &id, const std::vector<int> &nodeTags);

XMLElement *QueryElementByAttribute(XMLElement *root, const std::string &Attri_Name, const std::string &value);

int main(int argc, char **argv)
{
  /* code */

  XMLDocument *doc = new XMLDocument;
  XMLCreater(doc);
  std::cout << "Step 0\n";

  gmsh::initialize(argc, argv);
  gmsh::open("t3.msh");
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
    XMLElement *nodes = doc->FirstChildElement("problem")
                            ->FirstChildElement("mesh")
                            ->FirstChildElement("nodes");
    for (int inode : nodeTags)
    {
      int nodeNum = nodes->IntAttribute("size") + 1;
      nodes->SetAttribute("size", nodeNum);
      XMLElement *newNode = doc->NewElement("node");
      nodes->InsertEndChild(newNode);
      newNode->SetAttribute("id", nodeNum);
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
    XMLElement *elements = doc->FirstChildElement("problem")
                               ->FirstChildElement("mesh")
                               ->FirstChildElement("elements");
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
        int mat_pos = PhyName.find("mat");
        if (mat_pos != std::string::npos)
          mat_ID = std::stoi(PhyName.substr(mat_pos + 3, 1));
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
  //get all boundary condition information
  {
    XMLElement *phases = doc->FirstChildElement("problem")
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
          initialPhase(doc, phase, std::stoi(phase_ID));
          phase->SetAttribute("id", phase_ID.c_str());
        }
        // Dirichlet detected
        if (tokens[1].find("Ux") != std::string::npos)
        {
          // Ux boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Ux_id, nodeTags);
          break;
        }
        if (tokens[1].find("Uy") != std::string::npos)
        {
          // Uy boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Uy_id, nodeTags);
          break;
        }
        if (tokens[1].find("Uz") != std::string::npos)
        {
          // Uz boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Uz_id, nodeTags);
          break;
        }
        if (tokens[1].find("pore") != std::string::npos)
        {
          // pore boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Dirichlet");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++pore_id, nodeTags);
          break;
        }
        // Neumann detected
        if (tokens[1].find("Fx") != std::string::npos)
        {
          // Fx boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Fx_id, nodeTags);
          break;
        }
        if (tokens[1].find("Fy") != std::string::npos)
        {
          // Fy boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Fy_id, nodeTags);
          break;
        }
        if (tokens[1].find("Fz") != std::string::npos)
        {
          // Fz boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++Fz_id, nodeTags);
          break;
        }
        if (tokens[1].find("flow") != std::string::npos)
        {
          // flow boundary condition
          XMLElement *entry = phase->FirstChildElement("BndryCondition")
                                  ->FirstChildElement("Neumann");
          gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
          insertBC(entry, tokens, ++flow_id, nodeTags);
          break;
        }
      }
    }
  }
  std::cout << "done\n";
  doc->SaveFile("output.xml");
  delete doc;
  return 0;
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
    XMLText *describeText =
        doc->NewText("This is a example of input with xml formate");
    describe->InsertEndChild(describeText);
  }

  // root node Problem
  XMLElement *problem = doc->NewElement("problem");
  doc->InsertEndChild(problem);
  {
    problem->SetAttribute("name", "inputTest");
    // output node in problem
    XMLElement *output = doc->NewElement("output");
    problem->InsertEndChild(output);
    {
      // deformation amplify node in output
      XMLElement *deformAmplify = doc->NewElement("deformAmplify");
      output->InsertEndChild(deformAmplify);
      deformAmplify->SetAttribute("value", 200);
    } // end output
    // modelConfig node in problem
    XMLElement *modelConfig = doc->NewElement("modelConfig");
    problem->InsertEndChild(modelConfig);
    {
      // iterLimit node in modelConfig
      XMLElement *iterLimit = doc->NewElement("iterLimit");
      modelConfig->InsertEndChild(iterLimit);
      iterLimit->SetAttribute("value", 100);
      // gravity node in modelConfig
      XMLElement *gravity = doc->NewElement("gravity");
      modelConfig->InsertEndChild(gravity);
      gravity->SetAttribute("value", false);
      // fluid_solid_coupling node in modelConfig
      XMLElement *fluid_solid_coupling =
          doc->NewElement("fluid_solid_coupling");
      modelConfig->InsertEndChild(fluid_solid_coupling);
      fluid_solid_coupling->SetAttribute("value", false);
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

    // phases node in problem
    XMLElement *phases = doc->NewElement("phases");
    problem->InsertEndChild(phases);
    phases->SetAttribute("size", 0);

    // materials node in problem
    XMLElement *materials = doc->NewElement("materials");
    problem->InsertEndChild(materials);
    {
      materials->SetAttribute("size", 1);
      // mat1 node in materials
      XMLElement *mat1 = doc->NewElement("mat1");
      materials->InsertEndChild(mat1);
      mat1->SetAttribute("id", 1);
      mat1->SetAttribute("type", "elastic");
    }
    // mesh node in problem
    XMLElement *mesh = doc->NewElement("mesh");
    problem->InsertEndChild(mesh);
    {
      mesh->SetAttribute("dim", 2);
      // nodes node in mesh
      XMLElement *nodes = doc->NewElement("nodes");
      mesh->InsertEndChild(nodes);
      nodes->SetAttribute("size", 0);
      // elements node in mesh
      XMLElement *elements = doc->NewElement("elements");
      mesh->InsertEndChild(elements);
      elements->SetAttribute("size", 0);
    }
  } // end problem
  return 0;
}

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id)
{
  phase->SetAttribute("id", id);
  {
    // totalTime node in phase
    XMLElement *totalTime = doc->NewElement("totalTime");
    phase->InsertEndChild(totalTime);
    totalTime->SetAttribute("value", 1000);
    // deltaTime node in phase
    XMLElement *deltaTime = doc->NewElement("deltaTime");
    phase->InsertEndChild(deltaTime);
    deltaTime->SetAttribute("value", 10);
    // elemNum node in phase
    XMLElement *elemNum = doc->NewElement("elemNum");
    phase->InsertEndChild(elemNum);
    elemNum->SetAttribute("value", 9);
    // BndryCondition node in phase
    XMLElement *BndryCondition = doc->NewElement("BndryCondition");
    phase->InsertEndChild(BndryCondition);
    {
      // Dirichlet node in BndryCondition
      XMLElement *Dirichlet = doc->NewElement("Dirichlet");
      BndryCondition->InsertEndChild(Dirichlet);
      XMLElement *Neumann = doc->NewElement("Neumann");
      BndryCondition->InsertEndChild(Neumann);
    }
  }
}

void insertBC(XMLElement *root, std::vector<std::string> &tokens, const int &id, const std::vector<int> &nodeTags)
{
  XMLDocument *doc = root->GetDocument();
  XMLElement *newBC = doc->NewElement(tokens[1].c_str());
  root->InsertEndChild(newBC);
  newBC->SetAttribute("id", id);
  newBC->SetAttribute("type", "progressive");
  for (int node : nodeTags)
  {
    XMLElement *newNode = doc->NewElement("BCnode");
    newBC->InsertEndChild(newNode);
    newNode->SetAttribute("id", node);
    newNode->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_start").c_str());
    newNode->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_end").c_str());
  }
}
