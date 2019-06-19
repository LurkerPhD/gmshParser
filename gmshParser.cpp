#include "gmshParser.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>

///----------------------------------------------------------------------------
void gmshParser::readFromGmshFile(const std::string &_filename)
{
  _file_name_with_path = _filename;
  _file_name = _file_name_with_path.substr(_file_name_with_path.rfind('/') + 1, _file_name_with_path.size());
  _file_name = _file_name.substr(0, _file_name.find('.'));
  ///----------------------------------------------------------------------------
  gmsh::initialize();
  gmsh::open(_file_name_with_path);
  std::cout << "------------------------------------\n"
            << "gmsh file imported.\n";
  dimension = gmsh::model::getDimension();
  gmsh::model::mesh::renumberNodes();
  gmsh::model::mesh::renumberElements();

  ///----------------------------------------------------------------------------
  buildMainFrame();
  std::cout << "xml file created.\n";

  ///----------------------------------------------------------------------------
  /// get all nodes
  addNodes();
  std::cout << "node information parsed.\n";

  ///----------------------------------------------------------------------------
  //get all elements
  addElements();
  std::cout << "element information parsed.\n";

  ///----------------------------------------------------------------------------
  ///get all boundary condition information
  addBoundaries();

  ///----------------------------------------------------------------------------
  const std::string raw_name = _file_name_with_path.substr(0, _file_name_with_path.find("."));
  doc->SaveFile((raw_name + ".xml").c_str());
  std::cout << "xml file has been written.\n";
}
///----------------------------------------------------------------------------

void gmshParser::addNodes()
{
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

  XMLElement *nodes = doc->FirstChildElement("topology")
                          ->FirstChildElement("nodes");
  for (int inode : nodeTags)
  {
    int node_num = nodes->IntAttribute("number") + 1;
    nodes->SetAttribute("number", node_num);
    XMLElement *new_node = addChildElement("node", nodes);
    new_node->SetAttribute("id", node_num);
    new_node->SetAttribute("x1", coord[(inode - 1) * 3 + 0]);
    new_node->SetAttribute("x2", coord[(inode - 1) * 3 + 1]);
    new_node->SetAttribute("x3", coord[(inode - 1) * 3 + 2]);
  }

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addElements()
{
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t>> elementTags;
  std::vector<std::vector<std::size_t>> nodeInElementTags;
  gmsh::vectorpair dimTags;
  gmsh::model::getEntities(dimTags, dimension);

  XMLElement *elements = doc->FirstChildElement("topology")
                             ->FirstChildElement("elements");

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
      std::cout << "Entity" << x.second << "has no s!";
      std::cin.get();
      break;
    }
    case 1:
    {
      elementTypes.clear();
      elementTags.clear();
      nodeInElementTags.clear();
      gmsh::model::mesh::getElements(elementTypes, elementTags, nodeInElementTags, x.first, x.second);
      gmsh::model::getPhysicalName(dimension, PhyTags.at(0), PhyName);
      /// parsing physical names
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
        addMaterial(tokens, mat_ID);
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
          XMLElement *elem = addChildElement("elem", elements);
          elem->SetAttribute("id", elemNum);
          /// std::string type = "C" + std::to_string(dimension) + "D" + std::to_string(numNodes) + "PE";
          std::string type;
          type = "Quadrangle"; //elemNameRule(dimension, numNodes, "Lin");
          elem->SetAttribute("type", type.c_str());
          // insertElemMetaData(type.c_str(), element_mega_data);
          elem->SetAttribute("mat", mat_ID);
          for (int k = 0; k < numNodes; k++)
            elem->SetAttribute(("v" + std::to_string(k + 1)).c_str(),
                               (int)(nodeInElementTags[i][j * numNodes + k]));
        }
      }
      break;
    }
    default:
      std::cout << "Entity" << x.second << "contains more than one materials!";
    }
  }

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addBoundaries()
{
  gmsh::vectorpair dimTags;
  gmsh::model::getPhysicalGroups(dimTags);
  XMLElement *phases = doc->FirstChildElement("phases");

  for (auto x : dimTags)
  {
    std::string PhyName;
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    gmsh::model::getPhysicalName(x.first, x.second, PhyName);
    /// parsing physical names
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
      /// instablish phase node in xml when phase number detected
      int phase_ID = std::stoi(PhyName.substr(P_pos + 1, 1));
      XMLElement *phase = nullptr;
      for (auto _phase = phases->FirstChildElement("phase"); _phase; _phase = phase->NextSiblingElement("phase"))
        if (_phase->IntAttribute("id") == phase_ID)
          phase = _phase;
      if (phase == nullptr)
      {
        phase = addChildElement("phase", phases);
        int phase_num = phases->IntAttribute("number");
        phases->SetAttribute("number", ++phase_num);
        initialPhase(phase, phase_ID);
        phase->SetAttribute("id", phase_ID);
      }
      /// Dirichlet detected
      if (tokens[1].find("U") != std::string::npos || tokens[1].find("pore") != std::string::npos)
      {
        /// Dirichlet condition
        //XMLElement *entry = phase->FirstChildElement("Dirichlet");
        gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
        insertDirichlet(phase->FirstChildElement("boundary_condition"), tokens, nodeTags);
        continue;
      }
      /// Neumann detected
      if (tokens[1].find("F") != std::string::npos || tokens[1].find("flow") != std::string::npos)
      {
        /// Neumann condition
        //XMLElement *entry = phase->FirstChildElement("Neumann");
        gmsh::model::mesh::getNodesForPhysicalGroup(x.first, x.second, nodeTags, coords);
        insertNeumann(phase->FirstChildElement("boundary_condition"), tokens, nodeTags);
        continue;
      }
    }
  };
}
///----------------------------------------------------------------------------

void gmshParser::initialPhase(XMLElement *root, const int &id)
{

  root->SetAttribute("id", id);
  root->SetAttribute("total_time", "1000");
  root->SetAttribute("delta_time", "10");
  root->SetAttribute("element_num", "9");
  XMLElement *newBC = addChildElement("boundary_condition", root);
  newBC->SetAttribute("Dirichlet", 0);
  newBC->SetAttribute("Neumann", 0);

  return;
}
///----------------------------------------------------------------------------

void gmshParser::insertDirichlet(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
{
  for (int bcNode : nodeTags)
  {
    XMLElement *newDirichlet = addChildElement("Dirichlet", root);
    int Dirichlet = root->IntAttribute("Dirichlet");
    root->SetAttribute("Dirichlet", ++Dirichlet);
    newDirichlet->SetAttribute("node", bcNode);
    newDirichlet->SetAttribute("dof", tokens[2].c_str());
    newDirichlet->SetAttribute("type", "linear");
    newDirichlet->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
    newDirichlet->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
  }

  return;
}
///----------------------------------------------------------------------------

void gmshParser::insertNeumann(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
{
  for (int bcNode : nodeTags)
  {
    XMLElement *newNeumann = addChildElement("Neumann", root);
    int Neumann = root->IntAttribute("Neumann");
    root->SetAttribute("Neumann", ++Neumann);
    newNeumann->SetAttribute("node", bcNode);
    newNeumann->SetAttribute("dof", tokens[2].c_str());
    newNeumann->SetAttribute("type", "linear");
    newNeumann->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
    newNeumann->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
  }

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addMaterial(std::vector<std::string> &tokens, const int &id)
{
  XMLElement *materials = doc->FirstChildElement("materials");
  for (auto _mat = materials->FirstChildElement("mat"); _mat; _mat = _mat->NextSiblingElement("mat"))
    if (id == _mat->IntAttribute("id"))
      return;

  XMLElement *_mat = addChildElement("mat", materials);
  int mat_id = std::stoi(materials->Attribute("number"));
  materials->SetAttribute("number", ++mat_id);
  std::string type = tokens[1];
  _mat->SetAttribute("constitutive", type.c_str());

  XMLElement *volume_weight = addChildElement("volume_weight", _mat);
  volume_weight->SetAttribute("value", "19");
  volume_weight->SetAttribute("unit", "kN/m^3");
  XMLElement *friction_angle = addChildElement("friction_angle", _mat);
  friction_angle->SetAttribute("value", "30");
  friction_angle->SetAttribute("unit", "degree");
  XMLElement *dilation_angle = addChildElement("dilation_angle", _mat);
  dilation_angle->SetAttribute("value", "30");
  dilation_angle->SetAttribute("unit", "degree");
  XMLElement *cohesion = addChildElement("cohesion", _mat);
  cohesion->SetAttribute("value", "10");
  cohesion->SetAttribute("unit", "kPa");
  XMLElement *perm_x1 = addChildElement("perm_x1", _mat);
  perm_x1->SetAttribute("value", "1e-5");
  perm_x1->SetAttribute("unit", "m/day");
  XMLElement *perm_x2 = addChildElement("perm_x2", _mat);
  perm_x2->SetAttribute("value", "1e-5");
  perm_x2->SetAttribute("unit", "m/day");
  XMLElement *perm_x3 = addChildElement("perm_x3", _mat);
  perm_x3->SetAttribute("value", "1e-5");
  perm_x3->SetAttribute("unit", "m/day");
  XMLElement *zeta = addChildElement("zeta", _mat);
  zeta->SetAttribute("value", "0.1");
  zeta->SetAttribute("unit", "-");
  XMLElement *K_IC = addChildElement("K_IC", _mat);
  K_IC->SetAttribute("value", "2.9");
  K_IC->SetAttribute("unit", "kN/m^0.5");
  XMLElement *elastic_modulus = addChildElement("elastic_modulus", _mat);
  elastic_modulus->SetAttribute("value", "10000");
  elastic_modulus->SetAttribute("unit", "kPa");
  XMLElement *elastic_poisson = addChildElement("elastic_poisson", _mat);
  elastic_poisson->SetAttribute("value", "0.3");
  elastic_poisson->SetAttribute("unit", "-");
  if (type.find("DuncanChang") != std::string::npos)
  {
    XMLElement *Rf = addChildElement("Rf", _mat);
    Rf->SetAttribute("value", "201");
    Rf->SetAttribute("unit", "-");
    XMLElement *k = addChildElement("k", _mat);
    k->SetAttribute("value", "202");
    k->SetAttribute("unit", "kPa");
    XMLElement *n = addChildElement("n", _mat);
    n->SetAttribute("value", "203");
    n->SetAttribute("unit", "-");
    XMLElement *G = addChildElement("G", _mat);
    G->SetAttribute("value", "204");
    G->SetAttribute("unit", "-");
    XMLElement *F = addChildElement("F", _mat);
    F->SetAttribute("value", "205");
    F->SetAttribute("unit", "-");
    XMLElement *D = addChildElement("D", _mat);
    D->SetAttribute("value", "206");
    D->SetAttribute("unit", "-");
  }
  else if (type.find("ModifiedCamClay") != std::string::npos)
  {
    XMLElement *swelling_slope = addChildElement("swelling_slope", _mat);
    swelling_slope->SetAttribute("value", "401");
    swelling_slope->SetAttribute("unit", "-");
    XMLElement *normal_compress_slope = addChildElement("normal_compress_slope", _mat);
    normal_compress_slope->SetAttribute("value", "402");
    normal_compress_slope->SetAttribute("unit", "-");
    XMLElement *initial_void = addChildElement("initial_void", _mat);
    initial_void->SetAttribute("value", "403");
    initial_void->SetAttribute("unit", "-");
  }
  else if (type.find("DruckPrager") != std::string::npos)
  {
    XMLElement *qf = addChildElement("qf", _mat);
    qf->SetAttribute("value", "501");
    qf->SetAttribute("unit", "-");
    XMLElement *kf = addChildElement("kf", _mat);
    kf->SetAttribute("value", "502");
    kf->SetAttribute("unit", "-");
    XMLElement *q_dilantion = addChildElement("q_dilantion", _mat);
    q_dilantion->SetAttribute("value", "503");
    q_dilantion->SetAttribute("unit", "degree");
  }

  return;
}

///----------------------------------------------------------------------------

void gmshParser::buildMainFrame()
{
  doc = new XMLDocument();
  doc->Parse("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");

  addChildComment("1. model type information", doc);
  addDescribe();

  addChildComment("2. input/output information", doc);
  addBreakPoint();
  addOutput();

  addChildComment("3. calculation configuration information", doc);
  addCalculationConfiguration();

  addChildComment("4. model's topology information", doc);
  addTopology();

  addChildComment("5. s information", doc);
  addMaterials();

  addChildComment("6. initial condition information", doc);
  addInitCondition();

  addChildComment(" 7. phase step information", doc);
  addPhases();
}
///----------------------------------------------------------------------------

XMLElement *gmshParser::addChildElement(const std::string &_name, XMLNode *root)
{
  XMLElement *child = doc->NewElement(_name.c_str());
  root->InsertEndChild(child);

  return child;
}
///----------------------------------------------------------------------------

XMLComment *gmshParser::addChildComment(const std::string &_content, XMLNode *root)
{
  XMLComment *child = doc->NewComment(_content.c_str());
  root->InsertEndChild(child);

  return child;
}
///----------------------------------------------------------------------------

XMLText *gmshParser::addChildText(const std::string &_content, XMLNode *root)
{
  XMLText *child = doc->NewText(_content.c_str());
  root->InsertEndChild(child);

  return child;
}
///----------------------------------------------------------------------------

void gmshParser::addDescribe()
{
  /// root node Describe
  XMLElement *describe = addChildElement("describe", doc);
  describe->SetAttribute("name", _file_name.c_str());
  describe->SetAttribute("method", "FEM");
  describe->SetAttribute("type", "Static");
  describe->SetAttribute("couple", "Solid");
  describe->SetAttribute("equation", "Implicit");
  describe->SetAttribute("solver", "default");
  describe->SetAttribute("dim", dimension);
  if (dimension == 2)
    describe->SetAttribute("plane_stress", false);

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addBreakPoint()
{
  XMLElement *breakpoint = addChildElement("breakpoint", doc);
  breakpoint->SetAttribute("file", "");

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addOutput()
{
  XMLElement *output = addChildElement("output", doc);
  output->SetAttribute("path", ("./result/" + _file_name).c_str());
  XMLElement *visualize = addChildElement("visualize", output);
  visualize->SetAttribute("format", "vtk");
  XMLElement *plot = addChildElement("plot", output);
  plot->SetAttribute("option", "mesh+gauss");
  XMLElement *displacement_amplifier = addChildElement("displacement_amplifier", output);
  displacement_amplifier->SetAttribute("value", 1);
  XMLElement *sequence_output = addChildElement("sequence_output", output);
  addChildComment("<node id=\"0\" attribute=\"x1\"/>", sequence_output);
  addChildComment("<node id=\"1\" attribute=\"x2\"/>", sequence_output);
  addChildComment("<elem id=\"0\" attribute=\"x1\"/>", sequence_output);
  addChildComment("<elem id=\"1\" attribute=\"sigma_x\"/>", sequence_output);

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addCalculationConfiguration()
{
  XMLElement *calConfig = addChildElement("calculation_configuration", doc);

  XMLElement *time_step = addChildElement("time_step", calConfig);
  {
    XMLElement *load_apply = addChildElement("load_apply", time_step);
    load_apply->SetAttribute("scheme", "incremental");
    XMLElement *stiff_matrix = addChildElement("stiff_matrix", time_step);
    stiff_matrix->SetAttribute("scheme", "tangent");
    XMLElement *stiff_update = addChildElement("stiff_update", time_step);
    stiff_update->SetAttribute("interval", 0);
    XMLElement *time_differential = addChildElement("time_differential", time_step);
    time_differential->SetAttribute("theta", .5);
    XMLElement *large_defomation = addChildElement("large_defomation", time_step);
    large_defomation->SetAttribute("on", false);
  }
  XMLElement *iteration = addChildElement("iteration", calConfig);
  {
    XMLElement *max_limit = addChildElement("max_limit", iteration);
    max_limit->SetAttribute("time", 100);
    XMLElement *nonlin_material = addChildElement("nonlin_material", iteration);
    nonlin_material->SetAttribute("scheme", "viscoplasticity");
    XMLElement *tolerance = addChildElement("tolerance", iteration);
    tolerance->SetAttribute("relative", 0.01);
    tolerance->SetAttribute("absolute", 0);
  }
  XMLElement *xfem = addChildElement("xfem", calConfig);
  {
    xfem->SetAttribute("active", true);
    XMLElement *crack_segment = addChildElement("crack_segment", xfem);
    {
      crack_segment->SetAttribute("default_delta_length", 0.03);
      crack_segment->SetAttribute("gauss_num", 4);
    }
    XMLElement *SIF = addChildElement("SIF", xfem);
    {
      SIF->SetAttribute("calculation_method", "Jint");
      XMLElement *integral_radius_factor = addChildElement("integral_radius_factor", SIF);
      integral_radius_factor->SetAttribute("factor", 9);
      XMLElement *q_function_shape_factor = addChildElement("q_function_shape_factor", SIF);
      q_function_shape_factor->SetAttribute("exponent", 1);
    }
  }
  XMLElement *gravity_calculation = addChildElement("gravity_calculation", calConfig);
  gravity_calculation->SetAttribute("scheme", "k0");
  XMLElement *water_density = addChildElement("water_density", calConfig);
  {
    water_density->SetAttribute("value", 10);
    water_density->SetAttribute("unit", "kN/m^3");
  }
  XMLElement *atm = addChildElement("atm", calConfig);
  {
    atm->SetAttribute("value", 100);
    atm->SetAttribute("unit", "kPa");
  }
  return;
}
///----------------------------------------------------------------------------

void gmshParser::addTopology()
{
  XMLElement *topology = addChildElement("topology", doc);
  XMLElement *nodes = addChildElement("nodes", topology);
  nodes->SetAttribute("number", 0);
  XMLElement *elements = addChildElement("elements", topology);
  elements->SetAttribute("number", 0);
  elements->SetAttribute("order", "clockwise");
  elements->SetAttribute("poly_degree", 1);
  XMLElement *element_mega_data = addChildElement("element_mega_data", topology);
  {
    XMLElement *Segment = addChildElement("Segment", element_mega_data);
    Segment->SetAttribute("gauss_num", 3);
    XMLElement *Triangle = addChildElement("Triangle", element_mega_data);
    Triangle->SetAttribute("gauss_num", 3);
    XMLElement *Quadrangle = addChildElement("Quadrangle", element_mega_data);
    Quadrangle->SetAttribute("gauss_num", 4);
    XMLElement *Tetrahedral = addChildElement("Tetrahedral", element_mega_data);
    Tetrahedral->SetAttribute("gauss_num", 4);
    XMLElement *TriPrism = addChildElement("TriPrism", element_mega_data);
    TriPrism->SetAttribute("gauss_num", 6);
    XMLElement *Octahedral = addChildElement("Octahedral", element_mega_data);
    Octahedral->SetAttribute("gauss_num", 8);
  }
  return;
}
///----------------------------------------------------------------------------

void gmshParser::addMaterials()
{
  XMLElement *materials = addChildElement("materials", doc);
  materials->SetAttribute("number", 0);

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addInitCondition()
{
  XMLElement *init_condition = addChildElement("init_condition", doc);

  XMLElement *init_crack = addChildElement("init_crack", init_condition);
  addChildComment(" <seg id=\"1\" x0=\"0\" y0=\"0\" x1=\"1\" y1=\"1\"/> ", init_crack);
  XMLElement *init_displacement = addChildElement("init_displacement", init_condition);
  addChildComment(" <node=\"1\" dof=\"x1\" value=\"some_value\"/> ", init_displacement);
  XMLElement *init_pore = addChildElement("init_pore", init_condition);
  addChildComment(" <node=\"1\" value=\"some_value\"/> ", init_pore);

  return;
}
///----------------------------------------------------------------------------

void gmshParser::addPhases()
{
  XMLElement *phases = addChildElement("phases", doc);
  phases->SetAttribute("number", 0);

  return;
}
///----------------------------------------------------------------------------
