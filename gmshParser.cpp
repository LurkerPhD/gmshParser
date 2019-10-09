#include "gmshParser.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>

class gmshParser::impl
{
private:
  ///
  XMLDocument *m_doc;

  std::string m_fileName;

  int m_dim;

public:
  impl() : m_doc(nullptr), m_fileName(""), m_dim(0) {}
  ~impl() {}

  std::string getFileName() const { return m_fileName; }
  void setFileName(const std::string &_name_with_path)
  {
    std::string temp = _name_with_path;
    const std::size_t start = temp.rfind('/') + 1;
    const std::size_t end = temp.rfind('.');
    m_fileName = temp.substr(start, end - start);
  }

  int getDimension() const { return m_dim; }
  void setDimension(const int &_value) { m_dim = _value; }

  XMLDocument *getDocument() const { return m_doc; }

  void addNodes(XMLNode *root)
  {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

    for (int inode : nodeTags)
    {
      int node_num = root->ToElement()->IntAttribute("number");
      XMLElement *new_node = addChildElement("node", root);
      new_node->SetAttribute("id", node_num);
      new_node->SetAttribute("x1", coord[(inode - 1) * 3 + 0]);
      new_node->SetAttribute("x2", coord[(inode - 1) * 3 + 1]);
      new_node->SetAttribute("x3", coord[(inode - 1) * 3 + 2]);
      root->ToElement()->SetAttribute("number", ++node_num);
    }

    return;
  }
  ///----------------------------------------------------------------------------

  void addElements(XMLNode *root, const int &_dim, const std::string &_phy_name = "mat")
  {
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeInElementTags;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, _dim);

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
        std::cout << "Entity" << x.second << "has no PhyTags!";
        std::cin.get();
        break;
      }
      case 1:
      {
        elementTypes.clear();
        elementTags.clear();
        nodeInElementTags.clear();
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeInElementTags, x.first, x.second);
        gmsh::model::getPhysicalName(_dim, PhyTags.at(0), PhyName);
        if (PhyName.find(_phy_name) == std::string::npos)
          return;
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
          addMaterial(root->GetDocument()->FirstChildElement("materials"), tokens, mat_ID);
        }
        int elemNum = 0;
        for (int i = 0; i < elementTypes.size(); i++)
        {
          std::string elemName;
          int order;
          int elem_dim;
          int numNodes;
          std::vector<double> parametricCoord;
          gmsh::model::mesh::getElementProperties(elementTypes[i], elemName, elem_dim, order, numNodes, parametricCoord);
          for (int j = 0; j < elementTags[i].size(); j++)
          {
            XMLElement *elem = addChildElement("elem", root);
            elem->SetAttribute("id", elemNum);
            std::string type;
            if (elem_dim == 2)
              switch (numNodes / order)
              {
              case 2:
                type = "Segment";
                break;
              case 3:
                type = "Triangle";
                break;
              case 4:
                type = "Quadrangle";
                break;
              default:
                break;
              }
            else if (elem_dim == 3)
              switch (numNodes / order)
              {
              case 4:
                type = "Tetrahedral";
                break;
              case 6:
                type = "TriPrism";
                break;
              case 8:
                type = "Octahedral";
                break;

              default:
                break;
              }

            elem->SetAttribute("type", type.c_str());
            // insertElemMetaData(type.c_str(), element_mega_data);
            elem->SetAttribute("mat", mat_ID);
            for (int k = 0; k < numNodes; k++)
              elem->SetAttribute(("v" + std::to_string(k)).c_str(), (int)(nodeInElementTags[i][j * numNodes + k] - 1));
            elemNum++;
            root->ToElement()->SetAttribute("number", elemNum);
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

  void addCrack(XMLNode *cracks)
  {
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeInElementTags;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, getDimension() - 1);

    for (auto x : dimTags)
    {
      std::vector<int> PhyTags;
      std::string PhyName;
      gmsh::model::getPhysicalGroupsForEntity(x.first, x.second, PhyTags);

      elementTypes.clear();
      elementTags.clear();
      nodeInElementTags.clear();
      gmsh::model::mesh::getElements(elementTypes, elementTags, nodeInElementTags, x.first, x.second);
      gmsh::model::getPhysicalName(getDimension() - 1, PhyTags.at(0), PhyName);
      if (PhyName.find("crack") == std::string::npos)
        continue;

      XMLElement *crack = addChildElement("crack", cracks);
      int crack_num = cracks->ToElement()->IntAttribute("number");
      crack->SetAttribute("id", crack_num);
      cracks->ToElement()->SetAttribute("number", crack_num + 1);
      int elemNum = 0;
      for (int i = 0; i < elementTypes.size(); i++)
      {
        std::string elemName;
        int order;
        int elem_dim;
        int numNodes;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getElementProperties(elementTypes[i], elemName, elem_dim, order, numNodes, parametricCoord);
        for (int j = 0; j < elementTags[i].size(); j++)
        {
          XMLElement *elem = addChildElement("elem", crack);
          elem->SetAttribute("id", elemNum);
          std::string type;
          if (elem_dim == 1)
            type = "Segment";
          else if (elem_dim == 2)
            switch (numNodes / order)
            {
            case 3:
              type = "Triangle";
              break;
            case 4:
              type = "Quadrangle";
              break;
            default:
              break;
            }
          else if (elem_dim == 3)
            switch (numNodes / order)
            {
            case 4:
              type = "Tetrahedral";
              break;
            case 6:
              type = "TriPrism";
              break;
            case 8:
              type = "Octahedral";
              break;

            default:
              break;
            }

          elem->SetAttribute("type", type.c_str());
          for (int k = 0; k < numNodes; k++)
            elem->SetAttribute(("v" + std::to_string(k)).c_str(), (int)(nodeInElementTags[i][j * numNodes + k] - 1));
          elemNum++;
        }
      }
    }

    return;
  }
  ///----------------------------------------------------------------------------

  void addBoundaries(XMLNode *root)
  {
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);

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
        int phase_ID = std::stoi(PhyName.substr(P_pos + 1, 1)) - 1;
        XMLElement *phase = nullptr;
        for (auto _phase = root->FirstChildElement("phase"); _phase; _phase = phase->NextSiblingElement("phase"))
          if (_phase->IntAttribute("id") == phase_ID)
            phase = _phase;
        if (phase == nullptr)
        {
          phase = addChildElement("phase", root);
          int phase_num = root->ToElement()->IntAttribute("number");
          root->ToElement()->SetAttribute("number", ++phase_num);
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

  void initialPhase(XMLNode *root, const int &id)
  {
    root->ToElement()->SetAttribute("id", id);
    root->ToElement()->SetAttribute("start_time", 0);
    root->ToElement()->SetAttribute("total_time", 1000);
    root->ToElement()->SetAttribute("delta_time", 10);
    root->ToElement()->SetAttribute("element_num", 9);
    XMLElement *newBC = addChildElement("boundary_condition", root);
    newBC->SetAttribute("Dirichlet", 0);
    newBC->SetAttribute("Neumann", 0);

    return;
  }
  ///----------------------------------------------------------------------------

  void insertDirichlet(XMLNode *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
  {
    for (int bcNode : nodeTags)
    {
      XMLElement *newDirichlet = addChildElement("Dirichlet", root);
      int Dirichlet = root->ToElement()->IntAttribute("Dirichlet");
      root->ToElement()->SetAttribute("Dirichlet", ++Dirichlet);
      newDirichlet->SetAttribute("node", bcNode - 1);
      newDirichlet->SetAttribute("dof", tokens[2].c_str());
      newDirichlet->SetAttribute("type", "linear");
      newDirichlet->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
      newDirichlet->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
    }

    return;
  }
  ///----------------------------------------------------------------------------

  void insertNeumann(XMLNode *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
  {
    for (int bcNode : nodeTags)
    {
      XMLElement *newNeumann = addChildElement("Neumann", root);
      int Neumann = root->ToElement()->IntAttribute("Neumann");
      root->ToElement()->SetAttribute("Neumann", ++Neumann);
      newNeumann->SetAttribute("node", bcNode - 1);
      newNeumann->SetAttribute("dof", tokens[2].c_str());
      newNeumann->SetAttribute("type", "linear");
      newNeumann->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
      newNeumann->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
    }

    return;
  }
  ///----------------------------------------------------------------------------

  void addMaterial(XMLNode *root, std::vector<std::string> &tokens, const int &id)
  {
    for (auto _mat = root->FirstChildElement("mat"); _mat; _mat = _mat->NextSiblingElement("mat"))
      if (id == _mat->IntAttribute("id"))
        return;

    XMLElement *_mat = addChildElement("mat", root);
    int mat_id = std::stoi(root->ToElement()->Attribute("number"));
    root->ToElement()->SetAttribute("number", ++mat_id);
    std::string type = tokens[1];
    _mat->SetAttribute("constitutive", type.c_str());

    XMLElement *volume_weight = addChildElement("volume_weight", _mat);
    volume_weight->SetAttribute("value", 19);
    volume_weight->SetAttribute("unit", "kN/m^3");
    XMLElement *friction_angle = addChildElement("friction_angle", _mat);
    friction_angle->SetAttribute("value", 30);
    friction_angle->SetAttribute("unit", "degree");
    XMLElement *dilation_angle = addChildElement("dilation_angle", _mat);
    dilation_angle->SetAttribute("value", 30);
    dilation_angle->SetAttribute("unit", "degree");
    XMLElement *cohesion = addChildElement("cohesion", _mat);
    cohesion->SetAttribute("value", 10);
    cohesion->SetAttribute("unit", "kPa");
    XMLElement *perm_x1 = addChildElement("perm_x1", _mat);
    perm_x1->SetAttribute("value", 1e-5);
    perm_x1->SetAttribute("unit", "m/day");
    XMLElement *perm_x2 = addChildElement("perm_x2", _mat);
    perm_x2->SetAttribute("value", 1e-5);
    perm_x2->SetAttribute("unit", "m/day");
    XMLElement *perm_x3 = addChildElement("perm_x3", _mat);
    perm_x3->SetAttribute("value", 1e-5);
    perm_x3->SetAttribute("unit", "m/day");
    XMLElement *zeta = addChildElement("zeta", _mat);
    zeta->SetAttribute("value", 0.1);
    zeta->SetAttribute("unit", "-");
    XMLElement *K_IC = addChildElement("K_IC", _mat);
    K_IC->SetAttribute("value", 2.9);
    K_IC->SetAttribute("unit", "kN/m^0.5");
    XMLElement *elastic_modulus = addChildElement("elastic_modulus", _mat);
    elastic_modulus->SetAttribute("value", 10000);
    elastic_modulus->SetAttribute("unit", "kPa");
    XMLElement *elastic_poisson = addChildElement("elastic_poisson", _mat);
    elastic_poisson->SetAttribute("value", 0.3);
    elastic_poisson->SetAttribute("unit", "-");
    if (type.find("DuncanChang") != std::string::npos)
    {
      XMLElement *Rf = addChildElement("Rf", _mat);
      Rf->SetAttribute("value", 201);
      Rf->SetAttribute("unit", "-");
      XMLElement *k = addChildElement("k", _mat);
      k->SetAttribute("value", 202);
      k->SetAttribute("unit", "kPa");
      XMLElement *n = addChildElement("n", _mat);
      n->SetAttribute("value", 203);
      n->SetAttribute("unit", "-");
      XMLElement *G = addChildElement("G", _mat);
      G->SetAttribute("value", 204);
      G->SetAttribute("unit", "-");
      XMLElement *F = addChildElement("F", _mat);
      F->SetAttribute("value", 205);
      F->SetAttribute("unit", "-");
      XMLElement *D = addChildElement("D", _mat);
      D->SetAttribute("value", 206);
      D->SetAttribute("unit", "-");
    }
    else if (type.find("ModifiedCamClay") != std::string::npos)
    {
      XMLElement *swelling_slope = addChildElement("swelling_slope", _mat);
      swelling_slope->SetAttribute("value", 401);
      swelling_slope->SetAttribute("unit", "-");
      XMLElement *normal_compress_slope = addChildElement("normal_compress_slope", _mat);
      normal_compress_slope->SetAttribute("value", 402);
      normal_compress_slope->SetAttribute("unit", "-");
      XMLElement *initial_void = addChildElement("initial_void", _mat);
      initial_void->SetAttribute("value", 403);
      initial_void->SetAttribute("unit", "-");
    }
    else if (type.find("DruckPrager") != std::string::npos)
    {
      XMLElement *qf = addChildElement("qf", _mat);
      qf->SetAttribute("value", 501);
      qf->SetAttribute("unit", "-");
      XMLElement *kf = addChildElement("kf", _mat);
      kf->SetAttribute("value", 502);
      kf->SetAttribute("unit", "-");
      XMLElement *q_dilantion = addChildElement("q_dilantion", _mat);
      q_dilantion->SetAttribute("value", 503);
      q_dilantion->SetAttribute("unit", "degree");
    }

    return;
  };

  ///----------------------------------------------------------------------------

  void buildMainFrame()
  {
    m_doc = new XMLDocument();
    m_doc->Parse("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");

    addChildComment("1. model type information", m_doc);
    addDescribe(m_doc);
    addChildComment("2. input information", m_doc);
    addBreakPoint(m_doc);

    addChildComment("3. calculation configuration information", m_doc);
    addCalculationConfiguration(m_doc);

    addChildComment("4. model's topology information", m_doc);
    addTopology(m_doc);

    addChildComment("5. s information", m_doc);
    addMaterials(m_doc);

    addChildComment("6. initial condition information", m_doc);
    addInitCondition(m_doc);

    addChildComment("7. phase step information", m_doc);
    addPhases(m_doc);

    addChildComment("8. output information", m_doc);
    addOutput(m_doc);
  }
  ///----------------------------------------------------------------------------

  XMLElement *addChildElement(const std::string &_name, XMLNode *root)
  {
    XMLElement *child = m_doc->NewElement(_name.c_str());
    root->InsertEndChild(child);

    return child;
  }
  ///----------------------------------------------------------------------------

  XMLComment *addChildComment(const std::string &_content, XMLNode *root)
  {
    XMLComment *child = m_doc->NewComment(_content.c_str());
    root->InsertEndChild(child);

    return child;
  }
  ///----------------------------------------------------------------------------

  XMLText *addChildText(const std::string &_content, XMLNode *root)
  {
    XMLText *child = m_doc->NewText(_content.c_str());
    root->InsertEndChild(child);

    return child;
  }
  ///----------------------------------------------------------------------------

  void addDescribe(XMLNode *root)
  {
    /// root node Describe
    XMLElement *describe = addChildElement("describe", root);
    describe->SetAttribute("name", m_fileName.c_str());
    describe->SetAttribute("method", "FEM");
    describe->SetAttribute("unknown", "Displace+Biot");
    describe->SetAttribute("displace_pattern", "Standard+Heaviside+Branch");
    describe->SetAttribute("equation", "Static");
    describe->SetAttribute("solver", "Implicit");
    describe->SetAttribute("load", "Incremental");
    if (m_dim == 2)
      describe->SetAttribute("plane_stress", false);

    return;
  }
  ///----------------------------------------------------------------------------

  void addBreakPoint(XMLNode *root)
  {
    XMLElement *breakpoint = addChildElement("breakpoint", root);
    breakpoint->SetAttribute("file", "");

    return;
  }
  ///----------------------------------------------------------------------------

  void addOutput(XMLNode *root)
  {
    XMLElement *output = addChildElement("output", root);
    output->SetAttribute("path", ("./result/" + m_fileName).c_str());
    output->SetAttribute("format", "vtk");
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

  void addCalculationConfiguration(XMLNode *root)
  {
    XMLElement *calConfig = addChildElement("calculation_configuration", root);

    XMLElement *time_step = addChildElement("time_step", calConfig);
    {
      XMLElement *stiff_update = addChildElement("stiff_update", time_step);
      stiff_update->SetAttribute("interval", 0);
      XMLElement *time_differential = addChildElement("time_differential", time_step);
      time_differential->SetAttribute("theta", .8);
      XMLElement *large_defomation = addChildElement("large_defomation", time_step);
      large_defomation->SetAttribute("on", false);
    }
    XMLElement *iteration = addChildElement("iteration", calConfig);
    {
      XMLElement *max_limit = addChildElement("max_limit", iteration);
      max_limit->SetAttribute("time", 100);
      XMLElement *nonlinear = addChildElement("nonlinear", iteration);
      nonlinear->SetAttribute("scheme", "initial_stress");
      XMLElement *tolerance = addChildElement("tolerance", iteration);
      tolerance->SetAttribute("relative", 0.01);
      tolerance->SetAttribute("absolute", 0);
    }
    XMLElement *xfem = addChildElement("xfem", calConfig);
    xfem->SetAttribute("active", "true");
    {
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

  void addTopology(XMLNode *root)
  {
    XMLElement *topology = addChildElement("topology", root);
    XMLElement *nodes = addChildElement("nodes", topology);
    nodes->SetAttribute("number", 0);
    XMLElement *elements = addChildElement("elements", topology);
    elements->SetAttribute("number", 0);
    elements->SetAttribute("order", "clockwise");
    elements->SetAttribute("poly_degree", 1);
    XMLElement *cracks = addChildElement("cracks", topology);
    cracks->SetAttribute("number", 0);
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

  void addMaterials(XMLNode *root)
  {
    XMLElement *materials = addChildElement("materials", root);
    materials->SetAttribute("number", 0);

    return;
  }
  ///----------------------------------------------------------------------------

  void addInitCondition(XMLNode *root)
  {
    XMLElement *init_condition = addChildElement("init_condition", root);

    XMLElement *init_displacement = addChildElement("init_displacement", init_condition);
    addChildComment(" <node=\"1\" dof=\"x1\" value=\"some_value\"/> ", init_displacement);
    XMLElement *init_pore = addChildElement("init_pore", init_condition);
    addChildComment(" <node=\"1\" value=\"some_value\"/> ", init_pore);

    return;
  }
  ///----------------------------------------------------------------------------

  void addPhases(XMLNode *root)
  {
    XMLElement *phases = addChildElement("phases", root);
    phases->SetAttribute("number", 0);

    return;
  }
};

gmshParser::gmshParser()
    : m_pImpl(std::make_unique<impl>())
{
}
gmshParser::~gmshParser() {}
///----------------------------------------------------------------------------

void gmshParser::readFromGmshFile(const std::string &_name_with_path)
{
  m_pImpl->setFileName(_name_with_path);
  ///----------------------------------------------------------------------------
  gmsh::initialize();
  gmsh::open(_name_with_path);
  std::cout << "------------------------------------\n"
            << "gmsh file \"" << _name_with_path << "\" imported.\n";
  m_pImpl->setDimension(gmsh::model::getDimension());
  gmsh::model::mesh::renumberNodes();
  gmsh::model::mesh::renumberElements();
  ///----------------------------------------------------------------------------

  m_pImpl->buildMainFrame();
  std::cout << "1. xml file created.\n";
  ///----------------------------------------------------------------------------

  /// get all nodes
  m_pImpl->addNodes(m_pImpl->getDocument()->FirstChildElement("topology")->FirstChildElement("nodes"));
  std::cout << "2. node information parsed.\n";
  ///----------------------------------------------------------------------------

  //get all elements
  m_pImpl->addElements(m_pImpl->getDocument()->FirstChildElement("topology")->FirstChildElement("elements"), m_pImpl->getDimension());
  std::cout << "3. element information parsed.\n";
  ///----------------------------------------------------------------------------

  /// get crack mesh
  m_pImpl->addCrack(m_pImpl->getDocument()->FirstChildElement("topology")->FirstChildElement("cracks"));
  std::cout << "4. crack information parsed.\n";
  ///----------------------------------------------------------------------------

  ///get all boundary condition information
  m_pImpl->addBoundaries(m_pImpl->getDocument()->FirstChildElement("phases"));
  std::cout << "5. boundary information parsed.\n";
  ///----------------------------------------------------------------------------

  const std::string raw_name = _name_with_path.substr(0, _name_with_path.find("."));
  m_pImpl->getDocument()->SaveFile((raw_name + ".xml").c_str());
  std::cout << "6. xml file has been written.\n";
}
///----------------------------------------------------------------------------

///----------------------------------------------------------------------------
