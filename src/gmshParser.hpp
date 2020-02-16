/*=============================================================================
*
*   Copyright (C) 2020 All rights reserved.
*
*   Filename:		gmshParser.h
*
*   Author:		pt2.cpp
*
*   Date:		pt2.cpp
*
*   Last Editors:		pt2.cpp
*
*   Last modified:	2020-02-16 19:41
*
*   Description:		pt2.cpp
*
=============================================================================*/
#ifndef GMSHPARSER_H
#define GMSHPARSER_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <gmsh.h>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

inline std::pair<bool, std::string> parseName(const std::string &_src,
                                              const std::string &_patternStr,
                                              const std::size_t &_index = 1) {
  std::smatch result;
  std::regex _pattern(_patternStr.c_str(), std::regex::icase);

  if (std::regex_match(_src, result, _pattern))
    return {true, std::string(result[_index])};
  else
    return {false, std::string()};
}
inline bool ifMatch(const std::string &_src, const std::string &_patternStr) {
  std::smatch result;
  std::regex _pattern(_patternStr.c_str(), std::regex::icase);
  return std::regex_match(_src, result, _pattern);
}

template <typename T> class gmshParser {
  using pt = boost::property_tree::ptree;

private:
  pt m_ptree;
  std::string m_fileName;
  int m_dim;
  std::size_t m_elemnum;
  std::set<int> m_mat_list;
  std::set<int> m_node_list;

  template <typename P>
  void arrayAppend(pt &_root, const std::string &_name, const P &_value) {
    pt child;
    child.put(_name, _value);
    _root.push_back(std::make_pair("", child));
  }

  void buildMainFrame() {
    /* add description node */
    m_ptree.put("description.name", m_fileName.c_str());
    m_ptree.put("description.method", "FEM");
    m_ptree.put("description.unknown", "Displace+Biot");
    m_ptree.put("description.displace_pattern", "Standard+Signeded+Branch");
    m_ptree.put("description.equation", "Static");
    m_ptree.put("description.solver", "Implicit");
    m_ptree.put("description.load", "Incremental");
    if (m_dim == 2)
      m_ptree.put("description.plane_stress", false);

    /* add Describe(root); */
    m_ptree.put("io.log.log_level", "info");
    m_ptree.put("io.log.detail_file_name", "");
    m_ptree.put("io.breakpoint", "");
    m_ptree.put("io.output.file_name", "");

    /* addCalculationConfiguration(root); */
    m_ptree.put("calculation.paralleled", true);
    m_ptree.put("calculation.time_step.theta", 0.8);
    m_ptree.put("calculation.time_step.large_defomation", false);
    m_ptree.put("calculation.iteration.max_time", 100);
    m_ptree.put("calculation.iteration.nonlinear_scheme", "initial_stress");
    m_ptree.put("calculation.iteration.tolerance.relative", 0.01);
    m_ptree.put("calculation.iteration.tolerance.absolute", 0);
    m_ptree.put("calculation.xfem.active", true);
    m_ptree.put("calculation.xfem.branch_enrich_radius", 0);
    m_ptree.put("calculation.xfem.crack_segment.default_delta_length", 0.03);
    m_ptree.put("calculation.xfem.crack_segment.gauss_num", 4);
    m_ptree.put("calculation.xfem.sif_integral.domain.shape", "circle");
    m_ptree.put("calculation.xfem.sif_integral.domain.radius", 9);
    m_ptree.put("calculation.xfem.sif_integral.qfunction.exponent", 1);
    m_ptree.put("calculation.gravity_scheme", "k0");
    m_ptree.put("calculation.water_density.value", 10);
    m_ptree.put("calculation.water_density.unit", "kN\\m^3");
    m_ptree.put("calculation.atm.value", 100);
    m_ptree.put("calculation.atm.unit", "kPa");

    /* addTopology(root); */
    m_ptree.put("topology.element_degree", 1);
    m_ptree.put("topology.element_order", "clockwise");
    pt gauss_num;
    arrayAppend(gauss_num, "segment", 3);
    arrayAppend(gauss_num, "triangle", 3);
    arrayAppend(gauss_num, "quadrangle", 4);
    arrayAppend(gauss_num, "tetrahedral", 4);
    arrayAppend(gauss_num, "tri_prism", 6);
    arrayAppend(gauss_num, "octahedral", 8);
    m_ptree.add_child("topology.gauss_num", gauss_num);
    m_ptree.put("topology.cracks", "");

    /* addMaterials(root); */
    m_ptree.put("materials", "");

    /* addInitCondition(root); */
    m_ptree.put("init_condition.init_displacement", "");
    m_ptree.put("init_condition.init_pore", "");

    /* addPhases(root); */
    m_ptree.put("phases", "");
  }
  void addMesh(const int &_dim = -1) {
    gmsh::vectorpair PhysicalTags;
    gmsh::model::getPhysicalGroups(PhysicalTags, _dim);
    std::vector<int> _listElement;

    for (auto PhysicalTag : PhysicalTags) {
      int mat_ID;
      std::string PhyName;
      gmsh::model::getPhysicalName(PhysicalTag.first, PhysicalTag.second,
                                   PhyName);

      auto [_isMesh, _matIDStr] = parseName(PhyName, ".*mat(\\d+)_.*");
      auto [_isCrack, _crackIDStr] = parseName(PhyName, ".*crack(\\d*).*");
      if (_isMesh) {
        mat_ID = std::stoi(_matIDStr);
        auto [_bool, _matName] =
            parseName(PhyName, ".*mat\\d+_([a-z|0-9|A-Z]+).*");
        addMaterial(m_ptree.get_child("materials"), _matName, mat_ID);
        addMeshForEntity(PhysicalTag, m_ptree.get_child("topology"), mat_ID);
      } else if (_isCrack) {
        /* XMLElement *crack = nullptr; */
        /* XMLElement *cracks = root->FirstChildElement("cracks"); */
        /* for (XMLElement *_crack = cracks->FirstChildElement("crack"); _crack;
         */
        /*      _crack = _crack->NextSiblingElement("crack")) */
        /*   if (_crack->IntAttribute("id") == std::stoi(_crackIDStr)) { */
        /*     crack = _crack; */
        /*     break; */
        /*   } */
        /* if (!crack) { */
        /*   crack = addChildElement("crack", cracks); */
        /*   crack->SetAttribute("id", std::stoi(_crackIDStr) - 1); */
        /*   addChildElement("nodes", crack); */
        /*   addChildElement("elements", crack); */
        /* } */
        /* addMeshForEntity(PhysicalTag, crack); */
      }
    }
  }
  void addBoundaries(pt &root) {
    gmsh::vectorpair PhysicalTags;
    gmsh::model::getPhysicalGroups(PhysicalTags);

    for (auto PhysicalTag : PhysicalTags) {
      std::string PhyName;
      std::vector<std::size_t> nodeTags;
      std::vector<double> coords;
      gmsh::model::getPhysicalName(PhysicalTag.first, PhysicalTag.second,
                                   PhyName);
      /// parsing physical names. sth like "Dirichlet_P1_x1_1"
      auto [_isBoundary1, _boundaryType] = parseName(PhyName, "(.*)_P\\d+_.*");
      auto [_isBoundary2, _phaseIDStr] = parseName(PhyName, ".*P(\\d+)_.*");
      auto [_isBoundary3, _dofPosStr] = parseName(PhyName, ".*P\\d+_(.*?)_.*");
      auto [_isBoundary4, _boundaryIDStr] =
          parseName(PhyName, ".*P\\d+_.*?_(\\d+?).*");
      if (_isBoundary1) {
        /// find phase
        int phase_ID = std::stoi(_phaseIDStr) - 1;
        pt::iterator it = root.begin();
        bool existed = false;
        for (; it != root.end(); ++it)
          if (it->second.get<int>("id") == phase_ID) {
            existed = true;
            break;
          }
        if (!existed) {
          pt phase;
          std::cout << "  (Phase " << phase_ID << " added!)\n";
          phase.put("id", phase_ID);
          phase.put("start_time", 0);
          phase.put("total_time", 1000);
          phase.put("delta_time", 10);
          phase.put("element_num", 9);
          phase.add_child("boundary_condition.Dirichlet", pt());
          phase.add_child("boundary_condition.Neumann", pt());
          it = root.push_back(std::make_pair("", phase));
        }
        std::cout << _boundaryType << std::endl;
        bool ifNeumann = ifMatch(_boundaryType, ".*neumann.*");
        bool ifDirichlet = ifMatch(_boundaryType, ".*dirichlet.*");
        gmsh::model::mesh::getNodesForPhysicalGroup(
            PhysicalTag.first, PhysicalTag.second, nodeTags, coords);
        for (std::size_t bcNode : nodeTags) {
          pt newBC;
          newBC.put("node", bcNode);
          newBC.put("dof", _dofPosStr);
          newBC.put("type", "linear");
          newBC.put("start", PhyName.substr(0, 1) + "_start");
          newBC.put("end", PhyName.substr(0, 1) + "_end");
          if (ifDirichlet) {
            it->second.get_child("boundary_condition.Dirichlet")
                .push_back(std::make_pair("", newBC));
          } else if (ifNeumann) {
            it->second.get_child("boundary_condition.Neumann")
                .push_back(std::make_pair("", newBC));
          } else
            std::cout << "error boundary condition type!\n";
        }
      }
    }
  }

  void addMaterial(pt &root, const std::string &type, const int &id) {
    if (m_mat_list.find(id) != m_mat_list.end())
      return;
    std::cout << "  (New " << type << " material added!)\n";
    m_mat_list.emplace(id);

    pt _new_mat;
    _new_mat.put("id", id);
    _new_mat.put("type", type);
    _new_mat.put("volume_weight", 19);
    _new_mat.put("friction_angle", 30);
    _new_mat.put("dilation_angle", 30);
    _new_mat.put("cohesion,", 10);
    _new_mat.put("perm_x1", 1e-5);
    _new_mat.put("perm_x2", 1e-5);
    _new_mat.put("perm_x3", 1e-5);
    _new_mat.put("zeta", 0.1);
    _new_mat.put("K_IC", 2.9);
    _new_mat.put("elastic_modulus", 1e4);
    _new_mat.put("elastic_poisson", 0.3);
    if (ifMatch(type, ".*Duncan.*")) {
      _new_mat.put("Rf", 201);
      _new_mat.put("k", 202);
      _new_mat.put("n", 203);
      _new_mat.put("G", 204);
      _new_mat.put("F", 205);
      _new_mat.put("D", 206);
    } else if (ifMatch(type, ".*CamClay.*")) {
      _new_mat.put("swelling_slope", 301);
      _new_mat.put("normal_compress_slope", 302);
      _new_mat.put("initial_void", 303);
    } else if (ifMatch(type, ".*DruckPrager.*")) {
      _new_mat.put("qf", 401);
      _new_mat.put("kf", 402);
      _new_mat.put("q_dilantion", 403);
    }
    root.push_back(std::make_pair("", _new_mat));

    return;
  }

  void addMeshForEntity(const std::pair<int, int> &PhysicalTag, pt &root,
                        const int &mat_ID = -1) {
    std::map<int, std::string> ElementTypeNameMap = {
        {1, "segment"},     {2, "triangle"},   {3, "quadrangle"},
        {4, "tetrahedron"}, {5, "hexahedron"}, {6, "prism"},
        {7, "pyramid"}};
    pt nodes, elements;
    int elemOrder = 0;
    std::vector<int> EntityTags;
    gmsh::model::getEntitiesForPhysicalGroup(PhysicalTag.first,
                                             PhysicalTag.second, EntityTags);
    for (auto EntityTag : EntityTags) {
      std::vector<std::size_t> nodeTags;
      std::vector<double> coord;
      std::vector<double> parametricCoord;
      gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord,
                                  PhysicalTag.first, EntityTag, true);
      for (std::size_t i = 0, nn = nodeTags.size(); i != nn; ++i) {
        int nodeTag = (int)nodeTags[i];
        if (m_node_list.find(nodeTag) != m_node_list.end())
          continue;
        m_node_list.emplace(nodeTag);
        pt node;
        node.put("id", nodeTag - 1);
        node.put("x1", coord[i * 3 + 0]);
        node.put("x2", coord[i * 3 + 1]);
        node.put("x3", coord[i * 3 + 2]);
        nodes.push_back(std::make_pair("", node));
      }
      std::vector<int> elementTypes;
      std::vector<std::vector<std::size_t>> elementTags;
      std::vector<std::vector<std::size_t>> nodeInElementTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags,
                                     nodeInElementTags, PhysicalTag.first,
                                     EntityTag);
      for (std::size_t i = 0, nTypes = elementTypes.size(); i < nTypes; ++i) {
        std::string elemName;
        int elemDim, numNodes;
        std::vector<double> parametricCoord1;
        std::string elementType = ElementTypeNameMap[elementTypes[i]];
        gmsh::model::mesh::getElementProperties(elementTypes[i], elemName,
                                                elemDim, elemOrder, numNodes,
                                                parametricCoord1);
        for (std::size_t nElement = elementTags[i].size();
             m_elemnum != nElement; ++m_elemnum) {
          pt elem;
          elem.put("id", m_elemnum);
          elem.put("type", elementType);
          elem.put("mat", mat_ID);
          pt nodelist;
          for (std::size_t j = 0; j != (std::size_t)numNodes; ++j) {
            int nodeTag =
                (int)nodeInElementTags[i]
                                      [m_elemnum * std::size_t(numNodes) + j] -
                1;
            pt tmp(std::to_string(nodeTag));
            nodelist.push_back(std::make_pair("", tmp));
          }
          elem.add_child("nodes", nodelist);
          elements.push_back(std::make_pair("", elem));
        }
      }
    }
    root.add_child("nodes", nodes);
    root.add_child("elements", elements);
  }

public:
  /// default constructor
  gmshParser()
      : m_ptree(), m_fileName(), m_dim(0), m_elemnum(0), m_mat_list(),
        m_node_list() {}
  /// constructor with file name
  gmshParser(const std::string &_fileName)
      : m_ptree(), m_fileName(_fileName), m_dim(0), m_elemnum(0), m_mat_list(),
        m_node_list() {}
  /// default destructor
  ~gmshParser() = default;

  void setFileName(const std::string &_fileName) { m_fileName = _fileName; }
  std::string getFileName() const { return m_fileName; }

  void setDimension(const int &_dim) { m_dim = _dim; }
  int getDimension() const { return m_dim; }

  void parseGmshFile(const std::string &_fileNameWithPath) {
    auto [_bool, _fileName] =
        parseName(_fileNameWithPath, ".*[/|\\\\](.*).msh");
    m_fileName =
        _bool ? _fileName
              : _fileNameWithPath.substr(0, _fileNameWithPath.size() - 4);
    gmsh::initialize();
    gmsh::open(_fileNameWithPath);
    std::cout << "------------------------------------\n"
              << "gmsh file \"" << _fileNameWithPath << "\" imported.\n";
    m_dim = gmsh::model::getDimension();
    gmsh::model::mesh::renumberNodes();
    gmsh::model::mesh::renumberElements();
    ////////////////////////////////////////////////////////////////

    buildMainFrame();
    std::cout << "1. " << getSuffix() << " file created.\n";
    ////////////////////////////////////////////////////////////////

    /// get crack mesh
    addMesh();
    std::cout << "2. mesh information parsed.\n";
    ////////////////////////////////////////////////////////////////

    /// get all boundary condition information
    addBoundaries(m_ptree.get_child("phases"));
    std::cout << "3. boundary information parsed.\n";
    ////////////////////////////////////////////////////////////////

    auto [_bool1, raw_name] = parseName(_fileNameWithPath, "(.*)\\..{3}$");
    writeFile(raw_name, m_ptree);
    std::cout << "4." << getSuffix() << " file has been written.\n";
  }

  std::string getSuffix() const { return std::string(); }
  void writeFile(const std::string &, const pt &) {}
};

struct gmsh2json {};
struct gmsh2xml {};

template <> std::string gmshParser<gmsh2json>::getSuffix() const {
  return "json";
}
template <> std::string gmshParser<gmsh2xml>::getSuffix() const {
  return "xml";
}

template <>
void gmshParser<gmsh2json>::writeFile(const std::string &_fileNameWithPath,
                                      const pt &_ptree) {
  boost::property_tree::write_json(_fileNameWithPath + "." + getSuffix(),
                                   _ptree);
}
template <>
void gmshParser<gmsh2xml>::writeFile(const std::string &_fileNameWithPath,
                                     const pt &_ptree) {
  pt root;
  root.add_child("geoxfem_input", _ptree);
  boost::property_tree::write_xml(_fileNameWithPath + "." + getSuffix(), root);
}

#endif /* GMSHPARSER_H */
