/*=============================================================================
*
*   Copyright (C) 2020 All rights reserved.
*
*   Filename:		gmshParser.hpp
*
*   Author:		pt2.cpp
*
*   Date:		pt2.cpp
*
*   Last Editors:		pt2.cpp
*
*   Last modified:	2020-03-27 02:15
*
*   Description:		pt2.cpp
*
=============================================================================*/
#ifndef GMSHPARSER_H
#define GMSHPARSER_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <filesystem>
#include <gmsh.h>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

inline std::pair<bool, std::string> parse_name(const std::string& _src, /*{{{*/
                                               const std::string& _patternStr,
                                               const std::size_t& _index = 1) {
  std::smatch result;
  std::regex  _pattern(_patternStr.c_str(), std::regex::icase);

  if(std::regex_match(_src, result, _pattern))
    return {true, std::string(result[_index])};
  else
    return {false, std::string()};
}
inline bool ifMatch(const std::string& _src, const std::string& _patternStr) {
  std::smatch result;
  std::regex  _pattern(_patternStr.c_str(), std::regex::icase);
  return std::regex_match(_src, result, _pattern);
} /*}}}*/

struct Element { /*{{{*/
  size_t              id;
  int                 type;
  std::string         info;
  std::vector<size_t> plist;
};            /*}}}*/
struct Node { /*{{{*/
  size_t      id;
  double      x;
  double      y;
  double      z;
  std::string info;
  std::string dof;
  double      start_value;
  double      end_value;
};                     /*}}}*/
struct PhysicalGroup { /*{{{*/
  std::string            m_name;
  std::map<int, Element> m_list_elem;
  std::map<int, Node>    m_list_node;
}; /*}}}*/

template <typename T>
class gmshParser {
  using ptree = boost::property_tree::iptree;

private:
  ptree                                m_ptree;
  std::filesystem::path                m_file_name;
  int                                  m_dim;
  std::size_t                          m_elemnum;
  std::set<int>                        m_mat_list;
  std::set<int>                        m_node_list;
  std::map<std::string, PhysicalGroup> m_list_physical;

  std::map<int, std::string> ElementTypeNameMap;

  ptree toElementPtree(const Element& elem) { /*{{{*/
    ptree root;

    root.put("id", m_elemnum++);
    root.put("type", ElementTypeNameMap.at(elem.type));
    auto [mat, id] = parse_name(elem.info, ".*mat(\\d+)_.*");
    if(mat)
      root.put("mat", std::stoi(id));
    ptree plist;
    for(auto& v : elem.plist) {
      ptree tmp(std::to_string(v));
      arrayAppend(plist, tmp, "vertex");
    }
    root.add_child("nodes", plist);
    return root;
  }                                           /*}}}*/
  ptree toNodePtree(const Node& node) const { /*{{{*/
    ptree root;

    root.put("id", node.id);
    root.put("x0", node.x);
    root.put("x1", node.y);
    root.put("x2", node.z);
    return root;
  }                                                   /*}}}*/
  ptree toBoundaryNodePtree(const Node& node) const { /*{{{*/
    ptree root;

    root.put("node", node.id);
    root.put("dof", node.dof);
    root.put("type", "Linear");
    root.put("start", node.start_value);
    root.put("end", node.end_value);

    return root;
  } /*}}}*/

  void parse() { /*{{{*/
    gmsh::vectorpair PhysicalTags;
    gmsh::model::getPhysicalGroups(PhysicalTags);
    for(auto PhysicalTag : PhysicalTags) {
      std::string name;
      gmsh::model::getPhysicalName(PhysicalTag.first, PhysicalTag.second, name);
      auto p = parsePhysicalGroup(PhysicalTag.first, PhysicalTag.second);
      m_list_physical.emplace(name, p);
    }
    return;
  }                                                /*}}}*/
  PhysicalGroup parsePhysicalGroup(const int& dim, /*{{{*/
                                   const int& physical_tag) {
    PhysicalGroup p;

    gmsh::model::getPhysicalName(dim, physical_tag, p.m_name);

    int              elemOrder = 0;
    std::vector<int> EntityTags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, physical_tag, EntityTags);

    for(auto EntityTag : EntityTags) {
      std::vector<double>                   parametricCoord;
      std::vector<int>                      elementTypes;
      std::vector<std::vector<std::size_t>> elementTags;
      std::vector<std::vector<std::size_t>> nodeInElementTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags,
                                     nodeInElementTags, dim, EntityTag);
      for(std::size_t i = 0, nTypes = elementTypes.size(); i < nTypes; ++i) {
        std::string elemName;
        int         elemDim, numNodes;
        std::string elementType = ElementTypeNameMap[elementTypes[i]];
        gmsh::model::mesh::getElementProperties(elementTypes[i], elemName,
                                                elemDim, elemOrder, numNodes,
                                                parametricCoord);
        size_t num_elem = elementTags[i].size();
        for(size_t elem_id = 0; elem_id != num_elem; ++elem_id) {
          Element elem;
          /* elem.id   = m_elemnum++; */
          elem.type = elementTypes[i];
          elem.info = p.m_name;
          for(size_t j = 0; j != static_cast<size_t>(numNodes); ++j) {
            int nodeTag =
              (int)nodeInElementTags[i][elem_id * std::size_t(numNodes) + j]
              - 1;
            elem.plist.emplace_back(nodeTag);
          }
          p.m_list_elem.emplace(p.m_list_elem.size(), elem);
        }
      }
      ////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////
      std::vector<size_t> nodeTags;
      std::vector<double> coord;
      gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim,
                                  EntityTag, true);
      for(std::size_t i = 0, nn = nodeTags.size(); i != nn; ++i) {
        Node node;
        node.id              = nodeTags[i] - 1;
        node.x               = coord[i * 3 + 0];
        node.y               = coord[i * 3 + 1];
        node.z               = coord[i * 3 + 2];
        node.info            = p.m_name;
        auto [if1, dof_name] = parse_name(p.m_name, ".*_P\\d+_(.*?)_.*");
        if(if1)
          node.dof = dof_name;
        p.m_list_node.emplace(node.id, node);
      }
    }
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    auto [if1, start_str] =
      parse_name(p.m_name, "dirichlet_P\\d+_.*_.*_(.*)_.*");
    auto [if2, end_str] = parse_name(p.m_name, "dirichlet_P\\d+_.*_.*_.*_(.*)");
    if(if1) {
      double start_value = std::stod(start_str);
      double end_value   = std::stod(end_str);
      for(auto& node : p.m_list_node) {
        node.second.start_value = start_value;
        node.second.end_value   = end_value;
      }
    }
    auto [if3, start_str1] =
      parse_name(p.m_name, "neumann_P\\d+_.*_.*_(.*)_.*");
    auto [if4, end_str1] = parse_name(p.m_name, "neumann_P\\d+_.*_.*_.*_(.*)");
    if(if3) {
      double start_value = std::stod(start_str1);
      double end_value   = std::stod(end_str1);
      for(auto& elem : p.m_list_elem) {
        switch(elem.second.type) {
        case 1: {
          Node&  p0     = p.m_list_node.at(elem.second.plist.at(0));
          Node&  p1     = p.m_list_node.at(elem.second.plist.at(1));
          double length = sqrt(pow(p0.x - p1.x, 2) + pow(p0.y - p1.y, 2));
          p0.start_value += length * 0.5 * start_value;
          p0.end_value += length * 0.5 * end_value;
          p1.start_value += length * 0.5 * start_value;
          p1.end_value += length * 0.5 * end_value;
        }
        case 8: {
          Node&  p0     = p.m_list_node.at(elem.second.plist.at(0));
          Node&  p1     = p.m_list_node.at(elem.second.plist.at(1));
          Node&  p2     = p.m_list_node.at(elem.second.plist.at(2));
          double length = sqrt(pow(p0.x - p1.x, 2) + pow(p0.y - p1.y, 2));
          p0.start_value += length / 6. * start_value;
          p0.end_value += length / 6. * end_value;
          p1.start_value += length / 6. * start_value;
          p1.end_value += length / 6. * end_value;
          p2.start_value += length * 2. / 3. * start_value;
          p2.end_value += length * 2. / 3. * end_value;
        }
        default:
          break;
        }
      }
    }
    return p;
  } /*}}}*/

  void buildMainFrame() { /*{{{*/
    /* add description node */
    m_ptree.put("description.name", m_file_name.string());
    m_ptree.put("description.method", "FEM");
    m_ptree.put("description.unknown", "Displace+Biot");
    m_ptree.put("description.displace_pattern", "Standard+Signeded+Branch");
    m_ptree.put("description.equation", "Static");
    m_ptree.put("description.solver", "Implicit");
    m_ptree.put("description.load", "Incremental");

    /* add Describe(root); */
    m_ptree.put("io.log.log_level", "info");
    m_ptree.put("io.log.detail_file_name", "");
    m_ptree.put("io.breakpoint", "");
    m_ptree.put("io.output.path", "");

    /* addCalculationConfiguration(root); */
    m_ptree.put("calculation.paralleled", true);
    m_ptree.put("calculation.time_step.theta", 1);
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
    m_ptree.put("topology.shape_function_degree", 1);
    m_ptree.put("topology.elem_nodes_order", "clockwise");
    ptree gauss_num, tmp;
    tmp.put("name", "Segment");
    tmp.put("value", 3);
    arrayAppend(gauss_num, tmp, "nGP");
    tmp.clear();
    tmp.put("name", "Triangle");
    tmp.put("value", 16);
    arrayAppend(gauss_num, tmp, "nGP");
    tmp.clear();
    tmp.put("name", "Quadrangle");
    tmp.put("value", 9);
    arrayAppend(gauss_num, tmp, "nGP");
    tmp.clear();
    tmp.put("name", "Tetrahedral");
    tmp.put("value", 4);
    arrayAppend(gauss_num, tmp, "nGP");
    tmp.clear();
    tmp.put("name", "Tri_prism");
    tmp.put("value", 6);
    arrayAppend(gauss_num, tmp, "nGP");
    tmp.clear();
    tmp.put("name", "Octahedral");
    tmp.put("value", 8);
    arrayAppend(gauss_num, tmp, "nGP");
    m_ptree.add_child("topology.elem_gauss_num", gauss_num);
    m_ptree.add_child("topology.cracks", ptree{});

    /* addMaterials(root); */
    m_ptree.add_child("materials", ptree{});

    /* addInitCondition(root); */
    m_ptree.add_child("init_condition.init_displacement", ptree{});
    m_ptree.add_child("init_condition.init_pore", ptree{});

    /* addPhases(root); */
    m_ptree.add_child("phases", ptree{});
  }                /*}}}*/
  void addMesh() { /*{{{*/
    ptree nodes;
    ptree elements;
    for(auto& p : m_list_physical) {
      std::string name      = p.second.m_name;
      auto [match, mat_str] = parse_name(name, ".*mat(\\d+)_.*");
      if(match) {
        int mat_id = std::stoi(mat_str);
        auto [not_use, mat_type] =
          parse_name(name, ".*mat\\d+_([a-z|0-9|A-Z]+).*");
        addMaterial(m_ptree.get_child("materials"), mat_type, mat_id);
        for(auto& node : p.second.m_list_node)
          arrayAppend(nodes, toNodePtree(node.second), "vertex");
        for(auto& elem : p.second.m_list_elem)
          arrayAppend(elements, toElementPtree(elem.second), "elem");
      }
    }
    m_ptree.add_child("topology.nodes", nodes);
    m_ptree.add_child("topology.elements", elements);
    return;
  }                 /*}}}*/
  void addCrack() { /*{{{*/
    for(auto& p : m_list_physical) {
      std::string name        = p.second.m_name;
      auto [match, crack_str] = parse_name(name, ".*crack(\\d+)_.*");

      if(match) {
        ptree cracks = m_ptree.get_child("topology.cracks");
        ptree crack;
        ptree nodes;
        ptree elements;
        crack.put("id", std::stoi(crack_str));
        for(auto& node : p.second.m_list_node)
          arrayAppend(nodes, toNodePtree(node.second), "node");
        for(auto& elem : p.second.m_list_elem)
          arrayAppend(elements, toElementPtree(elem.second), "elem");
        crack.add_child("nodes", nodes);
        crack.add_child("elements", elements);
        arrayAppend(cracks, crack, "crack");
      }
    }
    return;
  }                    /*}}}*/
  void addBoundary() { /*{{{*/
    ptree& root = m_ptree.get_child("phases");
    for(auto& p : m_list_physical) {
      std::string name           = p.second.m_name;
      auto [match, boundary_str] = parse_name(name, ".*_P\\d+_.*");
      if(match) {
        auto [isBoundary1, _phaseIDStr] = parse_name(name, ".*P(\\d+)_.*");
        auto [isBoundary2, bc_type]     = parse_name(name, "(.*)_P\\d+_.*");
        if(isBoundary1) {
          /// find phase
          int             phase_ID = std::stoi(_phaseIDStr) - 1;
          ptree::iterator it       = root.begin();
          bool            existed  = false;
          for(; it != root.end(); ++it)
            if(it->second.get<int>("id") == phase_ID) {
              existed = true;
              break;
            }
          // creat new phase node if not exsited
          if(!existed)
            it = addPhase(phase_ID);
          ptree& boundary =
            it->second.get_child("boundary_condition." + bc_type);
          for(auto& node : p.second.m_list_node)
            arrayAppend(boundary, toBoundaryNodePtree(node.second), "bcnode");
        }
      }
    }
    return;
  } /*}}}*/

  ptree::iterator addPhase(const size_t& id) { /*{{{*/
    ptree& phases = m_ptree.get_child("phases");
    ptree  phase;
    phase.put("id", id);
    phase.put("start_time", 0);
    phase.put("total_time", 10);
    phase.put("delta_time", 10);
    phase.add_child("boundary_condition.Dirichlet", ptree());
    phase.add_child("boundary_condition.Neumann", ptree());
    std::cout << "  (Phase #" << id << " added!)\n";
    return phases.push_back(std::make_pair("", phase));
  } /*}}}*/

  void addMaterial(/*{{{*/
                   ptree&             root,
                   const std::string& type,
                   const int&         id) {
    if(m_mat_list.find(id) != m_mat_list.end())
      return;
    std::cout << "  (New " << type << " material added!)\n";
    m_mat_list.emplace(id);

    ptree _new_mat;
    _new_mat.put("id", id);
    _new_mat.put("type", type);
    _new_mat.put("volume_weight", 19);
    _new_mat.put("friction_angle", 30);
    _new_mat.put("dilation_angle", 30);
    _new_mat.put("cohesion,", 10);
    _new_mat.put("perm_x0", 1e-5);
    _new_mat.put("perm_x1", 1e-5);
    _new_mat.put("perm_x2", 1e-5);
    _new_mat.put("zeta", 0.1);
    _new_mat.put("K_IC", 2.9);
    _new_mat.put("elastic_modulus", 1e4);
    _new_mat.put("elastic_poisson", 0.3);
    if(ifMatch(type, ".*Duncan.*")) {
      _new_mat.put("Rf", 201);
      _new_mat.put("k", 202);
      _new_mat.put("n", 203);
      _new_mat.put("G", 204);
      _new_mat.put("F", 205);
      _new_mat.put("D", 206);
    }
    else if(ifMatch(type, ".*CamClay.*")) {
      _new_mat.put("swelling_slope", 301);
      _new_mat.put("normal_compress_slope", 302);
      _new_mat.put("initial_void", 303);
    }
    else if(ifMatch(type, ".*DruckPrager.*")) {
      _new_mat.put("qf", 401);
      _new_mat.put("kf", 402);
      _new_mat.put("q_dilantion", 403);
    }
    arrayAppend(root, _new_mat, "mat");

    return;
  } /*}}}*/

  std::string getSuffix() const { /*{{{*/
    return std::string();
  }                                                            /*}}}*/
  void writeFile(const std::filesystem::path&, const ptree&) { /*{{{*/
    throw "not implemented yet";
  } /*}}}*/

public:
  /// default constructor
  gmshParser()
      : m_ptree(),
        m_file_name(),
        m_dim(0),
        m_elemnum(0),
        m_mat_list(),
        m_node_list() {
    ElementTypeNameMap = {{1, "segment"},     {8, "segment"},      //
                          {2, "triangle"},    {9, "triangle"},     //
                          {3, "quadrangle"},  {16, "quadrangle"},  //
                          {4, "tetrahedron"},                      //
                          {5, "hexahedron"},                       //
                          {6, "prism"},                            //
                          {7, "pyramid"}};
  }
  /// constructor with file name
  gmshParser(const std::filesystem::path& _fileName)
      : m_ptree(),
        m_file_name(_fileName),
        m_dim(0),
        m_elemnum(0),
        m_mat_list(),
        m_node_list() {}
  /// default destructor
  ~gmshParser() = default;

  void setFileName(const std::filesystem::path& _fileName) {
    m_file_name = _fileName;
  }

  void arrayAppend(ptree&, const ptree&, const std::string& = "item") const {
    throw "not implemented yet";
  }

  void parseGmshFile(const std::filesystem::path& file_name_with_path) { /*{{{*/
    m_file_name = file_name_with_path.stem();
    gmsh::initialize();
    gmsh::open(file_name_with_path);
    std::cout << "------------------------------------\n"
              << "gmsh file \"" << file_name_with_path << "\" imported.\n";
    m_dim = gmsh::model::getDimension();
    gmsh::model::mesh::renumberNodes();
    gmsh::model::mesh::renumberElements();
    parse();
    ////////////////////////////////////////////////////////////////

    buildMainFrame();
    std::cout << "1. " << getSuffix() << " file created.\n";
    addMesh();
    std::cout << "2. mesh information parsed.\n";
    addBoundary();
    std::cout << "3. boundary information parsed.\n";
    // addCrack();

    auto [_bool1, raw_name] = parse_name(file_name_with_path, "(.*)\\..{3}$");
    writeFile(raw_name, m_ptree);
    std::cout << "4." << getSuffix() << " file has been written.\n";
  } /*}}}*/
};
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

struct gmsh2json { /*{{{*/
  template <typename P>
  void addAttribute(boost::property_tree::ptree& _root,
                    const std::string&           _name,
                    const P&                     _value) const {
    _root.put(_name, _value);
  }
};                /*}}}*/
struct gmsh2xml { /*{{{*/
  template <typename P>
  void addAttribute(boost::property_tree::ptree& _root,
                    const std::string&           _name,
                    const P&                     _value) const {
    _root.add("<xmlattr>." + _name, _value);
  }
}; /*}}}*/
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <>
void gmshParser<gmsh2json>::arrayAppend(ptree&       _root, /*{{{*/
                                        const ptree& _child,
                                        const std::string&) const {
  _root.push_back(std::make_pair("", _child));
} /*}}}*/
template <>
void gmshParser<gmsh2xml>::arrayAppend(ptree&             _root, /*{{{*/
                                       const ptree&       _child,
                                       const std::string& _name) const {
  _root.add_child(_name, _child);
} /*}}}*/
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <>
std::string gmshParser<gmsh2json>::getSuffix() const { /*{{{*/
  return "json";
} /*}}}*/
template <>
std::string gmshParser<gmsh2xml>::getSuffix() const { /*{{{*/
  return "xml";
} /*}}}*/
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

template <>
void gmshParser<gmsh2json>::writeFile(/*{{{*/
                                      const std::filesystem::path&
                                                   file_name_with_path,
                                      const ptree& _ptree) {
  boost::property_tree::write_json(
    file_name_with_path.string() + "." + getSuffix(), _ptree);
} /*}}}*/
template <>
void gmshParser<gmsh2xml>::writeFile(/*{{{*/
                                     const std::filesystem::path&
                                                  file_name_with_path,
                                     const ptree& _ptree) {
  ptree root;
  root.add_child("geoxfem_input", _ptree);
  boost::property_tree::write_xml(
    file_name_with_path.string() + "." + getSuffix(), root);
} /*}}}*/

#endif /* GMSHPARSER_H */
