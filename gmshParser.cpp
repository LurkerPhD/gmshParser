#include <algorithm>
#include <gmsh.h>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>

#include "tinyxml2.h"

using namespace tinyxml2;

int createMainFrame(XMLDocument *doc);

void insertElemMetaData(const char *typeID, XMLElement *element_mega_data);

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id);

void insertDirichlet(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags);

void insertNeumann(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags);

void insertMat(XMLElement *root, std::vector<std::string> &tokens, const int &id);

const char *elemType(const int &dim, const int &nSize);

const std::string elemNameRule(const int &dim, const int &nEdge, const char *interpolation);

XMLElement *QueryElementByAttribute(XMLElement *root, const std::string &Attri_Name, const std::string &value);

static int dimension;

int main(int argc, char *argv[])
{
  /* code */
  char *filename;
  if (argc > 1)
    filename = argv[1];
  else
    filename = "demo/test1.msh";

  gmsh::initialize();
  if (filename == nullptr)
    filename = (char *)"demo/test1.msh";
  gmsh::open(filename);
  std::cout << "gmsh file imported.\n";
  dimension = gmsh::model::getDimension();
  gmsh::model::mesh::renumberNodes();
  gmsh::model::mesh::renumberElements();

  XMLDocument *doc = new XMLDocument;
  createMainFrame(doc);
  std::cout << "xml file created.\n";
  /// get all nodes
  {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);
    /// std::cout << "number of nodes: " << nodeTags.size() << "\n";
    XMLElement *nodes = doc->FirstChildElement("topology")
                            ->FirstChildElement("nodes");
    for (int inode : nodeTags)
    {
      int nodeNum = nodes->IntAttribute("number") + 1;
      nodes->SetAttribute("number", nodeNum);
      XMLElement *newNode = doc->NewElement("node");
      nodes->InsertEndChild(newNode);
      newNode->SetAttribute("id", nodeNum);
      /// newNode->SetAttribute("type", ("femNode" + std::to_string(dimension) + "D").c_str());
      std::vector<int> temp = {(inode - 1) * 3 + 0, (inode - 1) * 3 + 1, (inode - 1) * 3 + 2};
      std::vector<double> temp2 = {coord[(inode - 1) * 3 + 0], coord[(inode - 1) * 3 + 1], coord[(inode - 1) * 3 + 2]};
      newNode->SetAttribute("x1", coord[(inode - 1) * 3 + 0]);
      newNode->SetAttribute("x2", coord[(inode - 1) * 3 + 1]);
      newNode->SetAttribute("x3", coord[(inode - 1) * 3 + 2]);
    }
  }
  std::cout << "node information parsed.\n";
  //get all elements
  {
    XMLElement *elements = doc->FirstChildElement("topology")
                               ->FirstChildElement("elements");
    /// std::cout << "step 0.\n";

    XMLElement *materials = doc->FirstChildElement("materials");
    /// std::cout << "step 1.\n";
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
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
            /// std::string type = "C" + std::to_string(dimension) + "D" + std::to_string(numNodes) + "PE";
            std::string type;
            type = "Quadrangle"; //elemNameRule(dimension, numNodes, "Lin");
            newElem->SetAttribute("type", type.c_str());
            // insertElemMetaData(type.c_str(), element_mega_data);
            newElem->SetAttribute("mat", mat_ID);
            for (int k = 0; k < numNodes; k++)
              newElem->SetAttribute(("v" + std::to_string(k + 1)).c_str(),
                                    (int)(nodeTags[i][j * numNodes + k]));
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
    XMLElement *phases = doc->FirstChildElement("phases");
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
    }
  }
  std::cout << "boundary condition information parsed.\n";
  doc->SaveFile(strcat(filename, ".xml"));
  //doc->SaveFile("demo.xml");
  std::cout << "xml file has been written.\n";
  delete doc;
  return 0;
}

void insertElemMetaData(const char *MetaDataID, XMLElement *element_mega_data)
{
  if (!element_mega_data->FirstChildElement(MetaDataID))
  {
    XMLElement *newElemMetaData = element_mega_data->GetDocument()->NewElement(MetaDataID);
    element_mega_data->InsertEndChild(newElemMetaData);
    newElemMetaData->SetAttribute("gauss_num", 4);
    /// newElemMetaData->SetAttribute("polyDegree", 2);
    /// newElemMetaData->SetAttribute("plane_stress", false);
    int n = (element_mega_data->IntAttribute("number"));
    element_mega_data->SetAttribute("number", ++n);
  }
  return;
}

int createMainFrame(XMLDocument *doc)
{
  const char *declaration =
      "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
  doc->Parse(declaration);
  XMLComment *comment1 = doc->NewComment("1. model type information");
  doc->InsertEndChild(comment1);
  /// root node Describe
  XMLElement *describe = doc->NewElement("describe");
  doc->InsertEndChild(describe);
  {
    /// XMLText *describeText =
    ///     doc->NewText("This is a example of input with xml formate");
    /// describe->InsertEndChild(describeText);
    describe->SetAttribute("name", "test1");
    describe->SetAttribute("method", "FEM");
    describe->SetAttribute("problem", "Static");
    describe->SetAttribute("couple", "Solid");
    describe->SetAttribute("equation", "Implicit");
    describe->SetAttribute("solver", "default");
    describe->SetAttribute("dim", dimension);
    if (dimension == 2)
      describe->SetAttribute("plane_stress", false);
  }

  XMLComment *comment2 = doc->NewComment("2. input/output information");
  doc->InsertEndChild(comment2);
  XMLElement *breakpoint = doc->NewElement("breakpoint");
  doc->InsertEndChild(breakpoint);
  breakpoint->SetAttribute("file", "");
  /// output node in doc
  XMLElement *output = doc->NewElement("output");
  doc->InsertEndChild(output);
  output->SetAttribute("path", "");
  output->SetAttribute("format", "vtk");
  output->SetAttribute("plot", "mesh+gauss+bc");
  output->SetAttribute("displacement_amplifier", 200);

  XMLComment *comment3 = doc->NewComment("3. calculation configuration information");
  doc->InsertEndChild(comment3);
  /// calConfig node in doc
  XMLElement *calConfig = doc->NewElement("calculation_configuration");
  doc->InsertEndChild(calConfig);
  {
    /// time_step node in calConfig
    XMLElement *time_step = doc->NewElement("time_step");
    calConfig->InsertEndChild(time_step);
    {
      /// load_apply node in time_step
      XMLElement *load_apply = doc->NewElement("load_apply");
      time_step->InsertEndChild(load_apply);
      load_apply->SetAttribute("scheme", "incremental");
      /// stiff_matrix node in time_step
      XMLElement *stiff_matrix = doc->NewElement("stiff_matrix");
      time_step->InsertEndChild(stiff_matrix);
      stiff_matrix->SetAttribute("scheme", "tangent");
      /// stiff_update node in time_step
      XMLElement *stiff_update = doc->NewElement("stiff_update");
      time_step->InsertEndChild(stiff_update);
      stiff_update->SetAttribute("interval", 0);
      /// memory_for_stiff node in time_step
      XMLElement *memory_for_stiff = doc->NewElement("memory_for_stiff");
      time_step->InsertEndChild(memory_for_stiff);
      memory_for_stiff->SetAttribute("sufficient", false);
      /// time_differential node in time_step
      XMLElement *time_differential = doc->NewElement("time_differential");
      time_step->InsertEndChild(time_differential);
      time_differential->SetAttribute("theta", .5);
      /// large_defomation node in time_step
      XMLElement *large_defomation = doc->NewElement("large_defomation");
      time_step->InsertEndChild(large_defomation);
      large_defomation->SetAttribute("on", false);
    }

    /// iteration node in calConfig
    XMLElement *iteration = doc->NewElement("iteration");
    calConfig->InsertEndChild(iteration);
    {
      /// max_limit node in iteration
      XMLElement *max_limit = doc->NewElement("max_limit");
      iteration->InsertEndChild(max_limit);
      max_limit->SetAttribute("time", 100);
      /// nonlin_material node in iteration
      XMLElement *nonlin_material = doc->NewElement("nonlin_material");
      iteration->InsertEndChild(nonlin_material);
      nonlin_material->SetAttribute("scheme", "viscoplasticity");
      /// tolerance node in iteration
      XMLElement *tolerance = doc->NewElement("tolerance");
      iteration->InsertEndChild(tolerance);
      tolerance->SetAttribute("relative", 0.01);
      tolerance->SetAttribute("absolute", 0);
    }

    /// xfem node in calConfig
    XMLElement *xfem = doc->NewElement("xfem");
    calConfig->InsertEndChild(xfem);
    {
      xfem->SetAttribute("active", true);
      /// crack_segment node in xfem
      XMLElement *crack_segment = doc->NewElement("crack_segment");
      xfem->InsertEndChild(crack_segment);
      crack_segment->SetAttribute("default_delta_length", 0.03);
      crack_segment->SetAttribute("gauss_num", 4);
      /// SIF node in xfem
      XMLElement *SIF = doc->NewElement("SIF");
      xfem->InsertEndChild(SIF);
      SIF->SetAttribute("calculation_method", "Jint");
      {
        /// integral_radius_factor node in SIF
        XMLElement *integral_radius_factor = doc->NewElement("integral_radius_factor");
        SIF->InsertEndChild(integral_radius_factor);
        integral_radius_factor->SetAttribute("factor", 9);
        /// q_function_shape_factor node in SIF
        XMLElement *q_function_shape_factor = doc->NewElement("q_function_shape_factor");
        SIF->InsertEndChild(q_function_shape_factor);
        q_function_shape_factor->SetAttribute("exponent", 1);
      }
    } /// end xfem

    /// gravity_calculation node in calConfig
    XMLElement *gravity_calculation = doc->NewElement("gravity_calculation");
    calConfig->InsertEndChild(gravity_calculation);
    gravity_calculation->SetAttribute("scheme", "k0");
    /// water_density node in calConfig
    XMLElement *water_density = doc->NewElement("water_density");
    calConfig->InsertEndChild(water_density);
    water_density->SetAttribute("value", 10);
    /// atm node in calConfig
    XMLElement *atm = doc->NewElement("atm");
    calConfig->InsertEndChild(atm);
    atm->SetAttribute("value", 100);

  } /// end calConfig

  XMLComment *comment4 = doc->NewComment("4. model's topology information");
  doc->InsertEndChild(comment4);
  /// mesh node in doc
  XMLElement *mesh = doc->NewElement("topology");
  doc->InsertEndChild(mesh);
  {
    /// nodes node in mesh
    XMLElement *nodes = doc->NewElement("nodes");
    mesh->InsertEndChild(nodes);
    nodes->SetAttribute("number", 0);
    /// elements node in mesh
    XMLElement *elements = doc->NewElement("elements");
    mesh->InsertEndChild(elements);
    elements->SetAttribute("number", 0);
    elements->SetAttribute("order", "clockwise");
    elements->SetAttribute("poly_degree", 1);
    /// element_mega_data node in doc
    XMLElement *element_mega_data = doc->NewElement("element_mega_data");
    mesh->InsertEndChild(element_mega_data);
    {
      XMLElement *Segment = doc->NewElement("Segment");
      element_mega_data->InsertEndChild(Segment);
      Segment->SetAttribute("gauss_num", 3);
      XMLElement *Triangle = doc->NewElement("Triangle");
      element_mega_data->InsertEndChild(Triangle);
      Triangle->SetAttribute("gauss_num", 3);
      XMLElement *Quadrangle = doc->NewElement("Quadrangle");
      element_mega_data->InsertEndChild(Quadrangle);
      Quadrangle->SetAttribute("gauss_num", 4);
      XMLElement *Tetrahedral = doc->NewElement("Tetrahedral");
      element_mega_data->InsertEndChild(Tetrahedral);
      Tetrahedral->SetAttribute("gauss_num", 4);
      XMLElement *TriPrism = doc->NewElement("TriPrism");
      element_mega_data->InsertEndChild(TriPrism);
      TriPrism->SetAttribute("gauss_num", 6);
      XMLElement *Octahedral = doc->NewElement("Octahedral");
      element_mega_data->InsertEndChild(Octahedral);
      Octahedral->SetAttribute("gauss_num", 8);
    }
  }

  XMLComment *comment5 = doc->NewComment("5. material information");
  doc->InsertEndChild(comment5);
  /// materials node in doc
  XMLElement *materials = doc->NewElement("materials");
  doc->InsertEndChild(materials);
  materials->SetAttribute("number", 0);

  XMLComment *comment6 = doc->NewComment("6. initial condition information");
  doc->InsertEndChild(comment6);
  /// initial condition in doc
  XMLElement *init_condition = doc->NewElement("init_condition");
  doc->InsertEndChild(init_condition);
  /// init_crack node in init_condition
  XMLElement *init_crack = doc->NewElement("init_crack");
  init_condition->InsertEndChild(init_crack);
  init_crack->SetAttribute("number", 0);
  XMLComment *init_crack_example = doc->NewComment(" <seg id=\"1\" x0=\"0\" y0=\"0\" x1=\"1\" y1=\"1\"/> ");
  init_crack->InsertEndChild(init_crack_example);
  /// init_displacement node in init_condition
  XMLElement *init_displacement = doc->NewElement("init_displacement");
  init_condition->InsertEndChild(init_displacement);
  init_displacement->SetAttribute("number", 0);
  XMLComment *init_disp_example = doc->NewComment(" <node=\"1\" dof=\"x1\" value=\"some_value\"/> ");
  init_displacement->InsertEndChild(init_disp_example);
  /// init_pore node in init_condition
  XMLElement *init_pore = doc->NewElement("init_pore");
  init_condition->InsertEndChild(init_pore);
  init_pore->SetAttribute("number", 0);
  XMLComment *init_pore_example = doc->NewComment(" <node=\"1\" value=\"some_value\"/> ");
  init_pore->InsertEndChild(init_pore_example);

  XMLComment *comment7 = doc->NewComment("7. phase step information");
  doc->InsertEndChild(comment7);
  /// phases node in doc
  XMLElement *phases = doc->NewElement("phases");
  doc->InsertEndChild(phases);
  phases->SetAttribute("number", 0);

  return 0;
}

void initialPhase(XMLDocument *doc, XMLElement *phase, const int &id)
{
  phase->SetAttribute("id", id);
  phase->SetAttribute("total_time", "1000");
  phase->SetAttribute("delta_time", "10");
  phase->SetAttribute("element_num", "9");
  XMLElement *newBC = doc->NewElement("boundary_condition");
  phase->InsertEndChild(newBC);
  newBC->SetAttribute("Dirichlet", 0);
  newBC->SetAttribute("Neumann", 0);
}

void insertDirichlet(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
{
  XMLDocument *doc = root->GetDocument();
  for (int bcNode : nodeTags)
  {
    XMLElement *newDirichlet = doc->NewElement("Dirichlet");
    root->InsertEndChild(newDirichlet);
    int Dirichlet = root->IntAttribute("Dirichlet");
    root->SetAttribute("Dirichlet", ++Dirichlet);
    newDirichlet->SetAttribute("node", bcNode);
    newDirichlet->SetAttribute("dof", tokens[2].c_str());
    newDirichlet->SetAttribute("type", "linear");
    newDirichlet->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
    newDirichlet->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
  }
}

void insertNeumann(XMLElement *root, std::vector<std::string> &tokens, const std::vector<std::size_t> &nodeTags)
{
  XMLDocument *doc = root->GetDocument();
  for (int bcNode : nodeTags)
  {
    XMLElement *newNeumann = doc->NewElement("Neumann");
    root->InsertEndChild(newNeumann);
    int Neumann = root->IntAttribute("Neumann");
    root->SetAttribute("Neumann", ++Neumann);
    newNeumann->SetAttribute("node", bcNode);
    newNeumann->SetAttribute("dof", tokens[2].c_str());
    newNeumann->SetAttribute("type", "linear");
    newNeumann->SetAttribute("start", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_start").c_str());
    newNeumann->SetAttribute("end", (tokens[0] + "_" + tokens[1] + "_" + tokens[2] + "_" + tokens[3] + "_end").c_str());
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
  XMLElement *volume_weight = doc->NewElement("volume_weight");
  mat->InsertEndChild(volume_weight);
  volume_weight->SetAttribute("value", "19");
  volume_weight->SetAttribute("unit", "kN/m^3");
  XMLElement *friction_angle = doc->NewElement("friction_angle");
  mat->InsertEndChild(friction_angle);
  friction_angle->SetAttribute("value", "30");
  friction_angle->SetAttribute("unit", "degree");
  XMLElement *dilation_angle = doc->NewElement("dilation_angle");
  mat->InsertEndChild(dilation_angle);
  dilation_angle->SetAttribute("value", "30");
  dilation_angle->SetAttribute("unit", "degree");
  XMLElement *cohesion = doc->NewElement("cohesion");
  mat->InsertEndChild(cohesion);
  cohesion->SetAttribute("value", "10");
  cohesion->SetAttribute("unit", "kPa");
  XMLElement *perm_x1 = doc->NewElement("perm_x1");
  mat->InsertEndChild(perm_x1);
  perm_x1->SetAttribute("value", "1e-5");
  perm_x1->SetAttribute("unit", "m/day");
  XMLElement *perm_x2 = doc->NewElement("perm_x2");
  mat->InsertEndChild(perm_x2);
  perm_x2->SetAttribute("value", "1e-5");
  perm_x2->SetAttribute("unit", "m/day");
  XMLElement *perm_x3 = doc->NewElement("perm_x3");
  mat->InsertEndChild(perm_x3);
  perm_x3->SetAttribute("value", "1e-5");
  perm_x3->SetAttribute("unit", "m/day");
  /// XMLElement *ft = doc->NewElement("ft");
  /// mat->InsertEndChild(ft);
  /// ft->SetAttribute("value", "15");
  XMLElement *zeta = doc->NewElement("zeta");
  mat->InsertEndChild(zeta);
  zeta->SetAttribute("value", "0.1");
  XMLElement *K_IC = doc->NewElement("K_IC");
  mat->InsertEndChild(K_IC);
  K_IC->SetAttribute("value", "2.9");
  K_IC->SetAttribute("unit", "kPa*m^0.5");
  XMLElement *elastic_modulus = doc->NewElement("elastic_modulus");
  mat->InsertEndChild(elastic_modulus);
  elastic_modulus->SetAttribute("value", "10000");
  elastic_modulus->SetAttribute("unit", "kPa");
  XMLElement *elastic_poisson = doc->NewElement("elastic_poisson");
  mat->InsertEndChild(elastic_poisson);
  elastic_poisson->SetAttribute("value", "0.3");
  if (type == "elastic")
  {
    /// XMLElement *elastic_modulus = doc->NewElement("elastic_modulus");
    /// mat->InsertEndChild(elastic_modulus);
    /// elastic_modulus->SetAttribute("value", "101");
    /// XMLElement *elastic_poisson = doc->NewElement("elastic_poisson");
    /// mat->InsertEndChild(elastic_poisson);
    /// elastic_poisson->SetAttribute("value", "102");
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
    XMLElement *D = doc->NewElement("D");
    mat->InsertEndChild(D);
    D->SetAttribute("value", "206");
  }
  else if (type == "MohrColumb")
  {
    ;
  }
  else if (type == "ModifiedCamClay")
  {
    XMLElement *swelling_slope = doc->NewElement("swelling_slope");
    mat->InsertEndChild(swelling_slope);
    swelling_slope->SetAttribute("value", "401");
    XMLElement *normal_compress_slope = doc->NewElement("normal_compress_slope");
    mat->InsertEndChild(normal_compress_slope);
    normal_compress_slope->SetAttribute("value", "402");
    XMLElement *initial_void = doc->NewElement("initial_void");
    mat->InsertEndChild(initial_void);
    initial_void->SetAttribute("value", "403");
  }
  else if (type == "DruckPrager")
  {
    XMLElement *qf = doc->NewElement("qf");
    mat->InsertEndChild(qf);
    qf->SetAttribute("value", "501");
    XMLElement *kf = doc->NewElement("kf");
    mat->InsertEndChild(kf);
    kf->SetAttribute("value", "502");
    XMLElement *q_dilantion = doc->NewElement("q_dilantion");
    mat->InsertEndChild(q_dilantion);
    q_dilantion->SetAttribute("value", "503");
  }
  mat->SetAttribute("constitutive", type.c_str());
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
  std::string dimension;
  switch (dim)
  {
  case 1:
    dimension = "T";
    break;
  case 2:
    dimension = "P";
    break;
  case 3:
    dimension = "C";
    break;
  default:
    std::cout << "Error dimension!\n";
    break;
  }
  return dimension + std::to_string(dim) + "D" + std::to_string(nEdge) + "B" + interpolation;
}
