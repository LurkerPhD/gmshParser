
#include <gmsh.h>
#include <iostream>
#include <vector>

int main()
{

    // step 1 of 2 -> import data for plotting (need the mesh and the simulation data)

    gmsh::initialize();

    // open mesh file to use as reference for plotting
    // we will write our results to another mesh file for output
    gmsh::open("demo/output/try.msh");

    // tell gmsh what elements of the mesh we want to plot data for
    std::vector<std::size_t> elt_numbers;
    // only plot data for element 1 for simplicity
    elt_numbers.emplace_back(1);
    elt_numbers.emplace_back(2);
    elt_numbers.emplace_back(3);
    elt_numbers.emplace_back(4);
    elt_numbers.emplace_back(5);
    elt_numbers.emplace_back(6);
    elt_numbers.emplace_back(7);
    elt_numbers.emplace_back(8);

    // fill a data vector corresponding to the element vector
    // the data vector's elements are vectors of length 9 (3x3 tensor)
    std::vector<std::vector<double>> elt_data;

    ///////////////////////////////////////////////////////////////////////////////

    // this is the tricky part -> how to import data in an intelligent way?
    // here we make a single zero vector of length 9 for demonstration
    elt_data.emplace_back(std::vector<double>({2, 3, 2, 1, 7, 8, 5, 6, 1}));
    elt_data.emplace_back(std::vector<double>({2, 6, 8, 4, 2, 3, 9, 0, 1}));
    elt_data.emplace_back(std::vector<double>({9, 7, 0, 2, 3, 5, 7, 5, 3}));
    elt_data.emplace_back(std::vector<double>({7, 9, 4, 7, 5, 6, 1, 2, 3}));
    elt_data.emplace_back(std::vector<double>({8, 3, 4, 6, 7, 8, 0, 7, 8}));
    elt_data.emplace_back(std::vector<double>({3, 5, 7, 8, 3, 5, 1, 6, 5}));
    elt_data.emplace_back(std::vector<double>({7, 7, 9, 2, 3, 4, 1, 1, 2}));
    elt_data.emplace_back(std::vector<double>({0, 8, 9, 6, 6, 8, 4, 5, 4}));

    ///////////////////////////////////////////////////////////////////////////////

    // step 2 of 2 -> ask Gmsh API to plot the data (see gmsh.h for reference)

    //get model name
    std::vector<std::string> modelNames;
    gmsh::model::list(modelNames);
    auto modelName = modelNames[0];

    auto viewTag0 = gmsh::view::add("stress-data");

    gmsh::view::addModelData(viewTag0, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 0);

    elt_data.clear();
    elt_data.emplace_back(std::vector<double>({2 * 0.9, 3 * 0.9, 2 * 0.9, 1 * 0.9, 7 * 0.9, 8 * 0.9, 5 * 0.9, 6 * 0.9, 1 * 0.9}));
    elt_data.emplace_back(std::vector<double>({2 * 0.9, 6 * 0.9, 8 * 0.9, 4 * 0.9, 2 * 0.9, 3 * 0.9, 9 * 0.9, 0 * 0.9, 1 * 0.9}));
    elt_data.emplace_back(std::vector<double>({9 * 0.9, 7 * 0.9, 0 * 0.9, 2 * 0.9, 3 * 0.9, 5 * 0.9, 7 * 0.9, 5 * 0.9, 3 * 0.9}));
    elt_data.emplace_back(std::vector<double>({7 * 0.9, 9 * 0.9, 4 * 0.9, 7 * 0.9, 5 * 0.9, 6 * 0.9, 1 * 0.9, 2 * 0.9, 3 * 0.9}));
    elt_data.emplace_back(std::vector<double>({8 * 0.9, 3 * 0.9, 4 * 0.9, 6 * 0.9, 7 * 0.9, 8 * 0.9, 0 * 0.9, 7 * 0.9, 8 * 0.9}));
    elt_data.emplace_back(std::vector<double>({3 * 0.9, 5 * 0.9, 7 * 0.9, 8 * 0.9, 3 * 0.9, 5 * 0.9, 1 * 0.9, 6 * 0.9, 5 * 0.9}));
    elt_data.emplace_back(std::vector<double>({7 * 0.9, 7 * 0.9, 9 * 0.9, 2 * 0.9, 3 * 0.9, 4 * 0.9, 1 * 0.9, 1 * 0.9, 2 * 0.9}));
    elt_data.emplace_back(std::vector<double>({0 * 0.9, 8 * 0.9, 9 * 0.9, 6 * 0.9, 6 * 0.9, 8 * 0.9, 4 * 0.9, 5 * 0.9, 4 * 0.9}));

    auto viewTag1 = gmsh::view::add("stress-data");
    gmsh::view::addModelData(viewTag1, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 1);

    elt_data.clear();
    elt_data.emplace_back(std::vector<double>({2 * 0.8, 3 * 0.8, 2 * 0.8, 1 * 0.8, 7 * 0.8, 8 * 0.8, 5 * 0.8, 6 * 0.8, 1 * 0.8}));
    elt_data.emplace_back(std::vector<double>({2 * 0.8, 6 * 0.8, 8 * 0.8, 4 * 0.8, 2 * 0.8, 3 * 0.8, 9 * 0.8, 0 * 0.8, 1 * 0.8}));
    elt_data.emplace_back(std::vector<double>({9 * 0.8, 7 * 0.8, 0 * 0.8, 2 * 0.8, 3 * 0.8, 5 * 0.8, 7 * 0.8, 5 * 0.8, 3 * 0.8}));
    elt_data.emplace_back(std::vector<double>({7 * 0.8, 9 * 0.8, 4 * 0.8, 7 * 0.8, 5 * 0.8, 6 * 0.8, 1 * 0.8, 2 * 0.8, 3 * 0.8}));
    elt_data.emplace_back(std::vector<double>({8 * 0.8, 3 * 0.8, 4 * 0.8, 6 * 0.8, 7 * 0.8, 8 * 0.8, 0 * 0.8, 7 * 0.8, 8 * 0.8}));
    elt_data.emplace_back(std::vector<double>({3 * 0.8, 5 * 0.8, 7 * 0.8, 8 * 0.8, 3 * 0.8, 5 * 0.8, 1 * 0.8, 6 * 0.8, 5 * 0.8}));
    elt_data.emplace_back(std::vector<double>({7 * 0.8, 7 * 0.8, 9 * 0.8, 2 * 0.8, 3 * 0.8, 4 * 0.8, 1 * 0.8, 1 * 0.8, 2 * 0.8}));
    elt_data.emplace_back(std::vector<double>({0 * 0.8, 8 * 0.8, 9 * 0.8, 6 * 0.8, 6 * 0.8, 8 * 0.8, 4 * 0.8, 5 * 0.8, 4 * 0.8}));

    auto viewTag2 = gmsh::view::add("stress-data");
    gmsh::view::addModelData(viewTag2, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 2);

    elt_data.clear();
    elt_data.emplace_back(std::vector<double>({2 * 0.7, 3 * 0.7, 2 * 0.7, 1 * 0.7, 7 * 0.7, 8 * 0.7, 5 * 0.7, 6 * 0.7, 1 * 0.7}));
    elt_data.emplace_back(std::vector<double>({2 * 0.7, 6 * 0.7, 8 * 0.7, 4 * 0.7, 2 * 0.7, 3 * 0.7, 9 * 0.7, 0 * 0.7, 1 * 0.7}));
    elt_data.emplace_back(std::vector<double>({9 * 0.7, 7 * 0.7, 0 * 0.7, 2 * 0.7, 3 * 0.7, 5 * 0.7, 7 * 0.7, 5 * 0.7, 3 * 0.7}));
    elt_data.emplace_back(std::vector<double>({7 * 0.7, 9 * 0.7, 4 * 0.7, 7 * 0.7, 5 * 0.7, 6 * 0.7, 1 * 0.7, 2 * 0.7, 3 * 0.7}));
    elt_data.emplace_back(std::vector<double>({8 * 0.7, 3 * 0.7, 4 * 0.7, 6 * 0.7, 7 * 0.7, 8 * 0.7, 0 * 0.7, 7 * 0.7, 8 * 0.7}));
    elt_data.emplace_back(std::vector<double>({3 * 0.7, 5 * 0.7, 7 * 0.7, 8 * 0.7, 3 * 0.7, 5 * 0.7, 1 * 0.7, 6 * 0.7, 5 * 0.7}));
    elt_data.emplace_back(std::vector<double>({7 * 0.7, 7 * 0.7, 9 * 0.7, 2 * 0.7, 3 * 0.7, 4 * 0.7, 1 * 0.7, 1 * 0.7, 2 * 0.7}));
    elt_data.emplace_back(std::vector<double>({0 * 0.7, 8 * 0.7, 9 * 0.7, 6 * 0.7, 6 * 0.7, 8 * 0.7, 4 * 0.7, 5 * 0.7, 4 * 0.7}));

    auto viewTag3 = gmsh::view::add("stress-data");
    gmsh::view::addModelData(viewTag3, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 3);

    elt_data.clear();
    elt_data.emplace_back(std::vector<double>({2 * 0.6, 3 * 0.6, 2 * 0.6, 1 * 0.6, 7 * 0.6, 8 * 0.6, 5 * 0.6, 6 * 0.6, 1 * 0.6}));
    elt_data.emplace_back(std::vector<double>({2 * 0.6, 6 * 0.6, 8 * 0.6, 4 * 0.6, 2 * 0.6, 3 * 0.6, 9 * 0.6, 0 * 0.6, 1 * 0.6}));
    elt_data.emplace_back(std::vector<double>({9 * 0.6, 7 * 0.6, 0 * 0.6, 2 * 0.6, 3 * 0.6, 5 * 0.6, 7 * 0.6, 5 * 0.6, 3 * 0.6}));
    elt_data.emplace_back(std::vector<double>({7 * 0.6, 9 * 0.6, 4 * 0.6, 7 * 0.6, 5 * 0.6, 6 * 0.6, 1 * 0.6, 2 * 0.6, 3 * 0.6}));
    elt_data.emplace_back(std::vector<double>({8 * 0.6, 3 * 0.6, 4 * 0.6, 6 * 0.6, 7 * 0.6, 8 * 0.6, 0 * 0.6, 7 * 0.6, 8 * 0.6}));
    elt_data.emplace_back(std::vector<double>({3 * 0.6, 5 * 0.6, 7 * 0.6, 8 * 0.6, 3 * 0.6, 5 * 0.6, 1 * 0.6, 6 * 0.6, 5 * 0.6}));
    elt_data.emplace_back(std::vector<double>({7 * 0.6, 7 * 0.6, 9 * 0.6, 2 * 0.6, 3 * 0.6, 4 * 0.6, 1 * 0.6, 1 * 0.6, 2 * 0.6}));
    elt_data.emplace_back(std::vector<double>({0 * 0.6, 8 * 0.6, 9 * 0.6, 6 * 0.6, 6 * 0.6, 8 * 0.6, 4 * 0.6, 5 * 0.6, 4 * 0.6}));

    auto viewTag4 = gmsh::view::add("stress-data");
    gmsh::view::addModelData(viewTag4, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 4);

    elt_data.clear();
    elt_data.emplace_back(std::vector<double>({2 * 0.1, 3 * 0.1, 2 * 0.1, 1 * 0.1, 7 * 0.1, 8 * 0.1, 5 * 0.1, 6 * 0.1, 1 * 0.1}));
    elt_data.emplace_back(std::vector<double>({2 * 0.1, 6 * 0.1, 8 * 0.1, 4 * 0.1, 2 * 0.1, 3 * 0.1, 9 * 0.1, 0 * 0.1, 1 * 0.1}));
    elt_data.emplace_back(std::vector<double>({9 * 0.1, 7 * 0.1, 0 * 0.1, 2 * 0.1, 3 * 0.1, 5 * 0.1, 7 * 0.1, 5 * 0.1, 3 * 0.1}));
    elt_data.emplace_back(std::vector<double>({7 * 0.1, 9 * 0.1, 4 * 0.1, 7 * 0.1, 5 * 0.1, 6 * 0.1, 1 * 0.1, 2 * 0.1, 3 * 0.1}));
    elt_data.emplace_back(std::vector<double>({8 * 0.1, 3 * 0.1, 4 * 0.1, 6 * 0.1, 7 * 0.1, 8 * 0.1, 0 * 0.1, 7 * 0.1, 8 * 0.1}));
    elt_data.emplace_back(std::vector<double>({3 * 0.1, 5 * 0.1, 7 * 0.1, 8 * 0.1, 3 * 0.1, 5 * 0.1, 1 * 0.1, 6 * 0.1, 5 * 0.1}));
    elt_data.emplace_back(std::vector<double>({7 * 0.1, 7 * 0.1, 9 * 0.1, 2 * 0.1, 3 * 0.1, 4 * 0.1, 1 * 0.1, 1 * 0.1, 2 * 0.1}));
    elt_data.emplace_back(std::vector<double>({0 * 0.1, 8 * 0.1, 9 * 0.1, 6 * 0.1, 6 * 0.1, 8 * 0.1, 4 * 0.1, 5 * 0.1, 4 * 0.1}));
    auto viewTag5 = gmsh::view::add("stress-data");
    gmsh::view::addModelData(viewTag5, 0, modelName,
                             "ElementData", elt_numbers, elt_data, 5);

    gmsh::view::combine("steps", "stress-data", true);

    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);

    //toggle here to see behaviour
    gmsh::view::write(viewTags[0], "stress-data.pos");

    //end gmsh run
    gmsh::finalize();
    return 0;
}
