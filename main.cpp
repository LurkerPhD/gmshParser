#include "gmshParser.h"

int main(int argc, char const *argv[])
{
    gmshParser parser;
    std::string filename;

    if (argc > 1)
        filename = argv[1];
    else
        filename = "/Users/lurker/Desktop/geoxfem/demo/test1.msh";

    parser.readFromGmshFile(filename);

    return 0;
}
