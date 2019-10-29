#include "gmshParser.h"
#include <dirent.h>
#include <iostream>
#include <regex>

void parseDirectory(const std::string &_path)
{
    std::string path = _path;
    if (path.back() != 47)
        path.append("/");

    struct dirent *next_file = NULL;
    DIR *dir = opendir(path.c_str());

    if (!dir)
        return;

    while ((next_file = readdir(dir)))
    {
        DIR *sub_dir = nullptr;
        FILE *file = nullptr;
        std::string abs_path;
        if (*(next_file->d_name) != '.')
        {
            abs_path = path + next_file->d_name;
            if ((sub_dir = opendir(abs_path.c_str())))
            {
                closedir(sub_dir);
                parseDirectory(abs_path);
            }
            else
            {
                auto [_ifXML, postfix] = parseName(abs_path, ".*[/|\\\\].+?\\.(.{3})$");
                if (_ifXML)
                {
                    if ((file = fopen(abs_path.c_str(), "r")) && postfix.compare("msh") == 0)
                    {
                        fclose(file);
                        gmshParser parser;
                        parser.readFromGmshFile(abs_path);
                    }
                }
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    std::string filename;

    if (argc > 1)
        filename = argv[1];
    else
        filename = "/Users/lurker.phd/Desktop/geoxfem/gmshParser";
    // filename = "/Users/lurker/Desktop/geoxfem/demo";

    parseDirectory(filename.c_str());

    return 0;
}
