#include "argumentparser.h"

#include <algorithm>

ArgumentParser::ArgumentParser(int &argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
    {
        this->tokens.push_back(std::string(argv[i]));
    }
}

const std::string &ArgumentParser::getCmdOption(const std::string &option) const
{
    std::vector<std::string>::const_iterator itr;

    itr = std::find(this->tokens.begin(), this->tokens.end(), option);

    if (itr != this->tokens.end() && ++itr != this->tokens.end())
    {
        return *itr;
    }

    static const std::string empty_string("");

    return empty_string;
}

bool ArgumentParser::cmdOptionExists(const std::string &option) const
{
    return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
}
