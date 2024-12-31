#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include <vector>
#include <string>

//From https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
class ArgumentParser
{
public:
    ArgumentParser (int &argc, char **argv);

    /// @author iain
    const std::string& getCmdOption(const std::string &option) const;

    /// @author iain
    bool cmdOptionExists(const std::string &option) const;

private:
    std::vector<std::string> tokens;
};

#endif // ARGUMENTPARSER_H
