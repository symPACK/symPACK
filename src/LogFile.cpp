#include "sympack/LogFile.hpp"
#include <iostream>

namespace symPACK{

LogFile * logfileptr = nullptr;
LogFile * profileptr = nullptr;
LogFile * progressptr = nullptr;
LogFile * statfileptr = nullptr;
std::stringstream * progstr = nullptr;
std::ostream symPACKOS( std::cout.rdbuf() );

}

