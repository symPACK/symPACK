find ../../include/sympack/*.hpp -type f | xargs perl -p -i -ne 's/(?<!#include <)(std::){0,1}(?<!SYMPACK::)vector/SYMPACK::vector/g'
