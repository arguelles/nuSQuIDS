
#ifndef NUSQUIDS_RESOURCES_H
#define NUSQUIDS_RESOURCES_H

#include <string>

namespace nusquids{
    
///Get the path to where data tables distributed with the library are stored. 
///This may be overridden by setting the NUSQUIDS_DATA_PATH environment variable. 
std::string getResourcePath();
    
}

#endif //NUSQUIDS_RESOURCES_H 
