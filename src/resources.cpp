#include <nuSQuIDS/resources.h>

#include <cstdlib>
#include "resource_paths.h"

namespace nusquids{
    
std::string getResourcePath(){
    char* pathFromEnv=getenv("NUSQUIDS_DATA_PATH");
    if(pathFromEnv)
        return pathFromEnv;
    if(getInstallBit())
        return INSTALL_DATA_PATH;
    return SOURCE_DATA_PATH;
}
    
}
