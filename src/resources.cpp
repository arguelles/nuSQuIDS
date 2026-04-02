#include <nuSQuIDS/resources.h>
#include <nuSQuIDS/version.h>

#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include "resource_paths.h"

namespace nusquids{

namespace{
    bool dirExists(const std::string& path){
        struct stat info;
        return stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR);
    }
}

std::string getResourcePath(){
    // 1. Explicit environment variable (highest priority)
    char* pathFromEnv=getenv("NUSQUIDS_DATA_PATH");
    if(pathFromEnv)
        return pathFromEnv;

    // 2. Standard paths from build configuration, but only if they exist
    if(getInstallBit() && dirExists(INSTALL_DATA_PATH))
        return INSTALL_DATA_PATH;
    if(dirExists(SOURCE_DATA_PATH))
        return SOURCE_DATA_PATH;

    // 3. Pooch-style cache directory (used by pip installations)
    //    Check NUSQUIDS_DATA_HOME first, then platform default
    char* dataHome=getenv("NUSQUIDS_DATA_HOME");
    std::string cacheBase;
    if(dataHome){
        cacheBase=dataHome;
    } else {
        char* home=getenv("HOME");
        if(home){
            char* xdgData=getenv("XDG_DATA_HOME");
            if(xdgData)
                cacheBase=std::string(xdgData)+"/nuSQuIDS";
            else
                cacheBase=std::string(home)+"/.local/share/nuSQuIDS";
        }
    }
    if(!cacheBase.empty() && dirExists(cacheBase)){
        // Look for versioned subdirectories (e.g., v1.13.3/)
        // containing the expected data layout
        std::string versionedPath=cacheBase+"/v"+NUSQUIDS_VERSION_STR;
        if(dirExists(versionedPath))
            return versionedPath;
    }

    // 4. Fall back to the build-time paths even if they don't exist
    //    (preserves original behavior; downstream code will report the actual error)
    if(getInstallBit())
        return INSTALL_DATA_PATH;
    return SOURCE_DATA_PATH;
}

}
