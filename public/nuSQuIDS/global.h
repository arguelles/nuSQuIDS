
#ifndef __GLOBAL_H
#define __GLOBAL_H

#include <cstdlib>
#include <string>

#define XSECTION_LOCATION "/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/data/xsections/"
#define SUN_MODEL_LOCATION  "/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/data/astro/bs05_agsop.dat"
#define SUN_MODEL_NELECTRON_LOCATION "/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/data/astro/nele_bs05op.dat"
// #define EARTH_MODEL_LOCATION "resources/data/astro/EARTH_MODEL_PREM.dat"
const std::string i3_testdata = getenv("I3_TESTDATA");
const std::string EARTH_MODEL_LOCATION = i3_testdata + "/nuSQuIDS/astro/EARTH_MODEL_PREM.dat";

#endif

