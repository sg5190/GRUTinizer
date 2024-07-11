#include <iostream>
#include <limits.h>

#include "TEnv.h"
#include "TPython.h"
#include "TString.h"

#include "TGRUTOptions.h"

int main(int argc, char** argv){
  // Set the GRUTSYS variable based on the executable path.
  // If GRUTSYS has already been defined, don't overwrite.
  if(!getenv("GRUTSYS")){
    char buff[PATH_MAX+1];
    size_t len = readlink("/proc/self/exe", buff, sizeof(buff)-1);
    buff[len] = '\0';
    std::string exe_path = buff;
    exe_path = exe_path.substr(0, exe_path.find_last_of('/')-4);
    setenv("GRUTSYS", exe_path.c_str(), 0);
  }
  std::string grut_path = Form("%s/.grutrc",getenv("GRUTSYS")); // + "/../.grutrc";
  gEnv->ReadFile(grut_path.c_str(),kEnvChange);


  TGRUTOptions::Get(argc-1, argv+1);

  if(argc > 1){
    TPython::ExecScript(argv[1], argc-1, const_cast<const char**>(argv+1));
  } else {
    TPython::Prompt();
  }
}
