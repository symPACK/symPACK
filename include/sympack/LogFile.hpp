#ifndef _LOGFILE_HEADER_
#define _LOGFILE_HEADER_

#include <sstream>
#include <fstream>

namespace symPACK{

class LogFile{
protected:

  void init_(const char * prefix, const char * suffix){
    std::stringstream  ss;
    ss << prefix << suffix;
    if(verbose){
      myOFS_.open( ss.str().c_str() );
    }
  }
public:
  bool verbose;
  std::ofstream myOFS_;

  LogFile(const char * prefix, const char * suffix,bool averbose = true){
    verbose = averbose;
    init_(prefix,suffix);
  }

  LogFile(int iam,bool averbose = true){
    verbose = averbose;
    std::stringstream ss;
    ss<<iam;
    init_("logTest",ss.str().c_str());
  }


  LogFile(){
  }




  ~LogFile(){
    if(myOFS_.is_open()){
      myOFS_.close();
    }
  }


  template <typename T>
  const std::ofstream& operator<<(const T& obj) 
  {
    // write obj to stream
    std::stringstream ss;
    ss<<obj;
    myOFS_<<ss.str();
    return myOFS_;
  }

  std::ofstream & OFS(){
    return myOFS_;
  }

};

extern LogFile * logfileptr;
extern LogFile * profileptr;
extern LogFile * progressptr;
extern LogFile * statfileptr;
extern LogFile * dumpfileptr;
extern std::stringstream * progstr;
extern std::ostream symPACKOS;

}

#endif //_LOGFILE_HEADER_
