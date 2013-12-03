#ifndef _LOGFILE_HEADER_
#define _LOGFILE_HEADER_

#include <fstream>

class LogFile{
protected:
  std::ofstream myOFS_;

  void init_(const char * prefix, const char * suffix){
    stringstream  ss;
    ss << prefix << suffix;
    myOFS_.open( ss.str().c_str() );
  }
public:

  LogFile(const char * prefix, const char * suffix){
    init_(prefix,suffix);
  }

  LogFile(int iam){
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

//  template <typename T>
//  const LogFile& operator<<(const T& obj) 
//  {
//    // write obj to stream
//    mySS_<<obj;
//    myOFS_<<mySS_.str();
//    return *this;
//  }

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

#endif
