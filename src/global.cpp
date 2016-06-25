#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/SuperNode.hpp"

#include <upcxx.h>

namespace SYMPACK{



#ifdef _TRACK_MEMORY_
  std::map<char*,size_t> MemoryAllocator::cnt_;
  size_t MemoryAllocator::total_=0;
#endif


  upcxx::team * workteam;
// *********************************************************************
// IO
// *********************************************************************
  std::ofstream  statusOFS;

  SYMPACK::vector<std::ofstream>  statusOFSs;

// *********************************************************************
// Error handling
// *********************************************************************
	// If we are not in RELEASE mode, then implement wrappers for a
	// CallStack
#ifndef _RELEASE_
	std::stack<std::string> callStack;	

	void PushCallStack( std::string s )
	{ callStack.push(s); }

	void PopCallStack()
	{ callStack.pop(); }

	void DumpCallStack()
	{
		std::ostringstream msg;
		while( ! callStack.empty() )
		{
			msg << "Stack[" << callStack.size() << "]: " << callStack.top() << "\n";
			callStack.pop();
		}
		std::cerr << msg.str() << std::endl;
	}

#endif // ifndef _RELEASE_
} // namespace SYMPACK
