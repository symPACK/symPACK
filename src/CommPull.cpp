#include "ngchol/CommPull.hpp"

namespace LIBCHOLESKY{
  std::list< IncomingMessage * > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
  int gMaxIrecv = 0;
}
