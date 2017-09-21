

#this is required to disable ASLR (randomized address space) on OSX
set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_pie" CACHE BOOL "...")

#optional Ordering libraries
#set(SCOTCH_DIR /path/to/scotch)
#set(ENABLE_SCOTCH ON CACHE BOOL "...")

#set(METIS_DIR /path/to/metis)
#set(ENABLE_METIS ON CACHE BOOL "...")

#set(PARMETIS_DIR /path/to/parmetis)
#set(ENABLE_PARMETIS ON CACHE BOOL "...")
