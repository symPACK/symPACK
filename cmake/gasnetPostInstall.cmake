file(INSTALL ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install DESTINATION ${CMAKE_INSTALL_PREFIX})



file(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include/*conduit/*.mak")
foreach(file ${files})
  #... calculate ${i} to get the test name
  #add_test(validate_${i}, "validator", ${file})
  file(READ ${file} file_content)
  string(REGEX REPLACE "\n" ";" ContentsAsList "${file_content}")
  unset(ModifiedContents)
  foreach(Line ${ContentsAsList})
    if("${Line}" MATCHES "GASNET_PREFIX = ")
      message(${Line})
      string(REGEX REPLACE "= ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install" "= ${CMAKE_INSTALL_PREFIX}/gasnet_install" Line ${Line})
      message(${Line})
    endif()

    set(ModifiedContents "${ModifiedContents}${Line}\n")

  endforeach()
  string(REGEX REPLACE "${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install" "${CMAKE_INSTALL_PREFIX}/gasnet_install" ofile ${file})
  message(${ofile})
  file(WRITE ${ofile} ${ModifiedContents})
endforeach()

