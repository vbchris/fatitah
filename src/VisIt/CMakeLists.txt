# Override this from configure if you want to install -public (to VisIt 
# install dir instead of user's .visit directory).
set(VISIT_PLUGIN_INSTALL_TYPE "-private")
set(VISIT_PLUGIN_SRC_DIR ${Uintah_SOURCE_DIR}/VisIt/udaReaderMTMD) 
set(VISIT_PLUGIN_BIN_DIR ${Uintah_BINARY_DIR}/VisIt/udaReaderMTMD) 
set(UDA2VIS_SRC_DIR ${Uintah_SOURCE_DIR}/StandAlone/tools/uda2vis) 
set(UDA2VIS_BIN_DIR ${Uintah_BINARY_DIR}/StandAlone/tools/uda2vis) 

if (WIN32)
  # Windows untested - BB
  set(LN_CMD mklink)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDFileFormat.h ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDFileFormat.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.C ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.C WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.h ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${UDA2VIS_BIN_DIR}/udaData.h ${UDA2VIS_SRC_DIR}/udaData.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
else ()
  set(LN_CMD ln -fs)
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDFileFormat.h ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDFileFormat.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.C ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.C WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.h ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND ln -fs ${UDA2VIS_SRC_DIR}/udaData.h ${UDA2VIS_BIN_DIR}/udaData.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
endif ()


# Execute the xml2cmake command to create CMake files
execute_process(COMMAND ${VISIT_DIR}/bin/xml2cmake ${VISIT_PLUGIN_INSTALL_TYPE} -clobber ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
execute_process(COMMAND ${VISIT_DIR}/bin/xml2info -clobber ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
#execute_process(COMMAND cmake -DVISIT_DISABLE_SETTING_COMPILER:BOOL=TRUE -DCMAKE_CXX_COMPILER:FILEPATH=${CXX} -DCMAKE_CXX_FLAGS:STRING="-I${Uintah_BINARY_DIR} -I${Uintah_SOURCE_DIR} ${CXXFLAGS}" ${VISIT_PLUGIN_BIN_DIR})

include(ExternalProject)
set(VISIT_CXX_FLAGS_C -DCMAKE_CXX_FLAGS:STRING=-I${Uintah_BINARY_DIR})
set(VISIT_INSTALL_FLAGS -DCMAKE_INSTALL_PREFIX=${VISIT_PLUGIN_BIN_DIR})
set(VISIT_CXX_FLAGS ${VISIT_CXX_FLAGS_C})
message("VISIT CXX Flags = ${VISIT_CXX_FLAGS}")
ExternalProject_Add (
  udaReaderMTMD 
  PREFIX ${VISIT_PLUGIN_BIN_DIR}
  DOWNLOAD_DIR ${VISIT_PLUGIN_SRC_DIR}
  SOURCE_DIR ${VISIT_PLUGIN_BIN_DIR}
  BINARY_DIR ${VISIT_PLUGIN_BIN_DIR}
  INSTALL_DIR ${VISIT_PLUGIN_BIN_DIR}
  INSTALL_COMMAND ""
  CMAKE_ARGS ${VISIT_CXX_FLAGS} ${VISIT_INSTALL_FLAGS}
)

