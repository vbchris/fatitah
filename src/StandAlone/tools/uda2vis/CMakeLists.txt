#
#  For more information, please see: http:__software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and_or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

# CMakeLists.txt for Packages_Uintah_StandAlone_tools_uda2vis

set(UDA2VIS_SRC
  uda2vis.cc
)

set(UDA2VIS_LIBS
  Packages_Uintah_Core_Exceptions    
  Packages_Uintah_Core_Grid          
  Packages_Uintah_Core_Util          
  Packages_Uintah_Core_Math          
  Packages_Uintah_Core_ProblemSpec   
  Packages_Uintah_Core_DataArchive   
  Packages_Uintah_Core_Parallel      
  Packages_Uintah_CCA_Ports          
  Packages_Uintah_CCA_Components_ProblemSpecification 
  Packages_Uintah_Core_Basis        
  Packages_Uintah_Core_Containers   
  Packages_Uintah_Core_Datatypes    
  Packages_Uintah_Core_Persistent   
  Packages_Uintah_Core_Thread       
  ${XML2_LIBRARY}
  ${Z_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
  ${F_LIBRARY}
)

add_library(Standalone_tools_uda2vis ${UDA2VIS_SRC})

target_link_libraries(Standalone_tools_uda2vis ${LIBS})
	
if (BUILD_SHARED_LIBS)  
  add_definitions(-DBUILD_Standalone_tools_uda2vis)
endif (BUILD_SHARED_LIBS)  
