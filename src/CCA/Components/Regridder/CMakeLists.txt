#
#  For more information, please see: http://software.sci.utah.edu
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
#  and/or sell copies of the Software, and to permit persons to whom the
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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Regridder

SET(Packages_Uintah_CCA_Components_Regridder_SRCS
  RegridderCommon.cc 
  RegridderFactory.cc 
  HierarchicalRegridder.cc 
  BNRTask.cc 
  BNRRegridder.cc 
  PatchFixer.cc  
  TiledRegridder.cc
)

ADD_LIBRARY(Packages_Uintah_CCA_Components_Regridder ${Packages_Uintah_CCA_Components_Regridder_SRCS})

TARGET_LINK_LIBRARIES(Packages_Uintah_CCA_Components_Regridder
  Packages_Uintah_CCA_Ports
  Packages_Uintah_Core_Disclosure
  Packages_Uintah_Core_ProblemSpec
  Packages_Uintah_Core_Parallel
  Packages_Uintah_Core_Grid
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_Thread
  ${XML2_LIBRARY}
  ${MPI_LIBS}
  ${M_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_CCA_Components_Regridder)
ENDIF(BUILD_SHARED_LIBS)
