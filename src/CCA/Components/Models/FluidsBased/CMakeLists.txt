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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Models/FluidsBased

SET(Packages_Uintah_CCA_Components_Models_SRCS
  ${Packages_Uintah_CCA_Components_Models_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/ArchesTable.cc    
  ${CMAKE_CURRENT_SOURCE_DIR}/TableInterface.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TableFactory.cc 
)
if(ICE)
  SET(Packages_Uintah_CCA_Components_Models_SRCS
    ${Packages_Uintah_CCA_Components_Models_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AdiabaticTable.cc     
    ${CMAKE_CURRENT_SOURCE_DIR}/flameSheet_rxn.cc     
    ${CMAKE_CURRENT_SOURCE_DIR}/MaterialProperties.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/Mixing.cc             
    ${CMAKE_CURRENT_SOURCE_DIR}/NonAdiabaticTable.cc  
    ${CMAKE_CURRENT_SOURCE_DIR}/PassiveScalar.cc      
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleRxn.cc          
    ${CMAKE_CURRENT_SOURCE_DIR}/TestModel.cc          
    ${CMAKE_CURRENT_SOURCE_DIR}/MassMomEng_src.cc
  )
endif(ICE)
SET(Packages_Uintah_CCA_Components_Models_SRCS
  ${Packages_Uintah_CCA_Components_Models_SRCS}
  PARENT_SCOPE
)

