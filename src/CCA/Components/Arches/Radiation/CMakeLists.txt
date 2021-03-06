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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Arches/Radiation

FORTRAN_DEPENDENCY(fortran/rordr.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rordrss.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rordrtn.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/radarray.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/radcal.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/radcoef.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/radwsgg.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rdomsolve.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rdomsrc.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rdomflux.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rdombmcalc.F CCA_Components_Arches_Radiation)
FORTRAN_DEPENDENCY(fortran/rdomvolq.F CCA_Components_Arches_Radiation)

SET(Packages_Uintah_CCA_Components_Arches_Radiation_SRCS
  ${CCA_Components_Arches_Radiation_FORTRAN_SRCS}
  RadiationSolver.cc
  DORadiationModel.cc

  #fortran files without fspec
  fortran/find.F
  fortran/fixradval.F
  fortran/eih2o.F
  fortran/eico2.F
  fortran/eco2.F
  fortran/eh2o.F
  fortran/soot.F
)

SET(LIBS 
  Packages_Uintah_Core_ProblemSpec
  Packages_Uintah_Core_Grid
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_Disclosure
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_Math
  Packages_Uintah_Core_Thread
)

IF(PETSC_LIBRARY)
  SET(Packages_Uintah_CCA_Components_Arches_Radiation_SRCS 
    ${Packages_Uintah_CCA_Components_Arches_Radiation_SRCS} RadPetscSolver.cc)
  SET(LIBS ${LIBS} ${PETSC_LIBRARIES})
ELSE(PETSC_LIBRARY)
  SET(Packages_Uintah_CCA_Components_Arches_Radiation_SRCS 
    ${Packages_Uintah_CCA_Components_Arches_Radiation_SRCS} FakeRadPetscSolver.cc)
ENDIF(PETSC_LIBRARY)

IF(HYPRE_LIBRARY)
  SET(Packages_Uintah_CCA_Components_Arches_Radiation_SRCS 
    ${Packages_Uintah_CCA_Components_Arches_Radiation_SRCS} RadHypreSolver.cc)
  SET(LIBS ${LIBS} ${HYPRE_LIBRARY})
ENDIF(HYPRE_LIBRARY)
  

ADD_LIBRARY(Packages_Uintah_CCA_Components_Arches_Radiation ${Packages_Uintah_CCA_Components_Arches_Radiation_SRCS})
SET_TARGET_PROPERTIES(Packages_Uintah_CCA_Components_Arches_Radiation PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(Packages_Uintah_CCA_Components_Arches_Radiation ${LIBS} ${F_LIBRARY})
IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_CCA_Components_Arches_Radiation)
ENDIF(BUILD_SHARED_LIBS)
