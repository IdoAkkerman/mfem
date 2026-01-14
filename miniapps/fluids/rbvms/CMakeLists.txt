# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.


##target_link_libraries(rbvms PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})

add_mfem_miniapp(rbvms
  MAIN rbvms.cpp 
  EXTRA_SOURCES integrator.cpp
  EXTRA_HEADERS integrator.hpp
  LIBRARIES mfem)

if (MFEM_ENABLE_TESTING)
  add_test(NAME rbvms_ser COMMAND rbvms -no-vis)
endif()

#if (MFEM_USE_MPI)
#  add_mfem_miniapp(prbvms
#  MAIN prbvms.cpp
#  EXTRA_SOURCES integrator.cpp
#  EXTRA_HEADERS integrator.hpp
#  LIBRARIES mfem)

  # Add the corresponding tests to the "test" target
#  if (MFEM_ENABLE_TESTING)
#    add_test(NAME prbvms_np${MFEM_MPI_NP}
#      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MFEM_MPI_NP}
#      ${MPIEXEC_PREFLAGS}
#      $<TARGET_FILE:prbvms> -no-vis
#      ${MPIEXEC_POSTFLAGS})
#  endif()
#endif()




