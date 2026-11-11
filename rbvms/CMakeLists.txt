# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.



add_mfem_miniapp(rbvms-tudship
  MAIN rbvms.cpp
  EXTRA_SOURCES integrator.cpp
  EXTRA_HEADERS integrator.hpp ../tudship-util.hpp
  LIBRARIES mfem)

set_target_properties(rbvms-tudship PROPERTIES OUTPUT_NAME rbvms)

add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/buildInfo.hpp"
                   COMMAND /bin/bash "${PROJECT_SOURCE_DIR}/buildInfo.sh" "${PROJECT_SOURCE_DIR}/buildInfo.hpp"
                   COMMENT "Getting git info")


foreach(target liddriven kovasznay vonkarman)
   add_library(${target} SHARED ${target}.cpp)
   target_include_directories(${target} PRIVATE
       ${PROJECT_SOURCE_DIR}
       ${PROJECT_BINARY_DIR}/config
   )
   add_dependencies(rbvms-tudship ${target})
endforeach()

set(rbvms_FILES
    README.md
    liddriven.cpp
    kovasznay.cpp
    vonkarman.cpp
    von-karman-nurbs.mesh
    von-karman-tri.geo
)

set(COPY_OUTPUTS)

foreach(FILE_ ${rbvms_FILES})
    set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILE_}")
    set(DST "${CMAKE_CURRENT_BINARY_DIR}/${FILE_}")

    add_custom_command(
        OUTPUT "${DST}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${SRC}" "${DST}"
        DEPENDS "${SRC}"
        COMMENT "Copying ${FILE_}"
        VERBATIM
    )

    list(APPEND COPY_OUTPUTS "${DST}")
endforeach()

add_custom_target(rbvms-files ALL
    DEPENDS ${COPY_OUTPUTS}
)

add_dependencies(rbvms-tudship rbvms-files)

if (MFEM_ENABLE_TESTING)
   add_test(NAME rbvms-ldc
         COMMAND  $<TARGET_FILE:rbvms-tudship> 
                  --mesh ${PROJECT_BINARY_DIR}/data/square-nurbs.mesh --strong-bdr 1
                  --lib $<TARGET_FILE:liddriven> -rho 1 -mu 0.001
                  -r 1 -o 2 -s 21 -dt 0.01 -tf 0.1)
   add_test(NAME rbvms-vkvs
         COMMAND  $<TARGET_FILE:rbvms-tudship>
                  --mesh von-karman-nurbs.mesh
                  --strong-bdr "1 2" --outflow-bdr 4 --weak-bdr 3
                  --lib $<TARGET_FILE:vonkarman> -rho 1 -mu 0.01
                  -r 1 -o 2 -s 45 -dt 0.1 -tf 1)
endif()

if (MFEM_USE_MPI)
  add_mfem_miniapp(prbvms-tudship
    MAIN prbvms.cpp
    EXTRA_SOURCES integrator.cpp
    EXTRA_HEADERS integrator.hpp ../tudship-util.hpp
    LIBRARIES mfem)

  set_target_properties(prbvms-tudship PROPERTIES OUTPUT_NAME prbvms)

  # Add the corresponding tests to the "test" target
  if (MFEM_ENABLE_TESTING)
    add_test(NAME prbvms-tudshipldc_np${MFEM_MPI_NP}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MFEM_MPI_NP}
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:prbvms-tudship>
               --mesh ${PROJECT_BINARY_DIR}/data/square-nurbs.mesh --strong-bdr 1
               --lib $<TARGET_FILE:liddriven>  -rho 1 -mu 0.001
               -o 2 -r 3 -s 21 -dt 0.01 -tf 0.1
      ${MPIEXEC_POSTFLAGS})

    add_test(NAME prbvms-tudshipvkvs_np${MFEM_MPI_NP}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MFEM_MPI_NP}
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:prbvms-tudship>
                  --mesh von-karman-nurbs.mesh
                  --strong-bdr "1 2" --outflow-bdr 4 --weak-bdr 3 
                  --lib $<TARGET_FILE:vonkarman> -rho 1 -mu 0.01
                  -r 1 -o 2 -s 45 -dt 0.1  -tf 1
      ${MPIEXEC_POSTFLAGS})
  endif()

endif()

