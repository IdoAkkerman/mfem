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

add_mfem_miniapp(rbvms-ls
  MAIN rbvms-ls.cpp
  EXTRA_SOURCES integrator.cpp dist_solver.cpp
  EXTRA_HEADERS integrator.hpp dist_solver.hpp ../tudship-util.hpp
  LIBRARIES mfem)

add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/buildInfo.hpp"
                   COMMAND /bin/bash "${PROJECT_SOURCE_DIR}/buildInfo.sh" "${PROJECT_BINARY_DIR}/buildInfo.hpp"
                   COMMENT "Getting git info")

add_custom_target(buildInfo
    DEPENDS "${PROJECT_BINARY_DIR}/buildInfo.hpp"
)

add_dependencies(rbvms-tudship buildInfo)

set(rbvms_ls_FILES
    README.md
    dambreak.geo
    dambreak.c
    dambreak.sh
)

set(COPY_OUTPUTS)

foreach(FILE_ ${rbvms_ls_FILES})
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

add_custom_target(rbvms-ls-files ALL
    DEPENDS ${COPY_OUTPUTS}
)

add_dependencies(rbvms-ls rbvms-ls-files)
