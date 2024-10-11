find_package(DataSpaces REQUIRED)

if(DATASPACES_FOUND)
    target_compile_definitions(lammps PRIVATE -DLMP_HAS_DSPACES)
    target_link_libraries(lammps PRIVATE DataSpaces::DataSpaces)
endif()