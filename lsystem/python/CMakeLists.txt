# python2 libraries
find_package(PkgConfig REQUIRED)
pkg_check_modules(PYTHON REQUIRED python2)

# boost python
find_package(Boost REQUIRED COMPONENTS python27)
message(STATUS "boost: ${Boost_FOUND}")
message(STATUS "boost::system: ${Boost_SYSTEM_FOUND}")
message(STATUS "boost::python: ${Boost_PYTHON_FOUND}")

if (NOT ${Boost_FOUND})

    #message(STATUS "boost: ${Boost_FOUND}")
    #message(STATUS "boost::system: ${Boost_SYSTEM_FOUND}")
    #message(STATUS "boost::python: ${Boost_PYTHON_FOUND}")

    if (${APPLE})
        message(FATAL_ERROR "boost python not found, try to run command\n"
                "$ brew install boost-python")
    else()
        message(FATAL_ERROR "boost python not found")
    endif ()

endif ()

if (${APPLE})
    set(PYTHON_LIBRARIES
            "/usr/local/Cellar/python@2/2.7.15_1/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib"
            ${PYTHON_LIBRARIES} )
endif ()

message("python includes: ${PYTHON_INCLUDE_DIRS}")
message("python libraries: ${PYTHON_LIBRARIES}")

# includes
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# add python wrapper target
add_library(lsystem_python SHARED
        EvolutionIndirect_python.cpp
        )
set_target_properties(lsystem_python PROPERTIES PREFIX "")
set_target_properties(lsystem_python PROPERTIES SUFFIX ".so")

# link libraries
target_link_libraries(lsystem_python
        lsystem_proto_lib
        ${Boost_LIBRARIES}
        ${PYTHON_LIBRARIES}
        )

# install
execute_process ( COMMAND python2 -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(OLD_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CMAKE_INSTALL_PREFIX "/")
install(TARGETS lsystem_python DESTINATION ${PYTHON_SITE_PACKAGES})
#install(TARGETS executable DESTINATION "${OLD_CMAKE_INSTALL_PREFIX}/bin")


message("Boost Libraries: ${Boost_LIBRARIES}")

