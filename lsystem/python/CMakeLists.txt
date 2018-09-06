# python3 libraries
find_package(PythonLibs REQUIRED)

# boost python
find_package(Boost REQUIRED COMPONENTS python)
message(STATUS "boost: ${Boost_FOUND}")
message(STATUS "boost::system: ${Boost_SYSTEM_FOUND}")
message(STATUS "boost::python: ${Boost_PYTHON_FOUND}")

if (NOT ${Boost_FOUND})
    if (${APPLE})
        message(FATAL_ERROR "boost python not found, try to run command\n"
                "$ brew install boost-python")
    endif ()

    if (NOT ${Boost_FOUND})
        message(FATAL_ERROR "boost python not found")
    endif ()
    message(STATUS "boost: ${Boost_FOUND}")
    message(STATUS "boost::system: ${Boost_SYSTEM_FOUND}")
    message(STATUS "boost::python: ${Boost_PYTHON_FOUND}")
endif ()

if (${APPLE})
    set(PYTHON_INCLUDE_DIRS
            "/usr/include/python2.7"
            ${PYTHON_INCLUDE_DIRS})
    set(Python_LIBRARIES
            "/usr/local/Cellar/python@2/2.7.15_1/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config/libpython2.7.dylib"
            ${Python_LIBRARIES} )
   # set(Python_LIBRARIES ${Python_LIBRARIES}
   #          "/usr/local/Cellar/python@2/2.7.15_1/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib")
endif ()

# includes
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# add python wrapper target
python_add_module(lsystem_python
        EvolutionIndirect_python.cpp
        )

# link libraries
target_link_libraries(lsystem_python
        lsystem_proto_lib
        ${Boost_LIBRARIES}
        ${Python_LIBRARIES}
        )

# install
#install(TARGETS revolve_brain_python DESTINATION lib)


message("Boost Libraries: ${Boost_LIBRARIES}")

