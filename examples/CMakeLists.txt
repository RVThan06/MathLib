# Find all .cpp files in this folder
file(GLOB EXAMPLE_SOURCES "*.cpp")

foreach(example_source ${EXAMPLE_SOURCES})
    # Extract the name (e.g., "LinearAlgebraDemo.cpp" -> "LinearAlgebraDemo")
    get_filename_component(example_name ${example_source} NAME_WE)

    # Create a standalone executable for this specific example
    add_executable(${example_name} ${example_source})

    # Link it to your library so it can find <mathlib/...> headers
    target_link_libraries(${example_name} PRIVATE MathLib)

    # Optional: Apply the same strict warnings we use for the library
    if(COMMAND set_project_warnings)
        set_project_warnings(${example_name})
    endif()

    message(STATUS "Registered Example: ${example_name}")
endforeach()
