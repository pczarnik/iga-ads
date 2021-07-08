if (NOT ADS_USE_SANITIZERS)
  return()
endif()

list(JOIN ADS_USE_SANITIZERS "," SANITIZER_LIST)
set(SANITIZER_FLAG "-fsanitize=${SANITIZER_LIST}")

target_compile_options(project-options INTERFACE ${SANITIZER_FLAG})
target_link_options(project-options INTERFACE ${SANITIZER_FLAG})

target_compile_options(project-options INTERFACE "-fno-omit-frame-pointer")

# Build types are case-insensitive
string(TOUPPER "${CMAKE_BUILD_TYPE}" UPPER_BUILD_TYPE)

# Enable symbol infor for build types that usually omit it
# Without it stack traces have no line numbers
if (NOT UPPER_BUILD_TYPE STREQUAL "DEBUG" AND
    NOT UPPER_BUILD_TYPE STREQUAL "RELWITHDEBINFO")
  target_compile_options(project-options INTERFACE -g)
endif()
