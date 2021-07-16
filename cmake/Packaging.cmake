# No need for packaging as a subproject
if (NOT ADS_IS_TOP_LEVEL)
  return()
endif()

# Metadata
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "ADS")
set(CPACK_PACKAGE_VENDOR "Marcin Łoś")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/README.md)

set(CPACK_VERBATIM_VARIABLES YES)

# Component archives settings
set(CPACK_ARCHIVE_COMPONENT_INSTALL YES)
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${PROJECT_VERSION}")
set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY YES)
set(CPACK_ARCHIVE_ADS-RUNTIME_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-runtime")
set(CPACK_ARCHIVE_ADS-DEVEL_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-devel")
set(CPACK_ARCHIVE_ADS-TOOLS_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-tools")

# Settings for source package
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git
  /CMakeUserPresets.json
  /.*build.*
  __pycache__/
  \\.swp
  /\\.idea
)

include(CPack)