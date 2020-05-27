
set(RADTSE_HEADER_FILES RadTseTrajGadget.h gadgetron_radtse_export.h)
set(RADTSE_SRC_FILES RadTseTrajGadget.cpp)

add_library(gadgetron_radtse SHARED ${RADTSE_HEADER_FILES} ${RADTSE_SRC_FILES})

set_target_properties(gadgetron_radtse PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_radtse
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_mri_core
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore)


install(TARGETS gadgetron_radtse
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install (FILES ${RADTSE_HEADER_FILES}
                     DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)