
configure_file(drivencavity.control drivencavity.control)

add_test(NAME DrivenCavity_explicit_stability WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_BINARY_DIR}/insac2ds
  ${CMAKE_CURRENT_BINARY_DIR}/drivencavity.control
  )
