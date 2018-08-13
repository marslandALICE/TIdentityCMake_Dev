
FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathMore )
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )
include_directories("./" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

add_executable(testIdentity_Sign  testIdentity_Sign.C)
add_executable(testIdentity_NetParticles testIdentity_NetParticles.C)
add_executable(testIdentity testIdentity.C)
add_executable(testIdentity_3D testIdentity_3D.C)
add_executable(identity_2D_ALICE identity_2D_ALICE.C)
add_executable(identity_3D_ALICE identity_3D_ALICE.C)

target_link_libraries(testIdentity_Sign TIdentity ${ROOT_LIBRARIES})
target_link_libraries(testIdentity_NetParticles TIdentity ${ROOT_LIBRARIES})
target_link_libraries(testIdentity TIdentity ${ROOT_LIBRARIES})
target_link_libraries(testIdentity_3D TIdentity ${ROOT_LIBRARIES})
target_link_libraries(identity_2D_ALICE TIdentity ${ROOT_LIBRARIES})
target_link_libraries(identity_3D_ALICE TIdentity ${ROOT_LIBRARIES})
