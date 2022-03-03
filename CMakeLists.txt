project(pythiaChargedDijet VERSION 1.0)

add_subdirectory(src)

list(APPEND EXTRA_LIBS ~/alice/PythiaMy/pythia8303)
list(APPEND EXTRA_INCL ~/alice/PythiaMy/pythia8303/include)

add_executable(pythiaChargedDijet pythiaChargedDijet.C)

target_link_directories(Tutorial PUBLIC ${EXTRA_LIBS})
target_link_libraries(Tutorial PUBLIC ${EXTRA_LIBS})
