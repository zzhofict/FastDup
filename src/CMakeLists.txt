# exe path
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

# source codes path
aux_source_directory(${PROJECT_SOURCE_DIR}/src MAIN_SRC)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/util UTIL_SRC)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/markdup MARKDUP_SRC)

set(KTHREAD_FILE ${PROJECT_SOURCE_DIR}/ext/klib/kthread.c)

# including path
include_directories("${PROJECT_SOURCE_DIR}/ext")
include_directories("${PROJECT_SOURCE_DIR}/ext/htslib")
include_directories("${PROJECT_SOURCE_DIR}/src")

# linking path
link_directories("${PROJECT_SOURCE_DIR}/ext/htslib")

# program name
set(PG_NAME "fastdup")

# dependency files
add_executable(${PG_NAME} ${MAIN_SRC} ${UTIL_SRC} ${MARKDUP_SRC} ${KTHREAD_FILE})

# link htslib
target_link_libraries(${PG_NAME} libhts.a)

# lib pthread
find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
    set_property(TARGET ${PG_NAME} PROPERTY COMPILE_OPTIONS "-pthread")
    set_property(TARGET ${PG_NAME} PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(${PG_NAME} "${CMAKE_THREAD_LIBS_INIT}")
endif()

# lib bzip2
find_package(BZip2 REQUIRED)
if(BZip2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIR})
    target_link_libraries(${PG_NAME} ${BZIP2_LIBRARIES})
else()
    message(FATAL_ERROR "BZIP2 is not found")
endif()

# lib deflate
find_library(DeflateLib deflate)
if(DeflateLib)
    target_link_libraries(${PG_NAME} ${DeflateLib})
else()
    message(WARNING "libdeflate is not found")
endif()

# lib lzma
find_package(LibLZMA MODULE)
if (LibLZMA_FOUND)
    include_directories(${LIBLZMA_INCLUDE_DIR})
    target_link_libraries(${PG_NAME} ${LIBLZMA_LIBRARY})
else()
    message(FATAL_ERROR "lzma is not found")
endif()

# lib curl
find_package(CURL REQUIRED)
if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIR})
    target_link_libraries(${PG_NAME} ${CURL_LIBRARY})
endif()

# lz
find_library(Z_LIBRARY z)
if(Z_LIBRARY)
    target_link_libraries(${PG_NAME} ${Z_LIBRARY})
else()
    message(FATAL_ERROR "lz is not found")
endif()

# lm
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(${PG_NAME} ${MATH_LIBRARY})
else()
    message(FATAL_ERROR "lm is not found")
endif()

# ssl
find_package(OpenSSL REQUIRED)
if(OpenSSL_FOUND)
    target_link_libraries(${PG_NAME} OpenSSL::SSL OpenSSL::Crypto)
else()
    message(FATAL_ERROR "OpenSSL is not found")
endif()

# install path
install(TARGETS ${PG_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION libstatic)