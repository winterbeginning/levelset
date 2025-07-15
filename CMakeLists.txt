cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(levelset LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG")

if(DEFINED ENV{WM_PROJECT_DIR})
	MESSAGE(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR "Cannot find OpenFOAM")
endif(DEFINED ENV{WM_PROJECT_DIR})

set(FOAMVERSION $ENV{WM_PROJECT_VERSION}) 
set(FOAMDIR $ENV{WM_PROJECT_DIR})
set(FOAMLIB $ENV{FOAM_LIBBIN})
set(FOAMSRC $ENV{FOAM_SRC})
set(FOAMMPI $ENV{FOAM_MPI})
set(FOAM_USER_APPBIN $ENV{FOAM_USER_APPBIN})
set(FOAM_USER_LIBBIN $ENV{FOAM_USER_LIBBIN})
set(WM_PROJECT_USER_DIR $ENV{WM_PROJECT_USER_DIR})

set(FOAMFLAG "-std=c++14 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wold-style-cast
-Wnon-virtual-dtor -Wno-unused-variable -Wno-invalid-offsetof -DNoRepository -ftemplate-depth-100 -fPIC")

# 验证必需的环境变量
foreach(var IN ITEMS FOAM_USER_APPBIN FOAM_USER_LIBBIN FOAMSRC WM_PROJECT_USER_DIR)
    if(NOT ${var})
        message(FATAL_ERROR "Environment variable not set: ${var}")
    endif()
endforeach()

set(FOAMINC
    ${FOAMSRC}/finiteVolume/lnInclude
    ${FOAMSRC}/dynamicFvMesh/lnInclude
    ${FOAMSRC}/meshTools/lnInclude
    ${FOAMSRC}/sampling/lnInclude
    ${FOAMSRC}/OpenFOAM/lnInclude
    ${FOAMSRC}/OSspecific/POSIX/lnInclude
)

# 添加链接目录
link_directories(
    ${FOAM_USER_LIBBIN}
    ${FOAMLIB}
    ${FOAMLIB}/${FOAMMPI}
)

# 添加头文件搜索路径
include_directories(${FOAMINC} .)

# 创建可执行文件
add_executable(${PROJECT_NAME} levelset.C)  #change source file name


add_definitions(${FOAMFLAG})
# 设置输出目录
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${FOAM_USER_APPBIN}
)

target_include_directories(${PROJECT_NAME} PUBLIC ${FOAMINC})

# 添加链接库
target_link_libraries(${PROJECT_NAME} PUBLIC
    finiteVolume
    dynamicFvMesh
    topoChangerFvMesh
    fvModels
    fvConstraints
    meshTools
    transportModels
    OpenFOAM
    Pstream
)
