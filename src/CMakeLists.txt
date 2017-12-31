
find_package( OpenCV 2.4 REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(CUDA_STATIC_CUDA_RUNTIME=OFF)
FIND_PACKAGE( PCL REQUIRED COMPONENTS common io visualization filters )
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4") 
ADD_DEFINITIONS( ${PCL_DEFINITIONS} )
INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS}  )
LINK_LIBRARIES( ${PCL_LIBRARY_DIRS} )


# 添加cmake模块路径
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

find_package( G2O REQUIRED)
find_package(CSparse REQUIRED) 
find_package(Cholmod REQUIRED) 
find_package(Eigen3 REQUIRED)
include_directories( ${G2O_INCLUDE_DIRS} ) 
INCLUDE_DIRECTORIES(${CSPARSE_INCLUDE_DIR})
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories(${CHOLMOD_INCLUDE_DIR}) 




set( G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse
)



ADD_LIBRARY( lm_badjust LM_badjust.cpp )
TARGET_LINK_LIBRARIES( lm_badjust
    ${OpenCV_LIBS}
    ${PCL_LIBRARIES}
    ${G2O_LIBS}
    ${CHOLMOD_LIBRARIES}
    )

ADD_LIBRARY( slamBase slambase.cpp )
TARGET_LINK_LIBRARIES( slamBase
    lm_badjust
    ${OpenCV_LIBS}
    ${PCL_LIBRARIES}
    ${G2O_LIBS}
    ${CHOLMOD_LIBRARIES}
    )
ADD_LIBRARY( reproject reproject.cpp )
TARGET_LINK_LIBRARIES( reproject
    lm_badjust
    ${OpenCV_LIBS}
    ${PCL_LIBRARIES}
    ${G2O_LIBS}
    ${CHOLMOD_LIBRARIES}
    )
	
ADD_EXECUTABLE( slamsystem slamsystem.cpp )
TARGET_LINK_LIBRARIES( slamsystem
    slamBase
	reproject
    ${OpenCV_LIBS}
    ${PCL_LIBRARIES}
    ${G2O_LIBS}
    ${CHOLMOD_LIBRARIES}
    )



