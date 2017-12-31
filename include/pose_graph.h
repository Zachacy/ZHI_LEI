
#include "slambase.h"
#include "LM_badjust.h"
#include <g2o/core/robust_kernel.h>
#include <g2o/core/robust_kernel_factory.h>
void setuppose_graph(g2o::SparseOptimizer* globalOptimizer,CAMERA_INTRINSIC_PARAMETERS camera);
bool detect_frame(FRAME* currframe, vector<FRAME>& map, int* index);
Eigen::Isometry3d  pose_graph(FRAME& lastFrame ,FRAME loopFrame,vector<FRAME>& map, Eigen::Isometry3d containere );
bool graph_container(FRAME& lastFrame,FRAME* currFrame,vector<FRAME>& frame_map,vector<DISPLAY>&
                     frame_group,CAMERA_INTRINSIC_PARAMETERS& camera,Eigen::Isometry3d totre);
