#include "slambase.h"
#include "LM_badjust.h"
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/types/icp/types_icp.h>
void setuprepro(g2o::SparseOptimizer* globalOptimizer,CAMERA_INTRINSIC_PARAMETERS cameras);//設定PnP解用的疊代器
bool project_point(vector<KEYFRAME> Map, FRAME &lastframe, FRAME &currframe, CAMERA_INTRINSIC_PARAMETERS camera,Eigen::Isometry3d& resual);//做深度更新與姿態校正
Eigen::Vector2d featrue_al(cv::Mat& orggray, FRAME& currgray, Eigen::Vector2d repro_pix, Eigen::Vector2d pxss_max_curr, Eigen::Vector2d pxss_min_curr
                           , double depth, double sigma, double scale, double& lenger, double& error_temp, bool &match_error);//從epipolar line 找匹配點
void depth_fit(vector<KEYFRAME>& Map, FRAME &lastframe, FRAME &currframe,Eigen::Isometry3d tcw,CAMERA_INTRINSIC_PARAMETERS camera);//更新深度值
bool compete_depth(double mean_1,double sigma_1,double mean_2,double sigma_2);//根據深度分布選擇合適的深度值
