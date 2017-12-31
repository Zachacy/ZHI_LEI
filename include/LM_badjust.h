#ifndef LM_BADJUST_H
#define LM_BADJUST_H


#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <chrono>
#include <ctime>
#include <climits>
#include <vector>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/core/optimizable_graph.h>
#include <g2o/types/sba/types_six_dof_expmap.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace g2o;

class EdgeSE3ProjectDirect: public BaseUnaryEdge< 1,double, VertexSE3Expmap>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectDirect() {}

    EdgeSE3ProjectDirect (
            Eigen::Vector3d point,
            float fx,
            float fy,
            float cx,
            float cy,
            vector<cv::Mat>& grayref,
            vector<cv::Mat>& graycur,
            int* Layer);

    virtual void computeError();

    // plus in manifold
    virtual void linearizeOplus( );


    // dummy read and write functions because we don't care...
    virtual bool read ( std::istream& in ) {}
    virtual bool write ( std::ostream& out ) const {}

protected:
    // get a gray scale value from reference image (bilinear interpolated)
    inline float getPixelValue ( cv::Mat img, float x, float y );



public:
    Eigen::Vector3d x_world_;   // 3D point in world frame
    float cx_=0, cy_=0, fx_=0, fy_=0; // Camera intrinsics
    vector<cv::Mat> graycur_ ; //當前幀的圖像金字塔
    vector<cv::Mat> grayref_ ;//當參考幀的圖像金字塔
    int* Layer_ =nullptr;
    int window_roi_ =0;

};
#endif // LM_BADJUST_H
