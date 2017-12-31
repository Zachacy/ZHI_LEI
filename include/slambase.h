#ifndef SLAMBASE_H
#define SLAMBASE_H

#include <fstream>
#include <vector>
#include <map>


using namespace std;

#include <boost/thread/thread.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <chrono>
#include <ctime>
#include <climits>
#include <math.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/video/tracking.hpp>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/pcd_io.h>


#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/solver.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/core/optimizable_graph.h>
#include <g2o/types/slam3d/se3quat.h>
#include <g2o/types/sba/types_six_dof_expmap.h>
#include "g2o/types/sba/types_sba.h"
#include "g2o/types/slam3d/parameter_camera.h"
#include "g2o/types/slam3d/edge_se3_pointxyz_disparity.h"
#include "g2o/solvers/structure_only/structure_only_solver.h"
#include "LM_badjust.h"
#include <g2o/types/icp/types_icp.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>



typedef pcl::PointXYZRGBA PointT;
typedef pcl::PointCloud<PointT> PointCloud;

typedef struct camera_intrinsic_parameters //相機參數
{
    double cx, cy, fx, fy, scale;
}CAMERA_INTRINSIC_PARAMETERS;

struct pointftr;


class Measurement //特徵點的內容
{
public:
   Measurement ( Eigen::Vector3d p, float g,int ID,PointT color) : pos_world ( p ), grayscale ( g ),frameID(ID) ,colorp(color){}
   cv::Mat sub_temp;
   int frameID;
   Eigen::Vector3d pos_world; //當前幀的三維資訊
   Eigen::Vector3d orig_pos; // 原始幀的三維資訊
   Eigen::Vector3d match_point; //投影點的三維資訊
   Eigen::Vector3d match_epip_point;//在做極線搜尋對應點的三維資訊
   double grayscale; //特徵點的灰階值
   PointT colorp; //顏色表示(用來display)
   double emtipro;  //計算特徵點與匹配點的灰階值誤差(用來計算權重)
   double pro_le;  //權重值
   double delta_hub; //huber的門檻值(已經沒使用了)
   double ser_mu; //當前幀的深度值 (基本上只要一更新,我就只會將資料處裡並儲存成當前幀的三維資訊)
   double ser_sigma; //當前幀的深度分佈
   double grad_def; //極線搜索結果與匹配點的灰階值誤差
   double grad_leng; //匹配點周圍的梯度量
   double epi_lenge; // 搜索到匹配點的總長度
   bool tack_point; //搜索過程不正確或是匹配結果有問題,就會false , 這樣這點就會在更新的時候被當成投影失敗的點
   bool estim_out; //還沒更新時,有些特徵點可能已經不符合判斷範圍,所以先用布林值表示此點是否能再做姿態估測
   int renew_thes;
};





class FRAME
{
public:
    cv::Mat rgb, depth,gray,grad; //该帧对应的彩色图与深度图
    int Level ; //圖像金字塔的層數
    vector<cv::Mat> Pry_img;  //圖像金字塔
    vector<Measurement*> kp; //特徵點
    map<int,int> pointmp;//每個幀內有保存多少其他幀所產生的特徵點數量
    map<int,bool> esmti_filter; //當特徵點數量為0的時候將會false 切除關系
    double chi2_e; //匹配的錯誤值
    int frameID; //此幀 ID
    float avg_itenity;//此幀的灰階圖的平均值
    double length_p;//所有特徵點在像素平面的位移量
    double angel_p;//所有特徵點在像素平面的轉動量
    double sigmaz;//權重變異數
    string timestep; //圖像名稱
    PointT color; 
    Eigen::Isometry3d g_t_r_w;// 幀在全局地圖的變換矩陣
    Eigen::Isometry3d l_t_r_w; // 幀與幀之間的變換矩陣


};

class KEYFRAME
{
public:
    cv::Mat rgb, depth,gray,grad; //该帧对应的彩色图与深度图
    int Level; //圖像金字塔的層數
    vector<cv::Mat> Pry_img; //圖像金字塔
    vector<Measurement> kp; //特徵點
    vector<int> esmti_filter; //
    string timestep;
    int frameID;
    float avg_itenity;
    PointT color;
    Eigen::Isometry3d g_t_r_w;// every keyframe in gobel
    KEYFRAME(FRAME tmp) //複製一個實際記憶體儲存的幀, 給只用指標的FRAME 做參考
    {
        tmp.rgb.copyTo(rgb);
        tmp.depth.copyTo(depth);
        tmp.gray.copyTo(gray);
        tmp.Pry_img.assign(Pry_img.begin(),Pry_img.end());
        g_t_r_w =Eigen::Isometry3d(tmp.l_t_r_w*tmp.g_t_r_w);
        color = tmp.color;
        timestep = tmp.timestep;
        frameID = tmp.frameID;
        avg_itenity = tmp.avg_itenity;
        for(std::map<int,bool >::iterator kc =tmp.esmti_filter.begin(); kc!=tmp.esmti_filter.end();kc++)
        {
            if(kc->second)
            {
                esmti_filter.push_back(kc->first);
            }
        }
    }
};



class ParameterReader //讀取相機參數用的方法
{
public:
    ParameterReader( string filename="./parameters.txt" )
    {
        ifstream fin( filename.c_str() );
        if (!fin)
        {
            cerr<<"parameter file does not exist."<<endl;
            return;
        }
        while(!fin.eof())
        {
            string str;
            getline( fin, str );
            if (str[0] == '#')
            {
                // 以‘＃’开头的是注释
                continue;
            }

            int pos = str.find("=");
            if (pos == -1)
                continue;
            string key = str.substr( 0, pos );
            string value = str.substr( pos+1, str.length() );
            data[key] = value;

            if ( !fin.good() )
                break;
        }
    }
    string getData( string key )
    {
        map<string, string>::iterator iter = data.find(key);
        if (iter == data.end())
        {
            cerr<<"Parameter name "<<key<<" not found!"<<endl;
            return string("NOT_FOUND");
        }
        return iter->second;
    }
public:
    map<string, string> data;
};

vector<cv::Mat> getGaussianPyramid(cv::Mat &img, int nLevels); //產生高斯金字塔
int getMaxLayer(cv::Mat &img);//計算最大階層數
Eigen::Vector3d project2Dto3D ( double x, double y, double d, float fx, float fy, float cx, float cy, float scale ); // 2D轉3D
 Eigen::Vector2d project3Dto2D ( float x, float y, float z, float fx, float fy, float cx, float cy );// 3D轉2D
bool depth_grade (cv::Mat& dpethw, int x,int y,float theshold); // 計算深度梯度 並自訂門檻,如果沒有超過門檻,會認為是平滑區域
double poseEstimationDirect (FRAME* grayref,FRAME* graycur, CAMERA_INTRINSIC_PARAMETERS& camera,Eigen::Isometry3d& Tcw );//計算相機姿態
void feature_arrange(vector<KEYFRAME>& ssss,FRAME* grayref,FRAME* graycur, CAMERA_INTRINSIC_PARAMETERS& camera,bool updata);//計算特徵點的權重
void  sobel_extract ( KEYFRAME& img,FRAME& REF_img,camera_intrinsic_parameters K ,int ID,bool check_point);//提取特徵點 (是用canny 不是sobel)
PointCloud::Ptr image2PointCloud( cv::Mat& rgb, cv::Mat& depth, CAMERA_INTRINSIC_PARAMETERS& camera  );//RGB-D轉成點雲
void joinPointCloud(PointCloud::Ptr original, cv::Mat rgb, cv::Mat depth, Eigen::Isometry3d T, CAMERA_INTRINSIC_PARAMETERS& camera );//加入RGB-D到某個點雲變數中
void copyframe_img(FRAME& refer ,const FRAME reft);//複製帧
PointCloud::Ptr tranPointCloud( PointCloud::Ptr original, Eigen::Isometry3d T, CAMERA_INTRINSIC_PARAMETERS& camera );//根據變換矩陣移動點雲
void frame_filter(KEYFRAME& image,FRAME& ref_img ,CAMERA_INTRINSIC_PARAMETERS pd,int ID );//sobel_extract的輸入函式
FRAME readFrame(string rgb_file ,string depth_file,string time_groud, ParameterReader& pd );//讀取圖片資料並做資料上的處理
FRAME readFrame(int index, ParameterReader& pd );
FRAME copyframe_value(const FRAME reft);//複製帧(用不道)
 CAMERA_INTRINSIC_PARAMETERS getDefaultCamera();//獲得相機參數
 PointCloud::Ptr connect_net(KEYFRAME *keyFrame , KEYFRAME accqu_frame, Eigen::Isometry3d totre); //顯示那些關鍵幀有相關(用不道)
PointCloud::Ptr drawcamera(float r ,float g ,float b);//建立相機模型
PointCloud::Ptr drawvertex(FRAME* keyFrame ,vector<KEYFRAME>& all_frame,Eigen::Isometry3d totre);//畫出當前幀中有相關的幀
void patch_patch(cv::Mat& chose_patch ,Eigen::Vector2d pixel_prev,int patch_filter );//設置mark到chose_patch中 ,防止重複點
bool baseline_check(vector<KEYFRAME>& Map, FRAME& lastframe , FRAME &currframe, camera_intrinsic_parameters camera);//確定特徵點在像素平面上的平移跟轉動
void repro_good(vector<KEYFRAME>& ssss, FRAME* grayref, FRAME* graycur, camera_intrinsic_parameters camera , bool update);//濾除不合格的特徵點
#endif // SLAMBASE_H
