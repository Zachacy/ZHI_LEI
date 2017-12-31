

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <chrono>
#include <ctime>
#include <climits>


#include "slambase.h"
#include "LM_badjust.h"
#include "reproject.h"
using namespace std;



boost::mutex io_mutex;//display result of visual odometry  in PCL model
boost::mutex ic_mutex;//display result of visual odometry  in opencv image

//set pcl visualization function
void display_fuc(FRAME* lastFrame, std::vector<KEYFRAME>* map_grop,PointCloud::Ptr dasdasd,Eigen::Isometry3d totre,CAMERA_INTRINSIC_PARAMETERS camera,bool* updaete )
{
     pcl::visualization::CloudViewer viewer("Cloud Viewer"); //set visual world
     PointCloud::Ptr camera_cloud=drawcamera(255.0,0.0,0.0); // build current frame camera model
     PointCloud::Ptr camera_cloud2=drawcamera(0.0,0.0,255.0); // build keyframe camera model
    PointCloud::Ptr outcloud(new PointCloud());
    PointCloud::Ptr pose_camera(new PointCloud());
    *outcloud+= *dasdasd; // add ground turth to pcl

    //if keyframe group increase , we will renew result of visualization.
    int map_size = map_grop->size();
    while(!viewer.wasStopped())
    {
        io_mutex.lock(); // common data need to lock!

        PointCloud::Ptr ouloud (new PointCloud());
         Eigen::Isometry3d dasd = lastFrame->l_t_r_w*lastFrame->g_t_r_w ;
         Eigen::Isometry3d abcd = totre*dasd.inverse();
         pcl::transformPointCloud( *camera_cloud, *pose_camera, abcd.matrix() );
         PointCloud::Ptr drawvex = drawvertex(lastFrame,*map_grop,totre);

         if(map_size!=map_grop->size())
         {
             for(int ind = map_size;ind<map_grop->size();ind++)
             {
                 Eigen::Isometry3d map_frame =  (*map_grop)[ind].g_t_r_w ; //get camera pose in global map
                 Eigen::Isometry3d data_local = totre*map_frame.inverse();//   tranform to ground turth world

                  if(map_grop->size()%4==0) //  this is setting how many keyframe do you show
                   {
                       PointCloud::Ptr pe_camera(new PointCloud());
                       pcl::transformPointCloud( *camera_cloud2, *pe_camera, data_local.matrix() );
                       *outcloud+= *pe_camera;

                       joinPointCloud( outcloud, (*map_grop)[ind].rgb,(*map_grop)[ind].depth, data_local, camera );
                  }
             }

             map_size = map_grop->size();
         }


         if(*updaete)//show camera movement
         {
             PointCloud::Ptr point_pose ( new PointCloud );
             PointT p;
             p.x=0.0;
             p.y=0.0;
             p.z=0.0;
             p.r=255.0;
             p.g=255.0;
             p.b=0.0;
             point_pose->points.push_back( p );
             pcl::transformPointCloud( *point_pose, *point_pose, abcd.matrix() );
             *outcloud +=*point_pose;
             *updaete = false;
         }

         io_mutex.unlock(); // common data need to lock!

         //add all model to variable of visualization
         *ouloud += *outcloud;
         *ouloud +=*pose_camera;
         *ouloud +=*drawvex;
        viewer.showCloud(ouloud);

    }

}
void display_picture(FRAME* lastFrame,FRAME* currFrame,bool* odo_finsh,bool* esti_cv,CAMERA_INTRINSIC_PARAMETERS camera)
{
    while(1)
    {
        ic_mutex.lock();// common data need to lock!
        if(*esti_cv)
        {

        cv::Mat img_show ( currFrame->rgb.rows, currFrame->rgb.cols, CV_8UC3 );
        currFrame->rgb.copyTo(img_show( cv::Rect ( 0,0,currFrame->rgb.cols, currFrame->rgb.rows ) ));

        // show which point project current image successfully.
        for(vector<Measurement*>::iterator pitd  = lastFrame->kp.begin();pitd != lastFrame->kp.end();++pitd)
        {
           if((*pitd)->estim_out)
           {
            Eigen::Vector3d p= lastFrame->l_t_r_w*(*pitd)->pos_world;//project point to current image
            Eigen::Vector3d p0= (*pitd)->pos_world;
            Eigen::Vector2d pixel_efes = project3Dto2D ( p0 ( 0,0 ),p0 ( 1,0 ), p0 ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy );
            Eigen::Vector2d pixel_prev = project3Dto2D ( p ( 0,0 ), p ( 1,0 ), p ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy );

            float b = (*pitd)->colorp.b;
            float g = (*pitd)->colorp.g;
            float r = (*pitd)->colorp.r;
             cv::circle ( img_show, cv::Point2d ( pixel_prev ( 0,0 ), pixel_prev ( 1,0 ) ), 1, cv::Scalar ( b,g,r ), -1 );
          }

        }
        cv::imshow ( "result", img_show );
        }
        ic_mutex.unlock();// common data need to lock!
        if(*esti_cv)//if camera is moving,we will set delay time
        {
           cv::waitKey(5);
           *esti_cv=false;
        }
         if(*odo_finsh)
         {
             break;
         }
    }
}


int main ( int argc, char** argv )
{

    ParameterReader pd; //導入相機參數與影像資校的變數
    int keyframe_index =0; //關鍵幀索引值

    string path_to_dataset = pd.getData("associate_dir");//讀取ground_turth資料的位置
    string associate_file = path_to_dataset + "/associate_with_groundtruth.txt";//讀取檔案
    ifstream fin ( associate_file ); // 載入檔案 
    string rgb_file, depth_file, time_rgb, time_depth,time_groud,tx,ty,tz,x,y,z,w; //根據ground_turth的格式給予對應的變數
    
    fin>>time_rgb>>rgb_file>>time_depth>>depth_file>>time_groud>>tx>>ty>>tz>>x>>y>>z>>w;
    cout<<time_rgb<<rgb_file<<time_depth<<depth_file<<time_groud<<tx<<ty<<tz<<x<<y<<z<<w<<endl;
    Eigen::Vector3d init_shift;
    init_shift[0]=atof(tx.c_str());
    init_shift[1]=atof(ty.c_str());
    init_shift[2]=atof(tz.c_str());
    Eigen::Quaterniond init_quar  =Eigen::Quaterniond(atof(w.c_str()),atof(x.c_str()),atof(y.c_str()),atof(z.c_str())) ;
    Eigen::Isometry3d totre(init_quar); //獲得第一個位置的變換矩陣
    totre(0,3)= init_shift[0];
    totre(1,3)= init_shift[1];
    totre(2,3)= init_shift[2];

    std::map<string,Eigen::Isometry3d > savedata;//儲存演算法的相機軌跡位置(ATE)

    FRAME lastFrame= readFrame( rgb_file,depth_file, time_groud, pd );//讀取RGB-D資料
    CAMERA_INTRINSIC_PARAMETERS camera = getDefaultCamera(); //讀取相機參數
    lastFrame.frameID = keyframe_index;


    std::vector<KEYFRAME> map; //儲存全局的關鍵幀
    //////////////////////////利用first_frame提取特徵後,存放到存放到lastFrame//////////////////////////////
	KEYFRAME first_frame(lastFrame);
    frame_filter(first_frame,lastFrame,camera,keyframe_index);//提取特徵
    lastFrame.pointmp.insert(pair<int,int>(first_frame.frameID,first_frame.kp.size()));//複製特徵給lastFrame
    lastFrame.esmti_filter.insert(pair<int,bool>(first_frame.frameID,true));
        map.push_back(first_frame);
	 map[first_frame.frameID].grad.copyTo(lastFrame.grad);//複製邊緣圖像給lastFrame
	///////////////////////////////////////////////////////////////////

	
	
	///////////////////////計算所有特徵點的總位移量+轉動角度///////////////////////////////////////////
    double leg_point =0;
    double leg_angle =0;
    Eigen::Vector2d center =Eigen::Vector2d(camera.cx,camera.cy) ;
    for(vector<Measurement>::iterator pt_tmp =map[first_frame.frameID].kp.begin();pt_tmp !=map[first_frame.frameID].kp.end();++pt_tmp)
    {
        Eigen::Vector2d ref_now = project3Dto2D ( (*pt_tmp).pos_world( 0,0 ), (*pt_tmp).pos_world( 1,0 ), (*pt_tmp).pos_world( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
        ref_now = ref_now-center;
        leg_point +=ref_now.norm();
        leg_angle +=abs(atan(ref_now[1]/ref_now[0])*180/M_PI);
        lastFrame.kp.push_back(&(*pt_tmp));
    }
    lastFrame.length_p = (double)leg_point/(double)lastFrame.kp.size();
    lastFrame.angel_p =  (double)leg_angle/(double)lastFrame.kp.size();
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	
	//////////////////////////////////////為了顯示ground turth的三維實際位置dasdasd//////////////////////////////////
    Eigen::Matrix3f K;
    K<<camera.fx,0.f,camera.cx,0.f,camera.fy,camera.cy,0.f,0.f,1.0f;
    string groundtruth_file = path_to_dataset + "/groundtruth.txt";  
    ifstream fgrow ( groundtruth_file.c_str() );
    PointCloud::Ptr pe_camera(new PointCloud());
   vector<Eigen::Isometry3d > czxczxczss;
    while(!fgrow.eof())
    {
        string str;
        getline( fgrow, str );
        string timestep,tx,ty,tz,x,y,z,w;
        Eigen::Vector3d tt;
        if (str[0] == '#')
        {
            continue;
        }
         fgrow>>timestep>>tx>>ty>>tz>>x>>y>>z>>w;
         tt[0]=atof(tx.c_str());
         tt[1]=atof(ty.c_str());
         tt[2]=atof(tz.c_str());
          Eigen::Quaterniond qq =Eigen::Quaterniond(atof(w.c_str()),atof(x.c_str()),atof(y.c_str()),atof(z.c_str())) ;
          Eigen::Isometry3d tmp(qq);
          tmp(0,3)= tt[0];
          tmp(1,3)= tt[1];
          tmp(2,3)= tt[2];
        czxczxczss.push_back(tmp);
    }
    PointCloud::Ptr dasdasd (new PointCloud());
    for(Eigen::Isometry3d group :czxczxczss)
    {
        PointCloud::Ptr point_pose ( new PointCloud );
        PointT p;
        p.x=0,0;
        p.y=0,0;
        p.z=0,0;
        p.r=0.0;
        p.g=255.0;
        p.b=0.0;
        point_pose->points.push_back( p );
         pcl::transformPointCloud( *point_pose, *point_pose, group.matrix() );
        *dasdasd += *point_pose;
    }
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool estim_upd =false;//有姿態計算時才顯示的布林值
     bool odo_fish =false;//當全局地圖有更新時才顯示的布林值
     bool esti_cv = false;//有姿態計算時才顯示的布林值
     FRAME currFrame;
   copyframe_img(currFrame,lastFrame);
   boost::thread m_Thread = boost::thread(display_fuc,&lastFrame,&map,dasdasd,totre,camera,&estim_upd);
   boost::thread mc_Thread = boost::thread(display_picture,&lastFrame,&currFrame,&odo_fish,&esti_cv,camera);
   
    while ( !fin.eof() )
    {
        ic_mutex.lock();
        fin>>time_rgb>>rgb_file>>time_depth>>depth_file>>time_groud>>tx>>ty>>tz>>x>>y>>z>>w;// 讀取測試資料庫的圖像位置與ground turth 的數值
        currFrame = readFrame( rgb_file,depth_file,time_groud,pd );
        ic_mutex.unlock();
        // 使用直接法计算相机运动
       lastFrame.chi2_e=poseEstimationDirect (&lastFrame,&currFrame, camera,lastFrame.l_t_r_w);
       estim_upd = true;
       esti_cv =true;
	   
        //feature aligen,refine!!
      // bool good_con = graph_container(lastFrame,&currFrame,map,frame_group,camera,totre);
     //  first_tast = good_con;

        bool ddd =baseline_check(map,lastFrame,currFrame, camera);
     //  feature_arrange(map,&lastFrame,&currFrame, camera,false);


        if(ddd)
        {
            io_mutex.lock();
            ic_mutex.lock();
			
            bool pecicion = project_point(map,lastFrame,currFrame,camera,lastFrame.l_t_r_w); //姿態校正
            keyframe_index++;
            repro_good(map,&lastFrame,&currFrame,camera,false); //濾除投影失敗點並更新深度數值

			////////////////////////////////////更新lastFrame的資訊///////////////////////////////
            copyframe_img(lastFrame,currFrame);//複製當前圖像資訊
            lastFrame.frameID=keyframe_index;
            int ref_count =0;
            KEYFRAME insertframe(lastFrame);
            frame_filter(insertframe,lastFrame,camera,lastFrame.frameID);
            ref_count = insertframe.kp.size();//紀錄最新的keyframe的特徵點數量
			//徹底清除lastFrame的特徵點在記憶體空間保存的數值//
            lastFrame.kp.clear();
            vector<Measurement*>(lastFrame.kp).swap(lastFrame.kp);
			//////////////////////////////////////////////////////////////
            map.push_back(insertframe); //this will shift address!!
            lastFrame.g_t_r_w =lastFrame.l_t_r_w*lastFrame.g_t_r_w; 
            lastFrame.l_t_r_w = Eigen::Isometry3d::Identity();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////在lastFrame標記最新的keyframe數量與連接確認////////////////////////////////////////////////
            if(lastFrame.pointmp.find(keyframe_index) == lastFrame.pointmp.end())
            {
                 lastFrame.pointmp.insert(pair<int,int>(map[lastFrame.frameID].frameID,ref_count));
            }
            if(lastFrame.esmti_filter.find(keyframe_index) == lastFrame.esmti_filter.end())
            {
                lastFrame.esmti_filter.insert(pair<int,bool>(map[lastFrame.frameID].frameID,true));
            }
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
            map[insertframe.frameID].grad.copyTo(lastFrame.grad);
			
            repro_good(map,&lastFrame,&currFrame,camera,true);//因為已經更新過了,所以只要提取確認過的結果
            io_mutex.unlock();
            ic_mutex.unlock();

        }

        Eigen::Isometry3d result =  lastFrame.l_t_r_w*lastFrame.g_t_r_w ;
        Eigen::Isometry3d groudset = totre*result.inverse();
        savedata.insert(std::pair<string,Eigen::Isometry3d>(time_groud,groudset));
    }
    odo_fish=true;
    m_Thread.join();
    mc_Thread.join();
    int i =0;
   /* for(KEYFRAME dra:map)
    {
        Eigen::Isometry3d curr_camera_pos = totre* dra.g_t_r_w.inverse();
        if(i%20==0)
        {
            joinPointCloud( outcloud, dra.rgb,dra.depth, curr_camera_pos, camera );
        }
    }
   pcl::io::savePCDFile( "result.pcd", *outcloud );*/

    string filename="HappyDay.txt";
    ofstream fp;
    /*  fp.open(filename.c_str());//開啟檔案
         for(KEYFRAME fram :map )
         {
             Eigen::Isometry3d experiment = totre*fram.g_t_r_w.inverse();
             Eigen::Vector3d tt =Eigen::Vector3d::Identity();
             tt[0]= experiment(0,3);
             tt[1]= experiment(1,3);
             tt[2]= experiment(2,3);
             Eigen::Quaterniond q_data(experiment.rotation());
             fp << fram.timestep  << " " << tt(0,0) << " " << tt(1,0)<< " " << tt(2,0)<< " " << q_data.x() << " " << q_data.y() << " " << q_data.z() << " " << q_data.w() << endl;
         }
          fp.close();*/
    fp.open(filename.c_str());//開啟檔案
    for(std::map<string,Eigen::Isometry3d>::iterator it_pid=savedata.begin();it_pid !=savedata.end();it_pid ++  ) //將演算法所計算的軌跡都存起來
    {

        Eigen::Vector3d tt =Eigen::Vector3d::Identity();
        tt[0]= it_pid->second(0,3);
        tt[1]= it_pid->second(1,3);
        tt[2]= it_pid->second(2,3);
        Eigen::Quaterniond q_data(it_pid->second.rotation());
        fp << it_pid->first  << " " << tt(0,0) << " " << tt(1,0)<< " " << tt(2,0)<< " " << q_data.x() << " " << q_data.y() << " " << q_data.z() << " " << q_data.w() << endl;
    }
     fp.close();



    return 0;
}
