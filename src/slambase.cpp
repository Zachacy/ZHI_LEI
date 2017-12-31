#include "slambase.h"
/* fr1_desk
 int point_cout_creat = 4000;
 float cannt_min = 1.0;
 float cannt_max = 1.5;
*/
int point_cout_creat = 3000;
float cannt_min = 2.0;
float cannt_max = 1.0;
Eigen::Vector3d project2Dto3D ( double x, double y, double d, float fx, float fy, float cx, float cy, float scale )
{
    float zz = float ( d ) /scale;
    float xx = zz* ( x-cx ) /fx;
    float yy = zz* ( y-cy ) /fy;
    return Eigen::Vector3d ( xx, yy, zz );
}
Eigen::Vector2d project3Dto2D ( float x, float y, float z, float fx, float fy, float cx, float cy )
{
   float z_abs = z;
   float u = fx*x/z_abs+cx;
   float v = fy*y/z_abs+cy;
   return Eigen::Vector2d ( u,v );
}
int getMaxLayer(cv::Mat &img){
    int width = img.cols;
    int height = img.rows;
    int res = 1;
    int p = 1;
    while(1){
        int tmp =  p*p;
        if(width % tmp == 0) ++ p;
        else break;
    }
    res = p;
    p = 1;
    while(1){
        int tmp =  p*p;
        if(height % tmp == 0) ++ p;
        else break;
    }
    res = res < p ? res : p;
    return res;
}

vector<cv::Mat> getGaussianPyramid(cv::Mat &img, int nLevels){
    vector<cv::Mat> pyr;
    pyr.push_back(img);
    for(int i = 0; i < nLevels - 1; i++){
        cv::Mat tmp;
        pyrDown(pyr[pyr.size() - 1], tmp);
        pyr.push_back(tmp);
    }
    return pyr;
}
void copyframe_img(FRAME& refer ,const FRAME reft)
{
    refer.avg_itenity=reft.avg_itenity;
    refer.color= PointT(reft.color);
    reft.depth.copyTo(refer.depth);
    refer.timestep=reft.timestep;
    refer.Pry_img.clear();
    vector<cv::Mat>(refer.Pry_img).swap(refer.Pry_img);
    refer.Pry_img.assign(reft.Pry_img.begin(),reft.Pry_img.end());
    reft.rgb.copyTo(refer.rgb);
    reft.gray.copyTo(refer.gray);
    refer.sigmaz = reft.sigmaz;
}
/*FRAME copyframe_value(const FRAME reft)
{
    FRAME f;
    f.avg_itenity=reft.avg_itenity;
    f.campos = reft.campos;
    f.chi2_e = reft.chi2_e;
    f.grad = reft.grad;
    f.color= PointT(reft.color);
    reft.depth.copyTo(f.depth);
    f.timestep=reft.timestep;
    f.esmti_filter.insert(reft.esmti_filter.begin(),reft.esmti_filter.end());
    f.pointmp.insert(reft.pointmp.begin(),reft.pointmp.end());
    f.esmti_size= reft.esmti_size;
    f.frameID= reft.frameID;
    f.g_t_r_w= Eigen::Isometry3d(reft.g_t_r_w);
    f.l_t_r_w= Eigen::Isometry3d(reft.l_t_r_w);
    f.kp.assign(reft.kp.begin(),reft.kp.end());
    f.Level= reft.Level;
    f.max_problie= reft.max_problie;
    f.min_mean= reft.min_mean;
    f.min_varience= reft.min_varience;
    f.Pry_img.assign(reft.Pry_img.begin(),reft.Pry_img.end());
    reft.rgb.copyTo(f.rgb);
    reft.gray.copyTo(f.gray);
    return f;
}*/
void patch_patch(cv::Mat& chose_patch ,Eigen::Vector2d pixel_prev,int patch_filter )
{
    int roi_filter_x = patch_filter;
    int roi_filter_y = patch_filter;

    if ( pixel_prev(0,0)<patch_filter/2 )
    {
        roi_filter_x = (int)pixel_prev(0,0)-1;
    }
    else if( pixel_prev(0,0)>=chose_patch.cols-patch_filter/2)
    {
        roi_filter_x = chose_patch.cols-pixel_prev(0,0)-1;
    }

    if (pixel_prev(1,0)-10<patch_filter/2)
    {
        roi_filter_y = (int)pixel_prev(1,0)-1;
    }
    else if( pixel_prev(1,0)+10>=chose_patch.rows-patch_filter/2)
    {
        roi_filter_y = chose_patch.rows-pixel_prev(1,0)-1;
    }

    int x =(int)pixel_prev(0,0);
    int y =(int)pixel_prev(1,0);
    for(int i=x-roi_filter_x/2;i<x+roi_filter_x/2;i++)
        for(int j=y-roi_filter_x/2;j<y+roi_filter_y/2;j++)
        {
            chose_patch.ptr<uchar>(j)[i] = 0;
        }
}

bool depth_grade (cv::Mat& dpethw, int x,int y,float theshold)
{
    bool check = false;
    if(dpethw.ptr<ushort>(y)[x]!=0 &&
        abs(dpethw.ptr<ushort>(y+1)[x]-dpethw.ptr<ushort>(y-1)[x])<theshold &&
        abs(dpethw.ptr<ushort>(y)[x+1]-dpethw.ptr<ushort>(y)[x-1])<theshold &&
        abs(dpethw.ptr<ushort>(y+1)[x+1]-dpethw.ptr<ushort>(y-1)[x-1])<theshold &&
        abs(dpethw.ptr<ushort>(y-1)[x+1]-dpethw.ptr<ushort>(y+1)[x-1])<theshold )

    {
        check = true;
    }

    return check;
}
void  sobel_extract ( KEYFRAME& img,FRAME& REF_img,camera_intrinsic_parameters K ,int ID,bool check_point)
{


    cv::Mat chose_patch =cv::Mat::ones(img.gray.rows,img.gray.cols,CV_8UC1);
    int point_threshold = point_cout_creat;
    int roi_filter=(int)sqrt(double(img.gray.rows*img.gray.cols)/(double)point_threshold);
    cv::Mat cannyer;
    cv::Canny(img.gray,cannyer,(int)(img.avg_itenity/cannt_min),(int)img.avg_itenity*cannt_max,3);
    cannyer.copyTo(img.grad);
	////////////////////////如果不是第一個關鍵幀,要考慮投影成功的點///////////////////////////////////
       if(check_point)
       {
           int patch_filter = 3;
           for(Measurement* m : REF_img.kp)
           {
               Eigen::Vector2d pixel_prev = project3Dto2D ( m->pos_world( 0,0 ), m->pos_world ( 1,0 ), m->pos_world ( 2,0 ), K.fx, K.fy, K.cx, K.cy );
               int prev_x= (int)pixel_prev(0,0);
               int prev_y= (int)pixel_prev(1,0);

               patch_patch(chose_patch,pixel_prev,patch_filter);

           }
       }
	   ////////////////////////////////////////////////////////////////////////////////////
	   //////////////////////////////////////////做深度的高斯分布//////////////////////////////////////////
       double gradep_thres = 0;
       double gradep_varien = 0;
       int depth_count=0;
         for ( int x=4; x<cannyer.cols; x++ )
         {
            for ( int y=10; y<cannyer.rows; y++ )
            {
                if(img.depth.ptr<ushort>(y)[x]>0&&cannyer.ptr<uchar>(y)[x]>0)
                {
                    if(check_point)
                    {
                        if(chose_patch.ptr<uchar>(y)[x]!= 0)
                        {
                            gradep_thres += (double)img.depth.ptr<ushort>(y)[x]/K.scale;
                            depth_count++;
                        }
                    }
                    else
                    {
                            gradep_thres += (double)img.depth.ptr<ushort>(y)[x]/K.scale;
                            depth_count++;
                    }
                }
            }
         }
         gradep_thres /=(double)depth_count;


         for ( int x=4; x<cannyer.cols; x++ )
         {
            for ( int y=10; y<cannyer.rows; y++ )
            {
                if(img.depth.ptr<ushort>(y)[x]>0&&cannyer.ptr<uchar>(y)[x]>0)
                {
                    if(check_point)
                    {
                        if(chose_patch.ptr<uchar>(y)[x]!= 0)
                        {
                            gradep_varien += pow((double)img.depth.ptr<ushort>(y)[x]/K.scale-gradep_thres,2);
                        }
                    }
                    else
                    {
                         gradep_varien += pow((double)img.depth.ptr<ushort>(y)[x]/K.scale-gradep_thres,2);
                    }
                }
            }
         }
         gradep_varien = sqrt(gradep_varien/(double)depth_count);
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //////////////////////////////////////////////計算window size 的大小///////////////////////////////////////////////////////
         int extract_co = 0;
         if(check_point)
         {
            extract_co=point_threshold;
         }
         else
         {
             if(REF_img.kp.size()<point_threshold)
             {
                extract_co = point_threshold-REF_img.kp.size();
             }
             else
             {
                 extract_co =(double)point_threshold/10.0;
             }
         }
        roi_filter = depth_count/extract_co;
        if(roi_filter <=0)
        {
            roi_filter =1;
        }
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
       int xx=10;
       int yy=10;
        int point_total =0;
         bool track_y = false;
         int idex =0;
           while(yy<img.gray.rows-10)
           {
               xx=10;
                while(xx<img.gray.cols-10)
                {
                    bool track_x = false;
                    bool check = depth_grade(img.depth,xx,yy,0.1*K.scale); //剔除深度梯度太平滑的區域
                    ushort d = img.depth.ptr<ushort> (yy)[xx];
                    if(cannyer.ptr<uchar>(yy)[xx]>0&&d>0 && check)
                    {

                        Eigen::Vector3d p3d = project2Dto3D ( xx, yy, d, K.fx, K.fy, K.cx, K.cy, K.scale );
                        Eigen::Vector2d checks2 = project3Dto2D ( p3d[0], p3d[1], p3d[2], K.fx, K.fy, K.cx, K.cy);
                        if(cannyer.ptr<uchar>((int)checks2[1])[(int)checks2[0]]>0&&img.depth.ptr<ushort> (yy)[xx]<(8*K.scale)
                          &&(double)img.depth.ptr<ushort>(yy)[xx]/K.scale<gradep_thres+gradep_varien
                          &&(double)img.depth.ptr<ushort>(yy)[xx]/K.scale>gradep_thres-gradep_varien)
                        {
                             double grayscale = (int)(img.gray.ptr<uchar> (yy)[xx]) ;
                            double grad_x = (double)(img.gray.ptr<uchar> (yy)[xx+1])-(double)(img.gray.ptr<uchar> (yy)[xx-1]) ;
                            double grad_y = (double)(img.gray.ptr<uchar> (yy+1)[xx])-(double)(img.gray.ptr<uchar> (yy-1)[xx]) ;
                            Measurement measurements( p3d, grayscale,ID,img.color ) ;
                            img.gray(cv::Rect((int)round(xx-2),(int)round(yy-2),5,5)).convertTo(measurements.sub_temp,CV_32FC1);

                            Eigen::Vector3d org_tmp = p3d;
                            org_tmp /=p3d[2];
                            measurements.match_point = Eigen::Vector3d(0,0,0);
                            measurements.match_epip_point =  Eigen::Vector3d(0,0,0);
                            measurements.estim_out=true;
                            measurements.tack_point = true;
                            measurements.orig_pos = org_tmp;
                            measurements.renew_thes =ID;
                            measurements.emtipro =1.0;
                            measurements.pro_le=1.0;
                            measurements.delta_hub=sqrt(grayscale);
                            measurements.ser_mu=p3d[2];
                            measurements.ser_sigma=0.01;
                            measurements.colorp = img.color;
                            measurements.grad_def = 0.0;
                            measurements.grad_leng = 0.0;
                            measurements.epi_lenge =0.0;
                            img.kp.push_back(measurements);
                        }
                       track_x = true;
                       point_total++;

                    }
                    if(track_x)
                    {

                          xx +=roi_filter;
                        track_y = true;
                    }
                    else
                    {
                        xx++;
                    }
                }
                if(track_y)
                {
                    yy +=roi_filter;
                }
                else{
                    yy++;
                }
           }
}

void frame_filter(KEYFRAME& image,FRAME& ref_img ,CAMERA_INTRINSIC_PARAMETERS pd,int ID )
{
    if(ID==0)
    {
         sobel_extract(image,ref_img,pd, ID,false);
    }
    else
    {
         sobel_extract(image,ref_img,pd, ID,true);
    }
}
bool  baseline_check(vector<KEYFRAME>& Map,FRAME& lastframe,FRAME& currframe ,camera_intrinsic_parameters camera )
{
    double leg_point =0;
    double leg_angle =0;
     Eigen::Isometry3d cur_tran = lastframe.l_t_r_w*lastframe.g_t_r_w;
     double epi_total=0.0;
     double out_false=0.0;
     bool dcheck_update =false;
     Eigen::Vector2d center =Eigen::Vector2d(camera.cx,camera.cy) ;
    for(vector<Measurement*> ::iterator it_pt = lastframe.kp.begin();it_pt!=lastframe.kp.end();++it_pt)
    {
        Eigen::Isometry3d w2c  = lastframe.l_t_r_w;
        Eigen::Vector3d ex_ref =  (*it_pt)->pos_world/(*it_pt)->pos_world[2];
        Eigen::Vector3d mean = w2c*(ex_ref*((*it_pt)->ser_mu));
        Eigen::Vector3d px_min_curr =  w2c*(ex_ref*((*it_pt)->ser_mu - sqrt((*it_pt)->ser_sigma))) ;	// 按最小深度投影的像素
        Eigen::Vector3d px_max_curr = w2c*(ex_ref*((*it_pt)->ser_mu + sqrt((*it_pt)->ser_sigma))) ;
        Eigen::Vector2d pxss_min_curr = project3Dto2D ( px_min_curr( 0,0 ), px_min_curr( 1,0 ), px_min_curr( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
        Eigen::Vector2d pxss_max_curr = project3Dto2D ( px_max_curr( 0,0 ), px_max_curr( 1,0 ), px_max_curr( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
        Eigen::Vector2d mean_pix = project3Dto2D ( mean( 0,0 ), mean( 1,0 ), mean( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
        Eigen::Vector2d threde =  mean_pix-center;
        leg_point +=threde.norm();
        leg_angle +=abs(atan(threde[1]/threde[0]))*180/M_PI;

        double epi_length = (pxss_min_curr-pxss_max_curr).norm();
        epi_total +=epi_length;


        if(mean_pix[0]<0||mean_pix[0]>currframe.gray.cols||mean_pix[1]<0||mean_pix[1]>currframe.gray.rows)
        {

             (*it_pt)->estim_out=false;
            continue;
        }
        else
        {
            (*it_pt)->estim_out=true;
        }
        double d = (double)currframe.depth.ptr<ushort>((int)mean_pix[1])[(int)mean_pix[0]]/camera.scale;
        if( d>px_max_curr[2]||d<px_min_curr[2])
        {
             (*it_pt)->estim_out=false;


        }
        else
        {
            (*it_pt)->estim_out=true;
              out_false+=1.0;
        }


    }
    epi_total /=(double)lastframe.kp.size();
    out_false /=(double)lastframe.kp.size();
    leg_point /=(double)lastframe.kp.size();
    leg_angle /=(double)lastframe.kp.size();
    double leng_th =abs(lastframe.length_p-leg_point)/ lastframe.length_p;
    double leng_age =abs(lastframe.angel_p-leg_angle)/ abs(lastframe.angel_p);

    if(epi_total>4)
   //if(epi_total>3.8)
   //if(epi_total>4||leng_th>0.01||leng_age>0.05)
    //if(epi_total>4||leng_th>0.085||leng_age>0.1)
    {
        dcheck_update=true;
    }
    return dcheck_update;
}


 CAMERA_INTRINSIC_PARAMETERS getDefaultCamera()
{
    ParameterReader pd;
    CAMERA_INTRINSIC_PARAMETERS camera;
    camera.fx = atof( pd.getData( "camera.fx" ).c_str());
    camera.fy = atof( pd.getData( "camera.fy" ).c_str());
    camera.cx = atof( pd.getData( "camera.cx" ).c_str());
    camera.cy = atof( pd.getData( "camera.cy" ).c_str());
    camera.scale = atof( pd.getData( "camera.scale" ).c_str() );
    return camera;
}



PointCloud::Ptr image2PointCloud( cv::Mat& rgb, cv::Mat& depth, CAMERA_INTRINSIC_PARAMETERS& camera  )
{
    PointCloud::Ptr cloud ( new PointCloud );
    ushort maxf =0;
    for (int m = 0; m < depth.rows; m++)
        for (int n=0; n < depth.cols; n++)
        {
            if(rgb.ptr<uchar>(m)[n*3]==0&&rgb.ptr<uchar>(m)[n*3+1]==0&&rgb.ptr<uchar>(m)[n*3+2]==0)
              continue;
            // 获取深度图中(m,n)处的值
            ushort d = depth.ptr<ushort>(m)[n];
            // d 可能没有值，若如此，跳过此点
            if (d < 300)
                continue;
            // d 存在值，则向点云增加一个点
            PointT p;

            // 计算这个点的空间坐标
            p.z = double(d) / camera.scale;
            p.x = (n - camera.cx) * p.z / camera.fx;
            p.y = (m - camera.cy) * p.z / camera.fy;

            // 从rgb图像中获取它的颜色
            // rgb是三通道的BGR格式图，所以按下面的顺序获取颜色
            p.b = rgb.ptr<uchar>(m)[n*3];
            p.g = rgb.ptr<uchar>(m)[n*3+1];
            p.r = rgb.ptr<uchar>(m)[n*3+2];


            // 把p加入到点云中
            cloud->points.push_back( p );
        }
    // 设置并保存点云
    cloud->height = 1;
    cloud->width = cloud->points.size();
    cloud->is_dense = false;
    printf("%u\n",maxf);
    return cloud;
}

void repro_good(vector<KEYFRAME>& ssss,FRAME* grayref,FRAME* graycur,camera_intrinsic_parameters camera ,bool update)
{
	
	//還沒確實更新前 將update = false 並 轉移所有特徵點的位置, 更新後設update = true 直接使用所有轉移好的特徵點
    if(grayref->kp.capacity()!=0||grayref->kp.size()!=0)
       {
        grayref->kp.clear();
        vector<Measurement*>(grayref->kp).swap(grayref->kp);
       }
	   
        double leg_point =0;
        double leg_angle =0;
        cv::Mat chose_patch =cv::Mat::ones(grayref->gray.rows,grayref->gray.cols,CV_8UC1);
        Eigen::Vector2d center =Eigen::Vector2d(camera.cx,camera.cy) ;
        for(std::map<int,bool >::iterator kc =grayref->esmti_filter.begin(); kc!=grayref->esmti_filter.end();kc++) //need to important!!!
        {

                grayref->pointmp[kc->first]=0;
                Eigen::Isometry3d c2w = grayref->l_t_r_w;//*grayref->g_t_r_w*ssss[kc->first].g_t_r_w.inverse();
                Eigen::Isometry3d c2c = grayref->g_t_r_w*ssss[kc->first].g_t_r_w.inverse();
                for(vector<Measurement>::iterator pt_tmp =ssss[kc->first].kp.begin();pt_tmp !=ssss[kc->first].kp.end();++pt_tmp)
                {


                         if(grayref->frameID!=(*pt_tmp).renew_thes)
                         {
                               continue;
                         }

                        if(!(*pt_tmp).tack_point)
                        {
                            continue;
                        }
                        Eigen::Vector3d p = ((*pt_tmp).pos_world/(*pt_tmp).pos_world[2]*(*pt_tmp).ser_mu);
                        Eigen::Vector3d p2 ;
                        if(!update)
                        {
                          p2 = c2w*p;
                        }
                        else
                        {
                            p2 = p;
                        }
                                // this is problem!! affect others point project error!!
                        Eigen::Vector2d pixel_now = project3Dto2D ( p2 ( 0,0 ), p2 ( 1,0 ), p2 ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );

                        if(!update)
                        {
                             if((int)pixel_now(0,0)<10||(int)pixel_now(0,0)>graycur->gray.cols-10||(int)pixel_now(1,0)<10
                                     ||(int)pixel_now(1,0)>graycur->gray.rows-10)
                             {
                                   continue;
                             }
                            /* ushort d = graycur->depth.ptr<ushort> ((int)pixel_now(1,0))[(int)pixel_now(0,0)];
                             if(d !=0)
                             {
                                   Eigen::Vector3d curr_depth = project2Dto3D ( pixel_now(0,0), pixel_now(1,0), (double)d, camera.fx, camera.fy, camera.cx, camera.cy, camera.scale );
                                    Eigen::Vector3d ref_cu_depth  = c2w.inverse()*curr_depth;
                                 if(ref_cu_depth[2]>((*pt_tmp).ser_mu+2.0*sqrt((*pt_tmp).ser_sigma))||ref_cu_depth[2]<((*pt_tmp).ser_mu-2.0*sqrt((*pt_tmp).ser_sigma)))
                                 {
                                    continue;
                                 }

                             }*/
                         }

                             Eigen::Vector3d p3 = p;
                             Eigen::Vector2d ref_now = project3Dto2D ( p3 ( 0,0 ), p3 ( 1,0 ), p3 ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );

                      int prev_x= (int)ref_now(0,0);
                      int prev_y= (int)ref_now(1,0);
                      if(grayref->grad.ptr<uchar>(prev_y)[prev_x]==0&&
                        grayref->grad.ptr<uchar>(prev_y+1)[prev_x]==0&&
                        grayref->grad.ptr<uchar>(prev_y+1)[prev_x-1]==0&&
                        grayref->grad.ptr<uchar>(prev_y+1)[prev_x+1]==0&&
                        grayref->grad.ptr<uchar>(prev_y-1)[prev_x]==0&&
                        grayref->grad.ptr<uchar>(prev_y-1)[prev_x-1]==0&&
                        grayref->grad.ptr<uchar>(prev_y-1)[prev_x+1]==0&&
                        grayref->grad.ptr<uchar>(prev_y)[prev_x+1]==0&&
                        grayref->grad.ptr<uchar>(prev_y)[prev_x-1]==0
                      )
                     {


                        continue;
                     }
                      if(chose_patch.ptr<uchar>((int)ref_now[1])[(int)ref_now[0]]== 0)
                      {

                          continue;
                      }
                      else
                      {
                           int patch_filter=5;
                           patch_patch(chose_patch,ref_now,patch_filter);
                      }
                      if(!update)
                      {

                          double d_0 =(*pt_tmp).ser_mu;
                          (*pt_tmp).ser_mu = p2[2];
                          (*pt_tmp).ser_sigma =pow((p2[2]/d_0),4)*(*pt_tmp).ser_sigma ;
                          (*pt_tmp).pos_world = p2;
                          (*pt_tmp).renew_thes =grayref->frameID+1;
                          (*pt_tmp).estim_out = true; //new line

                      }
                      if(update)
                      {
                          Eigen::Vector2d havegood= pixel_now-center;
                         leg_point +=havegood.norm();
                         leg_angle +=abs(atan(havegood[1]/havegood[0])*180/M_PI);
                      }

                    grayref->kp.push_back(&(*pt_tmp));
                    grayref->pointmp[kc->first]++;


                }


        }
        if(update)
        {
            grayref->length_p = (double)leg_point/(double)grayref->kp.size();
            grayref->angel_p =  (double)leg_angle/(double)grayref->kp.size();
        }

            for(std::map<int,int >::iterator kc =grayref->pointmp.begin(); kc!=grayref->pointmp.end();kc++)
            {

                      if(kc->first==grayref->frameID)
                      {
                          continue;
                      }


                Eigen::Isometry3d c2w =grayref->l_t_r_w*grayref->g_t_r_w*ssss[kc->first].g_t_r_w.inverse();
                double thres_leg = sqrt(camera.fx*camera.fx+camera.fy*camera.fy);
                double point_weight = (double)kc->second/(double)ssss[kc->first].kp.size();

                if(kc->second<=3)
                {
                    grayref->esmti_filter.erase(kc->first);
                    grayref->pointmp.erase(kc);
                }
            }


}
void feature_arrange(vector<KEYFRAME>& ssss,FRAME* grayref,FRAME* graycur, CAMERA_INTRINSIC_PARAMETERS& camera,bool updata)
{
    double g2o_avager=0;
    int g2oecount=0;
    double g2o_variance = 0;

    for(vector<Measurement*> ::iterator it_pt = grayref->kp.begin();it_pt!=grayref->kp.end();++it_pt)
    {

            Eigen::Vector3d p = (*it_pt)->pos_world;
            Eigen::Vector3d p2 = grayref->l_t_r_w*p;
            Eigen::Vector2d pixel_now = project3Dto2D ( p2 ( 0,0 ), p2 ( 1,0 ), p2 ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );

            if(pixel_now[0]<5||pixel_now[0]>grayref->gray.cols-5||pixel_now[1]<5||pixel_now[1]>grayref->gray.rows-5)
            {
                (*it_pt)->estim_out=false;
                continue;
            }
             (*it_pt)->emtipro =pow((double)(((*it_pt)->grayscale)-(double)graycur->gray.ptr<uchar>( (int)pixel_now(1,0) )[int(pixel_now(0,0))]),2);
              g2o_avager += (*it_pt)->emtipro*(6/(5+((*it_pt)->emtipro/grayref->sigmaz)));
              g2oecount++;
    }
    g2o_avager /=(double)g2oecount;
    grayref->sigmaz = g2o_avager;

    for(vector<Measurement*> ::iterator it_pt = grayref->kp.begin();it_pt!=grayref->kp.end();++it_pt)
    {
        if((*it_pt)->estim_out)
        {

            (*it_pt)->pro_le = (6/(5+((*it_pt)->emtipro/grayref->sigmaz)));
            if((*it_pt)->pro_le <=0.0)
            {
                (*it_pt)->pro_le =1.0;
            }
        }

    }

}

double poseEstimationDirect (FRAME* grayref,FRAME* graycur, CAMERA_INTRINSIC_PARAMETERS& camera,Eigen::Isometry3d& Tcw )
{
    // 初始化g2o
	////////////////////用直接法計算的前提///////////////////////////////////////////////
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6,1>> DirectBlock;  // 求解的向量是6＊1的
    DirectBlock::LinearSolverType* linearSolver = new g2o::LinearSolverDense< DirectBlock::PoseMatrixType > ();
    DirectBlock* solver_ptr = new DirectBlock ( linearSolver );
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg ( solver_ptr ); // L-M
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm ( solver );
    optimizer.setVerbose( false );
    g2o::VertexSE3Expmap* pose = new g2o::VertexSE3Expmap();


    pose->setEstimate ( g2o::SE3Quat ( Tcw.rotation(), Tcw.translation() ) );
    pose->setId ( 0 );
    optimizer.addVertex ( pose );
	//////////////////////////////////////////////////////////////////////////////////////////////
    int level = grayref->Level-1;
    int* Layer_ptr = &level;
    // 添加边
    int id=1;

    for(Measurement* it_pt:grayref->kp)
    {


        if( it_pt->estim_out)
        {
        EdgeSE3ProjectDirect*  edge = new EdgeSE3ProjectDirect (
            it_pt->pos_world, camera.fx, camera.fy, camera.cx, camera.cy,
            grayref->Pry_img,graycur->Pry_img,Layer_ptr );
        edge->setVertex ( 0, pose );
        float offset = graycur->avg_itenity / grayref->avg_itenity ;
        edge->setMeasurement ( it_pt->grayscale );

        Eigen::Matrix<double,1,1> imform = Eigen::Matrix<double,1,1>::Identity()*(it_pt->pro_le);
         edge->setInformation ( imform );
        /* g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber();
         rk->setDelta(it_pt->delta_hub);
        edge->setRobustKernel(rk);*/

        edge->setId ( id++ );

        optimizer.addEdge ( edge );
        }

    }
    optimizer.initializeOptimization();
    for(int n_level =grayref->Level-1;n_level>=0;n_level--) /////每當某一層疊代完畢後,會進入到下一層繼續疊代//////////
    {
        level = n_level ;

        optimizer.optimize ( 100 );
    }
    Tcw = pose->estimate();
    return optimizer.chi2();
}


void joinPointCloud( PointCloud::Ptr original, cv::Mat rgb,cv::Mat depth, Eigen::Isometry3d T, CAMERA_INTRINSIC_PARAMETERS& camera )
{
    PointCloud::Ptr newCloud = image2PointCloud( rgb, depth, camera );

    // 合并点云
    PointCloud::Ptr output (new PointCloud());
    pcl::transformPointCloud( *newCloud, *output, T.matrix() );

    // Voxel grid 滤波降采样

    static pcl::VoxelGrid<PointT> voxel;
    static ParameterReader pd;
    double gridsize = 0.01;
    voxel.setLeafSize( gridsize, gridsize, gridsize );
    voxel.setInputCloud( output );
    PointCloud::Ptr tmp( new PointCloud() );
    voxel.filter( *tmp );

    *original += *tmp;

}

PointCloud::Ptr tranPointCloud(PointCloud::Ptr original, Eigen::Isometry3d T, CAMERA_INTRINSIC_PARAMETERS& camera )
{
    PointCloud::Ptr output (new PointCloud());
    pcl::transformPointCloud( *original, *output, T.matrix() );
    return output;
}


FRAME readFrame(string rgb_file ,string depth_file,string time_groud, ParameterReader& pd )
{
    FRAME f;


    string path_to_dataset = pd.getData("associate_dir");
    f.rgb = cv::imread( path_to_dataset+"/"+rgb_file  );
    f.depth = cv::imread(  path_to_dataset+"/"+depth_file, -1 );
   /* stringstream ss;
    ss<<rgbDir<<index<<rgbExt;
    string filename;
    ss>>filename;

    f.rgb = cv::imread( filename );

    ss.clear();
    filename.clear();
    ss<<depthDir<<index<<depthExt;
    ss>>filename;

    f.depth = cv::imread( filename, -1 );
    */
    f.timestep=time_groud;
    f.frameID = -1;
     srand(time(NULL));
    f.color.r = rand() % 255 ;
    f.color.g = rand() % 255 ;
    f.color.b = rand() % 255 ;
    cv::cvtColor ( f.rgb, f.gray, cv::COLOR_BGR2GRAY );
    CvScalar scalar = mean(f.gray);
    f.avg_itenity = scalar.val[0];
    int n_layer = getMaxLayer( f.gray);
    f.Level = 5;
    f.Pry_img = getGaussianPyramid(f.gray,f.Level);
    f.g_t_r_w = Eigen::Isometry3d::Identity();
    f.l_t_r_w = Eigen::Isometry3d::Identity();
    f.sigmaz =f.avg_itenity*f.avg_itenity  ;
    return f;
}

FRAME readFrame(int index, ParameterReader& pd )
{
    FRAME f;

    string rgbDir   =   pd.getData("rgb_dir");
    string depthDir =   pd.getData("depth_dir");

    string rgbExt   =   pd.getData("rgb_extension");
    string depthExt =   pd.getData("depth_extension");

   stringstream ss;
    ss<<rgbDir<<index<<rgbExt;
    string filename;
    ss>>filename;

    f.rgb = cv::imread( filename );

    ss.clear();
    filename.clear();
    ss<<depthDir<<index<<depthExt;
    ss>>filename;

    f.depth = cv::imread( filename, -1 );

    f.frameID = index;
    f.color.r = rand() % 255 ;
    f.color.g = rand() % 255 ;
    f.color.b = rand() % 255 ;
    cv::cvtColor ( f.rgb, f.gray, cv::COLOR_BGR2GRAY );
    CvScalar scalar = mean(f.gray);
    f.avg_itenity = scalar.val[0];
    int n_layer = getMaxLayer( f.gray);
    f.Level = n_layer;
    f.Pry_img = getGaussianPyramid(f.gray,f.Level);
    f.g_t_r_w = Eigen::Isometry3d::Identity();
    f.l_t_r_w = Eigen::Isometry3d::Identity();

    f.chi2_e=0;
    return f;
}

PointCloud::Ptr drawcamera(float r ,float g ,float b)
{
    PointCloud::Ptr cloud ( new PointCloud );

        for (int m = 0; m < 50; m++)
        {
            PointT p;
            p.z = 0.0;
            p.x = -0.05+(float)(m*0.002);
            p.y = 0.025-(float)(m*0.001);
            p.b = b;
            p.g = g;
            p.r = r;

            PointT p2;
            p2.z = 0.0;
            p2.x = 0.05-(float)(m*0.002);
            p2.y = 0.025-(float)(m*0.001);
            p2.b = p.b;
            p2.g = p.g;
            p2.r = p.r;

            PointT p3;
            p3.z = 0.0;
            p3.x = 0.05-(float)(m*0.002);
            p3.y = 0.025;
            p3.b = p.b;
            p3.g = p.g;
            p3.r = p.r;

            PointT p4;
            p4.z = 0.0;
            p4.x = 0.05-(float)(m*0.002);
            p4.y = -0.025;
            p4.b = p.b;
            p4.g = p.g;
            p4.r = p.r;

            PointT p5;
            p5.z = 0.0;
            p5.x = 0.05;
            p5.y = 0.025-(float)(m*0.001);
            p5.b = p.b;
            p5.g = p.g;
            p5.r = p.r;

            PointT p6;
             p6.z = 0.0;
             p6.x = -0.05;
             p6.y = 0.025-(float)(m*0.001);
             p6.b = p.b;
             p6.g = p.g;
             p6.r = p.r;

             PointT p7;
             p7.z = 0.0-(float)(m*0.001);
             p7.x = 0.05-(float)(m*0.001);
             p7.y = 0.025-(float)(m*0.0005);
             p7.b = p.b;
             p7.g = p.g;
             p7.r = p.r;
             PointT p8;
             p8.z = 0.0-(float)(m*0.001);
             p8.x = -0.05+(float)(m*0.001);
             p8.y = -0.025+(float)(m*0.0005);
             p8.b = p.b;
             p8.g = p.g;
             p8.r = p.r;
             PointT p9;
             p9.z = 0.0-(float)(m*0.001);
             p9.x = 0.05-(float)(m*0.001);
             p9.y = -0.025+(float)(m*0.0005);
             p9.b = p.b;
             p9.g = p.g;
             p9.r = p.r;

             PointT p10;
             p10.z = 0.0-(float)(m*0.001);
             p10.x = -0.05+(float)(m*0.001);
             p10.y = 0.025-(float)(m*0.0005);
             p10.b = p.b;
             p10.g = p.g;
             p10.r = p.r;


            cloud->points.push_back( p );
            cloud->points.push_back( p2 );
            cloud->points.push_back( p3 );
            cloud->points.push_back( p4 );
            cloud->points.push_back( p5 );
            cloud->points.push_back( p6 );
            cloud->points.push_back( p7 );
            cloud->points.push_back( p8 );
            cloud->points.push_back( p9 );
            cloud->points.push_back( p10 );

        }

        cloud->height = 1;
        cloud->width = cloud->points.size();
        cloud->is_dense = false;

        return cloud;
}

PointCloud::Ptr connect_net(KEYFRAME* keyFrame ,KEYFRAME accqu_frame,Eigen::Isometry3d totre)
{
    PointCloud::Ptr cloud ( new PointCloud );
    Eigen::Vector3d origen =Eigen::Vector3d(0.0,0.0,0.0);

            Eigen::Isometry3d curr_pos = keyFrame->g_t_r_w ;
            Eigen::Isometry3d abcd = totre*curr_pos.inverse();
            Eigen::Isometry3d frame_pos =accqu_frame.g_t_r_w ;
            Eigen::Isometry3d ddcd = totre*frame_pos.inverse();

            Eigen::Vector3d step_line=(abcd*origen)-ddcd*origen;
            Eigen::Vector3d line_p0 =ddcd*origen;
            double step_x = step_line[0]/100.0;
            double step_y = step_line[1]/100.0;
            double step_z = step_line[2]/100.0;
            for (int m = 0; m < 100; m++)
            {
                PointT p_line;
                p_line.x = line_p0[0]+(float)(step_x*m);
                p_line.y = line_p0[1]+(float)(step_y*m);
                p_line.z = line_p0[2]+(float)(step_z*m);
                p_line.b = accqu_frame.color.b;
                p_line.g = accqu_frame.color.g;
                p_line.r = accqu_frame.color.r;
                cloud->points.push_back( p_line );
            }

    cloud->height = 1;
    cloud->width = cloud->points.size();
    cloud->is_dense = false;
    return cloud;
}

PointCloud::Ptr drawvertex(FRAME* keyFrame ,vector<KEYFRAME>& all_frame,Eigen::Isometry3d totre)
{
    PointCloud::Ptr cloud ( new PointCloud );
    Eigen::Vector3d origen =Eigen::Vector3d(0.0,0.0,0.0);
    for(std::map<int,bool>::iterator kc =keyFrame->esmti_filter.begin(); kc!=keyFrame->esmti_filter.end();kc++)
    {


            Eigen::Isometry3d curr_pos = keyFrame->l_t_r_w*keyFrame->g_t_r_w ;
            Eigen::Isometry3d abcd = totre*curr_pos.inverse();
            Eigen::Isometry3d frame_pos =all_frame[kc->first].g_t_r_w ;
            Eigen::Isometry3d ddcd = totre*frame_pos.inverse();

            Eigen::Vector3d step_line=(abcd*origen)-ddcd*origen;
            Eigen::Vector3d line_p0 =ddcd*origen;
            double step_x = step_line[0]/100.0;
            double step_y = step_line[1]/100.0;
            double step_z = step_line[2]/100.0;
            for (int m = 0; m < 100; m++)
            {
                PointT p_line;
                p_line.x = line_p0[0]+(float)(step_x*m);
                p_line.y = line_p0[1]+(float)(step_y*m);
                p_line.z = line_p0[2]+(float)(step_z*m);
                p_line.b = all_frame[kc->first].color.b;
                p_line.g = all_frame[kc->first].color.g;
                p_line.r = all_frame[kc->first].color.r;
                cloud->points.push_back( p_line );
            }


    }
    cloud->height = 1;
    cloud->width = cloud->points.size();
    cloud->is_dense = false;
    return cloud;
}
