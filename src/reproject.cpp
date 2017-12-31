#include "reproject.h"
#include <g2o/types/icp/types_icp.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
void setuprepro(g2o::SparseOptimizer* globalOptimizer,CAMERA_INTRINSIC_PARAMETERS cameras)
{
            g2o::BlockSolver_6_3::LinearSolverType* linearSolver = new  g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType> ();
                // 6*3 的参数
            g2o::BlockSolver_6_3* block_solver = new g2o::BlockSolver_6_3( linearSolver );
                // L-M 下降
            g2o::OptimizationAlgorithmLevenberg* algorithm = new g2o::OptimizationAlgorithmLevenberg( block_solver );
            globalOptimizer->setAlgorithm(algorithm);
             globalOptimizer->setVerbose( false );
             g2o::CameraParameters* camer = new g2o::CameraParameters(  cameras.fx, Eigen::Vector2d(cameras.cx, cameras.cy), 0 );
             camer->setId(0);

             globalOptimizer->addParameter( camer );

}


bool compete_depth(double mean_1,double sigma_1,double mean_2,double sigma_2)
{
    bool resule= false;
    if(mean_1+2.0*sigma_1>mean_2 &&mean_1-2.0*sigma_1<mean_2 )
    {
        return true;
    }

    if(mean_1+2.0*sigma_1>mean_2+2.0*sigma_2 &&mean_1-2.0*sigma_1<mean_2+2.0*sigma_2 )
    {
        return true;
    }
     if(mean_1+2.0*sigma_1>mean_2-2.0*sigma_2 &&mean_1-2.0*sigma_1<mean_2-2.0*sigma_2 )
     {
         return true;
     }
    return resule;
}


void depth_fit(vector<KEYFRAME>& Map, FRAME &lastframe, FRAME &currframe,Eigen::Isometry3d tcw,CAMERA_INTRINSIC_PARAMETERS camera)
{

    Eigen::Isometry3d cur_tran = tcw*lastframe.g_t_r_w;
    for(vector<Measurement*> ::iterator it_pt = lastframe.kp.begin();it_pt!=lastframe.kp.end();++it_pt)
    {
        if((*it_pt)->tack_point)
        {

             Eigen::Isometry3d w2c = tcw;

            Eigen::Vector3d orgipx_point =  w2c.inverse()*(*it_pt)->match_point;
            Eigen::Vector3d eps_pie_point =  w2c.inverse()*(*it_pt)->match_epip_point;
            Eigen::Vector2d new_pix = project3Dto2D ( eps_pie_point( 0,0 ), eps_pie_point( 1,0 ), eps_pie_point( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
            Eigen::Vector2d orgipx_pix = project3Dto2D ( orgipx_point( 0,0 ), orgipx_point( 1,0 ), orgipx_point( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
            Eigen::Vector3d p0=  (*it_pt)->pos_world/(*it_pt)->pos_world[2]*(*it_pt)->ser_mu; //(*it_pt)->orig_pos*(*it_pt)->ser_mu;
            Eigen::Vector2d pixel_efes = project3Dto2D ( p0 ( 0,0 ),p0 ( 1,0 ), p0 ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy );
            Eigen::Vector3d epipa_po = w2c.inverse().translation();
            Eigen::Vector2d epipa_pix = project3Dto2D ( epipa_po ( 0,0 ),epipa_po ( 1,0 ), epipa_po ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy );
            int inter_oxp_x = (int)pixel_efes(0,0);
            int inter_oxp_y = (int)pixel_efes(1,0);
            if((int)orgipx_pix(0,0)<5||(int)orgipx_pix(0,0)>currframe.gray.cols-5||(int)orgipx_pix(1,0)<5
                    ||(int)orgipx_pix(1,0)>currframe.gray.rows-5)
            {
                (*it_pt)->tack_point=false;
                    continue;
            }


            Eigen::Vector2d  epipa_direct = (epipa_pix-orgipx_pix);
            Eigen::Vector2d  epipa_direct_t = (epipa_pix-pixel_efes);
            epipa_direct.normalize();
            epipa_direct_t.normalize();

            double grad_x = (double)(lastframe.gray.ptr<uchar> (inter_oxp_y)[inter_oxp_x+1])-(double)(lastframe.gray.ptr<uchar> (inter_oxp_y)[inter_oxp_x-1]) ;
            double grad_y = (double)(lastframe.gray.ptr<uchar> (inter_oxp_y+1)[inter_oxp_x])-(double)(lastframe.gray.ptr<uchar> (inter_oxp_y-1)[inter_oxp_x]) ;
            if(grad_y==0.0&&grad_x==0)
            {
                (*it_pt)->tack_point=false;
                continue;
            }
            Eigen::Vector2d  grad_dirt = Eigen::Vector2d(grad_x,grad_y);
            grad_dirt.normalize();

             double reale = abs(epipa_direct_t.dot(grad_dirt));
            Eigen::Vector2d epsil_dir;
             Eigen::Vector2d  epi_lenge;
            double depth_sech =0.0;

                epsil_dir = pixel_efes-orgipx_pix;
                epi_lenge = pixel_efes-orgipx_pix;
                depth_sech =(*it_pt)->ser_mu-orgipx_point[2];



            double epsil =0.0;
            if(epsil_dir.norm()!=0 )
            {
                epsil = epsil_dir.norm();
            }
            epsil *= epsil;
            if(reale<0.02)
            {
                reale=0.02;
            }
            reale*=reale;
            double epie_error =0;
			
            if((*it_pt)->epi_lenge!=0)//如果匹配點不為投影點
            {

                    double depth_dov = 0.0;
                     Eigen::Vector2d  l0_go = (new_pix-pixel_efes);

                    if(l0_go.norm()==0.0 )
                    {
                      depth_dov =(*it_pt)->ser_mu -eps_pie_point[2] ;
                       depth_dov*=depth_dov;
                    }
                    else
                    {
                        double org_orgt =0.0;
                        double epi_lsa=l0_go.norm();
                        org_orgt = l0_go.norm();
                        org_orgt *= org_orgt;
                        l0_go.normalize();
                         double lo_tog = abs(epipa_direct_t.dot(grad_dirt));
                         if(lo_tog<0.02)
                         {
                             lo_tog=0.02;
                         }
                         lo_tog*=lo_tog;
                         double depth_diff = (*it_pt)->ser_mu-eps_pie_point[2];
                         depth_dov =pow(abs(depth_diff)/(epi_lsa) ,2)*((org_orgt/lo_tog));
                    }
                    if(epi_lenge.norm()!=0)
                    {
                        epie_error =pow(abs(depth_sech)/(epi_lenge.norm()) ,2)*((epsil/reale)+ 2.0*pow((*it_pt)->grad_def/(*it_pt)->grad_leng,2.0));
                    }
                    else
                    {
                        epie_error = orgipx_point[2]-(*it_pt)->ser_mu;
                        epie_error*=epie_error;
                    }
                    double epi_depth =orgipx_point[2];
                    double sensor_depth = eps_pie_point[2];
                    bool mat_fea =compete_depth((*it_pt)->ser_mu,sqrt((*it_pt)->ser_sigma),epi_depth,sqrt(epie_error));
                    bool fea_org =compete_depth((*it_pt)->ser_mu,sqrt((*it_pt)->ser_sigma),sensor_depth,sqrt(depth_dov));
                    if(mat_fea&&fea_org&&epie_error!=0.0&&depth_dov!=0.0)
                    {
                        double ser_fuse = ((sensor_depth/depth_dov)+(epi_depth/epie_error)+((*it_pt)->ser_mu/(*it_pt)->ser_sigma));
                         double sigma_ser2 =1.0/((1/depth_dov)+(1.0/ (*it_pt)->ser_sigma)+(1.0/epie_error));
                        ser_fuse *=sigma_ser2;
                        (*it_pt)->ser_sigma = sigma_ser2;
                        (*it_pt)->ser_mu = ser_fuse;
                    }
                    else if(mat_fea&&!fea_org&&epie_error!=0.0)
                    {
                        double ser_fuse = ((epi_depth/epie_error)+((*it_pt)->ser_mu/(*it_pt)->ser_sigma));
                         double sigma_ser2 =1.0/((1.0/ (*it_pt)->ser_sigma)+(1.0/epie_error));
                        ser_fuse *=sigma_ser2;
                        (*it_pt)->ser_sigma = sigma_ser2;
                        (*it_pt)->ser_mu = ser_fuse;
                    }
                    else
                    {
                        if(mat_fea&&fea_org&&(epie_error==0.0||depth_dov==0.0))
                        {
                            (*it_pt)->tack_point  = true;
                            continue;
                        }
                        else if (mat_fea&&!fea_org&&epie_error==0.0)
                        {
                            (*it_pt)->tack_point  = true;
                            continue;
                        }
                        else
                        {
                            (*it_pt)->tack_point  = false;
                            continue;
                        }
                    }


            }
            else if((*it_pt)->epi_lenge==0)//如果匹配點為投影點
            {
               Eigen::Vector2d  l0_go = (new_pix-pixel_efes);
                 if(l0_go.norm()!=0.0)
                 {
                     double org_orgt =0.0;
                     double epi_lsa=l0_go.norm();
                     org_orgt = l0_go.norm();

                     org_orgt *= org_orgt;
                     l0_go.normalize();
                      double lo_tog = abs(epipa_direct_t.dot(grad_dirt));
                      if(lo_tog<0.02)
                      {
                          lo_tog=0.02;
                      }
                      lo_tog*=lo_tog;
                      double depth_diff = (*it_pt)->ser_mu-eps_pie_point[2];
                      epie_error =pow(abs(depth_diff)/(epi_lsa) ,2)*((org_orgt/lo_tog)+ 2.0*pow((*it_pt)->grad_def/(*it_pt)->grad_leng,2.0));
                      double epi_depth =eps_pie_point[2];
                      bool mat_fea =compete_depth((*it_pt)->ser_mu,sqrt((*it_pt)->ser_sigma),epi_depth,sqrt(epie_error));
                      if(mat_fea&&epie_error!=0.0)
                      {
                          double ser_fuse = ((epi_depth/epie_error)+((*it_pt)->ser_mu/(*it_pt)->ser_sigma));
                           double sigma_ser2 =1.0/((1.0/ (*it_pt)->ser_sigma)+(1.0/epie_error));
                          ser_fuse *=sigma_ser2;
                          (*it_pt)->ser_sigma = sigma_ser2;
                          (*it_pt)->ser_mu = ser_fuse;
                      }
                      else
                      {
                          if(mat_fea&&epie_error==0.0)
                          {
                              (*it_pt)->tack_point  = true;
                              continue;
                          }
                          else
                          {
                              (*it_pt)->tack_point  = false;
                              continue;
                          }

                      }

                 }
                 else
                 {
                     double depth_dov = (orgipx_point[2])- (*it_pt)->ser_mu ;
                     depth_dov*=depth_dov;
                     bool mat_fea =compete_depth((*it_pt)->ser_mu,sqrt((*it_pt)->ser_sigma),orgipx_point[2],sqrt(depth_dov));
                     if(mat_fea&&depth_dov!=0.0)
                     {
                         epie_error =depth_dov;
                         double d_cove2 = epie_error;
                         double ser_fuse = (d_cove2*(*it_pt)->ser_mu+(*it_pt)->ser_sigma*(orgipx_point[2])) / ( (*it_pt)->ser_sigma+d_cove2);
                         double sigma_ser2 = ( (*it_pt)->ser_sigma * d_cove2 ) / ( (*it_pt)->ser_sigma + d_cove2 );
                         (*it_pt)->ser_sigma = sigma_ser2;
                         (*it_pt)->ser_mu = ser_fuse;
                     }
                     else
                     {
                         if(mat_fea&&depth_dov==0.0)
                         {
                             (*it_pt)->tack_point  = true;
                             continue;
                         }
                         else
                         {
                             (*it_pt)->tack_point  = false;
                             continue;
                         }
                     }
                 }

            }
            else
            {
                (*it_pt)->tack_point  = false;
                continue;
            }


             float b = (*it_pt)->colorp.b;
             float g = (*it_pt)->colorp.g;
             float r = (*it_pt)->colorp.r;

             (*it_pt)->pos_world = (*it_pt)->pos_world/(*it_pt)->pos_world[2]*(*it_pt)->ser_mu;


        }
    }


    //   cv::imshow ( "result_repro", img_show );



}


Eigen::Vector2d featrue_al(cv::Mat& orggray,FRAME& currgray,Eigen::Vector2d repro_pix,Eigen::Vector2d pxss_max_curr,Eigen::Vector2d pxss_min_curr
                           ,double depth,double sigma,double scale,double& lenger,double& error_temp)
{
    Eigen::Vector2d eip_pixlast =pxss_min_curr-pxss_max_curr;
    double scale_stp = 0.4;
    double eqline_fow =floor(eip_pixlast.norm());
    double step=0.0;
    eip_pixlast.normalize();
    Eigen::Vector2d epi_scal =eip_pixlast*scale_stp;
    cv::Mat initROI;
    currgray.gray(cv::Rect((int)round(repro_pix[0]-2),(int)round(repro_pix[1]-2),5,5)).convertTo(initROI,CV_16SC1);
    Eigen::Vector2d resualt=repro_pix;
    cv::Mat COMP_TMP =orggray-initROI;
     COMP_TMP=COMP_TMP.mul(COMP_TMP);
    double init_reds =cv::mean(COMP_TMP)[0];
    error_temp = sqrt(init_reds);
    Eigen::Vector2d step_epi =Eigen::Vector2d(0,0);
    double time_step =(repro_pix-pxss_max_curr).norm();
    while(step<eqline_fow) //在epipolar line上用ssd找相同的特徵點
    {

        double x_boud =pxss_max_curr[0]+step_epi[0];
        double y_boud =pxss_max_curr[1]+step_epi[1];
        if(x_boud<10||x_boud>currgray.gray.cols-10||y_boud<10||y_boud>currgray.gray.rows-10)
        {
            step_epi +=epi_scal;
            step +=epi_scal.norm();
            continue;
        }
        double d = (double)currgray.depth.ptr<ushort>((int)y_boud)[(int)x_boud]/scale;
        if(d<depth-sigma||d>depth+sigma)
        {
            step_epi +=epi_scal;
            step +=epi_scal.norm();
            continue;
        }
        cv::Mat curROI;
        currgray.gray(cv::Rect((int)round(x_boud-2),(int)round(y_boud-2),5,5)).convertTo(curROI,CV_16SC1);
        cv::Mat tmp =orggray-curROI;
        tmp =tmp.mul(tmp);
        double red =cv::mean(tmp)[0];
        if(red<init_reds)
        {
           init_reds= red;
           error_temp = sqrt(red);
           resualt=pxss_max_curr+step_epi;
           time_step = step_epi.norm();
        }
        step_epi +=epi_scal;
        step +=epi_scal.norm();
    }

        lenger =(repro_pix-pxss_max_curr).norm()-time_step; //this is problem ,because  any one not touch inite value, we need fix zero!!

    return resualt;
}
bool project_point(vector<KEYFRAME> Map, FRAME &lastframe, FRAME &currframe, CAMERA_INTRINSIC_PARAMETERS camera,Eigen::Isometry3d& resual)
{

	/////////////////////匹配點因該要在邊緣上,因此需做canny測試////////////////////
    cv::Mat cannyer;
    cv::Canny(currframe.gray,cannyer,(int)(currframe.avg_itenity/1.5),(int)currframe.avg_itenity,3);

	////////設定兩個相機的頂點/////////////////////
    g2o::SparseOptimizer globalOptimizer;
    setuprepro(&globalOptimizer,camera);
    int v_id=0;
    g2o::VertexSE3Expmap* v_frame1 = new g2o::VertexSE3Expmap();
     v_frame1->setId(v_id++);
     v_frame1->setFixed(true);
    v_frame1->setEstimate(  g2o::SE3Quat()   );
     globalOptimizer.addVertex( v_frame1 );

     g2o::VertexSE3Expmap* v_frame2 = new g2o::VertexSE3Expmap();
      v_frame2->setId(v_id++);
      v_frame2->setFixed(false);
     v_frame2->setEstimate(  g2o::SE3Quat()   );
      globalOptimizer.addVertex( v_frame2 );
		///////////////////////////////////////////////////////////////
       int match_p =0;
      for(vector<Measurement*> ::iterator it_pt = lastframe.kp.begin();it_pt!=lastframe.kp.end();++it_pt)
      {
          Eigen::Isometry3d w2c =  lastframe.l_t_r_w;
          Eigen::Vector3d ex_ref = (*it_pt)->pos_world/(*it_pt)->pos_world[2] ;
          Eigen::Vector3d extra =  w2c*(ex_ref*(*it_pt)->ser_mu);
          Eigen::Vector3d px_min_curr =  w2c*(ex_ref*((*it_pt)->ser_mu - 2.0*sqrt((*it_pt)->ser_sigma))) ;	// 根據深度誤差分布分別投影深度-sigma
          Eigen::Vector3d px_max_curr = w2c*(ex_ref*((*it_pt)->ser_mu + 2.0*sqrt((*it_pt)->ser_sigma))) ;// 根據深度誤差分布分別投影深度+sigma
          Eigen::Vector2d pxss_min_curr = project3Dto2D ( px_min_curr( 0,0 ), px_min_curr( 1,0 ), px_min_curr( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
          Eigen::Vector2d pxss_max_curr = project3Dto2D ( px_max_curr( 0,0 ), px_max_curr( 1,0 ), px_max_curr( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
          Eigen::Vector2d repro_pix = project3Dto2D ( extra( 0,0 ), extra ( 1,0 ), extra ( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );

          Eigen::Vector3d ref_point =(*it_pt)->pos_world;
          Eigen::Vector2d ref_pix = project3Dto2D ( ref_point( 0,0 ), ref_point( 1,0 ), ref_point( 2,0 ), camera.fx, camera.fy, camera.cx, camera.cy  );
          Eigen::Vector2d eip_pixlast = pxss_min_curr-pxss_max_curr;
            (*it_pt)->tack_point =false;

          if((int)repro_pix(0,0)<5||(int)repro_pix(0,0)>currframe.gray.cols-5||(int)repro_pix(1,0)<5
                  ||(int)repro_pix(1,0)>currframe.gray.rows-5)
          {
                  continue;
          }

          eip_pixlast.normalize();
          double inten_ratio = currframe.avg_itenity/Map[(*it_pt)->frameID].avg_itenity;
          double direct_tri =0;
          double errorer =0;
          cv::Mat orgROI;
          lastframe.gray(cv::Rect((int)round(ref_pix[0]-2),(int)round(ref_pix[1]-2),5,5)).convertTo(orgROI,CV_16SC1);
          orgROI *=inten_ratio;
		  //找匹配點//
          Eigen::Vector2d off_setmatch = featrue_al(orgROI,currframe,repro_pix,pxss_max_curr,pxss_min_curr,(*it_pt)->ser_mu , 2.0*sqrt((*it_pt)->ser_sigma),camera.scale,direct_tri,errorer);
          (*it_pt)->epi_lenge = (off_setmatch-repro_pix).norm();
		  ///////////
		  ///////////////////////確定匹配成功的條件///////////////
          if(errorer*(1.0+inten_ratio)<abs((*it_pt)->grayscale-(double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0))[(int)off_setmatch(0,0)])))
          {
              continue;
          }

          ushort ddd = currframe.depth.ptr<ushort> ((int)repro_pix(1,0))[(int)repro_pix(0,0)];
           if(ddd ==0)
           {
                 continue;
           }
           (*it_pt)->match_epip_point= project2Dto3D ( repro_pix(0,0), repro_pix(1,0), ddd, camera.fx, camera.fy, camera.cx, camera.cy, camera.scale );

          double grad_x = (double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0))[(int)off_setmatch(0,0)+1])-(double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0))[(int)off_setmatch(0,0)-1]) ;
          double grad_y = (double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0)+1)[(int)off_setmatch(0,0)])-(double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0)-1)[(int)off_setmatch(0,0)]) ;

          if(grad_y==0.0&&grad_x==0)
           {
              continue;
           }
           Eigen::Vector2d  grad_dirt = Eigen::Vector2d(grad_x,grad_y);
           (*it_pt)->grad_leng = grad_dirt.norm();


           ushort d = currframe.depth.ptr<ushort> ((int)off_setmatch(1,0))[(int)off_setmatch(0,0)];
            if(d ==0)
            {
                continue;
            }
            Eigen::Vector3d curr_depth = project2Dto3D ( off_setmatch(0,0), off_setmatch(1,0), d, camera.fx, camera.fy, camera.cx, camera.cy, camera.scale );

            (*it_pt)->grad_def = (*it_pt)->grayscale-(double)(currframe.gray.ptr<uchar> ((int)off_setmatch(1,0))[(int)off_setmatch(0,0)]);
            (*it_pt)->match_point = curr_depth;

            if(cannyer.ptr<uchar>((int)off_setmatch(1,0))[(int)off_setmatch(0,0)]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)+1)[(int)off_setmatch(0,0)]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)+1)[(int)off_setmatch(0,0)-1]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)+1)[(int)off_setmatch(0,0)+1]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)-1)[(int)off_setmatch(0,0)]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)-1)[(int)off_setmatch(0,0)-1]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0)-1)[(int)off_setmatch(0,0)+1]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0))[(int)off_setmatch(0,0)+1]==0&&
              cannyer.ptr<uchar>((int)off_setmatch(1,0))[(int)off_setmatch(0,0)-1]==0
            )
           {
              continue;
           }
			///////////////////////////////////////////////////////////////////////////
            (*it_pt)->tack_point =true;
           match_p++;




      }

	  //如果投影成功的點數量有超過一半 ,才會做深度更新與姿態調整
      if(match_p>lastframe.kp.size()*0.5)
      {
		  
        depth_fit(Map,lastframe,currframe,lastframe.l_t_r_w,camera);//深度更新
		
        for(vector<Measurement*> ::iterator it_pt = lastframe.kp.begin();it_pt!=lastframe.kp.end();++it_pt) //將更新後的特徵點在做三維校正
        {
            if((*it_pt)->tack_point)
            {
               Eigen::Vector2d ref_pie = project3Dto2D ( (*it_pt)->pos_world(0,0), (*it_pt)->pos_world(1,0),(*it_pt)->pos_world(2,0), camera.fx, camera.fy, camera.cx, camera.cy );
                Eigen::Vector2d curr_depth = project3Dto2D ( (*it_pt)->match_point(0,0), (*it_pt)->match_point(1,0), (*it_pt)->match_point(2,0), camera.fx, camera.fy, camera.cx, camera.cy );

                g2o::VertexSBAPointXYZ* v = new g2o::VertexSBAPointXYZ();
                  v->setId( v_id++ );
                  v->setFixed(false); //sparse to false , dense to true!!!!!
                  v->setMarginalized(true);
                  v->setEstimate( (*it_pt)->pos_world );
                  globalOptimizer.addVertex( v );
                EdgeProjectXYZ2UV* edge1 = new EdgeProjectXYZ2UV();
                edge1->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v));
                edge1->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v_frame1));
                edge1->setMeasurement(ref_pie);
                edge1->information() = Eigen::Matrix2d::Identity(2,2);
                edge1->setParameterId(0, 0);
                globalOptimizer.addEdge(edge1);

                EdgeProjectXYZ2UV* edge2 = new EdgeProjectXYZ2UV();
                edge2->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v));
                edge2->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v_frame2));
                edge2->setMeasurement(curr_depth);
                edge2->information() = Eigen::Matrix2d::Identity(2,2);
                edge2->setParameterId(0, 0);
                globalOptimizer.addEdge(edge2);
                //cv::circle ( img_show, cv::Point2d ( curr_depth ( 0,0 ), curr_depth ( 1,0 ) ), 1, cv::Scalar (0.0,0.0,255.0 ), -1 );
            }
        }


       globalOptimizer.setVerbose(false);
        globalOptimizer.initializeOptimization();
        globalOptimizer.optimize(100);



        resual=v_frame2->estimate();
        return true;
      }
      else
      {
          resual= lastframe.l_t_r_w;
          return false;
      }
}



