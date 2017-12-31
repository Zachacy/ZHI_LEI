#include "LM_badjust.h"



EdgeSE3ProjectDirect::EdgeSE3ProjectDirect ( Eigen::Vector3d point, float fx, float fy, float cx, float cy,vector<cv::Mat>& grayref,vector<cv::Mat>& graycur,int* Layer  )
    : x_world_ ( point ), fx_ ( fx ), fy_ ( fy ), cx_ ( cx ), cy_ ( cy ), grayref_(grayref),graycur_(graycur),Layer_ (Layer)
{}
 void EdgeSE3ProjectDirect::computeError()
    {
        const VertexSE3Expmap* v  =static_cast<const VertexSE3Expmap*> ( _vertices[0] );

        Eigen::Vector3d x_local = v->estimate().map ( x_world_ );

       float scale = 1.0f/(float)(1<< *Layer_);
        float x = (x_local[0]*fx_/x_local[2] + cx_)*scale;
        float y = (x_local[1]*fy_/x_local[2] + cy_)*scale;
        float xu = (x_world_[0]*fx_/x_world_[2] + cx_)*scale;
        float yv = (x_world_[1]*fy_/x_world_[2] + cy_)*scale;
        // check x,y is in the image
       if ( x-5<0 || ( x+5 ) >graycur_[*Layer_].cols || ( y-5) <0 || ( y+5 ) >graycur_[*Layer_].rows )
        {
            _error(0,0)  = 0.0;
            this->setLevel ( 1 );

        }
        else
        {
           this->setLevel ( 0 );
           _error ( 0,0 ) = (getPixelValue (graycur_[*Layer_] ,x,y ) - getPixelValue(grayref_[*Layer_],xu,yv)); //誤差值計算
          //  _error ( 0,0 ) = getPixelValue (graycur_[*Layer_] ,x,y ) - _measurement;
        }
    }

    // plus in manifold
 void EdgeSE3ProjectDirect::linearizeOplus( )
     {
         if ( level() == 1 )
         {
             _jacobianOplusXi = Eigen::Matrix<double, 1, 6>::Zero();
             return;
         }
         VertexSE3Expmap* vtx = static_cast<VertexSE3Expmap*> ( _vertices[0] );
         Eigen::Vector3d xyz_trans = vtx->estimate().map ( x_world_ );   // q in book


         double x = xyz_trans[0];
         double y = xyz_trans[1];
         double invz = 1.0/xyz_trans[2];
         double invz_2 = invz*invz;
         float scale = 1.0f/(float)(1<<*Layer_);
         float u = (x*fx_*invz + cx_)*scale;
         float v = (y*fy_*invz + cy_)*scale;
         // jacobian from se3 to u,v
         // NOTE that in g2o the Lie algebra is (\omega, \epsilon), where \omega is so(3) and \epsilon the translation
        /* if ( (u-5)<0 || ( u+5 ) >graycur_[*Layer_].cols || ( v-5) <0 || ( v+5 ) >graycur_[*Layer_].rows )
           {
               _error(0,0)  = 0.0;
               this->setLevel ( 1 );
               _jacobianOplusXi = Eigen::Matrix<double, 1, 6>::Zero();
               return;
          }*/
         Eigen::Matrix<double, 2, 6> jacobian_uv_ksai; //2D-3D的矩陣

         jacobian_uv_ksai ( 0,0 ) = - x*y*invz_2 *fx_;
         jacobian_uv_ksai ( 0,1 ) = ( 1+ ( x*x*invz_2 ) ) *fx_;
         jacobian_uv_ksai ( 0,2 ) = - y*invz *fx_;
         jacobian_uv_ksai ( 0,3 ) = invz *fx_;
         jacobian_uv_ksai ( 0,4 ) = 0;
         jacobian_uv_ksai ( 0,5 ) = -x*invz_2 *fx_;

         jacobian_uv_ksai ( 1,0 ) = - ( 1+y*y*invz_2 ) *fy_;
         jacobian_uv_ksai ( 1,1 ) = x*y*invz_2 *fy_;
         jacobian_uv_ksai ( 1,2 ) = x*invz *fy_;
         jacobian_uv_ksai ( 1,3 ) = 0;
         jacobian_uv_ksai ( 1,4 ) = invz *fy_;
         jacobian_uv_ksai ( 1,5 ) = -y*invz_2 *fy_;

         Eigen::Matrix<double, 1, 2> jacobian_pixel_uv;

		 //計算梯度//
            float scale_inv =1.0f/scale;
            float warp_u = ( getPixelValue (graycur_[*Layer_] ,u+1,v )-getPixelValue (graycur_[*Layer_], u-1,v ) ) /2; //x軸方向
            float warp_v =( getPixelValue (graycur_[*Layer_], u,v+1 )-getPixelValue (graycur_[*Layer_], u,v-1 ) ) /2;  //y 軸方向
			///////////////////////////
	
          jacobian_pixel_uv ( 0,0 ) = (warp_u)*scale_inv;
          jacobian_pixel_uv ( 0,1 ) = (warp_v)*scale_inv;


        _jacobianOplusXi = jacobian_pixel_uv*jacobian_uv_ksai*scale;

    }



    // get a gray scale value from reference image (bilinear interpolated)//change pyhrimd secalce!!
inline float EdgeSE3ProjectDirect::getPixelValue ( cv::Mat img, float x, float y )
    {
        uchar* data = & img.data[ int ( y ) * img.step + int ( x ) ];
        float xx = x - floor ( x );
        float yy = y - floor ( y );
        return float (
                   ( 1-xx ) * ( 1-yy ) * data[0] +
                   xx* ( 1-yy ) * data[1] +
                   ( 1-xx ) *yy*data[ img.step ] +
                   xx*yy*data[img.step+1]
               );
    }



