
#include "pose_graph.h"

bool graph_container(FRAME& lastFrame,FRAME* currFrame,vector<FRAME>& frame_map,vector<DISPLAY>& frame_group,CAMERA_INTRINSIC_PARAMETERS& camera,Eigen::Isometry3d totre)
{
    bool result=false;
    int loop_index=-1;
    Eigen::Isometry3d containere =Eigen::Isometry3d::Identity();
    Eigen::Isometry3d containe_new =Eigen::Isometry3d::Identity();
    if(detect_frame(&lastFrame,frame_map, &loop_index))
    {
        if(frame_group[lastFrame.frameID].conect.find(loop_index) == frame_group[lastFrame.frameID].conect.end())
        {
        double error =  poseEstimationDirect (&frame_map[loop_index],currFrame, camera,containere);
        bool can_contain = (error<lastFrame.chi2_e&&containere.translation().norm()<0.1)?true:false;
        if(can_contain)
        {
            PointCloud::Ptr line_net = connect_net(&lastFrame,frame_map[loop_index],totre);
            frame_group[lastFrame.frameID].conect.insert(pair<int,PointCloud::Ptr>(loop_index,line_net));
            frame_group[lastFrame.frameID].conecteigen.insert(pair<int,Eigen::Isometry3d>(loop_index,lastFrame.l_t_r_w*containere.inverse()));
            if(frame_group[loop_index].conect.find(lastFrame.frameID) == frame_group[loop_index].conect.end())
            {
                 frame_group[loop_index].conect.insert(pair<int,PointCloud::Ptr>(lastFrame.frameID,line_net));
                  frame_group[loop_index].conecteigen.insert(pair<int,Eigen::Isometry3d>(lastFrame.frameID,containere*lastFrame.l_t_r_w.inverse()));
            }
            if(abs(lastFrame.frameID-loop_index)>=2)
            {
              containe_new=pose_graph(lastFrame,frame_map[loop_index],frame_map,containere);
            }
            else
            {
                containe_new = containere;
            }
            cout<<"number:"<<lastFrame.frameID<<"detect!!"<<loop_index<<endl;
            lastFrame=copyframe_value(frame_map[loop_index]);
           // lastFrame.l_t_r_w = containere;
            lastFrame.l_t_r_w =containe_new;
            result=true;
        }
        }
    }


    return result;
}
bool detect_frame(FRAME* currframe,vector<FRAME>& map, int* index)
{
    double theslod = (currframe->l_t_r_w.translation().norm())*1.5;
    Eigen::Vector3d currpos =(currframe->l_t_r_w*currframe->g_t_r_w).translation();
    bool tmp =false;
    for(std::map<int,bool>::iterator kc =currframe->esmti_filter.begin(); kc!=currframe->esmti_filter.end();kc++)
    {
        if(kc->second&&abs(kc->first-currframe->frameID)>1)
        {
            double distan =  (currpos-map[kc->first].campos).norm();
            if(distan<theslod)
            {
                theslod=distan;
                *index = kc->first;
                tmp=true;
            }
        }

    }
    return tmp;
}
void setuppose_graph(g2o::SparseOptimizer* globalOptimizer)
{
             g2o::BlockSolver_6_3::LinearSolverType* linearSolver = new  g2o::LinearSolverCSparse<g2o::BlockSolver_6_3::PoseMatrixType> ();
             g2o::BlockSolver_6_3* block_solver = new g2o::BlockSolver_6_3( linearSolver );
             g2o::OptimizationAlgorithmLevenberg* algorithm = new g2o::OptimizationAlgorithmLevenberg( block_solver );
             globalOptimizer->setAlgorithm( algorithm );
             globalOptimizer->setVerbose( false );
}

/*
Eigen::Isometry3d  pose_graph(FRAME& lastFrame ,FRAME loopFrame,vector<FRAME>& map, Eigen::Isometry3d containere )
{
    g2o::SparseOptimizer globalOptimizer;
    setuppose_graph(&globalOptimizer);
    g2o::VertexSE3* v = new g2o::VertexSE3();
    v->setId( loopFrame.frameID );
    v->setEstimate(  loopFrame.g_t_r_w ); //估计为单位矩阵
    v->setFixed( true ); //第一个顶点固定，不用优化
    globalOptimizer.addVertex( v );

    for(int famid = loopFrame.frameID+1;famid<=lastFrame.frameID;famid++)
    {

        g2o::VertexSE3 *v = new g2o::VertexSE3();
        v->setId( map[famid].frameID );
        v->setEstimate( map[famid].g_t_r_w  );
        globalOptimizer.addVertex(v);

        g2o::EdgeSE3* edge = new g2o::EdgeSE3();

        edge->vertices() [0] = globalOptimizer.vertex( famid-1);
        edge->vertices() [1] = globalOptimizer.vertex( famid );

        Eigen::Matrix<double, 6, 6> information = Eigen::Matrix< double, 6,6 >::Identity();
        information(0,0) = information(1,1) = information(2,2) = 100;
        information(3,3) = information(4,4) = information(5,5) = 100;
        edge->setInformation( information );
        Eigen::Isometry3d T;
        T =  map[famid].g_t_r_w*map[famid-1].g_t_r_w.inverse();
        edge->setMeasurement( T );
        globalOptimizer.addEdge(edge);
    }
   g2o::EdgeSE3* edgeend = new g2o::EdgeSE3();

    edgeend->vertices() [0] = globalOptimizer.vertex( lastFrame.frameID );
    edgeend->vertices() [1] = globalOptimizer.vertex( loopFrame.frameID  );

    Eigen::Matrix<double, 6, 6> informationend = Eigen::Matrix< double, 6,6 >::Identity();
    informationend(0,0) = informationend(1,1) = informationend(2,2) = 100;
    informationend(3,3) = informationend(4,4) = informationend(5,5) = 100;
    edgeend->setInformation( informationend );
    edgeend->setMeasurement( lastFrame.l_t_r_w*containere.inverse() );
    globalOptimizer.addEdge(edgeend);

    globalOptimizer.initializeOptimization();
    globalOptimizer.optimize( 100 );
    for(int famid = loopFrame.frameID;famid<=lastFrame.frameID;famid++)
    {
        g2o::VertexSE3* vertex = dynamic_cast<g2o::VertexSE3*>(globalOptimizer.vertex( famid ));
       //  Eigen::Isometry3d pose = vertex->estimate()*map[famid].g_t_r_w;
        map[famid].g_t_r_w = vertex->estimate();
    }

}
*/

Eigen::Isometry3d  pose_graph(FRAME& lastFrame ,FRAME loopFrame,vector<FRAME>& map, Eigen::Isometry3d containere )
{
    g2o::SparseOptimizer globalOptimizer;
    setuppose_graph(&globalOptimizer);
    g2o::VertexSE3* v1 = new g2o::VertexSE3();
    v1->setId( loopFrame.frameID );
    v1->setEstimate(  loopFrame.g_t_r_w ); //估计为单位矩阵
    v1->setFixed( true ); //第一个顶点固定，不用优化
    globalOptimizer.addVertex( v1 );

    g2o::VertexSE3* v2 = new g2o::VertexSE3();
    v2->setId( lastFrame.frameID );
    v2->setEstimate(  lastFrame.g_t_r_w ); //估计为单位矩阵
    //v2->setFixed( true ); //第一个顶点固定，不用优化
    globalOptimizer.addVertex( v2 );

    g2o::VertexSE3* v3 = new g2o::VertexSE3();
    v3->setId( lastFrame.frameID+1 );
    v3->setEstimate(  loopFrame.g_t_r_w*containere ); //估计为单位矩阵
    v3->setFixed( true );
    globalOptimizer.addVertex( v3 );



    g2o::EdgeSE3* edge1 = new g2o::EdgeSE3();
    edge1->vertices() [0] = globalOptimizer.vertex( loopFrame.frameID );
    edge1->vertices() [1] = globalOptimizer.vertex( lastFrame.frameID  );
    Eigen::Matrix<double, 6, 6> information1 = Eigen::Matrix< double, 6,6 >::Identity();
    information1(0,0) = information1(1,1) = information1(2,2) = 100;
    information1(3,3) = information1(4,4) = information1(5,5) = 100;
    edge1->setInformation( information1 );
    edge1->setMeasurement(lastFrame.g_t_r_w*loopFrame.g_t_r_w.inverse());
    globalOptimizer.addEdge(edge1);

    g2o::EdgeSE3* edge2 = new g2o::EdgeSE3();
    edge2->vertices() [0] = globalOptimizer.vertex( lastFrame.frameID );
    edge2->vertices() [1] = globalOptimizer.vertex( lastFrame.frameID+1);
    Eigen::Matrix<double, 6, 6> information2 = Eigen::Matrix< double, 6,6 >::Identity();
    information2(0,0) = information2(1,1) = information2(2,2) = 100;
    information2(3,3) = information2(4,4) = information2(5,5) = 100;
    edge2->setInformation( information2 );
    edge2->setMeasurement(lastFrame.l_t_r_w);
    globalOptimizer.addEdge(edge2);

    g2o::EdgeSE3* edge3 = new g2o::EdgeSE3();
    edge3->vertices() [1] = globalOptimizer.vertex( loopFrame.frameID );
    edge3->vertices() [0] = globalOptimizer.vertex( lastFrame.frameID+1);
    Eigen::Matrix<double, 6, 6> information3 = Eigen::Matrix< double, 6,6 >::Identity();
    information3(0,0) = information3(1,1) = information3(2,2) = 100;
    information3(3,3) = information3(4,4) = information3(5,5) = 100;
    edge3->setInformation( information3 );
    edge3->setMeasurement(containere.inverse());
    globalOptimizer.addEdge(edge3);

    globalOptimizer.initializeOptimization();
    globalOptimizer.optimize( 100 );
    map[lastFrame.frameID].g_t_r_w = v2->estimate();
   return v3->estimate()*loopFrame.g_t_r_w.inverse();
}
