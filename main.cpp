#include <igl/viewer/Viewer.h>
#include <igl/upsample.h>
#include <igl/avg_edge_length.h>

#include <iostream>
#include <fstream>
// For making directories
#include <sys/stat.h>
// #include <direct.h>

#include "DataStructures.h"

void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename)
{
#ifndef WIN32
    char folderpath[50];
    sprintf(folderpath, "log/%s", foldername.c_str());
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char logpath[50];
    sprintf(logpath, "%s/%s.txt", folderpath, filename.c_str());
    std::ofstream myfile (logpath);
    for(int i = 0; i < W.rows(); i++)
    {
        if (myfile.is_open())
        {
            myfile << W.row(i) << "\n";
        }

        else
        {
            std::cout << "Unable to open file";
            break;
        }
    }
    myfile.close();
#endif
}

void computeCentroids(const Eigen::MatrixXi &F,const Eigen::MatrixXd &V, Eigen::MatrixXd &centroids)
{
    int nfaces = F.rows();
    int nverts = V.rows();

    centroids.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d pos(0,0,0);
        for (int j = 0; j < 3; j++)
        {
            pos += V.row(F(i,j));
        }
        centroids.row(i) = pos/3;
    }
}

Eigen::Vector3d faceNormal(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int faceidx)
{
    Eigen::Vector3d p0 = V.row(F(faceidx, 0));
    Eigen::Vector3d p1 = V.row(F(faceidx, 1));
    Eigen::Vector3d p2 = V.row(F(faceidx, 2));
    Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
    n /= n.norm();
    return n;
}


int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  Eigen::MatrixXi E;
  Eigen::MatrixXi F_edges;
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);

  Eigen::MatrixXd V_div;
  Eigen::MatrixXi F_div;
 
  igl::upsample(V,F,V_div, F_div,0);

  Eigen::MatrixXd F_centroids;
  computeCentroids(F_div, V_div, F_centroids); 
 
  Eigen::MatrixXd field = Eigen::MatrixXd::Zero(F_div.rows(), 3);
  field.row(0) = Eigen::Vector3d(0, 1, 0).transpose();
  const Eigen::RowVector3d red(0.9,.1,.1),green(0.1,0.9,0.2),blue(0.1,0.2,0.8);

  // Plot the mesh
  igl::viewer::Viewer *viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(V_div, F_div);
  viewer->data.set_face_based(true);
  viewer->data.add_edges( F_centroids  + field, F_centroids, green);

  viewer->callback_init = [&](igl::viewer::Viewer& viewer) 
  { 
      viewer.ngui->window()->setVisible(false); 
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
