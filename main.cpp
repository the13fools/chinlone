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
    sprintf(folderpath, "ex/%s", foldername.c_str());
    mkdir("ex", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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

void computeCentroids(const Eigen::MatrixXi &_F,const Eigen::MatrixXd &_V, Eigen::MatrixXd &centroids)
{
    int nfaces = _F.rows();
    int nverts = _V.rows();

    centroids.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d pos(0,0,0);
        for (int j = 0; j < 3; j++)
        {
            pos += _V.row(_F(i,j));
        }
        centroids.row(i) = pos/3;
    }
}

Eigen::Vector3d faceNormal(const Eigen::MatrixXi &_F, const Eigen::MatrixXd &_V, int faceidx)
{
    Eigen::Vector3d p0 = _V.row(_F(faceidx, 0));
    Eigen::Vector3d p1 = _V.row(_F(faceidx, 1));
    Eigen::Vector3d p2 = _V.row(_F(faceidx, 2));
    Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
    n /= n.norm();
    return n;
}


void propogateField(const Eigen::MatrixXi &F_div, const Eigen::MatrixXd &V_div, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, Eigen::MatrixXd &field)
{
  bool done = false;
  while (!done) 
  {
      for (int i = 0; i < F_div.rows(); i++)
      {
	  for (int e = 0; e < 3; e++) 
	  {
	      int edgeIdx  = F_edges(i, e);
	      int neighbor = E( edgeIdx, 2 );
	      if (neighbor == i) { neighbor = E( edgeIdx, 3 ); }
	      if (field.row(neighbor).norm() < .01)
	      {
		  Eigen::Vector3d n1 = faceNormal(F_div, V_div, i);
		  Eigen::Vector3d n2 = faceNormal(F_div, V_div, neighbor);
		  Eigen::Vector3d commone = V_div.row( E(edgeIdx, 0) ) - V_div.row( E(edgeIdx, 1) );
		  commone.normalize();

		  Eigen::Vector3d t1 = n1.cross(commone);
		  Eigen::Vector3d t2 = n2.cross(commone);	  
		  
		  double alpha = commone.dot( field.row(i) );
		  double beta  = t1.dot( field.row(i) );

		  field.row(neighbor) = alpha * commone + beta * t2;
	      }   	  
	  }	  
      }
      done = true;
      for (int i = 0; i < F_div.rows(); i++)
      {
          if (field.row(i).norm() < .01)
	  {
              done = false;
	  }
      }
  } 
}

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
/*  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
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
*/
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ("../tet.obj", V, F);

  Eigen::MatrixXd V_div;
  Eigen::MatrixXi F_div;

  int divFactor = 5;

  igl::upsample(V,F,V_div, F_div, divFactor);

  
  Eigen::MatrixXi E;
  Eigen::MatrixXi F_edges;
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges); 


  
  // Plot the mesh
  igl::viewer::Viewer *viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(V_div, F_div);
  viewer->data.set_face_based(true);

  const Eigen::RowVector3d red(0.9,.1,.1),green(0.1,0.9,0.2),blue(0.1,0.2,0.8),black(0,0,0);
  Eigen::RowVector3d colors[] = {red, green, blue,black};

  Eigen::MatrixXd F_centroids;
  computeCentroids(F_div, V_div, F_centroids); 

  for (int ax = 0; ax < 3; ax ++) 
  {

      Eigen::MatrixXd field = Eigen::MatrixXd::Zero(F.rows(), 3);

      Eigen::Matrix3d t(Eigen::AngleAxisd(M_PI *(1. / 12. +  2./3. *ax), faceNormal(F,V,0) ));
      Eigen::Vector3d r = Eigen::Vector3d(0,-1, 0).dot(faceNormal(F,V,0)) * faceNormal(F,V,0);;
      r = (Eigen::Vector3d(0, -1, 0) - r);
      r.normalize();
      field.row(0) = t * r;

      propogateField(F, V, E, F_edges, field);
      Eigen::MatrixXd field_div =  Eigen::MatrixXd::Zero(F_div.rows(), 3);

      for (int i = 0; i < field.rows(); i++) 
      {
	  for (int j = 0; j < pow(4, divFactor); j++)
	  {
	      field_div.row( pow(4, divFactor) * i + j) = field.row(i);
	  } 
      }
      logToFile(field_div, "tet", std::to_string(ax));
      igl::writeOBJ("tet.obj",V_div,F_div);

      viewer->data.add_edges( F_centroids  + field_div * .1 / 5., F_centroids, colors[ax]);
  }


  viewer->callback_init = [&](igl::viewer::Viewer& viewer) 
  { 
      viewer.ngui->window()->setVisible(false); 
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
