#include <igl/viewer/Viewer.h>
#include <igl/upsample.h>
#include <igl/avg_edge_length.h>
#include <igl/unproject_onto_mesh.h>


#include <igl/remove_unreferenced.h>

#include <igl/loop.h>

#include <iostream>
#include <fstream>
// For making directories
#include <sys/stat.h>
// #include <direct.h>

using namespace std;

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


#define MAXBUFSIZE  ((int) 1e6)
Eigen::MatrixXd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

   //  Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
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



int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  Eigen::MatrixXd V_orig;
  Eigen::MatrixXi F_orig;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ("../sphere.obj", V_orig, F_orig);


  Eigen::MatrixXd V_unref;
  Eigen::MatrixXi F_unref;
  Eigen::MatrixXi unref;
  igl::remove_unreferenced(V_orig, F_orig, V_unref, F_unref, unref);
  
  std::cout << unref;

  igl::upsample(V_unref, F_unref, V, F, 0);


  igl::writeOBJ("../sphere_clean.obj",V,F);

  Eigen::MatrixXd vf = Eigen::MatrixXd::Zero(F.rows(), 3);
  
  for (int i = 0; i < F.rows(); i++) 
  {
      Eigen::Vector3d n = faceNormal(F, V, i); 
      Eigen::Vector3d z = Eigen::Vector3d(0, 0, 1); 
      vf.row(i) = n.cross(z);
      vf.row(i) = Eigen::AngleAxisd( M_PI / 2., n) * vf.row(i).transpose();
      vf.row(i).normalize();
  } 

  Eigen::MatrixXd F_centroids;
  computeCentroids(F, V, F_centroids);

  // Initialize white
  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(),3,1);
  igl::viewer::Viewer viewer;
  viewer.callback_mouse_down = 
    [&V,&F,&C](igl::viewer::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
      viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
      // paint hit red
      C.row(fid)<<1,0,0;
      std::cout << fid << " " << F(fid, 0) + 1 << " "<< F(fid, 1) + 1 << " " << F(fid, 2) + 1 << " " << std::endl;
      viewer.data.set_colors(C);
      return true;
    }
    return false;
  };
  std::cout<<R"(Usage:
  [click]  Pick face on shape

)";



  viewer.data.set_mesh(V, F);
  viewer.data.set_face_based(true);
  viewer.data.set_colors(C); 

  const Eigen::RowVector3d red(0.9,.1,.1);
  viewer.data.add_edges( F_centroids, F_centroids + vf * .1, red);
  logToFile(vf, "a", "sphere_field.txt"); 
 /* viewer->callback_init = [&](igl::viewer::Viewer& viewer) 
  { 
      viewer.ngui->window()->setVisible(false); 
      viewer.screen->performLayout();
      return false;
  }; */

  viewer.launch();
}
