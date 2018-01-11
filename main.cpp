#include <igl/viewer/Viewer.h>
#include <igl/upsample.h>
#include <igl/avg_edge_length.h>

#include <iostream>
#include <fstream>
// For making directories
#include <sys/stat.h>
// #include <direct.h>

using namespace std;

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
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ("../tet.obj", V, F);

  Eigen::MatrixXd inp = readMatrix("../test.txt");
  Eigen::MatrixXd line_starts = inp.block(0, 0, inp.rows() - 1, 3);
  Eigen::MatrixXd line_ends  = inp.block(1, 0, inp.rows() - 1, 3);

  igl::viewer::Viewer *viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(V, F);
  viewer->data.set_face_based(true);
  
  const Eigen::RowVector3d red(0.9,.1,.1);
  viewer->data.add_edges( line_starts, line_ends, red);
  
  viewer->callback_init = [&](igl::viewer::Viewer& viewer) 
  { 
      viewer.ngui->window()->setVisible(false); 
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
