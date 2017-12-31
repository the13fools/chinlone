#include <map>
#include <iostream>

#include <Eigen/Geometry>

#include "DataStructures.h"

using namespace std;
using namespace Eigen;

void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E)
{
    map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
        Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>
        edgemap;

    int nfaces = (int)F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int idx1 = F(i, j);
            int idx2 = F(i, (j + 1) % 3);
            int slot = 0;

            if (idx1 > idx2)
            {
                swap(idx1, idx2);
                slot = 1;
            }
            map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
                Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>::iterator it = edgemap.find(pair<int, int>(idx1, idx2));
            if (it == edgemap.end())
            {
                Vector4i newedge;
                newedge[0] = idx1;
                newedge[1] = idx2;
                newedge[2] = newedge[3] = -1;
                newedge[2 + slot] = i;
                edgemap[pair<int, int>(idx1, idx2)] = newedge;
            }
            else
            {
                edgemap[pair<int, int>(idx1, idx2)][2 + slot] = i;
            }
        }
    }

    int nedges = (int)edgemap.size();
    E.resize(nedges, 4);
    int idx = 0;
    for (map<pair<int, int>, Vector4i>::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        E.row(idx) = it->second.transpose();
        idx++;
    }
}

bool consistencyCheckEdges(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges)
{
    int nfaces = (int)F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        for (int e = 0; e < 3; e++)
        {
            //each edge had better point back to the face
            int edgeidx = F_edges(i, e);
            if (E(edgeidx, 2) != i && E(edgeidx, 3) != i)
            {
                std::cerr << "Edge assigned to triangle " << i << " side " << e << " does not have triangle " << i << " as an adjacent face!" << std::endl;
                return false;
            }
            //the edge endpoints need to be opposite vertex e of triangle i
            for (int vtex = 0; vtex < 2; vtex++)
            {
                int facevertidx = (e + 1 + vtex) % 3;
                int vertidx = F(i, facevertidx);
                if (E(edgeidx, 0) != vertidx && E(edgeidx, 1) != vertidx)
                {
                    std::cerr << "Vertex " << vertidx << " is opposite vertex " << e << " on triangle " << i << " but is not part of the edge assigned opposite vertex " << e << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

void buildEdgesPerFace(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, Eigen::MatrixXi &F_edges)
{
    int nfaces = (int)F.rows();
    int nedges = (int)E.rows();

    F_edges.resize(nfaces, 3);
    F_edges = Eigen::MatrixXi::Constant(nfaces,3,-1);
    for (int i = 0; i < nedges; i++)
    {
        int f1 = E(i, 2);
        int f2 = E(i, 3);

        if (f1 > -1)
        {
            int insIdx = -1;
            for (int j = 0; j < 3; j++)
            {
                if (F(f1, j) == E(i, 1))
                {
                    insIdx = (j+1)%3;
                }
            }
            F_edges(f1, insIdx) = i;
        }

        if (f2 > -1)
        {
            int insIdx = -1;
            for (int j = 0; j < 3; j++)
            {
                if (F(f2, j) == E(i, 0))
                {
                    insIdx = (j+1)%3;
                }
            }
            F_edges(f2, insIdx) = i;
        }
    }
    if (!consistencyCheckEdges(F, E, F_edges))
    {
        assert(false);
        exit(-1);
    }
}




