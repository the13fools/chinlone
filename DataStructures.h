#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H 

#include <Eigen/Core>

void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E);
void buildEdgesPerFace(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, Eigen::MatrixXi &F_edges);

#endif
