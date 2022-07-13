#ifndef SPATIALDATA_H // include guard
#define SPATIALDATA_H

class SpatialData
    {
    public:
        SpatialData(Eigen::MatrixXd meshVer, Eigen::MatrixXi meshFace);
        std::set<int> getAdj(int vertex);
        Eigen::MatrixXd meshV;
        Eigen::MatrixXi meshF;
        int V;
        int F;
        std::vector<std::set<int>> adjVerts;
    };

#endif