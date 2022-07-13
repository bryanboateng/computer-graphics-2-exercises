#ifndef SPATIALDATA_H // include guard
#define SPATIALDATA_H

class SpatialData
    {
    public:
        SpatialData(Eigen::MatrixXd meshVer, Eigen::MatrixXi meshFace);
        Eigen::MatrixXd meshV;
        Eigen::MatrixXi meshF;
    };

#endif