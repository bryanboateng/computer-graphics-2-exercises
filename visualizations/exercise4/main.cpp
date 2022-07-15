#include <Eigen/Dense>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "portable-file-dialogs.h"
#include "SpatialData.h"
#include <queue>
#include <algorithm>
#include <vector>
#include <iterator>
#include <math.h> 

bool laplace_smoothing_enabled = false;
bool laplace_cotan_smoothing_enabled = false;

std::set<int> SpatialData::getAdj(int vertex)
{
    if(!adjVerts[vertex].empty()){
        return adjVerts[vertex]; 
    }
    else{
        std::set<int> adjVertex;
        for(int i = 0; i < F; ++i){
            if(vertex == meshF(i,0)){
                adjVertex.insert(meshF(i,1));
                adjVertex.insert(meshF(i,2));
            }
            if(vertex == meshF(i,1)){
                adjVertex.insert(meshF(i,0));
                adjVertex.insert(meshF(i,2));
            }
            if(vertex == meshF(i,2)){
                adjVertex.insert(meshF(i,0));
                adjVertex.insert(meshF(i,1));
            }
        }
        adjVerts[vertex] = adjVertex;
        return adjVertex;
    }
    
}

std::set<int> SpatialData::getAdjFaces(int vertex)
{
    if(!adjFaces[vertex].empty()){
        return adjFaces[vertex]; 
    }
    else{
        std::set<int> adjFace;
        for(int i = 0; i < F; ++i){
            if(vertex == meshF(i,0) || vertex == meshF(i,1) || vertex == meshF(i,2)){
                adjFace.insert(i);
            }
        }
        adjFaces[vertex] = adjFace;
        return adjFace;
    }
    
}

SpatialData::SpatialData(Eigen::MatrixXd meshVer, Eigen::MatrixXi meshFace)
{
    meshV = meshVer;
    meshF = meshFace;
    V = meshV.rows();
    F = meshF.rows();
    adjVerts.resize(V);
    adjFaces.resize(V);
    barycentric = Eigen::MatrixXd::Zero(V,V);
    cotanLaplacian = Eigen::MatrixXd::Zero(V,V);
}

void SpatialData::calculateBarycentric() //broken
{
    for(int i = 0; i < V; ++i)
    {
        std::set<int> faces = getAdjFaces(i);
        double barycentricArea = 0;
        for(auto face : faces){
            double a = (meshV.row(meshF(face,1)) - meshV.row(meshF(face,0))).norm();
            double b = (meshV.row(meshF(face,2)) - meshV.row(meshF(face,0))).norm();
            double c = (meshV.row(meshF(face,2)) - meshV.row(meshF(face,1))).norm();
            double s = (a+b+c)/2.0;
            barycentricArea += (1.0/3.0) * sqrt(s*(s-a)*(s-b)*(s-c));
        }


        barycentric(i,i) = barycentricArea;

    }
}

float cotan(const Eigen::Vector3d &a, const Eigen::Vector3d &b){
    return a.dot(b) / (a.cross(b)).norm();
}

void SpatialData::computeLaplacian()
{
    for(int i = 0; i < V; ++i)
    {
        std::set<int> neighs = getAdj(i);
        for(auto j: neighs)
        {
            std::set<int> facesi = getAdjFaces(i);
            std::set<int> facesj = getAdjFaces(j);
            std::set<int> jointFaces;
            std::set_intersection(facesi.begin(), facesi.end(), facesj.begin(), facesj.end(), std::inserter(jointFaces, jointFaces.begin()));
            std::vector<int> adjVertices;

            for(int f : jointFaces){
                for(int k = 0; k < 3; ++k){
                    int vs = meshF(f,k);
                    if(vs != i && vs != j){
                        adjVertices.push_back(vs);
                    }
                }
            }
            int v1 = adjVertices[0];
            int v2 = adjVertices[1];
            Eigen::Vector3d e1 = meshV.row(i) - meshV.row(v1);
            Eigen::Vector3d e2 = meshV.row(j) - meshV.row(v1);
            Eigen::Vector3d e3 = meshV.row(i) - meshV.row(v2);
            Eigen::Vector3d e4 = meshV.row(j) - meshV.row(v2);

            float cotan_alpha = cotan(e1, e2);
            float cotan_beta  = cotan(e3, e4);
            cotanLaplacian(i,j) = 0.5*(cotan_alpha + cotan_beta);
            cotanLaplacian(i,i) -= 0.5*(cotan_alpha + cotan_beta);

        }
    }

}


std::unique_ptr<SpatialData> spatial_data;


void GraphLaplace()
{
    Eigen::MatrixXd newMeshV(spatial_data->V,3);
    Eigen::MatrixXd laplaceOperator = Eigen::MatrixXd::Zero(spatial_data->V,spatial_data->V);
    Eigen::MatrixXd barycentric = Eigen::MatrixXd::Zero(spatial_data->V,spatial_data->V);

    if(!laplace_smoothing_enabled){
        polyscope::registerSurfaceMesh("input mesh", spatial_data->meshV, spatial_data->meshF);
        return;
    }


    for(int i = 0; i < spatial_data->V; ++i)
    {
        std::set<int> neighs = spatial_data->getAdj(i);
        barycentric(i,i) = 1.0/(double) neighs.size();
        laplaceOperator(i,i) = -1 * (double)neighs.size();
        for(auto adjecent: neighs)
        {
            laplaceOperator(i, adjecent) = 1;
        }
    }
    newMeshV = spatial_data->meshV + (barycentric * laplaceOperator * spatial_data->meshV);
 

    polyscope::registerSurfaceMesh("input mesh", newMeshV, spatial_data->meshF);
    
}

void GraphLaplaceCotan()
{
    Eigen::MatrixXd newMeshV(spatial_data->V,3);

    if(!laplace_cotan_smoothing_enabled){
        polyscope::registerSurfaceMesh("input mesh", spatial_data->meshV, spatial_data->meshF);
        return;
    }

    newMeshV = spatial_data->meshV + (spatial_data->barycentric * spatial_data->cotanLaplacian * spatial_data->meshV);

    polyscope::registerSurfaceMesh("input mesh", newMeshV, spatial_data->meshF);
    int i = 1;
}

void callback()
{
    if (ImGui::Button("Load Points"))
    {
        auto paths = pfd::open_file("Load Points", "", std::vector<std::string>{"point data (*.obj)", "*.obj"}, pfd::opt::none).result();
        if (!paths.empty())
        {
            std::filesystem::path path(paths[0]);

            // Read the point cloud
            try
            {
                Eigen::MatrixXd meshV;
                Eigen::MatrixXi meshF;

                igl::readOBJ(path.string(), meshV, meshF);
                spatial_data = std::make_unique<SpatialData>(meshV, meshF);

                // Register the mesh with Polyscope
                laplace_smoothing_enabled = false;
                polyscope::registerSurfaceMesh("input mesh", meshV, meshF);
            }
            catch (const std::invalid_argument &e)
            {
                polyscope::error(e.what());
                return;
            }
        }
    }
    if (ImGui::Checkbox("Laplace Smoothing Enabled", &laplace_smoothing_enabled))
        GraphLaplace();
    if (ImGui::Checkbox("Cotan Laplace Smoothing Enabled", &laplace_cotan_smoothing_enabled))
        GraphLaplaceCotan();

    if (ImGui::Button("Calculate Barycentric")){
        spatial_data->calculateBarycentric();
    }
    if (ImGui::Button("Compute Cotan Laplacian")){
        spatial_data->computeLaplacian();
    }

}

int main(int argc, char **argv)
{
    polyscope::init();


    // Show the GUI

    // Add a few gui elements
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();
}