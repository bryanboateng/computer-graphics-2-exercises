#include <Eigen/Dense>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "portable-file-dialogs.h"
#include "SpatialData.h"
#include <queue>

bool laplace_smoothing_enabled = false;

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
    barycentric.resize(V,V);
}

void SpatialData::calculateBarycentric()
{
    return;
}




std::unique_ptr<SpatialData> spatial_data;


void GraphLaplace()
{
    Eigen::MatrixXd newMeshV(spatial_data->V,3);
    int start = spatial_data->meshF(0,0);
    std::vector<bool> visited;
    visited.resize(spatial_data->V,false);
    std::queue<int> queue;
    visited[start] = true;
    queue.push(start);

    if(!laplace_smoothing_enabled){
        polyscope::registerSurfaceMesh("input mesh", spatial_data->meshV, spatial_data->meshF);
        return;
    }


    while(!queue.empty())
    {
        start = queue.front();
        queue.pop();
        
        newMeshV.row(start) = spatial_data->meshV.row(start);
        Eigen::MatrixXd laplaceOperator {{0,0,0}};
        for(auto adjecent: spatial_data->getAdj(start)){
            laplaceOperator += (spatial_data->meshV.row(adjecent) - spatial_data->meshV.row(start));
            if(!visited[adjecent]){
                visited[adjecent] = true;
                queue.push(adjecent);
            }
        }
        newMeshV.row(start) += (laplaceOperator/(double) spatial_data->getAdj(start).size());
    }

    polyscope::registerSurfaceMesh("input mesh", newMeshV, spatial_data->meshF);
    
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