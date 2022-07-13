#include <Eigen/Dense>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "portable-file-dialogs.h"
#include "SpatialData.h"

SpatialData::SpatialData(Eigen::MatrixXd meshVer, Eigen::MatrixXi meshFace)
{
    meshV = meshVer;
    meshF = meshFace;
}

std::unique_ptr<SpatialData> spatial_data;

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
                spatial_data = std::make_unique<SpatialData>(meshV, meshF);

                igl::readOBJ(path.string(), meshV, meshF);

                // Register the mesh with Polyscope
                polyscope::registerSurfaceMesh("input mesh", meshV, meshF);
            }
            catch (const std::invalid_argument &e)
            {
                polyscope::error(e.what());
                return;
            }
        }
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