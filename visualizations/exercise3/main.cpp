#include <Eigen/Dense>

#include "../../shared/colors.cpp"
#include "../../shared/de_casteljau.cpp"
#include "../../shared/kd_tree_2.cpp"
#include "../../shared/off_file_parser.cpp"
#include "../../shared/typealiases.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "portable-file-dialogs.h"

class SpatialData
{
public:
    std::unique_ptr<KdTree2> kd_tree_2;
    std::vector<Eigen::Vector3f> normals;
    std::vector<float> minima;
    std::vector<float> maxima;
    SpatialData(std::vector<Eigen::Vector3f> const &points, std::vector<Eigen::Vector3f> const &norms, std::vector<float> mins, std::vector<float> maxs)
    {
        kd_tree_2 = std::make_unique<KdTree2>(points);
        normals = norms;
        minima = std::move(mins);
        maxima = std::move(maxs);
    }
    std::vector<Eigen::Vector3f> control_mesh_nodes;
    std::vector<std::array<int, 2>> control_mesh_edges;
};

std::unique_ptr<SpatialData> spatial_data;

void callback()
{
    if (ImGui::Button("Load Points"))
    {
        auto paths = pfd::open_file("Load Points", "", std::vector<std::string>{"point data (*.off)", "*.off"}, pfd::opt::none).result();
        if (!paths.empty())
        {
            std::filesystem::path path(paths[0]);

            // Read the point cloud
            try
            {
                auto parse_result = parseOff(path.string());
                std::vector<Eigen::Vector3f> points = std::get<0>(parse_result);
                std::vector<Eigen::Vector3f> normals = std::get<1>(parse_result);
                std::vector<float> minima = std::get<2>(parse_result);
                std::vector<float> maxima = std::get<3>(parse_result);

                // Build spatial data structure
                spatial_data = std::make_unique<SpatialData>(points, normals, minima, maxima);

                // Create the polyscope geometry
                polyscope::registerPointCloud("Points", points)
                        ->setPointRadius(0.0025)
                        ->setPointColor(kOrange)
                        ->addVectorQuantity("normals", normals)
                        ->setVectorColor(kBlue)
                        ->setEnabled(false);
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
    // Configure the argument parser
    args::ArgumentParser parser("Computer Graphics 2 Sample Code.");

    // Parse args
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Help &)
    {
        std::cout << parser;
        return 0;
    }
    catch (const args::ParseError &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Options
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::shadowBlurIters = 6;
    polyscope::view::upDir = polyscope::UpDir::YUp;
    // Initialize polyscope
    polyscope::init();

    // Add a few gui elements
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();

    return 0;
}