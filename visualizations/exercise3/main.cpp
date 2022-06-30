#include <Eigen/Dense>

#include "../../shared/colors.cpp"
#include "../../shared/de_casteljau.cpp"
#include "../../shared/kd_tree.cpp"
#include "../../shared/off_file_parser.cpp"
#include "../../shared/typealiases.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "portable-file-dialogs.h"

class SpatialData
{
public:
    std::unique_ptr<KdTree> kd_tree;
    std::vector<Eigen::Vector3f> points;
    std::vector<Eigen::Vector3f> normals;
    std::vector<float> minima;
    std::vector<float> maxima;
    SpatialData(std::vector<Eigen::Vector3f> const &pts, std::vector<Eigen::Vector3f> const &norms, std::vector<float> mins, std::vector<float> maxs)
    {
        kd_tree = std::make_unique<KdTree>(pts);
        points = pts;
        normals = norms;
        minima = std::move(mins);
        maxima = std::move(maxs);
    }
    std::vector<Eigen::Vector3f> control_mesh_nodes;
    std::vector<std::array<int, 2>> control_mesh_edges;
};

std::unique_ptr<SpatialData> spatial_data;

bool pointIsNearestToOtherPoint(const Eigen::Vector3f &point1, const Eigen::Vector3f &point2)
{
    Eigen::Vector3f nearestPoint = spatial_data->kd_tree->collectKNearest(point1, 1)[0];
    return nearestPoint == point2;
}

Eigen::Vector3f calculateOffsetPoint(
    const Eigen::Vector3f &point,
    const Eigen::Vector3f &normal,
    float startingAlpha,
    bool isPositiveDirection)
{
    float alpha = startingAlpha;
    float factor = isPositiveDirection ? 1 : -1;

    Eigen::Vector3f offsetPoint = point + (factor * alpha * normal);
    while (!pointIsNearestToOtherPoint(offsetPoint, point))
    {
        alpha = alpha / 2;
        offsetPoint = point + (factor * alpha * normal);
    }
    return offsetPoint;
}

void createOffsetPoints()
{
    std::vector<Eigen::Vector3f> positive;
    std::vector<Eigen::Vector3f> negative;
    std::vector<float> p1 = spatial_data->minima;
    std::vector<float> p2 = spatial_data->maxima;
    float dx = p1[0] - p2[0];
    float dy = p1[1] - p2[1];
    float dz = p1[2] - p2[2];
    float starting_alpha = 0.01 * std::sqrt(dx * dx + dy * dy + dz * dz) * 1.2;
    for (size_t i = 0; i < spatial_data->points.size(); i++)
    {
        Eigen::Vector3f point = spatial_data->points[i];
        Eigen::Vector3f normal = spatial_data->normals[i].normalized();
        Eigen::Vector3f positive_offset_point = calculateOffsetPoint(
            point,
            normal,
            starting_alpha,
            true);
        Eigen::Vector3f negative_offset_point = calculateOffsetPoint(
            point,
            normal,
            starting_alpha,
            false);
        positive.push_back(positive_offset_point);
        negative.push_back(negative_offset_point);
    }

    polyscope::registerPointCloud("Offset points - positive", positive)
        ->setPointRadius(0.0025)
        ->setPointColor(kGreen)
        ->setEnabled(false);

    polyscope::registerPointCloud("Offset points - negative", negative)
        ->setPointRadius(0.0025)
        ->setPointColor(kRed)
        ->setEnabled(false);
}

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
                createOffsetPoints();
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