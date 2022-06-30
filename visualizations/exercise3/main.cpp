#include <Eigen/Dense>

#include "../../shared/colors.cpp"
#include "../../shared/de_casteljau.cpp"
#include "../../shared/kd_tree.cpp"
#include "../../shared/kd_tree_pair.cpp"
#include "../../shared/off_file_parser.cpp"
#include "../../shared/typealiases.cpp"
#include "../../shared/wendland.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "portable-file-dialogs.h"

int grid_x_count = 10;
int grid_y_count = 10;
int grid_z_count = 10;
bool show_grid = false;
float radius = 75;


std::unique_ptr<KdTreePair> kd_tree_pair;

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
};

std::unique_ptr<SpatialData> spatial_data;

bool pointIsNearestToOtherPoint(const Eigen::Vector3f &point1, const Eigen::Vector3f &point2)
{
    Eigen::Vector3f nearestPoint = spatial_data->kd_tree->collectKNearest(point1, 1)[0];
    return nearestPoint == point2;
}

std::pair<Eigen::Vector3f, float> calculateOffsetPoint(
    const Eigen::Vector3f &point,
    const Eigen::Vector3f &normal,
    float startingAlpha)
{
    float alpha = startingAlpha;

    Eigen::Vector3f offsetPoint = point + (alpha * normal);
    while (!pointIsNearestToOtherPoint(offsetPoint, point))
    {
        alpha = alpha / 2;
        offsetPoint = point + (alpha * normal);
    }
    return std::make_pair(offsetPoint, alpha);
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
    std::vector<std::pair<Eigen::Vector3f, float> > scored_points;
    for (size_t i = 0; i < spatial_data->points.size(); i++)
    {
        Eigen::Vector3f point = spatial_data->points[i];
        Eigen::Vector3f normal = spatial_data->normals[i].normalized();
        std::pair<Eigen::Vector3f, float> positive_offset_point = calculateOffsetPoint(
            point,
            normal,
            starting_alpha);
        std::pair<Eigen::Vector3f, float> negative_offset_point = calculateOffsetPoint(
            point,
            normal,
            -starting_alpha);
        scored_points.push_back(positive_offset_point);
        scored_points.push_back(negative_offset_point);
        positive.push_back(positive_offset_point.first);
        negative.push_back(negative_offset_point.first);
    }
    for (const Eigen::Vector3f &point : spatial_data->points)
    {
        scored_points.emplace_back(point, 0);
    }

    kd_tree_pair = std::make_unique<KdTreePair>(scored_points);

    polyscope::registerPointCloud("Offset points - positive", positive)
        ->setPointRadius(0.0025)
        ->setPointColor(kGreen)
        ->setEnabled(false);

    polyscope::registerPointCloud("Offset points - negative", negative)
        ->setPointRadius(0.0025)
        ->setPointColor(kRed)
        ->setEnabled(false);
}

void createGrid()
{
    std::vector<Eigen::Vector3f> insideGridPoints;
    std::vector<Eigen::Vector3f> outsideGridPoints;

    float min_x = spatial_data->minima[0];
    float max_x = spatial_data->maxima[0];
    float min_y = spatial_data->minima[1];
    float max_y = spatial_data->maxima[1];
    float min_z = spatial_data->minima[2];
    float max_z = spatial_data->maxima[2];

    Eigen::VectorXf xs = Eigen::VectorXf::LinSpaced(grid_x_count + 1, min_x, max_x);
    Eigen::VectorXf ys = Eigen::VectorXf::LinSpaced(grid_y_count + 1, min_y, max_y);
    Eigen::VectorXf zs = Eigen::VectorXf::LinSpaced(grid_z_count + 1, min_z, max_z);
    for (float x : xs)
    {
        for (float y : ys)
        {
            for (float z : zs)
            {
                Eigen::Vector3f p{x, y, z};

                std::vector<std::pair<Eigen::Vector3f, float> > collected_points = kd_tree_pair->collectInRadius(p, radius);
                if (collected_points.empty())
                {
                    collected_points = kd_tree_pair->collectKNearest(p, 1);
                }

                double A = 0;
                double b = 0;
                Eigen::Vector3d d_p = p.cast<double>();
                for (const std::pair<Eigen::Vector3f, float> &collected_point : collected_points)
                {
                    Eigen::Vector3d d_collected_point = collected_point.first.cast<double>();
                    double wendland_value = wendland((d_collected_point - d_p).norm());

                    auto A_i = wendland_value;
                    A += A_i;
                    auto b_i = collected_point.second * wendland_value;
                    b += b_i;
                }

                if (b / A < 0)
                {
                    insideGridPoints.push_back(p);
                }
                else
                {
                    outsideGridPoints.push_back(p);
                }
            }
        }
    }

    polyscope::registerPointCloud("Grid points inside", insideGridPoints)
        ->setPointRadius(0.00125)
        ->setPointColor(kGreen);

    polyscope::registerPointCloud("Grid points outside", outsideGridPoints)
        ->setPointRadius(0.00125)
        ->setPointColor(kPink);
}

void updateGrid()
{
    if (spatial_data == nullptr)
        return;

    if (show_grid)
    {
        createGrid();
    }
    else
    {
        if (polyscope::hasPointCloud("Grid points inside"))
        {
            polyscope::removePointCloud("Grid points inside");
        }
        if (polyscope::hasPointCloud("Grid points outside"))
        {
            polyscope::removePointCloud("Grid points outside");
        }
    }
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
                updateGrid();
            }
            catch (const std::invalid_argument &e)
            {
                polyscope::error(e.what());
                return;
            }
        }
    }
    ImGui::Text("Control Points");
    ImGui::Separator();
    ImGui::Text("Grid");
    ImGui::SameLine();
    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x * 0.25f);
    if (ImGui::SliderInt("##grid_x_count", &grid_x_count, 1, 30))
    {
        updateGrid();
    }
    ImGui::SameLine();
    if (ImGui::SliderInt("##grid_y_count", &grid_y_count, 1, 30))
    {
        updateGrid();
    }
    ImGui::SameLine();
    if (ImGui::SliderInt("##grid_z_count", &grid_z_count, 1, 30))
    {
        updateGrid();
    }
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Text("#X, #Y, #Z");
    if (ImGui::SliderFloat("Radius", &radius, 1.0f, 200.0f, "%.3f"))
    {
        updateGrid();
    }
    if (ImGui::Checkbox("Show Grid", &show_grid))
        updateGrid();
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