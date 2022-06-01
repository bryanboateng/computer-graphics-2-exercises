#include <Eigen/Dense>

#include "../../colors.cpp"
#include "../../kd_tree.cpp"
#include "../../kd_tree_2.cpp"
#include "../../off_file_parser.cpp"
#include "../../typealiases.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "portable-file-dialogs.h"

int grid_x_count = 10;
int grid_y_count = 10;
float control_points_radius = 0.1;
bool show_grid = false;
bool show_control_mesh = false;

int is_bezier = 1; // 0 -> mls surface
bool show_surface_mesh = false;
bool show_normals = false;
int subdivision_count;
float mls_radius = 0.1;

class SpatialData
{
public:
    std::unique_ptr<KdTree> kd_tree;
    std::unique_ptr<KdTree2> kd_tree_2;
    std::vector<float> minima;
    std::vector<float> maxima;
    SpatialData(std::vector<Point> points, std::vector<float> mins, std::vector<float> maxs)
    {
        kd_tree = std::make_unique<KdTree>(points);
        kd_tree_2 = std::make_unique<KdTree2>(points);
        minima = mins;
        maxima = maxs;
    }
};

std::unique_ptr<SpatialData> spatial_data;

float wendland(float d)
{
    return std::pow((1 - d), 4) * (4 * d + 1);
}

float weightedLeastSquares(float p_u, float p_v)
{
    auto collected_points = spatial_data->kd_tree_2->collectInRadius(Point{p_u, p_v, 0}, control_points_radius);

    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(6, 6);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(6);
    for (Point collected_point : collected_points)
    {
        float wendland_value = wendland((Eigen::Vector2f(collected_point[0], collected_point[1]) - Eigen::Vector2f(p_u, p_v)).norm());
        float u_i = collected_point[0];
        float v_i = collected_point[1];

        Eigen::VectorXf base(6);
        base << 1, u_i, v_i, u_i * u_i, u_i * v_i, v_i * v_i;

        Eigen::MatrixXf A_i = (base * base.transpose()) * wendland_value;
        A += A_i;
        Eigen::VectorXf b_i = base * collected_point[2] * wendland_value;
        b += b_i;
    }

    Eigen::VectorXf coefficients = A.llt().solve(b);
    return coefficients[0] + coefficients[1] * p_u + coefficients[2] * p_v + coefficients[3] * p_u * p_u + coefficients[4] * p_u * p_v + coefficients[5] * p_v * p_v;
}

void createControlMesh()
{
    PointList mesh_nodes;
    std::vector<std::array<int, 2>> edges;
    float min_x = spatial_data->minima[0];
    float max_x = spatial_data->maxima[0];
    float min_y = spatial_data->minima[1];
    float max_y = spatial_data->maxima[1];

    float x_spacing = (max_x - min_x) / grid_x_count;
    float y_spacing = (max_y - min_y) / grid_y_count;
    for (int i_y = 0; i_y <= grid_y_count; i_y++)
    {
        for (int i_x = 0; i_x <= grid_x_count; i_x++)
        {
            float x_i = x_spacing * i_x;
            float y_i = y_spacing * i_y;
            mesh_nodes.push_back(Point{x_i, y_i, weightedLeastSquares(x_i, y_i)});
        }
    }

    int mesh_nodes_size = mesh_nodes.size();
    for (int i = 0; i < mesh_nodes_size; i++)
    {
        if (i % (grid_x_count + 1) != grid_x_count)
        {
            edges.push_back({i, i + 1});
        }

        int vertical_neighbor_index = i + grid_x_count + 1;
        if (vertical_neighbor_index < mesh_nodes_size)
        {
            edges.push_back({i, vertical_neighbor_index});
        }
    }

    polyscope::registerCurveNetwork("controlMesh", mesh_nodes, edges)
        ->setColor(kRed)
        ->setRadius(0.00064);

    polyscope::registerPointCloud("controlMeshPoints", mesh_nodes)
        ->setPointRadius(0.004)
        ->setPointColor(kRed);
}

void createGrid()
{
    PointList meshNodes;
    std::vector<std::array<int, 2>> edges;

    float min_x = spatial_data->minima[0];
    float max_x = spatial_data->maxima[0];
    float min_y = spatial_data->minima[1];
    float max_y = spatial_data->maxima[1];

    for (float i = 0; i <= grid_x_count; i++)
    {
        auto frac = 1 - (i / grid_x_count);
        meshNodes.push_back(Point{frac * min_x + (1 - frac) * max_x, min_y, 0});
        meshNodes.push_back(Point{frac * min_x + (1 - frac) * max_x, max_y, 0});
    }

    for (float i = 0; i <= grid_y_count; i++)
    {
        auto frac = 1 - (i / grid_y_count);
        meshNodes.push_back(Point{min_x, frac * min_y + (1 - frac) * max_y, 0});
        meshNodes.push_back(Point{max_x, frac * min_y + (1 - frac) * max_y, 0});
    }

    int edge_count = meshNodes.size() / 2;
    for (int i = 0; i < edge_count; i++)
    {
        edges.push_back({i * 2, i * 2 + 1});
    }

    polyscope::registerCurveNetwork("grid", meshNodes, edges)
        ->setColor(kPink)
        ->setRadius(0.00064);
}

void updateGrid()
{
    if (spatial_data == nullptr)
        return;
    if (polyscope::hasCurveNetwork("grid"))
    {
        polyscope::removeCurveNetwork("grid");
    }
    if (show_grid)
        createGrid();
}

void updateControlMesh()
{
    if (spatial_data == nullptr)
        return;
    if (polyscope::hasCurveNetwork("controlMesh"))
    {
        polyscope::removeCurveNetwork("controlMesh");
    }
    if (polyscope::hasPointCloud("controlMeshPoints"))
    {
        polyscope::removePointCloud("controlMeshPoints");
    }
    if (show_control_mesh)
        createControlMesh();
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
                std::vector<Point> points = std::get<0>(parse_result);
                std::vector<float> minima = std::get<2>(parse_result);
                std::vector<float> maxima = std::get<3>(parse_result);

                // Build spatial data structure
                spatial_data = std::make_unique<SpatialData>(points, minima, maxima);

                // Create the polyscope geometry
                polyscope::registerPointCloud("Points", points)
                    ->setPointRadius(0.0025)
                    ->setPointColor(kOrange);
                updateGrid();
                updateControlMesh();
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
    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x * 0.35f);
    if (ImGui::SliderInt("##grid_x_count", &grid_x_count, 1, 30))
    {
        updateGrid();
        updateControlMesh();
    }
    ImGui::SameLine();
    if (ImGui::SliderInt("##grid_y_count", &grid_y_count, 1, 30))
    {
        updateGrid();
        updateControlMesh();
    }
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Text("Width, Height");
    if (ImGui::SliderFloat("Radius", &control_points_radius, 0.0f, 1.0f, "%.3f"))
        updateControlMesh();
    if (ImGui::Checkbox("Show Grid", &show_grid))
        updateGrid();
    if (ImGui::Checkbox("Show Control Mesh", &show_control_mesh))
        updateControlMesh();
    ImGui::Text("Surface");
    ImGui::Separator();
    ImGui::RadioButton("Beziér Surface", &is_bezier, 1);
    ImGui::SameLine();
    ImGui::RadioButton("MLS Surface", &is_bezier, 0);
    ImGui::Checkbox("Show Surface Mesh", &show_surface_mesh);
    ImGui::Checkbox("Show Normals", &show_normals);
    ImGui::SliderInt("Subdivisions", &subdivision_count, 1, 10);
    ImGui::SliderFloat("Radius (MLS only)", &mls_radius, 0.0f, 1.0f, "%.3f");
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
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    // Initialize polyscope
    polyscope::init();

    // Add a few gui elements
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();

    return 0;
}
