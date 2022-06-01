#include "../../kd_tree.cpp"
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
bool show_tangents = false;
int subdivision_count;
float mls_radius = 0.1;

std::unique_ptr<KdTree> kd_tree;

void createGrid(float min_x, float max_x, float min_y, float max_y)
{
    PointList meshNodes;
    std::vector<std::array<int, 2>> edges;

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
        ->setColor({255, 0, 0})
        ->setRadius(0.001);
}

void updateGrid()
{
    if (kd_tree == nullptr)
        return;
    if (polyscope::hasCurveNetwork("grid"))
    {
        polyscope::removeCurveNetwork("grid");
    }
    if (show_grid)
    {
        createGrid(kd_tree->minima[0], kd_tree->maxima[0], kd_tree->minima[1], kd_tree->maxima[1]);
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
                std::vector<Point> points = std::get<0>(parse_result);
                std::vector<float> minima = std::get<2>(parse_result);
                std::vector<float> maxima = std::get<3>(parse_result);

                // Build spatial data structure
                kd_tree = std::make_unique<KdTree>(points, minima, maxima);

                // Create the polyscope geometry
                polyscope::registerPointCloud("Points", points);
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
    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x * 0.35f);
    if (ImGui::SliderInt("##grid_x_count", &grid_x_count, 1, 30))
        updateGrid();
    ImGui::SameLine();
    if (ImGui::SliderInt("##grid_y_count", &grid_y_count, 1, 30))
        updateGrid();
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Text("Width, Height");
    ImGui::SliderFloat("Radius", &control_points_radius, 0.0f, 1.0f, "%.3f");
    if (ImGui::Checkbox("Show Grid", &show_grid))
        updateGrid();
    ImGui::Checkbox("Show Control Mesh", &show_control_mesh);
    ImGui::Text("Surface");
    ImGui::Separator();
    ImGui::RadioButton("Beziér Surface", &is_bezier, 1);
    ImGui::SameLine();
    ImGui::RadioButton("MLS Surface", &is_bezier, 0);
    ImGui::Checkbox("Show Surface Mesh", &show_surface_mesh);
    ImGui::Checkbox("Show Normals", &show_normals);
    ImGui::Checkbox("Show Tangents (Bézier only)", &show_tangents);
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
