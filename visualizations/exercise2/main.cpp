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

int is_bezier = 1;  // 0 -> mls surface
bool show_surface_mesh = false;
bool show_normals = false;
bool show_tangents = false;
int subdivision_count;
float mls_radius = 0.1;

std::vector<std::string> structures;
polyscope::PointCloud* point_cloud = nullptr;
std::unique_ptr<KdTree> kd_tree;

void callback() {
    if (ImGui::Button("Load Points")) {
        auto paths = pfd::open_file("Load Points", "", std::vector<std::string>{"point data (*.off)", "*.off"}, pfd::opt::none).result();
        if (!paths.empty()) {
            std::filesystem::path path(paths[0]);

            // Read the point cloud
            try {
                auto parse_result = parseOff(path.string());
                std::vector<Point> points = std::get<0>(parse_result);
                std::vector<float> minima = std::get<2>(parse_result);
                std::vector<float> maxima = std::get<3>(parse_result);
                while (!structures.empty()) {
                    polyscope::removeStructure(structures.back());
                    structures.pop_back();
                }

                // Create the polyscope geometry
                point_cloud = polyscope::registerPointCloud("Points", points);

                // Build spatial data structure
                kd_tree = std::make_unique<KdTree>(points, minima, maxima);
            } catch (const std::invalid_argument& e) {
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
    ImGui::SliderInt("##grid_x_count", &grid_x_count, 1, 30);
    ImGui::SameLine();
    ImGui::SliderInt("##grid_y_count", &grid_y_count, 1, 30);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Text("Width, Height");
    ImGui::SliderFloat("Radius", &control_points_radius, 0.0f, 1.0f, "%.3f");
    ImGui::Checkbox("Show Grid", &show_grid);
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

int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("Computer Graphics 2 Sample Code.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);

    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
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
