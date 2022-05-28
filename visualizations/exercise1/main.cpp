#include "../../kd_tree.cpp"
#include "../../off_file_parser.cpp"
#include "../../typealiases.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "portable-file-dialogs.h"

int point_index = 0;
int k = 0;
float radius = 0.0314;
std::vector<std::string> structures;
polyscope::PointCloud* point_cloud = nullptr;
std::unique_ptr<KdTree> kd_tree;

void createHyperplane(int axis, float median, std::vector<float> minima, std::vector<float> maxima, int id) {
    std::vector<std::array<int, 2>> edges{{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    std::string uid = "kd grid_" + std::to_string(id);

    if (axis == 0) {
        PointList meshNodes{Point{median, minima[1], minima[2]}, Point{median, maxima[1], minima[2]}, Point{median, maxima[1], maxima[2]}, Point{median, minima[1], maxima[2]}};
        polyscope::registerCurveNetwork(uid, meshNodes, edges);
    }
    if (axis == 1) {
        PointList meshNodes{Point{minima[0], median, minima[2]}, Point{maxima[0], median, minima[2]}, Point{maxima[0], median, maxima[2]}, Point{minima[0], median, maxima[2]}};
        polyscope::registerCurveNetwork(uid, meshNodes, edges);
    }
    if (axis == 2) {
        PointList meshNodes{Point{minima[0], minima[1], median}, Point{maxima[0], minima[1], median}, Point{maxima[0], maxima[1], median}, Point{minima[0], maxima[1], median}};
        polyscope::registerCurveNetwork(uid, meshNodes, edges);
    }
    structures.push_back(uid);
    polyscope::getCurveNetwork(uid)->setColor({255, 0, 0});
    polyscope::getCurveNetwork(uid)->setRadius(0.001);
}

void createMainBox(std::vector<float> minima, std::vector<float> maxima) {
    PointList meshNodes{
        Point{minima[0], minima[1], minima[2]}, Point{minima[0], maxima[1], minima[2]},
        Point{minima[0], maxima[1], maxima[2]}, Point{minima[0], minima[1], maxima[2]},
        Point{maxima[0], minima[1], minima[2]}, Point{maxima[0], maxima[1], minima[2]},
        Point{maxima[0], maxima[1], maxima[2]}, Point{maxima[0], minima[1], maxima[2]}};
    std::vector<std::array<int, 2>> edges{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {4, 5}, {2, 6}, {3, 7}, {6, 7}, {4, 7}, {5, 6}};
    structures.push_back("mainbb");
    polyscope::registerCurveNetwork("mainbb", meshNodes, edges);
}

void renderKDTree(KdTreeNode* node, std::vector<float> mins, std::vector<float> maxs, int& id) {
    int axis = node->depth % 3;
    if (node->left != nullptr && node->right != nullptr) createHyperplane(axis, node->median, mins, maxs, id);
    if (node->left != nullptr) {
        std::vector<float> maxs_new = maxs;
        maxs_new[axis] = node->median;
        id = id + 1;
        renderKDTree(node->left, mins, maxs_new, id);
    }
    if (node->right != nullptr) {
        std::vector<float> mins_new = mins;
        mins_new[axis] = node->median;
        id = id + 1;
        renderKDTree(node->right, mins_new, maxs, id);
    }
}

void callback() {
    if (ImGui::Button("Load Off file")) {
        auto paths = pfd::open_file("Load Off", "", std::vector<std::string>{"point data (*.off)", "*.off"}, pfd::opt::none).result();
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

    ImGui::InputInt("Point Index", &point_index);
    ImGui::InputFloat("radius", &radius);
    ImGui::InputInt("k", &k);
    if (ImGui::Button("Collect in Radius")) {
        PointList storedPoints = kd_tree->getPoints();
        std::vector<Point> radiusPoints = kd_tree->collectInRadius(storedPoints[point_index], radius);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Collect K Nearest")) {
        PointList storedPoints = kd_tree->getPoints();
        std::vector<Point> radiusPoints = kd_tree->collectKNearest(storedPoints[point_index], k);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Render KD Tree")) {
        int id = 0;
        createMainBox(kd_tree->minima, kd_tree->maxima);
        polyscope::getCurveNetwork("mainbb")->setColor({255, 0, 0});
        polyscope::getCurveNetwork("mainbb")->setRadius(0.001);
        renderKDTree(kd_tree->root, kd_tree->minima, kd_tree->maxima, id);
    }
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
    // Initialize polyscope
    polyscope::init();

    // Add a few gui elements
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();

    return 0;
}
