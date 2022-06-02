#include <Eigen/Dense>

#include "../../colors.cpp"
#include "../../kd_tree.cpp"
#include "../../kd_tree_2.cpp"
#include "../../off_file_parser.cpp"
#include "../../typealiases.cpp"
#include "args/args.hxx"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "portable-file-dialogs.h"

int grid_x_count = 10;
int grid_y_count = 10;
float control_points_radius = 0.1;
bool show_grid = false;
bool show_control_mesh = false;

int is_bezier = 1; // 0 -> mls surface
bool show_surface_mesh = false;
bool show_normals = false;
int subdivision_count = 0;
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
    PointList control_mesh_nodes;
    std::vector<std::array<int, 2>> control_mesh_edges;
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

std::pair<Point, Point> deCasteljau(PointList const &points, int i, int r, float u)
{
    if (r == 0)
    {
        return std::pair<Point, Point>(points[i], {0, 0, 0});
    }
    else
    {
        Point a = deCasteljau(points, i + 1, r - 1, u).first;
        Point b = deCasteljau(points, i, r - 1, u).first;
        Eigen::Vector3f vec_a;
        vec_a << a[0], a[1], a[2];
        Eigen::Vector3f vec_b;
        vec_b << b[0], b[1], b[2];
        Eigen::Vector3f result = (u * vec_a) + ((1 - u) * vec_b);
        Eigen::Vector3f tangent = r * (vec_a - vec_b);
        return std::pair<Point, Point>({result[0], result[1], result[2]}, {tangent[0], tangent[1], tangent[2]});
    }
}

std::pair<Point, Point> get_bezier_surface_point(float p_u, float p_v)
{
    PointList points0;
    for (int i_y = 0; i_y <= grid_y_count; i_y++)
    {
        std::vector<Point> v1(spatial_data->control_mesh_nodes.begin() + i_y * (grid_x_count + 1), spatial_data->control_mesh_nodes.begin() + (i_y + 1) * (grid_x_count + 1));
        points0.push_back(deCasteljau(v1, 0, grid_x_count, p_u).first);
    }
    Point point = deCasteljau(points0, 0, grid_y_count, p_v).first;
    Point tangent_0 = deCasteljau(points0, 0, grid_y_count, p_v).second;

    PointList points1;
    for (int i_x = 0; i_x <= grid_x_count; i_x++)
    {
        std::vector<Point> v2;
        for (int j = i_x; j < (grid_x_count + 1) * (grid_y_count + 1); j += (grid_x_count + 1))
        {
            v2.push_back(spatial_data->control_mesh_nodes[j]);
        }

        points1.push_back(deCasteljau(v2, 0, grid_y_count, p_v).first);
    }
    Point tangent_1 = deCasteljau(points1, 0, grid_x_count, p_u).second;
    Eigen::Vector3f t0;
    t0 << tangent_0[0], tangent_0[1], tangent_0[2];
    Eigen::Vector3f t1;
    t1 << tangent_1[0], tangent_1[1], tangent_1[2];
    Eigen::Vector3f p;
    p << point[0], point[1], point[2];

    Eigen::Vector3f ef = (t1).cross((t0));
    Eigen::Vector3f normal = (ef / ef.norm());
    return std::pair<Point, Point>(point, {normal[0], normal[1], normal[2]});
}

void updateControlMeshData()
{
    if (spatial_data == nullptr)
        return;

    PointList nodes;
    std::vector<std::array<int, 2>> edges;
    float min_x = spatial_data->minima[0];
    float max_x = spatial_data->maxima[0];
    float min_y = spatial_data->minima[1];
    float max_y = spatial_data->maxima[1];

    Eigen::VectorXf wleofi = Eigen::VectorXf::LinSpaced(grid_x_count + 1, min_x, max_x);
    Eigen::VectorXf ghwerhe = Eigen::VectorXf::LinSpaced(grid_y_count + 1, min_y, max_y);
    int i = 0;
    for (float y_i : ghwerhe)
    {
        for (float x_i : wleofi)
        {
            nodes.push_back(Point{x_i, y_i, weightedLeastSquares(x_i, y_i)});

            if (i % (grid_x_count + 1) != grid_x_count)
            {
                edges.push_back({i, i + 1});
            }

            int vertical_neighbor_index = i + grid_x_count + 1;
            if (vertical_neighbor_index < ((grid_x_count + 1) * (grid_y_count + 1)))
            {
                edges.push_back({i, vertical_neighbor_index});
            }
            i++;
        }
    }
    spatial_data->control_mesh_nodes = nodes;
    spatial_data->control_mesh_edges = edges;
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

void createControlMesh()
{
    polyscope::registerCurveNetwork("controlMesh", spatial_data->control_mesh_nodes, spatial_data->control_mesh_edges)
        ->setColor(kRed)
        ->setRadius(0.00064);

    polyscope::registerPointCloud("controlMeshPoints", spatial_data->control_mesh_nodes)
        ->setPointRadius(0.004)
        ->setPointColor(kRed);
}

void createSurfaceMesh()
{
    PointList nodes;
    int subdivided_grid_x_count = grid_x_count * (subdivision_count + 1) + 1;
    int subdivided_grid_y_count = grid_y_count * (subdivision_count + 1) + 1;
    Eigen::VectorXf xs = Eigen::VectorXf::LinSpaced(subdivided_grid_x_count, 0, 1);
    Eigen::VectorXf ys = Eigen::VectorXf::LinSpaced(subdivided_grid_y_count, 0, 1);

    std::vector<std::array<int, 3>> faces;
    std::vector<std::array<float, 3>> vector_quantity;

    for (float y : ys)
    {
        for (float x : xs)
        {
            auto iwn = get_bezier_surface_point(x, y);
            nodes.push_back(iwn.first);
            vector_quantity.push_back(iwn.second);
        }
    }

    int nodes_size = nodes.size();
    for (int i = 0; i < nodes_size; i++)
    {
        if (i % subdivided_grid_x_count != 0)
        {
            int diagonal_partner = i + subdivided_grid_x_count - 1;
            if (diagonal_partner < nodes_size)
            {
                faces.push_back({i - 1, i, diagonal_partner});
                faces.push_back({diagonal_partner + 1, i, diagonal_partner});
            }
        }
    }

    polyscope::registerSurfaceMesh("surfaceMesh", nodes, faces)
        ->setSurfaceColor(kGreen)
        ->addVertexVectorQuantity("normals", vector_quantity)
        ->setVectorColor(kBlue)
        ->setEnabled(true);
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
        if (polyscope::hasCurveNetwork("grid"))
        {
            polyscope::removeCurveNetwork("grid");
        }
    }
}

void updateControlMesh()
{
    if (spatial_data == nullptr)
        return;

    if (show_control_mesh)
    {
        createControlMesh();
    }
    else
    {
        if (polyscope::hasCurveNetwork("controlMesh"))
        {
            polyscope::removeCurveNetwork("controlMesh");
        }
        if (polyscope::hasPointCloud("controlMeshPoints"))
        {
            polyscope::removePointCloud("controlMeshPoints");
        }
    }
}

void updateSurfaceMesh()
{
    if (spatial_data == nullptr)
        return;

    if (show_surface_mesh)
    {
        createSurfaceMesh();
    }
    else
    {
        if (polyscope::hasSurfaceMesh("surfaceMesh"))
        {
            polyscope::removeSurfaceMesh("surfaceMesh");
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
                updateControlMeshData();
                updateControlMesh();
                updateSurfaceMesh();
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
        updateControlMeshData();
        updateControlMesh();
        updateSurfaceMesh();
    }
    ImGui::SameLine();
    if (ImGui::SliderInt("##grid_y_count", &grid_y_count, 1, 30))
    {
        updateGrid();
        updateControlMeshData();
        updateControlMesh();
        updateSurfaceMesh();
    }
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Text("Width, Height");
    if (ImGui::SliderFloat("Radius", &control_points_radius, 0.0f, 1.0f, "%.3f"))
    {
        updateControlMeshData();
        updateControlMesh();
        updateSurfaceMesh();
    }
    if (ImGui::Checkbox("Show Grid", &show_grid))
        updateGrid();
    if (ImGui::Checkbox("Show Control Mesh", &show_control_mesh))
    {
        updateControlMesh();
    }
    ImGui::Text("Surface");
    ImGui::Separator();
    if (ImGui::RadioButton("Beziér Surface", &is_bezier, 1))
        updateSurfaceMesh();
    ImGui::SameLine();
    if (ImGui::RadioButton("MLS Surface", &is_bezier, 0))
        updateSurfaceMesh();
    if (ImGui::Checkbox("Show Surface Mesh", &show_surface_mesh))
        updateSurfaceMesh();
    ImGui::Checkbox("Show Normals", &show_normals);
    if (ImGui::SliderInt("Subdivisions", &subdivision_count, 0, 9))
        updateSurfaceMesh();
    if (ImGui::SliderFloat("Radius (MLS only)", &mls_radius, 0.0f, 1.0f, "%.3f"))
        updateSurfaceMesh();
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
