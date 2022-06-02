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

int grid_x_count = 10;
int grid_y_count = 10;
float control_points_radius = 0.1;
bool show_grid = false;
bool show_control_mesh = false;

int is_bezier = 1; // 0 -> mls surface
bool show_surface_mesh = false;
int subdivision_count = 0;
float mls_radius = 0.1;

class SpatialData
{
public:
    std::unique_ptr<KdTree2> kd_tree_2;
    std::vector<float> minima;
    std::vector<float> maxima;
    SpatialData(std::vector<Eigen::Vector3f> const &points, std::vector<float> mins, std::vector<float> maxs)
    {
        kd_tree_2 = std::make_unique<KdTree2>(points);
        minima = std::move(mins);
        maxima = std::move(maxs);
    }
    std::vector<Eigen::Vector3f> control_mesh_nodes;
    std::vector<std::array<int, 2>> control_mesh_edges;
};

std::unique_ptr<SpatialData> spatial_data;

double wendland(double d)
{
    return std::pow((1 - d), 4) * (4 * d + 1);
}

float weightedLeastSquares(float u, float v, Eigen::Matrix<double, 6, 1> coefficients)
{
    return coefficients[0] + coefficients[1] * u + coefficients[2] * v + coefficients[3] * u * u + coefficients[4] * u * v + coefficients[5] * v * v;
}

Eigen::Vector3f weightedLeastSquaresNormal(float u, float v, Eigen::Matrix<double, 6, 1> coefficients)
{
    Eigen::Vector3f derivative_cross_product;
    derivative_cross_product << (-1) * (coefficients[1] + 2 * coefficients[3] * u + coefficients[4] * v), (-1) * (coefficients[2] + coefficients[4] * u + 2 * coefficients[5] * v), 1;
    return derivative_cross_product / derivative_cross_product.norm();
}

Eigen::Matrix<double, 6, 1> weightedLeastSquaresCoefficients(float u, float v, float radius)
{
    Eigen::Vector2f p;
    p << u, v;
    std::vector<Eigen::Vector3f> collected_points = spatial_data->kd_tree_2->collectInRadius(p, radius);

    Eigen::Matrix<double, 6, 6> A = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 6, 1> b = Eigen::Matrix<double, 6, 1>::Zero();
    Eigen::Vector2d d_p = p.cast<double>();
    for (Eigen::Vector3f collected_point : collected_points)
    {
        Eigen::Vector3d d_collected_point = collected_point.cast<double>();
        double wendland_value = wendland((d_collected_point.head<2>() - d_p).norm());
        double u_i = d_collected_point[0];
        double v_i = d_collected_point[1];

        Eigen::Matrix<double, 6, 1> base;
        base << 1, u_i, v_i, u_i * u_i, u_i * v_i, v_i * v_i;

        Eigen::Matrix<double, 6, 6> A_i = (base * base.transpose()) * wendland_value;
        A += A_i;
        Eigen::Matrix<double, 6, 1> b_i = base * collected_point[2] * wendland_value;
        b += b_i;
    }

    return A.llt().solve(b);
}

std::pair<Eigen::Vector3f, Eigen::Vector3f> get_bezier_surface_point_and_normal(float p_u, float p_v)
{
    std::vector<Eigen::Vector3f> bezier_points_x;
    for (int i_y = 0; i_y <= grid_y_count; i_y++)
    {
        std::vector<Eigen::Vector3f> points(spatial_data->control_mesh_nodes.begin() + i_y * (grid_x_count + 1), spatial_data->control_mesh_nodes.begin() + (i_y + 1) * (grid_x_count + 1));
        DeCasteljau de_casteljau = DeCasteljau(points, p_u);
        bezier_points_x.push_back(de_casteljau.calculate(0, grid_x_count));
    }
    DeCasteljau de_casteljau_y = DeCasteljau(bezier_points_x, p_v);
    Eigen::Vector3f b0_y = de_casteljau_y.calculate(0, grid_y_count - 1);
    Eigen::Vector3f b1_y = de_casteljau_y.calculate(1, grid_y_count - 1);
    Eigen::Vector3f bezier_point = (p_v * b1_y) + ((1 - p_v) * b0_y);
    Eigen::Vector3f tangent_xy = grid_y_count * (b1_y - b0_y);

    std::vector<Eigen::Vector3f> bezier_points_y;
    for (int i_x = 0; i_x <= grid_x_count; i_x++)
    {
        std::vector<Eigen::Vector3f> points;
        for (int j = i_x; j < (grid_x_count + 1) * (grid_y_count + 1); j += (grid_x_count + 1))
        {
            points.push_back(spatial_data->control_mesh_nodes[j]);
        }

        DeCasteljau de_casteljau = DeCasteljau(points, p_v);
        bezier_points_y.push_back(de_casteljau.calculate(0, grid_y_count));
    }
    DeCasteljau de_casteljau_x = DeCasteljau(bezier_points_y, p_u);
    Eigen::Vector3f b0_x = de_casteljau_x.calculate(0, grid_x_count - 1);
    Eigen::Vector3f b1_x = de_casteljau_x.calculate(1, grid_x_count - 1);
    Eigen::Vector3f tangent_yx = grid_x_count * (b1_x - b0_x);

    Eigen::Vector3f cross_product = tangent_yx.cross(tangent_xy);
    Eigen::Vector3f normal = (cross_product / cross_product.norm());
    return {bezier_point, normal};
}

void updateControlMeshData()
{
    if (spatial_data == nullptr)
        return;

    std::vector<Eigen::Vector3f> nodes;
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
            Eigen::Matrix<double, 6, 1> coefficients = weightedLeastSquaresCoefficients(x_i, y_i, control_points_radius);
            Eigen::Vector3f three_d_point;
            three_d_point << x_i, y_i, weightedLeastSquares(x_i, y_i, coefficients);
            nodes.push_back(three_d_point);

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

    Eigen::VectorXf xs = Eigen::VectorXf::LinSpaced(grid_x_count + 1, min_x, max_x);
    for (float x : xs)
    {
        meshNodes.push_back(Point{x, min_y, 0});
        meshNodes.push_back(Point{x, max_y, 0});
    }

    Eigen::VectorXf ys = Eigen::VectorXf::LinSpaced(grid_y_count + 1, min_y, max_y);
    for (float y : ys)
    {
        meshNodes.push_back(Point{min_x, y, 0});
        meshNodes.push_back(Point{max_x, y, 0});
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
    std::vector<Eigen::Vector3f> nodes;
    int subdivided_grid_x_count = grid_x_count * (subdivision_count + 1) + 1;
    int subdivided_grid_y_count = grid_y_count * (subdivision_count + 1) + 1;
    Eigen::VectorXf xs = Eigen::VectorXf::LinSpaced(subdivided_grid_x_count, 0, 1);
    Eigen::VectorXf ys = Eigen::VectorXf::LinSpaced(subdivided_grid_y_count, 0, 1);

    std::vector<std::array<int, 3>> faces;
    std::vector<Eigen::Vector3f> vector_quantity;

    for (float y : ys)
    {
        for (float x : xs)
        {
            if (is_bezier)
            {
                std::pair<Eigen::Vector3f, Eigen::Vector3f> bezier_result = get_bezier_surface_point_and_normal(x, y);
                nodes.push_back(bezier_result.first);
                vector_quantity.push_back(bezier_result.second);
            }
            else
            {
                Eigen::Matrix<double, 6, 1> coefficients = weightedLeastSquaresCoefficients(x, y, mls_radius);

                Eigen::Vector3f point;
                point << x, y, weightedLeastSquares(x, y, coefficients);
                nodes.push_back(point);
                vector_quantity.push_back(weightedLeastSquaresNormal(x, y, coefficients));
            }
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
        ->setEnabled(false);
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
                std::vector<Eigen::Vector3f> points = std::get<0>(parse_result);
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
