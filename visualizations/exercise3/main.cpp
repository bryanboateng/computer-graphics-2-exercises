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
int grid_z_count = 10;
float control_points_radius = 0.1;
bool show_grid = false;
bool show_control_mesh = false;
bool show_grid_points = false;
bool show_grid_points_fv = false;

int is_bezier = 1; // 0 -> mls surface
bool show_surface_mesh = false;
int subdivision_count = 0;
float mls_radius = 0.1;

template<typename T>
struct matrix_hash : std::unary_function<T, size_t> {
    std::size_t operator()(T const& matrix) const {
        // Note that it is oblivious to the storage order of Eigen matrix (column- or
        // row-major). It will give you the same hash value for two different matrices if they
        // are the transpose of each other in different storage order.
        size_t seed = 0;
        for (size_t i = 0; i < (size_t) matrix.size(); ++i) {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

class SpatialData
{
public:
    void BoundingBoxDiagonal();
    std::unique_ptr<KdTree2> kd_tree_2;
    Eigen::Vector3f minima;
    Eigen::Vector3f maxima;
    std::vector<Eigen::Vector3f> points;
    std::vector<Eigen::Vector3f> normals;
    std::vector<Eigen::Vector3f> grid_points;
    std::vector<Eigen::Vector3f> points_offset_pos;
    std::vector<Eigen::Vector3f> points_offset_neg;
    std::unordered_map<Eigen::Vector3f, double, matrix_hash<Eigen::Vector3f>> function_map;
    double boundingbox_diagonal;
    double alpha;
    SpatialData(std::vector<Eigen::Vector3f> const &pts, Eigen::Vector3f mins, Eigen::Vector3f maxs, std::vector<Eigen::Vector3f> const &norms)
    {
        points = pts;
        normals = norms;
        kd_tree_2 = std::make_unique<KdTree2>(points);
        minima << std::move(mins);
        maxima << std::move(maxs);
        boundingbox_diagonal = (maxima - minima).norm();
        alpha = 0.01 *boundingbox_diagonal;
        CreatePointGrid();
    }
    std::vector<Eigen::Vector3f> control_mesh_nodes;
    std::vector<std::array<int, 2>> control_mesh_edges;
    void computeOffsetPoints();
    void CreatePointGrid();
    void ComputeFunctionValues();
};

void SpatialData::computeOffsetPoints()
{
    for(int i = 0; i < (int) points.size(); ++i){
        Eigen::Vector3f vec_pos(points[i] + normals[i].normalized() * alpha);
        Eigen::Vector3f vec_neg(points[i] - normals[i].normalized() * alpha);
        points_offset_pos.push_back(vec_pos);
        points_offset_neg.push_back(vec_neg);
        function_map[points[i]] = 0;
        function_map[vec_pos] = alpha;
        function_map[vec_neg] = -alpha;
        polyscope::registerPointCloud("Points2", points_offset_pos)
                ->setPointRadius(0.0025)
                ->setPointColor(kGreen);
        polyscope::registerPointCloud("Points3", points_offset_neg)
                ->setPointRadius(0.0025)
                ->setPointColor(kRed);
    }

}

void SpatialData::CreatePointGrid(){
    float min_x = minima[0];
    float max_x = maxima[0];
    float min_y = minima[1];
    float max_y = maxima[1];
    float min_z = minima[2];
    float max_z = maxima[2];

    Eigen::VectorXf grid_x = Eigen::VectorXf::LinSpaced(grid_x_count + 1, min_x, max_x);
    Eigen::VectorXf grid_y = Eigen::VectorXf::LinSpaced(grid_y_count + 1, min_y, max_y);
    Eigen::VectorXf grid_z = Eigen::VectorXf::LinSpaced(grid_z_count + 1, min_z, max_z);

    for(float z : grid_z)
    {
        for (float y : grid_y)
        {
            for (float x : grid_x)
            {
                Eigen::Vector3f grid_point(x,y,z);
                grid_points.push_back(grid_point);
            }
        }
    }

}



std::unique_ptr<SpatialData> spatial_data;





double wendland(double d)
{
    return std::pow((1 - d), 4) * (4 * d + 1);
}

float weightedLeastSquares(Eigen::Vector3f point, Eigen::Matrix<double, 1, 1> coefficients)
{
    return coefficients[0];
}

Eigen::Vector3f weightedLeastSquaresNormal(float u, float v, Eigen::Matrix<double, 1, 1> coefficients)
{
    Eigen::Vector3f derivative_cross_product;
    derivative_cross_product <<  1,1,1;
    return derivative_cross_product / derivative_cross_product.norm();
}

Eigen::Matrix<double, 1, 1> weightedLeastSquaresCoefficients(Eigen::Vector3f point, float radius)
{
    std::vector<Eigen::Vector3f> collected_points = spatial_data->kd_tree_2->collectInRadius(point, radius);
    Eigen::Matrix<double, 1, 1> A = Eigen::Matrix<double, 1, 1>::Zero();
    Eigen::Matrix<double, 1, 1> b = Eigen::Matrix<double, 1, 1>::Zero();
    for (Eigen::Vector3f collected_point : collected_points)
    {
        Eigen::Vector3d d_collected_point = collected_point.cast<double>();
        Eigen::Vector3d d_point = point.cast<double>();
        double wendland_value = wendland((d_collected_point - d_point).norm());

        Eigen::Matrix<double, 1, 1> base;
        base << 1;

        Eigen::Matrix<double, 1, 1> A_i = (base * base.transpose()) * wendland_value;
        A += A_i;
        Eigen::Matrix<double, 1, 1> b_i = base * spatial_data->function_map[collected_point] * wendland_value;
        b += b_i;
    }

    return A.llt().solve(b);
}

void ComputeFunctionValues(){
    std::vector<std::array<double, 3>> colors;

    for(auto point : spatial_data->grid_points){
        double func_val = weightedLeastSquaresCoefficients(point, mls_radius)[0];
        spatial_data->function_map[point] = func_val;
        if(func_val <= 0){
            colors.push_back({255.,0.,0.});
        }
        else{
            colors.push_back({0.,255.,0.});;
        }
    }
    polyscope::getPointCloud("PointGrid")->addColorQuantity("Color Values", colors);
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
    float min_z = spatial_data->minima[2];
    float max_z = spatial_data->maxima[2];

    Eigen::VectorXf grid_x = Eigen::VectorXf::LinSpaced(grid_x_count + 1, min_x, max_x);
    Eigen::VectorXf grid_y = Eigen::VectorXf::LinSpaced(grid_y_count + 1, min_y, max_y);
    Eigen::VectorXf grid_z = Eigen::VectorXf::LinSpaced(grid_z_count + 1, min_z, max_z);
    int i = 0;
    for(float z_i : grid_z){
        for (float y_i : grid_y)
        {
            for (float x_i : grid_x)
            {
                Eigen::Vector3f grid_point(x_i, y_i, z_i);
                Eigen::Matrix<double, 1, 1> coefficients = weightedLeastSquaresCoefficients(grid_point, control_points_radius);
                Eigen::Vector3f three_d_point;
                three_d_point << grid_point;
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
    int subdivided_grid_z_count = grid_z_count * (subdivision_count + 1) + 1;
    Eigen::VectorXf xs = Eigen::VectorXf::LinSpaced(subdivided_grid_x_count, 0, 1);
    Eigen::VectorXf ys = Eigen::VectorXf::LinSpaced(subdivided_grid_y_count, 0, 1);
    Eigen::VectorXf zs = Eigen::VectorXf::LinSpaced(subdivided_grid_z_count, 0, 1);

    std::vector<std::array<int, 3>> faces;
    std::vector<Eigen::Vector3f> vector_quantity;

    for(float z : zs)
    {
        for (float y : ys)
        {
            for (float x : xs)
            {
                Eigen::Vector3f grid_point(x,y,z);
                Eigen::Matrix<double, 1, 1> coefficients = weightedLeastSquaresCoefficients(grid_point, mls_radius);

                Eigen::Vector3f point;
                point << weightedLeastSquares(grid_point, coefficients);
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

void updatePointGrid(){
    if (spatial_data == nullptr)
        return;
    if (show_grid_points){
        polyscope::registerPointCloud("PointGrid", spatial_data->grid_points)
                ->setPointRadius(0.0025)
                ->setPointColor(kBlue);
    }
    else{
        if (polyscope::hasPointCloud("PointGrid"))
        {
            polyscope::removePointCloud("PointGrid");
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
                Eigen::Vector3f minima = std::get<2>(parse_result);
                Eigen::Vector3f maxima = std::get<3>(parse_result);

                // Build spatial data structure
                spatial_data = std::make_unique<SpatialData>(points, minima, maxima, normals);

                // Create the polyscope geometry
                polyscope::registerPointCloud("Points", points)
                        ->setPointRadius(0.0025)
                        ->setPointColor(kOrange)->addVectorQuantity("normals", normals)->setEnabled(true);
                updateGrid();
                updateControlMeshData();
                updateControlMesh();
                updateSurfaceMesh();
                spatial_data->computeOffsetPoints();

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
    // if (ImGui::Checkbox("Show Grid", &show_grid))
    //     updateGrid();
    // if (ImGui::Checkbox("Show Control Mesh", &show_control_mesh))
    // {
    //     updateControlMesh();
    // }
    // ImGui::Text("Surface");
    // ImGui::Separator();
    // if (ImGui::RadioButton("Beziér Surface", &is_bezier, 1))
    //     updateSurfaceMesh();
    // ImGui::SameLine();
    // if (ImGui::RadioButton("MLS Surface", &is_bezier, 0))
    //     updateSurfaceMesh();
    // if (ImGui::Checkbox("Show Surface Mesh", &show_surface_mesh))
    //     updateSurfaceMesh();
    // if (ImGui::SliderInt("Subdivisions", &subdivision_count, 0, 9))
    //     updateSurfaceMesh();
    // if (ImGui::SliderFloat("Radius (MLS only)", &mls_radius, 0.0f, 1.0f, "%.3f"))
    //     updateSurfaceMesh();
    if (ImGui::Checkbox("Show Grid Points", &show_grid_points))
        updatePointGrid();
    if (ImGui::Checkbox("Show Grid Points Function Values", &show_grid_points_fv))
        ComputeFunctionValues();
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