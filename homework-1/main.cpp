#include "polyscope/polyscope.h"

#include "polyscope/combining_hash_functions.h"
#include "polyscope/messages.h"

#include "polyscope/file_helpers.h"
#include "polyscope/point_cloud.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <fstream>

#include "args/args.hxx"
#include "portable-file-dialogs.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

using Point = std::array<float, 3>;
using Normal = std::array<float, 3>;

using PointList = std::vector<Point>;

/**
 * Teach polyscope how to handle our datatype
 */
float adaptorF_custom_accessVector3Value(const Point& v, unsigned int ind)
{
    return v[ind];
}

void readOff(std::string const& filename, std::vector<Point>* points, std::vector<Normal>* normals = nullptr)
{
    points->clear();
    if (normals) normals->clear();
    std::string line_fileformat;
    int vertices, faces, normals_amt;
    float x, y, z;

    std::ifstream offile(filename);

    if (offile.is_open()){
        if(!(offile >> line_fileformat)){
            polyscope::warning("Unable to read file."); 
        }
        if(line_fileformat == "OFF"){
            offile >> vertices >> faces >> normals_amt;
            for(int i = 0; i < vertices; ++i){
                offile >> x >> y >> z;
                points->push_back(Point{x,y,z});
            }
        }
        else polyscope::warning("This only works for OFF files."); 
    offile.close();
  }

  else polyscope::warning("Unable to read file."); 
    
}

struct EuclideanDistance
{
    static float measure(Point const& p1, Point const& p2)
    {
        float dx = p1[0] - p2[0];
        float dy = p1[1] - p2[1];
        float dz = p1[2] - p2[2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};

/*
 * This is not yet a spatial data structure :)
 */
class SpatialDataStructure
{
public:
    SpatialDataStructure(PointList const& points)
        : m_points(points)
    {
    }

    virtual ~SpatialDataStructure() = default;

    PointList const& getPoints() const
    {
        return m_points;
    }

    virtual std::vector<std::size_t> collectInRadius(const Point& p, float radius) const
    {
        std::vector<std::size_t> result;

        // Dummy brute-force implementation
        // TODO: Use spatial data structure for sub-linear search
        for (std::size_t i = 0; i < m_points.size(); ++i)
        {
            float distance = EuclideanDistance::measure(p, m_points[i]);
            if (distance <= radius)
                result.push_back(i);
        }

        return result;
    }

    virtual std::vector<std::size_t> collectKNearest(const Point& p, unsigned int k) const
    {
        std::vector<std::size_t> result;

        // Bogus knn implementation, giving you the first k points!
        // TODO: Use spatial data structure for sub-linear search
        for (std::size_t i = 0; (i < k) && (i < m_points.size()); ++i)
        {
            result.push_back(i);
        }

        return result;
    }

private:
    PointList m_points;
};

// Application variables
polyscope::PointCloud* pc = nullptr;
std::unique_ptr<SpatialDataStructure> sds;

void callback() {
    if (ImGui::Button("Load Off file")) {
        auto paths = pfd::open_file("Load Off", "", std::vector<std::string>{ "point data (*.off)", "*.off" }, pfd::opt::none).result();
        if (!paths.empty())
        {
            std::filesystem::path path(paths[0]);

            if (path.extension() == ".off")
            {
                // Read the point cloud
                std::vector<Point> points;
                readOff(path.string(), &points);

                // Create the polyscope geometry
                pc = polyscope::registerPointCloud("Points", points);

                // Build spatial data structure
                sds = std::make_unique<SpatialDataStructure>(points);
            }
        }
    }

    // TODO: Implement radius search
    // TODO: Implement visualizations
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
