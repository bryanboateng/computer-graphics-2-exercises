#include "polyscope/polyscope.h"

#include "polyscope/combining_hash_functions.h"
#include "polyscope/messages.h"
#include "polyscope/curve_network.h"

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
#include <stdlib.h>
#include <queue> 


#include "args/args.hxx"
#include "portable-file-dialogs.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

using Point = std::array<float, 3>;
using Normal = std::array<float, 3>;

using PointList = std::vector<Point>;
const int BUCKET_SIZE = 1000;
int nPts = 0;
int k = 0;
float radius = 0.0314;
float max_x, max_y, max_z = -10000;
float min_x, min_y, min_z = 10000;

/**
 * Teach polyscope how to handle our datatype
 */

class kd_tree_node{
    public:
    float median;
    kd_tree_node *left;
    kd_tree_node *right;
    PointList bucket;
    int depth;
    kd_tree_node(float p, kd_tree_node *l, kd_tree_node *r, PointList b, int d){
        median = p;
        left = l;
        right = r;
        bucket = b;
        depth = d;
    }
};


float adaptorF_custom_accessVector3Value(const Point& v, unsigned int ind)
{
    return v[ind];
}

void readOff(std::string const& filename, std::vector<Point>* points, std::vector<Normal>* normals = nullptr)
{
    points->clear();
    if (normals) normals->clear();
    max_x = -10000;
    max_y = -10000;
    max_z = -10000;
    min_x = 10000;
    min_y = 10000;
    min_z = 10000;

    std::ifstream off_file(filename);
    if (off_file.is_open()) {
        std::string header_keyword;
        if(!(off_file >> header_keyword)){
            off_file.close();
            throw std::invalid_argument("Unable to read file.");
        } else if (header_keyword == "OFF") {
            int vertex_count, face_count, edge_count;
            off_file >> vertex_count >> face_count >> edge_count;
            for (int i = 0; i < vertex_count; ++i) {
                float x, y, z;
                off_file >> x >> y >> z;
                if(x > max_x) max_x = x;
                if(y > max_y) max_y = y;
                if(z > max_z) max_z = z;
                if(x < min_x) min_x = x;
                if(y < min_y) min_y = y;
                if(z < min_z) min_z = z;
                
                points->push_back(Point{x,y,z});
            }
            off_file.close();
        } else {
            off_file.close();
            throw std::invalid_argument("Incorrect file format. This program only supports point data from OFF-files that contain the header keyword OFF. More about the specifications of OFF-Files can be found at: http://paulbourke.net/dataformats/oogl/#OFF");
        }
    } else {
        throw std::invalid_argument("Unable to read file.");
    }
}

void createHyperplane(int axis, float median, std::vector<float> mins, std::vector<float> maxs, int id){

        std::vector<std::array<int, 2>> edges{{0,1}, {1,2}, {2,3}, {3,0}};
        std::string uid = "kd grid_" + std::to_string(id);

        if(axis == 0){
             PointList meshNodes{Point{median,mins[1],mins[2]},Point{median,maxs[1],mins[2]},Point{median,maxs[1],maxs[2]},Point{median,mins[1],maxs[2]}};
             polyscope::registerCurveNetwork(uid, meshNodes, edges);
        }
        if(axis == 1){
             PointList meshNodes{Point{mins[0],median,mins[2]},Point{maxs[0],median,mins[2]},Point{maxs[0],median,maxs[2]},Point{mins[0],median,maxs[2]}};
             polyscope::registerCurveNetwork(uid, meshNodes, edges);
        }
        if(axis==2){
            PointList meshNodes{Point{mins[0],mins[1],median},Point{maxs[0],mins[1],median},Point{maxs[0],maxs[1],median},Point{mins[0],maxs[1],median}};
            polyscope::registerCurveNetwork(uid, meshNodes, edges);
        }
        polyscope::getCurveNetwork(uid)->setColor({255,0,0});
        polyscope::getCurveNetwork(uid)->setRadius(0.001);
        
}
void createMainBox(std::vector<float> mins, std::vector<float> maxs){
    PointList meshNodes{
        Point{mins[0],mins[1],mins[2]},Point{mins[0],maxs[1],mins[2]},
        Point{mins[0],maxs[1],maxs[2]},Point{mins[0],mins[1],maxs[2]},
        Point{maxs[0],mins[1],mins[2]},Point{maxs[0],maxs[1],mins[2]},
        Point{maxs[0],maxs[1],maxs[2]},Point{maxs[0],mins[1],maxs[2]}
        };
    std::vector<std::array<int, 2>> edges{{0,1}, {1,2}, {2,3}, {3,0},{0,4},{1,5},{4,5},{2,6},{3,7},{6,7},{4,7},{5,6}};
    polyscope::registerCurveNetwork("mainbb", meshNodes, edges);
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
 * This is a KD tree
 */
class SpatialDataStructure
{
public:
    kd_tree_node* root;
    SpatialDataStructure(PointList const& points)
        : m_points(points)
    {
        root = build_tree_using_sort(m_points, 0);

        float median = medianSearch(y);

    }
    kd_tree_node* build_tree_using_sort(PointList &pts, int depth){
            if (pts.size() == 0){
                return NULL;
            }

            std::sort(pts.begin(), pts.end(), [&](Point a, Point b) {
                return a[depth%3] < b[depth%3];
             });
             int median = pts.size()/2;
             kd_tree_node *leftChild = nullptr;
             kd_tree_node *rightChild = nullptr;
             PointList bucket;

             if(pts.size() > BUCKET_SIZE){
                PointList left,right;
                left.assign(pts.begin(), pts.begin()+median);
                right.assign(pts.begin() + median, pts.end());
                leftChild = build_tree_using_sort(left, depth+1);
                rightChild = build_tree_using_sort(right, depth+1);
             }
             else{
                 bucket = pts;
             }
             float median_val = pts[median][depth%3];
             kd_tree_node *root = new kd_tree_node(median_val, leftChild, rightChild, bucket, depth);
            return root;

        }


    virtual ~SpatialDataStructure() = default;

    PointList const& getPoints() const
    {
        return m_points;
    }
    
    enum axis {x=0, y=1, z=2};


    int medianSearch(axis ax) const
    {   
        float median = 0;
        std::vector<PointList> shortLists;
        std::vector<float> medianList;
        for(size_t i=0; i<m_points.size(); i=i+5){
            if(m_points.begin() + i + 5 > m_points.end()){
                shortLists.push_back(PointList(m_points.begin() + i, m_points.end()));
            }
            else{
                shortLists.push_back(PointList(m_points.begin() + i, m_points.begin()+i+5));
            }
        }
        for(auto shortList : shortLists){
            std::sort(shortList.begin(), shortList.end(),
                       [ax](const Point& a, const Point& b) {
                         return a[ax] < b[ax];});
            int shortListLength = shortList.size();
            medianList.push_back(shortList[shortListLength/2][ax]);
        }
        std::sort(medianList.begin(), medianList.end());
        median = medianList[medianList.size()/2];
        return median;
    }


    virtual std::vector<Point> collectInRadiusDUMMY(const Point& p, float radius) const
    {
        std::vector<Point> result;

        // Dummy brute-force implementation
        // TODO: Use spatial data structure for sub-linear search
        for (size_t i = 0; i < m_points.size(); ++i)
        {
            float distance = EuclideanDistance::measure(p, m_points[i]);
            if (distance <= radius)
                result.push_back(m_points[i]);
        }

        return result;
    }

    virtual std::vector<Point> collectInRadius(const Point& p, float radius) const
    {
        PointList list;
        collectInRadiusKnn(&list, root, p, radius, 0);
        std::cout << "Found: " << list.size() << " Points" << std::endl;
        return list;
    }

     static void collectInRadiusKnn(PointList *list, kd_tree_node * cursor, const Point& p, float radius, int axis){
        if (cursor != nullptr){
            if(cursor->bucket.size() > 0){
                for(size_t i=0; i<cursor->bucket.size();i++){
                    float distance = EuclideanDistance::measure(p, cursor->bucket[i]);
                    if (distance <= radius){
                        list->push_back(cursor->bucket[i]);

                    }
                }
            }
            kd_tree_node *nonMatchingSide = nullptr;
            kd_tree_node *matchingSide = nullptr;
            //Left child exists
            if (cursor->median > p[axis]){
                matchingSide = cursor->left;
                nonMatchingSide = cursor->right;
            }else{
                matchingSide = cursor->right;
                nonMatchingSide = cursor->left;
            }
            collectInRadiusKnn(list,matchingSide, p, radius, (axis+1)%3);
            if (nonMatchingSide != nullptr){
                float x = (axis == 0) ? cursor->median: p[0];
                float y = (axis == 1) ? cursor->median: p[1];
                float z = (axis == 2) ? cursor->median: p[2];
                if(EuclideanDistance::measure(p, Point{x,y,z}) <= radius){
                    collectInRadiusKnn(list,nonMatchingSide,p, radius, (axis+1)%3);
                }
            } 
        }
        return;
    }


    virtual std::vector<Point> collectKNearest(const Point& p, unsigned int k) const
    {
        std::vector<Point> result;

        // Bogus knn implementation, giving you the first k points!
        // TODO: Use spatial data structure for sub-linear search
        for (std::size_t i = 0; (i < k) && (i < m_points.size()); ++i)
        {
            result.push_back(m_points[i]);
        }

        return result;
    }

    void renderKDTree(kd_tree_node* node, std::vector<float> mins, std::vector<float> maxs, int& id) const
    {   
        int axis= node->depth%3;
        if(node->left != nullptr && node->right != nullptr) createHyperplane(axis, node->median, mins, maxs, id);
        if(node->left != nullptr){
            std::vector<float> mins_new = mins;
            mins_new[axis] = node->median;
            id = id + 1;
            renderKDTree(node->left, mins_new, maxs, id);
        }
        if(node->right != nullptr){
            std::vector<float>maxs_new = maxs;
            maxs_new[axis] = node->median;
            id = id + 1;
            renderKDTree(node->right, mins, maxs_new,id);
        }
        

    }
    void print_points(PointList pts)
    {
        std::cout<<"Printing points"<<std::endl;
        for(size_t i=0;i<pts.size();i++){
            printf("(%f %f %f)\n",pts[i][0],pts[i][1],pts[i][2]);
        }
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
                try {
                    std::vector<Point> points;
                    readOff(path.string(), &points);

                    // Create the polyscope geometry
                    pc = polyscope::registerPointCloud("Points", points);

                    // Build spatial data structure
                    sds = std::make_unique<SpatialDataStructure>(points);
                    PointList storedPoints = sds->getPoints();
                }
                catch( const std::invalid_argument& e ) {
                    polyscope::error(e.what());
                    return;
                }
            }
        }
    }


    ImGui::InputInt("Point Number", &nPts);            
    ImGui::InputFloat("radius", &radius);
    ImGui::InputInt("k", &k); 
    if (ImGui::Button("Collect in Radius")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectInRadius(storedPoints[nPts], radius);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Collect K Nearest")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectKNearest(storedPoints[nPts], k);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }   
    if (ImGui::Button("Render KD Tree")){
        int id = 0;
        std::vector<float> mins{min_x,min_y,min_z};
        std::vector<float> maxs{max_x,max_y,max_z};
        createMainBox(mins, maxs);
        polyscope::getCurveNetwork("mainbb")->setColor({255,0,0});
        polyscope::getCurveNetwork("mainbb")->setRadius(0.001);
        sds->renderKDTree(sds->root,mins,maxs, id);
        

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
