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
#include <unordered_map>
#include <limits>


#include "args/args.hxx"
#include "portable-file-dialogs.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

#include "median_search.cpp"

using Point = std::array<float, 3>;
using Normal = std::array<float, 3>;

using PointList = std::vector<Point>;
const int BUCKET_SIZE = 400;
int nPts = 0;
int k = 0;
float radius = 0.0314;
float max_x, max_y, max_z = -10000;
float min_x, min_y, min_z = 10000;
std::vector<std::string> structures;

/**
 * Teach polyscope how to handle our datatype
 */

class kd_tree_node{
    public:
    float median;
    kd_tree_node *left;
    kd_tree_node *right;
    kd_tree_node *parent = nullptr;
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
        while(!structures.empty()){
            polyscope::removeStructure(structures.back());
            structures.pop_back();
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
        structures.push_back(uid);
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
    structures.push_back("mainbb");
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
        root = build_tree_linear_median_search(&m_points,0);
    }

    kd_tree_node *build_tree_linear_median_search(std::vector<Point> *pts, int depth){
            if (pts->size() == 0){
            return NULL;
        }
        int axis = depth % 3;
        float median = median_search::search(*pts,axis,pts->size()/2);
        std::vector<Point> left,right;
        kd_tree_node *leftChild = nullptr;
        kd_tree_node *rightChild = nullptr;
        PointList bucket;
        if(pts->size() > BUCKET_SIZE){
            for(std::size_t i = 0; i < pts->size(); ++i) {
                if (pts->at(i)[axis] < median){
                    left.push_back(pts->at(i));
                }else{
                    right.push_back(pts->at(i));
                }
            }
            leftChild = build_tree_linear_median_search(&left, depth+1);
            rightChild = build_tree_linear_median_search(&right, depth+1);
        }else{
            bucket = *pts;
        }
        kd_tree_node *retVal = new kd_tree_node(median,leftChild,rightChild, bucket, depth);
        if (leftChild != nullptr) leftChild->parent = retVal;
        if (rightChild != nullptr) rightChild->parent = retVal;
        return retVal;
    }

    virtual ~SpatialDataStructure() = default;

    PointList const& getPoints() const
    {
        return m_points;
    }
    
    enum axis {x=0, y=1, z=2};

    virtual std::vector<Point> collectInRadiusBruteForce(const Point& p, float radius) const
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
                float distance = EuclideanDistance::measure(p, Point{x,y,z});
                bool compareValue = (nonMatchingSide == cursor->left) ? distance <= radius : distance < radius;
                if(compareValue){
                    collectInRadiusKnn(list,nonMatchingSide,p, radius, (axis+1)%3);
                }
            } 
        }
        return;
    }

    void collectDistanceToBuckets(kd_tree_node * cursor, const Point & p, std::vector<std::pair<kd_tree_node *, float>> *leafDistances) const {
        if (cursor != nullptr){
            if(cursor->bucket.size() > 0){
                if (cursor == root){
                    leafDistances->push_back(std::pair(cursor, 0));
                }else{
                    float splitMedian = cursor->parent->median;
                    int splitAxis = (cursor->parent->depth % 3);
                    float distance;
                    if (splitMedian > p[splitAxis] && cursor==cursor->parent->left){
                        distance = 0;
                    }else if(splitMedian <= p[splitAxis] && cursor==cursor->parent->right){
                        distance = 0;
                    }else{
                        float x = (splitAxis == 0) ? splitMedian: p[0];
                        float y = (splitAxis == 1) ? splitMedian: p[1];
                        float z = (splitAxis == 2) ? splitMedian: p[2];
                        float distance = EuclideanDistance::measure(p,Point{x,y,z});
                    }
                    leafDistances->push_back(std::pair(cursor, distance));
                }
            }
            collectDistanceToBuckets(cursor->left, p, leafDistances);
            collectDistanceToBuckets(cursor->right, p, leafDistances);
        }
    }

    virtual PointList collectNClosest(const PointList &list, const Point &p, size_t n) const{
        if (n >= list.size()){
            return list;
        }
        PointList result;
        std::vector<std::pair<int,float>> distances; //index+distance
        for (size_t i=0; i<list.size(); i++){
            float distance = EuclideanDistance::measure(p,list[i]);
            distances.push_back(std::pair<int,float>(i,distance));
        }
        std::sort(distances.begin(), distances.end(), [=](std::pair<int, float>& a, std::pair<int, float>& b)
        {
            return a.second < b.second;
        });
        for(size_t i=0;i<n;i++){
            result.push_back(list[distances[i].first]);
        }
        return result;
    }

    virtual PointList collectCloserThan(const PointList &list, const Point &p, float maxDistance, size_t maxPoints) const{
        PointList result;
        std::vector<std::pair<int,float>> distances;
        for (size_t i=0; i<list.size(); i++){
            float distance = EuclideanDistance::measure(p,list[i]);
            if (distance <=maxDistance){
                distances.push_back(std::pair<int,float>(i,distance));
            }
        }
        std::sort(distances.begin(), distances.end(), [=](std::pair<int, float>& a, std::pair<int, float>& b)
        {
            return a.second < b.second;
        });
        for(size_t i=0;i<distances.size() && i<maxPoints;i++){
            result.push_back(list[distances[i].first]);
        }
        return result;
    }
    /**
     * Algorithm:
     *  -> traverse entire tree, calculate minimum distance from point to leaves
     *  -> fill up k nearest candidates candidates from closest bucket/s
     *  -> make sure to check all buckets where the minimum distance is less than
     *     the furthest point in the candidates list.
     * */
    virtual PointList collectKNearest(const Point& p, unsigned int k) const
    {
        kd_tree_node * cursor = root;
        bool found = false;
        std::vector<std::pair<kd_tree_node *, float>> bucketDistances;
        collectDistanceToBuckets(root,p,&bucketDistances);
        std::sort(bucketDistances.begin(), bucketDistances.end(), [=](std::pair<kd_tree_node *, float>& a, std::pair<kd_tree_node *, float>& b)
        {
            return a.second < b.second;
        });
        std::priority_queue<std::pair<float, Point>> queue;
        size_t position = 0;
        while(queue.size()<k && position<bucketDistances.size()){
            PointList bucket = bucketDistances[position].first->bucket;
            PointList nClosest = collectNClosest(bucket,p,k-queue.size());
            for(Point entry: nClosest){
                float distance = EuclideanDistance::measure(p,entry);
                queue.push(std::pair<float,Point>(distance,entry));
            }
            if(bucket.size() == nClosest.size()){
                position+=1;
            }
        }
        while(position<bucketDistances.size()){
            float maxDistance = queue.top().first;
            auto bucketPair = bucketDistances[position];
            if (bucketPair.second > maxDistance){
                break;
            }
            PointList bucket = bucketDistances[position].first->bucket;
            PointList nClosest = collectCloserThan(bucket,p,maxDistance,k);
            
            for(auto item: nClosest){
                float distance = EuclideanDistance::measure(p,item);
                auto maxValue = queue.top();
                if(distance<maxValue.first){
                    queue.pop();
                    queue.push(std::pair<float,Point>(distance,item));
                }else{
                    break;
                }
            }
            position += 1;
        }
        PointList result;
        while (!queue.empty())
        {
            result.push_back(queue.top().second);
            queue.pop();
        }
        return result;
    }

    virtual PointList collectKNearestBruteForce(const Point& p, unsigned int k) const {
        return collectNClosest(m_points, p, k);
    }


    void renderKDTree(kd_tree_node* node, std::vector<float> mins, std::vector<float> maxs, int& id) const
    {   
        int axis= node->depth%3;
        if(node->left != nullptr && node->right != nullptr) createHyperplane(axis, node->median, mins, maxs, id);
        if(node->left != nullptr){
            std::vector<float>maxs_new = maxs;
            maxs_new[axis] = node->median;
            id = id + 1;
            renderKDTree(node->left, mins, maxs_new,id);
        }
        if(node->right != nullptr){
            std::vector<float> mins_new = mins;
            mins_new[axis] = node->median;
            id = id + 1;
            renderKDTree(node->right, mins_new, maxs, id);
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
    if (ImGui::Button("Collect in Radius KD Tree")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectInRadius(storedPoints[nPts], radius);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Collect in Radius Brute Force")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectInRadiusBruteForce(storedPoints[nPts], radius);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Collect K Nearest KD Tree")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectKNearest(storedPoints[nPts], k);
        polyscope::PointCloud* rpc = polyscope::registerPointCloud("Points in Radius", radiusPoints);
        rpc->setPointRadius(0.0051);
    }
    if (ImGui::Button("Collect K Nearest Brute Force")){
        PointList storedPoints = sds->getPoints();
        std::vector<Point> radiusPoints =  sds->collectKNearestBruteForce(storedPoints[nPts], k);
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
