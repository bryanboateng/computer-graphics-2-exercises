#include <vector>
#include <tuple>
#include <algorithm>
#include <stdlib.h>
#include <time.h>   
using Point = std::array<float, 3>;
class kd_tree_node{
    public:
    Point data;
    kd_tree_node *parent = NULL;
    kd_tree_node *left;
    kd_tree_node *right;
    kd_tree_node(Point p, kd_tree_node *l, kd_tree_node *r ){
        data = p;
        left = l;
        right = r;
    }
    void setParent(kd_tree_node *p){
        parent = p;
    }
};

class kd_tree{
    static const bool DEBUG=false;
    private:
        std::vector<Point>* initPoints;
        kd_tree_node *root;
    public:
        kd_tree(std::vector<Point>* points){
            initPoints = points;
            root = build_tree_linear_median_search(initPoints, 0);
        }
        kd_tree_node *build_tree_using_sort(std::vector<Point> *pts, int depth){
            if (pts->size() == 0){
                return NULL;
            }
            //TODO: use linear median algorithm
            std::sort(pts->begin(), pts->end(), [&](Point a, Point b) {
                return a[depth%3] > b[depth%3];
             });
            
            int median = pts->size()/2;
            std::vector<Point> left,right;
            left.assign(pts->begin(), pts->begin()+median);
            right.assign(pts->begin() + median + 1, pts->end());
            kd_tree_node *leftChild = build_tree_using_sort(&left, depth+1);
            kd_tree_node *rightChild = build_tree_using_sort(&right, depth+1);
            kd_tree_node *retVal = new kd_tree_node(pts->at(median),leftChild,rightChild);
            if (leftChild){
                leftChild->setParent(retVal);
            }
            if (rightChild){
                rightChild->setParent(retVal);
            }
            return retVal;
        }
        kd_tree_node *build_tree_linear_median_search(std::vector<Point> *pts, int depth){
             if (pts->size() == 0){
                return NULL;
            }
            if (pts->size() == 1){
                return new kd_tree_node(pts->at(0),NULL,NULL);
            }
            int axis = depth % 3;
            Point median = median_search(*pts,axis,pts->size()/2,0);
            std::vector<Point> left,right;
            for(std::size_t i = 0; i < pts->size(); ++i) {
                if (pts->at(i)[axis] < median[axis]){
                    left.push_back(pts->at(i));
                }else{
                    right.push_back(pts->at(i));
                }
            }
            kd_tree_node *leftChild = build_tree_linear_median_search(&left, depth+1);
            kd_tree_node *rightChild = build_tree_linear_median_search(&right, depth+1);
            kd_tree_node *retVal = new kd_tree_node(median,leftChild,rightChild);
            if (leftChild){
                leftChild->setParent(retVal);
            }
            if (rightChild){
                rightChild->setParent(retVal);
            }
            return retVal;
        }
        std::vector<Point> collectInRadius(Point p, float radius){
            std::vector<Point> result;
            auto match = findMatch(p);
            return result;
        }
        std::vector<Point> collectKNearest(Point p, int knearest){
            std::vector<Point> result;
            //TODO:
            return result;
        }
        kd_tree_node * findMatch(kd_tree_node *cursor, Point query, int depth){
            if (cursor->data == query){
                return cursor;
            }
            if (query[depth%3] <= cursor->data[depth%3]){
                if (cursor->left == NULL){
                    return cursor;
                }
                return findMatch(cursor->left, query, depth+1);
            }else{
                if (cursor->right == NULL){
                    return cursor;
                }
                return findMatch(cursor->right, query, depth+1);
            }
        }
        kd_tree_node *findMatch(Point query){
            return findMatch(root,query,0);
        }
        kd_tree_node *findNearestNeighbor(Point query){
            auto match = findMatch(query);

            return NULL;
        }
        static void print_vec(std::vector<Point> points, std::string name){
            std::cout << name << ". " <<std::endl;
            for(size_t i=0; i<points.size();i++){
                std::cout << points[i][0] << " ";
            }
            std::cout<<std::endl;
        }

        static Point median_search(std::vector<Point> points, int axis, std::size_t medianPos, int iteration){
            if (DEBUG){
                print_vec(points,"Iteration: " + std::to_string(iteration));
            }
            if(points.size() == 1){
                return points[0];
            }
            srand (time(NULL));
            std::size_t randIndex = rand() % points.size();
            std::vector<Point> left;
            std::vector<Point> right;
            if (DEBUG){
                std::cout << "MEDIAN_POS:"<< medianPos << " ARR_SIZE:" << points.size()<< " RAND_IND:"<< randIndex << " " <<std::endl;
            }
            for(std::size_t i = 0; i < points.size(); ++i) {
                if (i == randIndex){
                    continue;
                }
                Point p = points[i];
                if (p[axis] < points[randIndex][axis]){
                    left.push_back(p);
                } else{
                    right.push_back(p);
                }
            }
            if(DEBUG){
                print_vec(left,"left:");
                print_vec(right,"right:");
            }
            if (left.size() == medianPos){
                return points[randIndex];
            }
            if (left.size() >= medianPos){
                return median_search(left,axis,medianPos, iteration+1);
            }else{
                return median_search(right,axis,medianPos-left.size()-1, iteration+1);
            }
        }
};

