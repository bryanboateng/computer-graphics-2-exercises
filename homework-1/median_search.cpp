#include <vector>
#include <tuple>
#include <algorithm>
#include <stdlib.h>
#include <time.h>   
using Point = std::array<float, 3>;
namespace median_search {
    static float search(std::vector<Point> points, int axis, std::size_t medianPos){
        if(points.size() == 1){
            return points[0][axis];
        }
        srand (time(NULL));
        std::size_t randIndex = rand() % points.size();
        std::vector<Point> left;
        std::vector<Point> right;
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
        if (left.size() == medianPos){
            return points[randIndex][axis];
        }
        if (left.size() >= medianPos){
            return search(left,axis,medianPos);
        }else{
            return search(right,axis,medianPos-left.size()-1);
        }
    }
}
