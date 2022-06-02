#include "Eigen/Dense"
#include <map>
#include <vector>

class DeCasteljau
{
public:
    std::vector<Eigen::Vector3f> points;
    float u;

    explicit DeCasteljau(std::vector<Eigen::Vector3f> p, float u_)
    {
        points = std::move(p);
        u = u_;
    }

    Eigen::Vector3f calculate(int i, int r)
    {
        if (r == 0)
        {
            return points[i];
        }
        else
        {
            std::pair<int, int> key{i, r};
            auto it = cache.find(key);
            if (it == cache.end())
            {
                Eigen::Vector3f a = calculate(i + 1, r - 1);
                Eigen::Vector3f b = calculate(i, r - 1);

                Eigen::Vector3f value = (u * a) + ((1 - u) * b);
                cache.emplace(key, value);
                return value;
            }
            else
            {
                return it->second;
            }
        }
    }

private:
    std::map<std::pair<int, int>, Eigen::Vector3f> cache;
};