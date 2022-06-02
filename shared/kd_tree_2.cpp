#include <cmath>
#include <queue>
#include <utility>

#include "median_search.cpp"
#include "typealiases.cpp"

const int kMaxBucketSize = 400;

static float euclideanDistance2(Eigen::Vector2f const &p1, Eigen::Vector2f const &p2)
{
    return (p1-p2).norm();
}

class KdTreeNode2
{
public:
    float median;
    KdTreeNode2 *left;
    KdTreeNode2 *right;
    KdTreeNode2 *parent = nullptr;
    std::vector<Eigen::Vector3f> bucket;
    int depth;
    KdTreeNode2(float p, KdTreeNode2 *l, KdTreeNode2 *r, std::vector<Eigen::Vector3f> const &b, int d)
    {
        median = p;
        left = l;
        right = r;
        bucket = b;
        depth = d;
    }
};

class KdTree2
{
public:
    KdTreeNode2 *root;
    explicit KdTree2(std::vector<Eigen::Vector3f> points)
            : m_points(std::move(points))
    {
        root = build_tree_linear_median_search(&m_points, 0);
    }

    virtual ~KdTree2() = default;

    [[nodiscard]] std::vector<Eigen::Vector3f> const &getPoints() const
    {
        return m_points;
    }

    [[nodiscard]] virtual std::vector<Eigen::Vector3f> collectInRadius(Eigen::Vector2f p, float radius) const
    {
        std::vector<Eigen::Vector3f> list;
        collectInRadiusKnn(&list, root, p, radius, 0);
        return list;
    }

private:
    KdTreeNode2 *build_tree_linear_median_search(std::vector<Eigen::Vector3f> const *pts, int depth)
    {
        if (pts->empty())
        {
            return nullptr;
        }
        int axis = depth % 2;

        std::vector<float> d;
        for(Eigen::Vector3f const &p : *pts) {
            d.push_back(p[axis]);
        }

        float median = median_search::search(d, pts->size() / 2);
        std::vector<Eigen::Vector3f> left, right;
        KdTreeNode2 *leftChild = nullptr;
        KdTreeNode2 *rightChild = nullptr;
        std::vector<Eigen::Vector3f> bucket;
        if (pts->size() > kMaxBucketSize)
        {
            for (Eigen::Vector3f const &pt : *pts)
            {
                if (pt[axis] < median)
                {
                    left.push_back(pt);
                }
                else
                {
                    right.push_back(pt);
                }
            }
            leftChild = build_tree_linear_median_search(&left, depth + 1);
            rightChild = build_tree_linear_median_search(&right, depth + 1);
        }
        else
        {
            bucket = *pts;
        }
        KdTreeNode2 *retVal = new KdTreeNode2(median, leftChild, rightChild, bucket, depth);
        if (leftChild != nullptr)
            leftChild->parent = retVal;
        if (rightChild != nullptr)
            rightChild->parent = retVal;
        return retVal;
    }

    static void collectInRadiusKnn(std::vector<Eigen::Vector3f> *list, KdTreeNode2 *cursor, Eigen::Vector2f const &p, float radius, int axis)
    {
        if (cursor != nullptr)
        {
            if (!cursor->bucket.empty())
            {
                for (Eigen::Vector3f & i : cursor->bucket)
                {
                    float distance = euclideanDistance2(p, i.head<2>());
                    if (distance <= radius)
                    {
                        list->push_back(i);
                    }
                }
            }
            KdTreeNode2 *nonMatchingSide = nullptr;
            KdTreeNode2 *matchingSide = nullptr;
            // Left child exists
            if (cursor->median > p[axis])
            {
                matchingSide = cursor->left;
                nonMatchingSide = cursor->right;
            }
            else
            {
                matchingSide = cursor->right;
                nonMatchingSide = cursor->left;
            }
            collectInRadiusKnn(list, matchingSide, p, radius, (axis + 1) % 2);
            if (nonMatchingSide != nullptr)
            {
                float x = (axis == 0) ? cursor->median : p[0];
                float y = (axis == 1) ? cursor->median : p[1];
                Eigen::Vector2f v;
                v << x, y;
                float distance = euclideanDistance2(p, v);
                bool compareValue = (nonMatchingSide == cursor->left) ? distance <= radius : distance < radius;
                if (compareValue)
                {
                    collectInRadiusKnn(list, nonMatchingSide, p, radius, (axis + 1) % 2);
                }
            }
        }
    }

    std::vector<Eigen::Vector3f> m_points;
};
