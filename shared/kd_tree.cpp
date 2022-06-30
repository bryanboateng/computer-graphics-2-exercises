#include <cmath>
#include <queue>

#include "median_search.cpp"
#include "typealiases.cpp"

const int kMaxBucketSize = 400;

static float euclideanDistance(Eigen::Vector3f const &p1, Eigen::Vector3f const &p2)
{
    return (p1-p2).norm();
}


struct Compare {
    bool operator()(std::pair<float, Eigen::Vector3f> const& p1, std::pair<float, Eigen::Vector3f> const& p2)
    {
        return p1.first < p2.first;
    }
};

class KdTreeNode
{
public:
    float median;
    KdTreeNode *left;
    KdTreeNode *right;
    KdTreeNode *parent = nullptr;
    std::vector<Eigen::Vector3f> bucket;
    int depth;
    KdTreeNode(float p, KdTreeNode *l, KdTreeNode *r, std::vector<Eigen::Vector3f> b, int d)
    {
        median = p;
        left = l;
        right = r;
        bucket = b;
        depth = d;
    }
};

class KdTree
{
public:
    KdTreeNode *root;
    KdTree(const std::vector<Eigen::Vector3f> &points)
        : m_points(std::move(points))
    {
        root = build_tree_linear_median_search(&m_points, 0);
    }

    virtual ~KdTree() = default;

    [[nodiscard]] std::vector<Eigen::Vector3f> const &getPoints() const
    {
        return m_points;
    }

    virtual std::vector<Eigen::Vector3f> collectInRadius(const Eigen::Vector3f &p, float radius) const
    {
        std::vector<Eigen::Vector3f> list;
        collectInRadiusKnn(&list, root, p, radius, 0);
        return list;
    }

    /**
     * Algorithm:
     *  -> traverse entire tree, calculate minimum distance from point to leaves
     *  -> fill up k nearest candidates candidates from closest bucket/s
     *  -> make sure to check all buckets where the minimum distance is less than
     *     the furthest point in the candidates list.
     * */
    virtual std::vector<Eigen::Vector3f> collectKNearest(const Eigen::Vector3f &p, unsigned int k) const
    {
        std::vector<Eigen::Vector3f> result;
        if (k == 0)
        {
            return result;
        }
        KdTreeNode *cursor = root;
        bool found = false;
        std::vector<std::pair<KdTreeNode *, float>> bucketDistances;
        collectDistanceToBuckets(root, p, &bucketDistances);
        std::sort(bucketDistances.begin(), bucketDistances.end(), [=](std::pair<KdTreeNode *, float> &a, std::pair<KdTreeNode *, float> &b)
                  { return a.second < b.second; });
        std::priority_queue<std::pair<float, Eigen::Vector3f>, std::vector<std::pair<float, Eigen::Vector3f>>, Compare> queue;
        size_t position = 0;
        while (queue.size() < k && position < bucketDistances.size())
        {
            std::vector<Eigen::Vector3f> bucket = bucketDistances[position].first->bucket;
            std::vector<Eigen::Vector3f> nClosest = collectNClosest(bucket, p, k - queue.size());
            for (Eigen::Vector3f entry : nClosest)
            {
                float distance = euclideanDistance(p, entry);
                queue.push(std::pair<float, Eigen::Vector3f>(distance, entry));
            }
            if (bucket.size() == nClosest.size())
            {
                position += 1;
            }
        }
        while (position < bucketDistances.size())
        {
            float maxDistance = queue.top().first;
            auto bucketPair = bucketDistances[position];
            if (bucketPair.second > maxDistance)
            {
                break;
            }
            std::vector<Eigen::Vector3f> bucket = bucketDistances[position].first->bucket;
            std::vector<Eigen::Vector3f> nClosest = collectCloserThan(bucket, p, maxDistance, k);

            for (Eigen::Vector3f item : nClosest)
            {
                float distance = euclideanDistance(p, item);
                std::pair<float, Eigen::Vector3f> maxValue = queue.top();
                if (distance < maxValue.first)
                {
                    queue.pop();
                    queue.push(std::pair<float, Eigen::Vector3f>(distance, item));
                }
                else
                {
                    break;
                }
            }
            position += 1;
        }
        while (!queue.empty())
        {
            result.push_back(queue.top().second);
            queue.pop();
        }
        return result;
    }

private:
    KdTreeNode *build_tree_linear_median_search(std::vector<Eigen::Vector3f> *pts, int depth)
    {
        if (pts->size() == 0)
        {
            return NULL;
        }
        int axis = depth % 3;

        std::vector<float> d;
        for(Eigen::Vector3f const &p : *pts) {
            d.push_back(p[axis]);
        }

        float median = median_search::search(d, pts->size() / 2);
        std::vector<Eigen::Vector3f> left, right;
        KdTreeNode *leftChild = nullptr;
        KdTreeNode *rightChild = nullptr;
        std::vector<Eigen::Vector3f> bucket;
        if (pts->size() > kMaxBucketSize)
        {
            for (std::size_t i = 0; i < pts->size(); ++i)
            {
                if (pts->at(i)[axis] < median)
                {
                    left.push_back(pts->at(i));
                }
                else
                {
                    right.push_back(pts->at(i));
                }
            }
            leftChild = build_tree_linear_median_search(&left, depth + 1);
            rightChild = build_tree_linear_median_search(&right, depth + 1);
        }
        else
        {
            bucket = *pts;
        }
        KdTreeNode *retVal = new KdTreeNode(median, leftChild, rightChild, bucket, depth);
        if (leftChild != nullptr)
            leftChild->parent = retVal;
        if (rightChild != nullptr)
            rightChild->parent = retVal;
        return retVal;
    }

    virtual std::vector<Eigen::Vector3f> collectNClosest(const std::vector<Eigen::Vector3f> &list, const Eigen::Vector3f &p, size_t n) const
    {
        if (n >= list.size())
        {
            return list;
        }
        std::vector<Eigen::Vector3f> result;
        std::vector<std::pair<int, float>> distances; // index+distance
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance(p, list[i]);
            distances.push_back(std::pair<int, float>(i, distance));
        }
        std::sort(distances.begin(), distances.end(), [=](std::pair<int, float> &a, std::pair<int, float> &b)
                  { return a.second < b.second; });
        for (size_t i = 0; i < n; i++)
        {
            result.push_back(list[distances[i].first]);
        }
        return result;
    }

    virtual std::vector<Eigen::Vector3f> collectCloserThan(const std::vector<Eigen::Vector3f> &list, const Eigen::Vector3f &p, float maxDistance, size_t maxPoints) const
    {
        std::vector<Eigen::Vector3f> result;
        std::vector<std::pair<int, float>> distances;
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance(p, list[i]);
            if (distance <= maxDistance)
            {
                distances.push_back(std::pair<int, float>(i, distance));
            }
        }
        std::sort(distances.begin(), distances.end(), [=](std::pair<int, float> &a, std::pair<int, float> &b)
                  { return a.second < b.second; });
        for (size_t i = 0; i < distances.size() && i < maxPoints; i++)
        {
            result.push_back(list[distances[i].first]);
        }
        return result;
    }

    void collectDistanceToBuckets(KdTreeNode *cursor, const Eigen::Vector3f &p, std::vector<std::pair<KdTreeNode *, float>> *leafDistances) const
    {
        if (cursor != nullptr)
        {
            if (cursor->bucket.size() > 0)
            {
                if (cursor == root)
                {
                    leafDistances->push_back(std::pair(cursor, 0));
                }
                else
                {
                    float splitMedian = cursor->parent->median;
                    int splitAxis = (cursor->parent->depth % 3);
                    float distance;
                    if (splitMedian > p[splitAxis] && cursor == cursor->parent->left)
                    {
                        distance = 0;
                    }
                    else if (splitMedian <= p[splitAxis] && cursor == cursor->parent->right)
                    {
                        distance = 0;
                    }
                    else
                    {
                        float x = (splitAxis == 0) ? splitMedian : p[0];
                        float y = (splitAxis == 1) ? splitMedian : p[1];
                        float z = (splitAxis == 2) ? splitMedian : p[2];
                        float distance = euclideanDistance(p, Eigen::Vector3f{x, y, z});
                    }
                    leafDistances->push_back(std::pair(cursor, distance));
                }
            }
            collectDistanceToBuckets(cursor->left, p, leafDistances);
            collectDistanceToBuckets(cursor->right, p, leafDistances);
        }
    }

    static void collectInRadiusKnn(std::vector<Eigen::Vector3f> *list, KdTreeNode *cursor, const Eigen::Vector3f &p, float radius, int axis)
    {
        if (cursor != nullptr)
        {
            if (cursor->bucket.size() > 0)
            {
                for (size_t i = 0; i < cursor->bucket.size(); i++)
                {
                    float distance = euclideanDistance(p, cursor->bucket[i]);
                    if (distance <= radius)
                    {
                        list->push_back(cursor->bucket[i]);
                    }
                }
            }
            KdTreeNode *nonMatchingSide = nullptr;
            KdTreeNode *matchingSide = nullptr;
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
            collectInRadiusKnn(list, matchingSide, p, radius, (axis + 1) % 3);
            if (nonMatchingSide != nullptr)
            {
                float x = (axis == 0) ? cursor->median : p[0];
                float y = (axis == 1) ? cursor->median : p[1];
                float z = (axis == 2) ? cursor->median : p[2];
                float distance = euclideanDistance(p, Eigen::Vector3f{x, y, z});
                bool compareValue = (nonMatchingSide == cursor->left) ? distance <= radius : distance < radius;
                if (compareValue)
                {
                    collectInRadiusKnn(list, nonMatchingSide, p, radius, (axis + 1) % 3);
                }
            }
        }
        return;
    }

    std::vector<Eigen::Vector3f> m_points;
};
