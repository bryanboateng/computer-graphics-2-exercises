#include <cmath>
#include <queue>
#include <utility>

#include "typealiases.cpp"

struct ComparePair {
    bool operator()(std::pair<float, std::pair<Eigen::Vector3f, float>> const& p1, std::pair<float, std::pair<Eigen::Vector3f, float>> const& p2)
    {
        return p1.first < p2.first;
    }
};

class KdTreeNodePair
{
public:
    float median;
    KdTreeNodePair *left;
    KdTreeNodePair *right;
    KdTreeNodePair *parent = nullptr;
    std::vector<std::pair<Eigen::Vector3f, float>> bucket;
    int depth;
    KdTreeNodePair(float p, KdTreeNodePair *l, KdTreeNodePair *r, std::vector<std::pair<Eigen::Vector3f, float>> b, int d)
    {
        median = p;
        left = l;
        right = r;
        bucket = std::move(b);
        depth = d;
    }
};

class KdTreePair
{
public:
    KdTreeNodePair *root;
    explicit KdTreePair(std::vector<std::pair<Eigen::Vector3f, float>> points)
        : m_points(std::move(points))
    {
        root = build_tree_linear_median_search(&m_points, 0);
    }

    virtual ~KdTreePair() = default;

    [[nodiscard]] std::vector<std::pair<Eigen::Vector3f, float>> const &getPoints() const
    {
        return m_points;
    }

    [[nodiscard]] virtual std::vector<std::pair<Eigen::Vector3f, float>> collectInRadius(const Eigen::Vector3f &p, float radius) const
    {
        std::vector<std::pair<Eigen::Vector3f, float>> list;
        collectInRadiusKnn(&list, root, p, radius, 0);
        return list;
    }

    [[nodiscard]] virtual std::vector<std::pair<Eigen::Vector3f, float>> collectKNearest(const Eigen::Vector3f &p, unsigned int k) const
    {
        std::vector<std::pair<Eigen::Vector3f, float>> result;
        if (k == 0)
        {
            return result;
        }
        KdTreeNodePair *cursor = root;
        bool found = false;
        std::vector<std::pair<KdTreeNodePair *, float>> bucketDistances;
        collectDistanceToBuckets(root, p, &bucketDistances);
        std::sort(bucketDistances.begin(), bucketDistances.end(), [=](std::pair<KdTreeNodePair *, float> &a, std::pair<KdTreeNodePair *, float> &b)
        { return a.second < b.second; });
        std::priority_queue<std::pair<float, std::pair<Eigen::Vector3f, float>>, std::vector<std::pair<float, std::pair<Eigen::Vector3f, float>>>, ComparePair> queue;
        size_t position = 0;
        while (queue.size() < k && position < bucketDistances.size())
        {
            std::vector<std::pair<Eigen::Vector3f, float>> bucket = bucketDistances[position].first->bucket;
            std::vector<std::pair<Eigen::Vector3f, float>> nClosest = collectNClosest(bucket, p, k - queue.size());
            for (const std::pair<Eigen::Vector3f, float>& entry : nClosest)
            {
                float distance = euclideanDistance(p, entry.first);
                queue.push(std::pair<float, std::pair<Eigen::Vector3f, float>>(distance, entry));
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
            std::vector<std::pair<Eigen::Vector3f, float>> bucket = bucketDistances[position].first->bucket;
            std::vector<std::pair<Eigen::Vector3f, float>> nClosest = collectCloserThan(bucket, p, maxDistance, k);

            for (const std::pair<Eigen::Vector3f, float>& item : nClosest)
            {
                float distance = euclideanDistance(p, item.first);
                std::pair<float, std::pair<Eigen::Vector3f, float>> maxValue = queue.top();
                if (distance < maxValue.first)
                {
                    queue.pop();
                    queue.push(std::pair<float, std::pair<Eigen::Vector3f, float>>(distance, item));
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
    KdTreeNodePair *build_tree_linear_median_search(std::vector<std::pair<Eigen::Vector3f, float>> *pts, int depth)
    {
        if (pts->size() == 0)
        {
            return NULL;
        }
        int axis = depth % 3;

        std::vector<float> d;
        for(std::pair<Eigen::Vector3f, float> const &p : *pts) {
            d.push_back(p.first[axis]);
        }

        float median = median_search::search(d, pts->size() / 2);
        std::vector<std::pair<Eigen::Vector3f, float>> left, right;
        KdTreeNodePair *leftChild = nullptr;
        KdTreeNodePair *rightChild = nullptr;
        std::vector<std::pair<Eigen::Vector3f, float>> bucket;
        if (pts->size() > kMaxBucketSize)
        {
            for (std::size_t i = 0; i < pts->size(); ++i)
            {
                if (pts->at(i).first[axis] < median)
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
        KdTreeNodePair *retVal = new KdTreeNodePair(median, leftChild, rightChild, bucket, depth);
        if (leftChild != nullptr)
            leftChild->parent = retVal;
        if (rightChild != nullptr)
            rightChild->parent = retVal;
        return retVal;
    }

    [[nodiscard]] virtual std::vector<std::pair<Eigen::Vector3f, float>> collectNClosest(const std::vector<std::pair<Eigen::Vector3f, float>> &list, const Eigen::Vector3f &p, size_t n) const
    {
        if (n >= list.size())
        {
            return list;
        }
        std::vector<std::pair<Eigen::Vector3f, float>> result;
        std::vector<std::pair<int, float>> distances; // index+distance
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance(p, list[i].first);
            distances.emplace_back(i, distance);
        }
        std::sort(distances.begin(), distances.end(), [=](std::pair<int, float> &a, std::pair<int, float> &b)
        { return a.second < b.second; });
        for (size_t i = 0; i < n; i++)
        {
            result.push_back(list[distances[i].first]);
        }
        return result;
    }

    [[nodiscard]] virtual std::vector<std::pair<Eigen::Vector3f, float>> collectCloserThan(const std::vector<std::pair<Eigen::Vector3f, float>> &list, const Eigen::Vector3f &p, float maxDistance, size_t maxPoints) const
    {
        std::vector<std::pair<Eigen::Vector3f, float>> result;
        std::vector<std::pair<int, float>> distances;
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance(p, list[i].first);
            if (distance <= maxDistance)
            {
                distances.emplace_back(i, distance);
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



    void collectDistanceToBuckets(KdTreeNodePair *cursor, const Eigen::Vector3f &p, std::vector<std::pair<KdTreeNodePair *, float>> *leafDistances) const
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

    static void collectInRadiusKnn(std::vector<std::pair<Eigen::Vector3f, float>> *list, KdTreeNodePair *cursor, const Eigen::Vector3f &p, float radius, int axis)
    {
        if (cursor != nullptr)
        {
            if (cursor->bucket.size() > 0)
            {
                for (size_t i = 0; i < cursor->bucket.size(); i++)
                {
                    float distance = euclideanDistance(p, cursor->bucket[i].first);
                    if (distance <= radius)
                    {
                        list->push_back(cursor->bucket[i]);
                    }
                }
            }
            KdTreeNodePair *nonMatchingSide = nullptr;
            KdTreeNodePair *matchingSide = nullptr;
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

    std::vector<std::pair<Eigen::Vector3f, float>> m_points;
};
