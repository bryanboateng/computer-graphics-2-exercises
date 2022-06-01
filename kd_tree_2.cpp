#include <cmath>
#include <queue>

#include "typealiases.cpp"


static float euclideanDistance2(Point const &p1, Point const &p2)
{
    float dx = p1[0] - p2[0];
    float dy = p1[1] - p2[1];
    return std::sqrt(dx * dx + dy * dy);
}

class KdTreeNode2
{
public:
    float median;
    KdTreeNode2 *left;
    KdTreeNode2 *right;
    KdTreeNode2 *parent = nullptr;
    PointList bucket;
    int depth;
    KdTreeNode2(float p, KdTreeNode2 *l, KdTreeNode2 *r, PointList b, int d)
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
    KdTree2(PointList const &points)
            : m_points(points)
    {
        root = build_tree_linear_median_search(&m_points, 0);
    }

    virtual ~KdTree2() = default;

    PointList const &getPoints() const
    {
        return m_points;
    }

    virtual std::vector<Point> collectInRadius(const Point &p, float radius) const
    {
        PointList list;
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
    virtual PointList collectKNearest(const Point &p, unsigned int k) const
    {
        PointList result;
        if (k == 0)
        {
            return result;
        }
        KdTreeNode2 *cursor = root;
        bool found = false;
        std::vector<std::pair<KdTreeNode2 *, float>> bucketDistances;
        collectDistanceToBuckets(root, p, &bucketDistances);
        std::sort(bucketDistances.begin(), bucketDistances.end(), [=](std::pair<KdTreeNode2 *, float> &a, std::pair<KdTreeNode2 *, float> &b)
        { return a.second < b.second; });
        std::priority_queue<std::pair<float, Point>> queue;
        size_t position = 0;
        while (queue.size() < k && position < bucketDistances.size())
        {
            PointList bucket = bucketDistances[position].first->bucket;
            PointList nClosest = collectNClosest(bucket, p, k - queue.size());
            for (Point entry : nClosest)
            {
                float distance = euclideanDistance2(p, entry);
                queue.push(std::pair<float, Point>(distance, entry));
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
            PointList bucket = bucketDistances[position].first->bucket;
            PointList nClosest = collectCloserThan(bucket, p, maxDistance, k);

            for (auto item : nClosest)
            {
                float distance = euclideanDistance2(p, item);
                auto maxValue = queue.top();
                if (distance < maxValue.first)
                {
                    queue.pop();
                    queue.push(std::pair<float, Point>(distance, item));
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
    KdTreeNode2 *build_tree_linear_median_search(std::vector<Point> *pts, int depth)
    {
        if (pts->size() == 0)
        {
            return NULL;
        }
        int axis = depth % 2;
        float median = median_search::search(*pts, axis, pts->size() / 2);
        std::vector<Point> left, right;
        KdTreeNode2 *leftChild = nullptr;
        KdTreeNode2 *rightChild = nullptr;
        PointList bucket;
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
        KdTreeNode2 *retVal = new KdTreeNode2(median, leftChild, rightChild, bucket, depth);
        if (leftChild != nullptr)
            leftChild->parent = retVal;
        if (rightChild != nullptr)
            rightChild->parent = retVal;
        return retVal;
    }

    virtual PointList collectNClosest(const PointList &list, const Point &p, size_t n) const
    {
        if (n >= list.size())
        {
            return list;
        }
        PointList result;
        std::vector<std::pair<int, float>> distances; // index+distance
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance2(p, list[i]);
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

    virtual PointList collectCloserThan(const PointList &list, const Point &p, float maxDistance, size_t maxPoints) const
    {
        PointList result;
        std::vector<std::pair<int, float>> distances;
        for (size_t i = 0; i < list.size(); i++)
        {
            float distance = euclideanDistance2(p, list[i]);
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

    void collectDistanceToBuckets(KdTreeNode2 *cursor, const Point &p, std::vector<std::pair<KdTreeNode2 *, float>> *leafDistances) const
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
                    int splitAxis = (cursor->parent->depth % 2);
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
                        float distance = euclideanDistance2(p, Point{x, y, z});
                    }
                    leafDistances->push_back(std::pair(cursor, distance));
                }
            }
            collectDistanceToBuckets(cursor->left, p, leafDistances);
            collectDistanceToBuckets(cursor->right, p, leafDistances);
        }
    }

    static void collectInRadiusKnn(PointList *list, KdTreeNode2 *cursor, const Point &p, float radius, int axis)
    {
        if (cursor != nullptr)
        {
            if (cursor->bucket.size() > 0)
            {
                for (size_t i = 0; i < cursor->bucket.size(); i++)
                {
                    float distance = euclideanDistance2(p, cursor->bucket[i]);
                    if (distance <= radius)
                    {
                        list->push_back(cursor->bucket[i]);
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
                float z = (axis == 2) ? cursor->median : p[2];
                float distance = euclideanDistance2(p, Point{x, y, z});
                bool compareValue = (nonMatchingSide == cursor->left) ? distance <= radius : distance < radius;
                if (compareValue)
                {
                    collectInRadiusKnn(list, nonMatchingSide, p, radius, (axis + 1) % 2);
                }
            }
        }
        return;
    }

    PointList m_points;
};
