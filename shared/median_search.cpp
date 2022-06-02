#include <random>

#include "typealiases.cpp"

namespace median_search
{
    static float search(std::vector<float> values, size_t position)
    {
        if (values.size() == 1)
        {
            return values[0];
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<unsigned long> distrib(0, values.size() - 1);

        size_t random_index = distrib(gen);
        std::vector<float> left;
        std::vector<float> right;
        for (size_t i = 0; i < values.size(); ++i)
        {
            if (i == random_index)
            {
                continue;
            }
            float p = values[i];
            if (p < values[random_index])
            {
                left.push_back(p);
            }
            else
            {
                right.push_back(p);
            }
        }
        if (left.size() == position)
        {
            return values[random_index];
        }
        if (left.size() >= position)
        {
            return search(left, position);
        }
        else
        {
            return search(right, position - left.size() - 1);
        }
    }
} // namespace median_search
