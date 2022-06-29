#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "typealiases.cpp"

std::tuple<std::vector<Eigen::Vector3f>, std::vector<Eigen::Vector3f>, Eigen::Vector3f, Eigen::Vector3f> parseOff(std::string const &filename)
{
    std::ifstream off_file(filename);
    if (off_file.is_open())
    {
        std::string header_keyword;
        if (!(off_file >> header_keyword))
        {
            off_file.close();
            throw std::invalid_argument("Unable to read file.");
        }
        else if (header_keyword == "OFF")
        {
            int vertex_count, face_count, edge_count;
            off_file >> vertex_count >> face_count >> edge_count;
            std::vector<Eigen::Vector3f> points;
            std::vector<Eigen::Vector3f> normals;
            float max_x, max_y, max_z = -10000;
            float min_x, min_y, min_z = 10000;
            for (int i = 0; i < vertex_count; ++i)
            {
                float x, y, z;
                off_file >> x >> y >> z;
                if (x > max_x)
                    max_x = x;
                if (y > max_y)
                    max_y = y;
                if (z > max_z)
                    max_z = z;
                if (x < min_x)
                    min_x = x;
                if (y < min_y)
                    min_y = y;
                if (z < min_z)
                    min_z = z;

                Eigen::Vector3f point;
                point << x, y, z;
                points.push_back(point);
            }
            off_file.close();
            return {points, normals, {min_x, min_y, min_z}, {max_x, max_y, max_z}};
        }
        else if (header_keyword == "NOFF")
        {
            int vertex_count, face_count, edge_count;
            off_file >> vertex_count >> face_count >> edge_count;
            std::vector<Eigen::Vector3f> points;
            std::vector<Eigen::Vector3f> normals;
            float max_x, max_y, max_z = -10000;
            float min_x, min_y, min_z = 10000;
            for (int i = 0; i < vertex_count; ++i)
            {
                float x, y, z;
                float nx, ny, nz;
                off_file >> x >> y >> z >> nx >> ny >> nz;
                if (x > max_x)
                    max_x = x;
                if (y > max_y)
                    max_y = y;
                if (z > max_z)
                    max_z = z;
                if (x < min_x)
                    min_x = x;
                if (y < min_y)
                    min_y = y;
                if (z < min_z)
                    min_z = z;

                Eigen::Vector3f point;
                Eigen::Vector3f normal;
                point << x, y, z;
                normal << nx, ny, nz;
                points.push_back(point);
                normals.push_back(normal);
            }
            off_file.close();
            return {points, normals, {min_x, min_y, min_z}, {max_x, max_y, max_z}};
        }
        else
        {
            off_file.close();
            throw std::invalid_argument("Incorrect file format. This program only supports point data from OFF-files that contain the header keyword OFF. More about the specifications of OFF-Files can be found at: http://paulbourke.net/dataformats/oogl/#OFF");
        }
    }
    else
    {
        throw std::invalid_argument("Unable to read file.");
    }
}