#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>

const double EPSILON = 1e-5;

struct Values {
    double x, y, vx, vy;
};

std::vector<Values> readFile(const std::string& filename) {
    std::vector<Values> values;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        if (line.find("x:") == std::string::npos) {
            continue;
        }

        Values v;
        std::istringstream iss(line);
        std::string token;
        
        // Skip "x:"
        iss >> token;
        iss >> v.x;
        
        // Skip "y:"
        iss >> token;
        iss >> v.y;
        
        // Skip "vx:"
        iss >> token;
        iss >> v.vx;
        
        // Skip "vy:"
        iss >> token;
        iss >> v.vy;

        values.push_back(v);
    }

    return values;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <file1.txt> <file2.txt>" << std::endl;
        return 1;
    }

    std::string file1 = argv[1];
    std::string file2 = argv[2];

    std::vector<Values> values1 = readFile(file1);
    std::vector<Values> values2 = readFile(file2);

    for (size_t i = 0; i < values1.size(); ++i) {
        double diff_x = std::abs(values1[i].x - values2[i].x);
        double diff_y = std::abs(values1[i].y - values2[i].y);
        double diff_vx = std::abs(values1[i].vx - values2[i].vx);
        double diff_vy = std::abs(values1[i].vy - values2[i].vy);

        if (diff_x > EPSILON || diff_y > EPSILON || diff_vx > EPSILON || diff_vy > EPSILON) {
            std::cout << "Test failed: Differences exceed epsilon (" << EPSILON << ")" << std::endl;
            std::cout << "Line " << i + 1 << " differences:" << std::endl;
            std::cout << "x: " << diff_x << std::endl;
            std::cout << "y: " << diff_y << std::endl;
            std::cout << "vx: " << diff_vx << std::endl;
            std::cout << "vy: " << diff_vy << std::endl;
            return 1;
        }
    }

    return 0;
} 