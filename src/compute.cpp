// TODO
// SOA with mass, velocity, position
// Calculate Forces
// Leapfrog integration scheme
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

template <typename T>
struct Vec3 {
    T x, y, z; 

    Vec3 operator+(const Vec3 &operand) {
        x += operand.x;
        y += operand.y;
        z += operand.z;
    }
    
    Vec3 operator*(int factor) {
        x *= factor;
        y *= factor;
        z *= factor;
    }

    T norm_square() {
        return pow(x, 2) + pow(y, 2) + pow(z, 2);
    }
};

template <typename T>
struct ParticleSystem {
    std::vector<Vec3<T>> position;
    std::vector<Vec3<T>> velocity;
    std::vector<T> mass;

    void initialize(int numOfPart) {
        std::random_device rand;
        std::mt19937 rng(rand());
        std::uniform_real_distribution<T> weightGen(6.0, 600.0); // Weights in 10^23 kg
        std::uniform_real_distribution<T> vectorGen(1.0, 1000.0);

        position.reserve(static_cast<typename std::vector<Vec3<T>>::size_type>(numOfPart));
        for (int i = 0; i < numOfPart; i++) {
           position.push_back(Vec3<T> {vectorGen(rng), vectorGen(rng), vectorGen(rng)});            
        }

        velocity.reserve(static_cast<typename std::vector<Vec3<T>>::size_type>(numOfPart));
        for (int i = 0; i < numOfPart; i++) {
           velocity.push_back(Vec3<T> {vectorGen(rng), vectorGen(rng), vectorGen(rng)});            
        }

        mass.reserve(static_cast<typename std::vector<T>::size_type>(numOfPart));
        for (int i = 0; i < numOfPart; i++) {
           mass.push_back(weightGen(rng));            
        }
    }
};

int main() {
    ParticleSystem<double> psys;
    psys.initialize(100);
    for (auto i: psys.velocity) {
        std::cout << i.x << ", " << i.y << ", " << i.z << "\n";
    } 
    return 0;
}