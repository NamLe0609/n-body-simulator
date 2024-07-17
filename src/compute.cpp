// TODO
// SOA with mass, velocity, position
// Calculate Forces
// Leapfrog integration scheme
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <random>

template <typename T>
struct Vec3 {
    T x, y, z; 

    Vec3 constexpr operator+(const Vec3 &operand) const {
        return {x + operand.x, y + operand.y, z + operand.z};
    }
    
    Vec3 constexpr operator-(const Vec3 &operand) const {
        return {x - operand.x, y - operand.y, z - operand.z};
    }
    
    Vec3 constexpr operator*(T factor) const {
        return {x * factor, y * factor, z * factor};
    }

    std::string print() {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
    }

    T constexpr norm_square() {
        return x*x + y*y + z*z;
    }
};


template <typename T>
void constexpr static resetForces(std::vector<Vec3<T>>& force) {
    for (Vec3<T> vector: force) {
        vector.x = 0;
        vector.y = 0;
        vector.z = 0; 
    }
}

template <typename T>
void constexpr static calculateForce(const std::vector<Vec3<T>>& position, const std::vector<T>& mass, std::vector<Vec3<T>>& force) {
    constexpr T gconst = 6.67e-11;
    constexpr T epsilonSquared = 1e-9;
    const std::size_t size = position.size();
    
    resetForces(force);

    // Last element has had forces from all other particles computed already
    // So we stop at size - 1
    for (std::size_t i = 0; i < size - 1; i++) {
        auto massi = mass[i];
        Vec3 posi = position[i];
        auto gm = gconst * massi;

        for (std::size_t j = i + 1; j < size; j++) {
            Vec3<T> rij = position[j] - posi;
            auto multiplier = gm * mass[j] / std::pow(rij.norm_square() + epsilonSquared, 1.5);
            Vec3 forceij = rij * multiplier;
            force[i] = force[i] + forceij;
            force[j] = force[j] - forceij;
        }
    }
}

template <typename T>
void constexpr static leapfrog(std::vector<Vec3<T>>& position, std::vector<Vec3<T>>& velocity, std::vector<Vec3<T>>& force, const std::vector<T>& mass) {
    constexpr T timestep = 1e6; 
    constexpr T timestepHalf = timestep / 2;
    const std::size_t size = position.size();

    // Velocity half step
    for (std::size_t i = 0; i < size; i++) {
       velocity[i] = velocity[i] + force[i] * timestepHalf * (1 / mass[i]); 
    }
    
    // Position next time step
    for (std::size_t i = 0; i < size; i++) {
        position[i] = position[i] + velocity[i] * timestep;
    }

    // Force at next time step
    calculateForce(position, mass, force);

    // Velocity next time step
    for (std::size_t i = 0; i < size; i++) {
       velocity[i] = velocity[i] + force[i] * timestepHalf * (1 / mass[i]); 
    }
}

template <typename T>
struct ParticleSystem {
    std::vector<Vec3<T>> position;
    std::vector<Vec3<T>> velocity;
    std::vector<T> mass;

    // Holds force on the particle to exploit symmetry 
    std::vector<Vec3<T>> force;

    void initialize(int numOfPart) {
        std::random_device rand;
        std::mt19937 rng(rand());
        std::uniform_real_distribution<T> weightGen(6.0e23, 6.0e25); 
        std::uniform_real_distribution<T> vectorGen(1.0, 1.0e2);

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
        
        force.reserve(static_cast<typename std::vector<Vec3<T>>::size_type>(numOfPart));
        for (int i = 0; i < numOfPart; i++) {
           force.push_back(Vec3<T> {0, 0, 0});            
        }
        calculateForce(position, mass, force);
    }
};

int main() {
    ParticleSystem<double> psys;
    psys.initialize(2);

    for (auto vec: psys.velocity) {
        std::cout << vec.print() << "\n";
    }
    for (auto vec: psys.position) {
        std::cout << vec.print() << "\n";
    }

    for (int i = 0; i < 3; i ++) {
        leapfrog(psys.position, psys.velocity, psys.force, psys.mass); 
        for (auto vec: psys.velocity) {
            std::cout << "vel: " << vec.print() << "\n";
        }
        for (auto vec: psys.position) {
            std::cout << "pos: " << vec.print() << "\n";
        }
    }
    return 0;
}