#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glew.h>
#include "Voxel.hpp"

typedef glm::vec<2, double> vec2;
/**
 * Voxel class that contains voxel properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this voxel.
 * @author Zugasti
 */
class Voxel;
class VoxelLink{
public:
    Voxel* v1;
    Voxel* v2;
    double threshold;
    VoxelLink(Voxel* v1, Voxel* v2, double threshold);
    void display(glm::mat2 rotation, vec2 cmPosition, glm::vec3 color, double radius);
};
