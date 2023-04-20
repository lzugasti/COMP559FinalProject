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
#include "VoxelLink.hpp"

typedef glm::vec<2, double> vec2;
/**
 * Voxel class that contains voxel properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this voxel.
 * @author Zugasti
 */
class VoxelLink;
class Voxel{
public:
    /** Identifies this particles position in the particle list */
    int index;

    glm::vec3 color{ 0.0f, 0.95f, 0.0f };

    double mass = 1;
    double radius = 25;

    vec2 pRelative;
    vec2 p0Relative;

    VoxelLink* top;
    VoxelLink* right;
    VoxelLink* bottom;
    VoxelLink* left;

    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Voxel(double x, double y, double radius);

    /**
     * Resets the position of this particle
     */
    void reset();

    /**
    * Display the current voxel on screen
    */
    void display(glm::mat2 rotation, vec2 cmPosition, glm::vec3 color);
};
