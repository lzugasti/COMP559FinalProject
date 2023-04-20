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
#include "VoxelLink.hpp"

typedef glm::vec<2, double> vec2;
/**
 * Voxel class that contains voxel properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this voxel.
 * @author Zugasti
 */
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Voxel::Voxel( double x, double y, double radius ) {
        this->p0Relative = vec2(x, y);
        this->radius = radius;
        this->mass = 1;
        this->top = NULL;
        this->right = NULL;
        this->bottom = NULL;
        this->left = NULL;
        reset();
    }

    /**
     * Resets the position of this particle
     */
    void Voxel::reset() {
        pRelative = p0Relative;
    }

    /**
    * Display the current voxel on screen
    */
    void Voxel::display(glm::mat2 rotation, vec2 cmPosition, glm::vec3 color) {
        glColor3f(color.r, color.g, color.b);
        glm::mat4x2 corners;
        corners[0][0] = pRelative.x - radius;
        corners[1][0] = pRelative.x + radius;
        corners[2][0] = pRelative.x + radius;
        corners[3][0] = pRelative.x - radius;
        corners[0][1] = pRelative.y - radius;
        corners[1][1] = pRelative.y - radius;
        corners[2][1] = pRelative.y + radius;
        corners[3][1] = pRelative.y + radius;

        vec2 p = cmPosition;

        corners = rotation * corners;
        corners[0][0] += p.x;
        corners[1][0] += p.x;
        corners[2][0] += p.x;
        corners[3][0] += p.x;
        corners[0][1] += p.y;
        corners[1][1] += p.y;
        corners[2][1] += p.y;
        corners[3][1] += p.y;
        glBegin(GL_POLYGON);
        glVertex2f(corners[0][0], corners[0][1]);
        glVertex2f(corners[1][0], corners[1][1]);
        glVertex2f(corners[2][0], corners[2][1]);
        glVertex2f(corners[3][0], corners[3][1]);
        glEnd();
    }
