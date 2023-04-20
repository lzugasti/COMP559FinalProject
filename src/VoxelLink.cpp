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
VoxelLink::VoxelLink(Voxel* v1, Voxel* v2, double threshold)
    {
        this->v1 = v1;
        this->v2 = v2;
        this->threshold = threshold;
    }

void VoxelLink::display(glm::mat2 rotation, vec2 cmPosition, glm::vec3 color, double radius) {
        glColor3f(color.r, color.g, color.b);
        glm::mat4x2 corners;
        vec2 pRelative = (v1->pRelative + v2->pRelative) / 2.0;
        corners[0][0] = pRelative.x - radius;
        corners[1][0] = pRelative.x + radius;
        corners[2][0] = pRelative.x + radius;
        corners[3][0] = pRelative.x - radius;
        corners[0][1] = pRelative.y - radius;
        corners[1][1] = pRelative.y - radius;
        corners[2][1] = pRelative.y + radius;
        corners[3][1] = pRelative.y + radius;

        vec2 p = cmPosition;//pRelative + cmPosition;
        //cout << "Displaying voxel at: " << pRelative.x << ", " << pRelative.y << endl;

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
