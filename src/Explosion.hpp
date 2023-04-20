#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#define GLM_FORCE_RADIANS
#define AMBIENT_PRESSURE 0;
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
class Explosion{
public:
    vec2 pos;
    double time;
    double speed; // speed in pixels per second;
    double peak;
    double T; // determines when the positive pressure ends (and when the negative phase starts);
    double b; 
    int decay;

    Explosion(double x, double y, double time, double speed, double peak, double T, double b, int decay) {
        this->pos = vec2(x, y);
        this->time = time;
        this->speed = speed;
        this->peak = peak;
        this->T = T;
        this->b = b;
        this->decay = decay;
    }

    double getPressureAt(double x, double y, double current_time)
    {
        vec2 target = vec2(x, y);
        double distance = glm::length(target - pos);
        double time_to_reach = distance / speed;
        double time_since_explosion = current_time - time;
        // otherwise no explosion
        if (time_to_reach <= time_since_explosion)
        {
            double t = time_since_explosion - time_to_reach;
            double e = exp(-b*t / T);
            double p = (1 - t / T);
            double distance_m = distance / 100;
            double decayAmount = 1;
            for (int i = 0; i < decay; i++)
            {
                decayAmount *= distance_m;
            }
            double shock = peak / (decayAmount);
            double pressure = shock * p * e;
            return pressure;
        }
        return 0;
    }

    void display(int numberOfSegments, double current_time)
    {
        double decayAmount = 1;
        for (int i = 0; i < decay; i++)
        {
            decayAmount *= (current_time - time);
        }
        glColor3f(1.0 * (1/ (decayAmount)), 0.0, 0.5 * (1 / (decayAmount)));
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < numberOfSegments; i++) {
            double theta = 2.0 * 3.1415926 * (double)i / ((double)numberOfSegments);
            double r = (current_time - time) * speed;
            double x = r * cosf(theta);
            double y = r * sinf(theta);
            glVertex2f(x + pos.x, y + pos.y);
        }
        glEnd();
    }
};
