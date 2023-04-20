#pragma once
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"

#include "Particle.hpp"
#include "Spring.hpp"
#include "BendingSpring.hpp"
#include "RobustCCD.hpp"
#include "Voxel.hpp"
#include "CollisonHandler.hpp"
#include "Explosion.hpp"

typedef glm::vec<2, double> vec2;

using namespace std;

/**
 * Implementation of a simple particle system
 * @author kry & Mercier-Aubin
 */
class VoxelSystem  {
    
public:
    string name = "";
    std::vector<Shape*> shapes;
    CollisionHandler collisionHandler;
    std::vector<Explosion*> explosions;


    /**
     * Creates an empty particle system
     */
    VoxelSystem(){
    }

    ~VoxelSystem() {
        for (int i = 1; i < shapes.size(); ++i) {
            delete shapes[i];
        }
    }
    
    /**
     * Deletes all voxels
     */
    void clearVoxels() {
        for (Shape* s : shapes) { delete s; }
        shapes.clear();
    }
    
    /** Elapsed simulation time */
    double time = 0;
        
    /**
     *  Compute nodal forces using gravity and spring contributions
     */
    void computeForce() {

        // reset forces and set gravity
        for (Shape* s : shapes) {
            if (s->pinned)
                continue;
            vec2 viscousDampingForce = -viscousDamping * s->mass * s->v;
            s->f = vec2(0, 1) * s->mass * (useGravity ? gravity : 0)  + viscousDampingForce;
            for (Explosion* e : explosions)
            {
                vec2 direction = glm::normalize(s->p - e->pos);
                s->f += e->getPressureAt(s->p.x, s->p.y, time) * s->numberOfExternalVoxels * glm::normalize(s->p - e->pos);
            }
        }
    }

    void checkForFractures() {
        for (Shape* s : shapes) {
            glm::mat2 rotation;
            rotation[0][0] = glm::cos(s->angle);
            rotation[0][1] = -glm::sin(s->angle);
            rotation[1][0] = glm::sin(s->angle);
            rotation[1][1] = glm::cos(s->angle);
            for (int i = s->links.size() - 1; i >= 0; i--)
            {
                VoxelLink* l = s->links[i];

                vec2 pRelative = (l->v1->pRelative + l->v2->pRelative) / 2.0;
                vec2 p = rotation * pRelative;
                p += s->p;
                for (Explosion* e : explosions)
                {
                    double pressure = e->getPressureAt(p.x, p.y, time);
                    l->threshold -= pressure;
                }
            }
        }
        std::vector<Shape*> toAdd;
        std::vector<Shape*> toKeep;
        for (int i = 0; i < shapes.size(); i++)
        {
            if(!shapes[i]->ExploreVoxelLinks(shapes[i]->voxels[0], &toAdd))
            {
                toKeep.push_back(shapes[i]);
            }
        }
        shapes.clear();
        for (Shape *s : toKeep)
        {
            shapes.push_back(s);
        }
        for (Shape* s : toAdd)
        {
            shapes.push_back(s);
        }
    }
    
    void stepVelocities(float h) {
        for (Shape* s : shapes) {
            if (s->pinned)
                continue;
            s->angle += h * s->angularV;
            s->v += (h / s->mass) * s->f;
        }
    }

    /**
     * Update the positions with current velocities given the time step
     * @param h time step
     */
    void stepPositions(double h) {
        for (Shape* s : shapes) {
            if (s->pinned)
                continue;
            s->p += h * s->v;
            s->f = vec2(0, 0);
        }
    }

    void removeOldExplosions()
    {
        std::vector<Explosion*> surviveTest;
        for each (Explosion * e in explosions)
        {
            if (e->speed * (time - e->time) < 1000)
            {
                surviveTest.push_back(e);
            }
        }
        explosions.clear();
        for each (Explosion * e in surviveTest)
        {
            explosions.push_back(e);
        }
    }

    /** Time in seconds that was necessary to advance the system */
    float computeTime = 0;
    
    /**
     * Advances the state of the system
     * @param elapsed
     */
    void advanceTime( double elapsed ) {
        //cout << "Time advanced by: " << elapsed << endl;
        double now = glfwGetTime();
        
        removeOldExplosions();
        checkForFractures();
        computeForce();
        stepVelocities(elapsed);
        //collisionHandler.applyImpulses(elapsed, &shapes, width, height);
        if(useCollisions)
            collisionHandler.applyImpulsesWithAngular(elapsed, &shapes, width, height);
        stepPositions(elapsed);
        
        time += elapsed;
        computeTime = (glfwGetTime() - now);
    }
    

    Shape* createShape(float x, float y, float vx, float vy, float angle, float angularVelocity, double r, double g, double b, double voxelRadius, int shapeId) {
        Shape* v = NULL;
        if (spawnType == 0)
            v = new Shape(x, y, vx, vy, angle, angularVelocity, threshold, r, g, b, voxelRadius, shapeId);
        else 
            v = new Shape(x, y, vx, vy, angle, angularVelocity, threshold, r, g, b, voxelRadius, rectangleWidth, rectangleHeight);

        v->index = shapes.size();
        shapes.push_back(v);
        return v;
    }

    Explosion* createExplosion(float x, float y) {
        Explosion* e = new Explosion(x, y, time, speed, peak, T, b, decay);
        explosions.push_back(e);
        return e;
    }
    
    void remove( Shape* v ) {

    	shapes.erase( std::remove(shapes.begin(), shapes.end(), v ) );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < shapes.size(); i++ ) {
            shapes[i]->index = i;
    	}
    }
    
    void reset()
    {
        explosions.clear();
        shapes.clear();
    }

    void init() {
    }

    int height = 0;
    int width = 0;

    void display() {

        glClear(GL_COLOR_BUFFER_BIT);
        for each (Shape *s in shapes)
        {
            s->display(displayLinks);
        }
        for each (Explosion * e in explosions)
        {
            e->display(50, time);
        }
    }
    
    bool useGravity = false;
    bool useCollisions = true;
    bool displayLinks = true;
    double peak = 2500;
    double T = 0.1;
    double b = 12;
    double speed = 200;
    int decay = 2;
    int spawnType = 1;
    int rectangleWidth = 4;
    int rectangleHeight = 3;
    double threshold = 2000;
    double gravity = 9.8;
    double viscousDamping = 0;
};