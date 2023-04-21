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
class Shape{
public:
    /** Identifies this particles position in the particle list */
    int index;

    bool pinned = false;

    std::vector<Voxel*> voxels;
    std::vector<VoxelLink*> links;
    double voxelRadius = 25;
    int numberOfExternalVoxels;

    glm::vec3 color{ 0.0f, 0.95f, 0.0f };

    double mass;
    double inertia;

    vec2 p;
    double angle;
    vec2 v;
    double angularV;
    vec2 vTemp;

    vec2 p0;
    double angle0;
    vec2 v0;
    double angularV0;
    vec2 f;

    double threshold;

    /** Default constructor */
    void createShapePrefabricated(int shapeId, int size)
    {
        if (size == 1)
        {
            switch (shapeId)
            {
            case 1:
                this->voxels.push_back(new Voxel(0, 0, voxelRadius));
                numberOfExternalVoxels = 1;
                break;

            case 2: {
                Voxel* left = new Voxel(0, 0, voxelRadius);
                Voxel* right = new Voxel(voxelRadius * 2, 0, voxelRadius);
                this->voxels.push_back(left);
                this->voxels.push_back(right);
                left->right = new VoxelLink(left, right, threshold);
                this->links.push_back(left->right);
                right->left = left->right;
                numberOfExternalVoxels = 2;

                break;
            }
            case 3: {

                Voxel* topLeft = new Voxel(0, 0, voxelRadius);
                Voxel* topRight = new Voxel(voxelRadius * 2, 0, voxelRadius);
                Voxel* bottomRight = new Voxel(voxelRadius * 2, voxelRadius * 2, voxelRadius);
                this->voxels.push_back(topLeft);
                this->voxels.push_back(topRight);
                this->voxels.push_back(bottomRight);
                topLeft->right = new VoxelLink(topLeft, topRight, threshold);
                topRight->left = topLeft->right;
                bottomRight->top = new VoxelLink(bottomRight, topRight, threshold);
                topRight->bottom = bottomRight->top;
                this->links.push_back(topLeft->right);
                this->links.push_back(bottomRight->top);
                numberOfExternalVoxels = 3;
                break;
            }
            case 4:
            {
                Voxel* middle = new Voxel(0, 0, voxelRadius);
                Voxel* middleRight = new Voxel(voxelRadius * 2, 0, voxelRadius);
                Voxel* bottomRight = new Voxel(voxelRadius * 2, voxelRadius * 2, voxelRadius);
                Voxel* middleLeft = new Voxel(-voxelRadius * 2, 0, voxelRadius);
                Voxel* topLeft = new Voxel(-voxelRadius * 2, -voxelRadius * 2, voxelRadius);
                this->voxels.push_back(middle);
                this->voxels.push_back(middleRight);
                this->voxels.push_back(bottomRight);
                this->voxels.push_back(middleLeft);
                this->voxels.push_back(topLeft);
                middle->right = new VoxelLink(middle, middleRight, threshold);
                this->links.push_back(middle->right);
                middleRight->left = middle->right;
                middle->left = new VoxelLink(middle, middleLeft, threshold);
                this->links.push_back(middle->left);
                middleLeft->right = middle->left;
                bottomRight->top = new VoxelLink(bottomRight, middleRight, threshold);
                this->links.push_back(bottomRight->top);
                middleRight->bottom = bottomRight->top;
                middleLeft->top = new VoxelLink(middleLeft, topLeft, threshold);
                this->links.push_back(middleLeft->top);
                topLeft->bottom = middleLeft->top;
                numberOfExternalVoxels = 5;
                break;
            }
            case 5:
            {
                Voxel* middle = new Voxel(0, 0, voxelRadius);
                Voxel* right = new Voxel(voxelRadius * 2, 0, voxelRadius);
                Voxel* left = new Voxel(-voxelRadius * 2, 0, voxelRadius);
                Voxel* top = new Voxel(0, -voxelRadius * 2, voxelRadius);
                Voxel* bottom = new Voxel(0, voxelRadius * 2, voxelRadius);
                this->voxels.push_back(middle);
                this->voxels.push_back(right);
                this->voxels.push_back(left);
                this->voxels.push_back(top);
                this->voxels.push_back(bottom);
                middle->bottom = new VoxelLink(middle, bottom, threshold);
                this->links.push_back(middle->bottom);
                bottom->top = middle->bottom;
                middle->top = new VoxelLink(middle, top, threshold);
                this->links.push_back(middle->top);
                top->bottom = middle->top;
                middle->right = new VoxelLink(middle, right, threshold);
                this->links.push_back(middle->right);
                right->left = middle->right;
                middle->left = new VoxelLink(middle, left, threshold);
                this->links.push_back(middle->left);
                left->right = middle->left;
                numberOfExternalVoxels = 4;
                break;
            }
            default:
                break;
            }
        }
        
    }

    void createShapeRectangle(int height, int width)
    {
        std::vector<std::vector<Voxel*>> rectangle;
        for (int y = 0; y < height; y++)
        {
            std::vector<Voxel*> row;
            for (int x = 0; x < width; x++)
            {
                Voxel* v = new Voxel(x * voxelRadius * 2, y * voxelRadius * 2, voxelRadius);
                this->voxels.push_back(v);
                row.push_back(v);
            }
            rectangle.push_back(row);
        }
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (y > 0)
                {
                    rectangle[y][x]->top = new VoxelLink(rectangle[y][x], rectangle[y - 1][x], threshold);
                    rectangle[y - 1][x]->bottom = rectangle[y][x]->top;
                    this->links.push_back(rectangle[y][x]->top);
                }
                if (x > 0)
                {
                    rectangle[y][x]->left = new VoxelLink(rectangle[y][x], rectangle[y][x - 1], threshold);
                    rectangle[y][x - 1]->right = rectangle[y][x]->left;
                    this->links.push_back(rectangle[y][x]->left);
                }
            }
        }
        this->numberOfExternalVoxels = width * 2 + height * 2 - 2;
    }

    bool ExploreVoxelLinks(Voxel* start, std::vector<Shape*>* newShapes)
    {
        std::vector<Voxel*> allVisited;
        std::vector<Voxel*> visited;
        ExploreVoxelLinksHelper(start, &visited, &allVisited);
        if (allVisited.size() != this->voxels.size())
        {
            allVisited.clear();
            do {
                visited.clear();
                // find another node
                for (Voxel* v : this->voxels)
                {
                    if (std::find(allVisited.begin(), allVisited.end(), v) == allVisited.end())
                    {
                        start = v;
                        break;
                    }
                }
                ExploreVoxelLinksHelper(start, &visited, &allVisited);

                glm::mat2 rotation;
                rotation[0][0] = glm::cos(this->angle);
                rotation[0][1] = -glm::sin(this->angle);
                rotation[1][0] = glm::sin(this->angle);
                rotation[1][1] = glm::cos(this->angle);
                vec2 cm(0.0, 0.0);
                for each (Voxel * v in visited)
                {
                    cm += v->p0Relative;
                }
                double newShapeMass = 0;
                for each (Voxel * v in visited)
                {
                    newShapeMass += v->mass;
                }

                cm = rotation * (cm / newShapeMass);
                Shape* newShape = new Shape(p.x + cm.x, p.y + cm.y, v.x, v.y, angle, angularV, threshold, color.r, color.g, color.b, voxelRadius);
                for (Voxel* v : visited)
                {
                    newShape->voxels.push_back(v);
                }
                for (VoxelLink* l : links)
                {
                    if (std::find(visited.begin(), visited.end(), l->v1) != visited.end()
                        &&
                        std::find(visited.begin(), visited.end(), l->v2) != visited.end())
                    {
                        newShape->links.push_back(l);
                    }
                    else if (std::find(visited.begin(), visited.end(), l->v1) != visited.end()
                        ||
                        std::find(visited.begin(), visited.end(), l->v2) != visited.end())
                    {
                        if (l->v1->bottom == l)
                        {
                            l->v1->bottom = NULL;
                        }
                        if (l->v1->top == l)
                        {
                            l->v1->top = NULL;
                        }
                        if (l->v1->right == l)
                        {
                            l->v1->right = NULL;
                        }
                        if (l->v1->left == l)
                        {
                            l->v1->left = NULL;
                        }
                        if (l->v2->bottom == l)
                        {
                            l->v2->bottom = NULL;
                        }
                        if (l->v2->top == l)
                        {
                            l->v2->top = NULL;
                        }
                        if (l->v2->right == l)
                        {
                            l->v2->right = NULL;
                        }
                        if (l->v2->left == l)
                        {
                            l->v2->left = NULL;
                        }
                    }
                }
                newShape->mass = newShape->getMass();
                newShape->adjustVoxelsToCenterOfMass();
                newShape->inertia = newShape->getInertia();
                newShape->numberOfExternalVoxels = visited.size();
                newShapes->push_back(newShape);
            } while (allVisited.size() != this->voxels.size());
            return true;
        }
        else
        {
            return false;
        }
    }

    void ExploreVoxelLinksHelper(Voxel* current, std::vector<Voxel*>* visited, std::vector<Voxel*>* allVisited)
    {
        if (std::find(allVisited->begin(), allVisited->end(), current) == allVisited->end())
        {
            visited->push_back(current);
            allVisited->push_back(current);
            if (current->left != NULL && current->left->threshold > 0)
            {
                if (current->left->v1 != current)
                    ExploreVoxelLinksHelper(current->left->v1, visited, allVisited);
                else
                    ExploreVoxelLinksHelper(current->left->v2, visited, allVisited);
            }
            if (current->right != NULL && current->right->threshold > 0)
            {
                if (current->right->v1 != current)
                    ExploreVoxelLinksHelper(current->right->v1, visited, allVisited);
                else
                    ExploreVoxelLinksHelper(current->right->v2, visited, allVisited);
            }
            if (current->top != NULL && current->top->threshold > 0)
            {
                if (current->top->v1 != current)
                    ExploreVoxelLinksHelper(current->top->v1, visited, allVisited);
                else
                    ExploreVoxelLinksHelper(current->top->v2, visited, allVisited);
            }
            if (current->bottom != NULL && current->bottom->threshold > 0)
            {
                if (current->bottom->v1 != current)
                    ExploreVoxelLinksHelper(current->bottom->v1, visited, allVisited);
                else
                    ExploreVoxelLinksHelper(current->bottom->v2, visited, allVisited);
            }
        }
    }
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Shape( double x, double y, double vx, double vy, double angle, double angularV, double threshold, double r, double g, double b, double voxelRadius, int shapeId) : vTemp(0, 0), v(0, 0) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);
        angle0 = angle;
        angularV0 = angularV;
        this->threshold = threshold;
        reset();
        this->angle = angle;
        this->angularV = angularV;
        this->color.r = r;
        this->color.g = g;
        this->color.b = b;
        this->voxelRadius = voxelRadius;
        createShapePrefabricated(shapeId, 1);
        this->mass = getMass();
        adjustVoxelsToCenterOfMass();
        this->inertia = getInertia();
     }

    Shape(double x, double y, double vx, double vy, double angle, double angularV, double threshold, double r, double g, double b, double voxelRadius, int width, int height) : vTemp(0, 0), v(0, 0) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);
        angle0 = angle;
        angularV0 = angularV;
        this->threshold = threshold;
        reset();
        this->angle = angle;
        this->angularV = angularV;
        this->color.r = r;
        this->color.g = g;
        this->color.b = b;
        this->voxelRadius = voxelRadius;
        createShapeRectangle(height, width);
        this->mass = getMass();
        adjustVoxelsToCenterOfMass();
        this->inertia = getInertia();
    }

    Shape(double x, double y, double vx, double vy, double angle, double angularV, double threshold, double r, double g, double b, double voxelRadius) : vTemp(0, 0), v(0, 0) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);
        angle0 = angle;
        angularV0 = angularV;
        this->threshold = threshold;
        reset();
        this->angle = angle;
        this->angularV = angularV;
        this->color.r = r;
        this->color.g = g;
        this->color.b = b;
        this->voxelRadius = voxelRadius;
    }

    double getMass()
    {
        double m = 0;
        for each (Voxel* v in voxels)
        {
            m += v->mass;
        }
        return m;
    }

    double getInertia()
    {
        double i = 0;
        for each (Voxel * v in voxels)
        {
            double dSquared = glm::dot(v->pRelative, v->pRelative);
            //i += (v->mass / 12.0) * (2*(2*voxelRadius) * (2*voxelRadius)) + v->mass * dSquared;
            //i += (v->mass / 12.0) * (2*(2*voxelRadius) * (2*voxelRadius)) + v->mass * dSquared;
            //i += 0.5 * v->mass * voxelRadius * voxelRadius + v->mass * dSquared;
            i += v->mass * (voxelRadius * voxelRadius / 6 + dSquared);
            //i += 0.5 * v->mass * voxelRadius * voxelRadius + v->mass * dSquared;
            //i += 0.5 * v->mass * voxelRadius * voxelRadius;
            //i += 0.5 * v->mass * voxelRadius * voxelRadius *glm::max(1.0, glm::sqrt(glm::dot(v->pRelative, v->pRelative)));
        }
        return i;
    }

    void adjustVoxelsToCenterOfMass()
    {
        vec2 cm = vec2(0.0, 0.0);
        for each (Voxel * v in voxels)
        {
            cm += v->mass * v->p0Relative;
        }
        cm = cm / this->mass;
        
        for each (Voxel * v in voxels)
        {
            v->p0Relative = v->p0Relative - cm;
            v->pRelative = v->pRelative - cm;
        }
    }

    /**
     * Resets the position of this particle
     */
    void reset() {
        p = p0;
        v = v0;
        f = vec2(0, 0);
    }

    /**
     * Clears all forces acting on this particle
     */
    void clearForce() {
        f = vec2(0, 0);
    }

    /**
     * Adds the given force to this particle
     * @param force
     */
    void addForce(vec2 force) {
        f += force;
    }

    /**
    * Display the current voxel on screen
    */
    void display(bool displayLinks) {
        //cout << "Displaing shape: " << index << endl;
        glm::mat2 rotation;
        rotation[0][0] = glm::cos(angle);
        rotation[0][1] = -glm::sin(angle);
        rotation[1][0] = glm::sin(angle);
        rotation[1][1] = glm::cos(angle);

        for each (Voxel* v in voxels)
        {
            v->display(rotation, p, color);
        }
        if (displayLinks)
        {
            for each (VoxelLink * l in links)
            {
                if (l->threshold > 0)
                    l->display(rotation, p, glm::vec3(1.0, 0.0, 0.0), voxelRadius / 3.0);
            }
        }
    }

    /**
     * Computes the distance of another voxel to this voxel
     * @param x
     * @param y
     * @return the distance
     */
    bool checkCollisionWithWalls(double height, double width, double eps, double e)
    {
        glm::mat2 rotation;
        rotation[0][0] = glm::cos(angle);
        rotation[0][1] = -glm::sin(angle);
        rotation[1][0] = glm::sin(angle);
        rotation[1][1] = glm::cos(angle);

        for each (Voxel * vox in voxels)
        {
            vec2 vpos = this->p + (vec2)(rotation*vox->pRelative);
            if (vpos.y + voxelRadius >= height - eps)
            {
                vec2 n = glm::vec2(0.0, 1.0);

                // moving away, so no collision
                if (glm::dot(this->v, -n) > eps)
                    continue;

                double J = -(1 + e) * glm::dot(v, n) / (glm::dot(n, n) * (1 / mass));

                // Update velocities
                v = v + (J / mass) * n;
            }

            if (vpos.x + voxelRadius >= width - eps)
            {
                vec2 n = glm::vec2(1.0, 0.0);

                // moving away, so no collision
                if (glm::dot(this->v, -n) > eps)
                    continue;

                double J = -(1 + e) * glm::dot(v, n) / (glm::dot(n, n) * (1 / mass));

                // Update velocities
                v = v + (J / mass) * n;
            }


            if (vpos.x - voxelRadius <= eps)
            {
                vec2 n = glm::vec2(-1.0, 0.0);

                // moving away, so no collision
                if (glm::dot(this->v, -n) > eps)
                    continue;

                double J = -(1 + e) * glm::dot(v, n) / (glm::dot(n, n) * (1 / mass));

                // Update velocities
                v = v + (J / mass) * n;
            }
            
        }
        return false;
    }

    bool checkCollisionWithWallsWithAngular(double height, double width, double eps, double e)
    {
        bool hadCollision = false;
        glm::mat2 rotation;
        rotation[0][0] = glm::cos(angle);
        rotation[0][1] = -glm::sin(angle);
        rotation[1][0] = glm::sin(angle);
        rotation[1][1] = glm::cos(angle);

        std::vector<vec2> impulses;
        std::vector<vec2> rPerps;
        for each (Voxel * vox in voxels)
        {
            vec2 vpos = this->p + (vec2)(rotation * vox->pRelative);
            if (vpos.y + voxelRadius >= height - eps)
            {
                vec2 n = vec2(0.0, 1.0);
                vec2 r = this->p - vec2(vpos.x, vpos.y+voxelRadius);

                vec2 rPerp = vec2(-r.y, r.x);

                vec2 angularV = rPerp * this->angularV;
                
                vec2 relativeVelocity = (this->v + angularV);

                // moving away, so no collision
                if (glm::dot(relativeVelocity, n) < eps)
                {
                    continue;
                }
                
                double rdotn = glm::dot(rPerp, n);
                double J = -(1 + e) * glm::dot(relativeVelocity, n);
                double denom = (glm::dot(n, n) * (1 / (mass)) + (rdotn * rdotn / (inertia)));
                J /= denom;
                vec2 impulse = J * n;
                impulses.push_back(impulse);
                rPerps.push_back(rPerp);
                hadCollision = true;
                this->v = this->v + impulse / mass;
                this->angularV = this->angularV + glm::dot(rPerp, impulse) / inertia;
            }

            if (vpos.x + voxelRadius >= width - eps)
            {
                //vec2 n = glm::vec2(1.0, 0.0);

                //// moving away, so no collision
                //if (glm::dot(this->v, -n) > eps)
                //    continue;

                //double J = -(1 + e) * glm::dot(v, n) / (glm::dot(n, n) * (1 / mass));

                //// Update velocities
                //v = v + (J / mass) * n;

                vec2 n = glm::vec2(1.0, 0.0);
                vec2 r = this->p - vec2(vpos.x+voxelRadius, vpos.y);

                vec2 rPerp = vec2(-r.y, r.x);

                vec2 angularV = rPerp * this->angularV;
                
                vec2 relativeVelocity = (this->v + angularV);

                // moving away, so no collision
                if (glm::dot(relativeVelocity, n) < eps)
                {
                    continue;
                }
                
                double rdotn = glm::dot(rPerp, n);
                double J = -(1 + e) * glm::dot(relativeVelocity, n);
                double denom = (glm::dot(n, n) * (1 / (mass)) + (rdotn * rdotn / (inertia)));
                J /= denom;
                vec2 impulse = J * n;
                impulses.push_back(impulse);
                rPerps.push_back(rPerp);
                hadCollision = true;
                this->v = this->v + impulse / mass;
                this->angularV = this->angularV + glm::dot(rPerp, impulse) / inertia;
            }


            if (vpos.x - voxelRadius <= eps)
            {
                //vec2 n = glm::vec2(-1.0, 0.0);

                //// moving away, so no collision
                //if (glm::dot(this->v, -n) > eps)
                //    continue;

                //double J = -(1 + e) * glm::dot(v, n) / (glm::dot(n, n) * (1 / mass));

                //// Update velocities
                //v = v + (J / mass) * n;

                vec2 n = glm::vec2(-1.0, 0.0);
                vec2 r = this->p - vec2(vpos.x - voxelRadius, vpos.y);

                vec2 rPerp = vec2(-r.y, r.x);

                vec2 angularV = rPerp * this->angularV;

                vec2 relativeVelocity = (this->v + angularV);

                // moving away, so no collision
                if (glm::dot(relativeVelocity, n) < eps)
                {
                    continue;
                }

                double rdotn = glm::dot(rPerp, n);
                double J = -(1 + e) * glm::dot(relativeVelocity, n);
                double denom = (glm::dot(n, n) * (1 / (mass)) + (rdotn * rdotn / (inertia)));
                J /= denom;
                vec2 impulse = J * n;
                impulses.push_back(impulse);
                rPerps.push_back(rPerp);
                hadCollision = true;
                this->v = this->v + impulse / mass;
                this->angularV = this->angularV + glm::dot(rPerp, impulse) / inertia;
            }

        }

        /*for (int i = 0; i < impulses.size(); i++)
        {
            vec2 impulse = impulses[i] / (double)impulses.size();
            v = v + impulse / mass;
            angularV = angularV + glm::dot(rPerps[i], impulse) / inertia;
        }*/ 
        return hadCollision;
        
    }

    double getEnergy(double height)
    {
        return 0.5 * mass * glm::dot(v, v) +mass * 9.8 * (height - p.y) + 0.5 * inertia * angularV*angularV;
    }

    bool checkCollision(Shape *s, std::vector<vec2> *collisions, std::vector<vec2> *normals, std::vector<vec2>* relativeVelocities, double eps) {
        bool hadCollision = false;
        glm::mat2 rotation1;
        rotation1[0][0] = glm::cos(angle);
        rotation1[0][1] = -glm::sin(angle);
        rotation1[1][0] = glm::sin(angle);
        rotation1[1][1] = glm::cos(angle);

        glm::mat2 rotation2;
        rotation2[0][0] = glm::cos(s->angle);
        rotation2[0][1] = -glm::sin(s->angle);
        rotation2[1][0] = glm::sin(s->angle);
        rotation2[1][1] = glm::cos(s->angle);

        for each (Voxel *v1 in voxels)
        {
            for each (Voxel * v2 in s->voxels)
            {
                vec2 v1pos = this->p + (vec2)(rotation1*v1->pRelative);
                vec2 v2pos = s->p + (vec2)(rotation2 * v2->pRelative);
                vec2 diff = v1pos - v2pos;
                double distance = (float)sqrt(diff.x * diff.x + diff.y * diff.y) - 2 * voxelRadius;
                if (distance < eps)
                {
                    vec2 position = (v1pos + v2pos) / 2.0;
                    vec2 normal = v1pos - v2pos;
                    vec2 relativeVelocity = this->v - s->v;
                    collisions->push_back(position);
                    normals->push_back(normal);
                    relativeVelocities->push_back(relativeVelocity);
                    hadCollision = true;
                }
            }
        }
        return hadCollision;
    }

    bool checkCollisionAngular(Shape* s, std::vector<vec2>* collisions, std::vector<vec2>* normals, std::vector<vec2>* relativeVelocities, std::vector<vec2>* angularV1s, std::vector<vec2>* angularV2s, std::vector<vec2>* r1perps, std::vector<vec2>* r2perps, std::vector<vec2>* r1s, std::vector<vec2>* r2s, double eps) {
        bool hadCollision = false;
        glm::mat2 rotation1;
        rotation1[0][0] = glm::cos(angle);
        rotation1[0][1] = -glm::sin(angle);
        rotation1[1][0] = glm::sin(angle);
        rotation1[1][1] = glm::cos(angle);

        glm::mat2 rotation2;
        rotation2[0][0] = glm::cos(s->angle);
        rotation2[0][1] = -glm::sin(s->angle);
        rotation2[1][0] = glm::sin(s->angle);
        rotation2[1][1] = glm::cos(s->angle);

        for each (Voxel * v1 in voxels)
        {
            for each (Voxel * v2 in s->voxels)
            {
                vec2 v1pos = this->p + (vec2)(rotation1 * v1->pRelative);
                vec2 v2pos = s->p + (vec2)(rotation2 * v2->pRelative);
                double distance = glm::distance(v1pos, v2pos);
                if (distance < 2 * voxelRadius + eps)
                {
                    vec2 position = v1pos + (v2pos - v1pos) / 2.0;
                    vec2 normal = glm::normalize(v1pos - v2pos);

                    vec2 r1 = this->p - position;
                    vec2 r2 = s->p - position;

                    vec2 r1Perp = vec2(-r1.y, r1.x);
                    vec2 r2Perp = vec2(-r2.y, r2.x);

                    vec2 angularV1 = r1Perp * this->angularV;
                    vec2 angularV2 = r2Perp * s->angularV;

                    vec2 relativeVelocity = (this->v + angularV1) - (s->v + angularV2);
                    collisions->push_back(position);
                    normals->push_back(normal);
                    angularV1s->push_back(angularV1);
                    angularV2s->push_back(angularV2);
                    r1perps->push_back(r1Perp);
                    r2perps->push_back(r2Perp);
                    r1s->push_back(r1);
                    r2s->push_back(r2);
                    relativeVelocities->push_back(relativeVelocity);
                    hadCollision = true;
                }
            }
        }
        return hadCollision;
    }
};
