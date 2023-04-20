#pragma once
#include "Particle.hpp"
#include "Spring.hpp"
#include "Shape.hpp"

class CollisionHandler {
public:
	
	bool collision;
	bool useRepulsion;
	int maxIterations;
	double restitutionValue;
	double thresholdDistance;
	double eps;
	double repulsionStiffness;


	std::vector<vec2> collisions;
	std::vector<vec2> normals;
	std::vector<vec2> relativeVelocities;
	std::vector<vec2> angularV1s;
	std::vector<vec2> angularV2s;
	std::vector<vec2> r1s;
	std::vector<vec2> r2s;
	std::vector<vec2> r1perps;
	std::vector<vec2> r2perps;
	std::vector<vec2> impulses;

	/** TODO: number of iterations in the last CCD processing loop, to keep an eye on how tricky the problem is */
	int iterations;
	/** TODO: keep track of highest count needed so far */
	int highestIterations = 0;

	CollisionHandler::CollisionHandler() :
		eps(0.000001),
		collision(true),
		useRepulsion(true),
		maxIterations(500),
		restitutionValue(0.0),
		thresholdDistance(2),
		iterations(0),
		repulsionStiffness(100) {}

	void applyImpulses(double h, std::vector<Shape*>* shapes, int width, int height)
	{
		double e = .9;
		collisions.clear();
		normals.clear();
		relativeVelocities.clear();
		for (int i = 0; i < shapes->size(); i++)
		{
			shapes->at(i)->checkCollisionWithWalls(height, width, eps, e);
		}
		for (int i = 0; i < shapes->size(); i++)
		{
			for (int j = i + 1; j < shapes->size(); j++)
			{
				Shape *shape1 = shapes->at(i);
				Shape *shape2 = shapes->at(j);
				if (shape1->checkCollision(shape2, &collisions, &normals, &relativeVelocities, eps))
				{
					for (int k = 0; k < collisions.size(); k++)
					{
						double m1 = shape1->mass;
						double m2 = shape2->mass;

						vec2 n = normals[k];
						vec2 relV = relativeVelocities[k];

						// moving away, so no collision
						if (glm::dot(relV, n) > eps)
							continue;

						double J = -(1 + e) * glm::dot(relV, n) / (glm::dot(n, n) * (1 / m1 + 1 / m2));

						// Update velocities
						shape1->v = shape1->v + (J / m1) * n;
						shape2->v = shape2->v - (J / m2) * n;
					}
					/*vec2 p1 = vox1->p;
					vec2 p2 = vox2->p;
					vec2 v1 = vox1->v;
					vec2 v2 = vox2->v;
					double m1 = vox1->mass;
					double m2 = vox2->mass;

					// Compute the normal
					vec2 n = p1 - p2;

					// compute relative velocity
					vec2 relV = v1 - v2;
					// if relative normal velocity is 0 then no collision
					if (glm::dot(relV, -n) < -eps)
					{
						continue;
					}
					// Compute impulse
					double e = 1;
					double J = -(1 + e) * glm::dot(relV, n) / (glm::dot(n, n) * (1 / m1 + 1 / m2));

					// Update velocities
					vox1->v = v1 + (J / m1) * n;
					vox2->v = v2 - (J / m2) * n;
					vec2 tempV1 = vox1->v;
					vec2 tempV2 = vox2->v;*/
				}
			}
		}
	}
	void applyImpulsesWithAngular(double h, std::vector<Shape*>* shapes, int width, int height)
	{
		double e = .3;
		bool hadCollisions = false;
		int iterations = 0;
		int maxIterations = 2;
		do {
			collisions.clear();
			normals.clear();
			relativeVelocities.clear();
			angularV1s.clear();
			angularV2s.clear();
			r1perps.clear();
			r2perps.clear();
			r1s.clear();
			r2s.clear();
			impulses.clear();
			hadCollisions = false;
			for (int i = 0; i < shapes->size(); i++)
			{
				if (shapes->at(i)->checkCollisionWithWallsWithAngular(height, width, eps, e))
				{
					hadCollisions = true;
				}
			}
			for (int i = 0; i < shapes->size(); i++)
			{
				for (int j = i + 1; j < shapes->size(); j++)
				{
					Shape* shape1 = shapes->at(i);
					Shape* shape2 = shapes->at(j);

					int totalCollisions = 0;
					if (shape1->checkCollisionAngular(shape2, &collisions, &normals, &relativeVelocities, &angularV1s, &angularV2s, &r1perps, &r2perps, &r1s, &r2s, eps))
					{
						double m1 = shape1->mass;
						double m2 = shape2->mass;
						for (int k = 0; k < collisions.size(); k++)
						{

							vec2 n = normals[k];
							vec2 relV = relativeVelocities[k];

							// moving away, so no collision
							if (glm::dot(relV, -n) < eps)
							{
								impulses.push_back(vec2(0, 0));
								continue;
							}
							hadCollisions = true;

							totalCollisions++;

							double r1dotn = glm::dot(r1perps[k], n);
							double r2dotn = glm::dot(r2perps[k], n);

							double J = -(1 + e) * glm::dot(relV, n) / (glm::dot(n, n) * (1 / m1 + 1 / m2) + ((r1dotn * r1dotn / shape1->inertia) + (r2dotn * r2dotn / shape2->inertia)));
							impulses.push_back(J * n);
							vec2 impulse = J * n;
							shape1->v = shape1->v + (impulse / m1);
							shape2->v = shape2->v - (impulse / m2);

							shape1->angularV = shape1->angularV + (glm::dot(r1perps[k], impulse) / shape1->inertia);
							shape2->angularV = shape2->angularV - (glm::dot(r2perps[k], impulse) / shape2->inertia);

						}

						/*for (int k = 0; k < collisions.size(); k++)
						{
							// Update velocities
							if (impulses[k].x == 0 && impulses[k].y == 0)
								continue;
							vec2 impulse = impulses[k];// / (double)totalCollisions;
							shape1->v = shape1->v + (impulse / m1);
							shape2->v = shape2->v - (impulse / m2);

							shape1->angularV = shape1->angularV + (glm::dot(r1perps[k], impulse) / shape1->inertia);
							shape2->angularV = shape2->angularV - (glm::dot(r2perps[k], impulse) / shape2->inertia);
						}*/
					}
				}
			}
			iterations++;
		} while (hadCollisions && iterations < maxIterations);

		
	
		/*for (int i = 0; i < shapes->size(); i++)
		{
			if (glm::dot(shapes->at(i)->v, shapes->at(i)->v) > 50000)
			{
				shapes->at(i)->v = shapes->at(i)->v / 1000.0;
			}

			if (shapes->at(i)->angularV > 1000)
			{
				shapes->at(i)->angularV = shapes->at(i)->angularV / 1000.0;
			}
				
		}*/
	}

	/**
	 * Try to deal with contacts before they happen
	 * @param h
	 * @param system
	 */
	/*void applyRepulsion(double h, std::vector<Voxel*>* particles) {
		//TODO: apply repulsion on all particles

		// - check for point-edge proximity between all particle-edge pairs
		/*for (int i = 0; i < particles->size(); i++)
		{
			for (int j = 0; j < springs->size(); j++)
			{
				Particle* p0 = particles->at(i);
				Particle* p1 = springs->at(j)->p1;
				Particle* p2 = springs->at(j)->p2;

				if (p1 != p0 && p2 != p0)
				{
					// - find the normal
					vec2 lineToPoint = p0->p - p1->p;
					vec2 u0 = p1->p;
					vec2 u1 = p2->p - p1->p;
					vec2 projectionLineToPoint = glm::dot(lineToPoint, u1) / glm::dot(u1, u1) * u1 +u0;
					vec2 n = projectionLineToPoint - (lineToPoint + u0);
					//double distance = glm::abs((p2->p.x - p1->p.x) * (p0->p.y - p1->p.y) - (p2->p.y - p1->p.y) * (p0->p.x - p1->p.x)) / glm::sqrt((p2->p.x - p1->p.x) * (p2->p.x - p1->p.x) + (p2->p.y - p1->p.y) * (p2->p.y - p1->p.y));
					double distance = glm::sqrt(glm::dot(n, n));
					n = n / glm::sqrt(glm::dot(n, n)); // normalize

					// - TODO: take care to deal with segment end points carefully
					double alpha = glm::dot((p2->p - p1->p), (projectionLineToPoint - p1->p)) / (glm::dot(p2->p - p1->p, p2->p - p1->p));

					// Note here that if the normal vector has length zero then they might be coming at eachother "parallel"
					// Which means that they will still impact but there is no movement relative to the normal.
					if (-eps <= distance && distance <= eps)
					{
						double relativeToNormal = 0;
						if (alpha < 0.5)
						{
							n = p1->p-p2->p;
							distance = - alpha * glm::length(n);
							n = -glm::normalize(n);
							relativeToNormal = glm::dot(p1->v - p0->v, n);
						}
						else
						{
							n = p2->p - p1->p;
							distance = -(1 - alpha) * glm::length(n);
							n = -glm::normalize(n);
							relativeToNormal = glm::dot(p2->v - p0->v, n);
						}
						// - compute an appropriate  impulse if the distance is less than H
						double d = thresholdDistance - distance;

						// - don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
						// if they are approaching eachother
						if (relativeToNormal < 0.1 * d / h + eps)
						{
							if (d >= eps)
							{
								double m0 = p0->pinned ? DBL_MAX : p0->mass;
								double m1 = p1->pinned ? DBL_MAX : p1->mass;
								double m2 = p2->pinned ? DBL_MAX : p2->mass;

								double J = -glm::min(h * repulsionStiffness * d, m0 * (0.1 * d / h - relativeToNormal));
								//cout << "1 --- d: " << d << " H: " << thresholdDistance << " actual distance: " << distance << " alpha: " << alpha << endl;

								p0->v = p0->v + (J / m0) * n;
								p1->v = p1->v - (J / m1) * n * (1 - alpha);
								p2->v = p2->v - (J / m2) * n * (alpha);
							}
						}
					}
					else if (-eps <= alpha && alpha <= 1 + eps)
					{
						double relativeToNormal = glm::dot((1 - alpha) * p1->v + alpha * p2->v - p0->v, n);

						// - compute an appropriate  impulse if the distance is less than H
						double d = thresholdDistance - distance;

						// - don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
						// if they are approaching eachother
						if (relativeToNormal < 0.1 * d / h + eps)
						{
							if (d >= eps)
							{
								double m0 = p0->pinned ? DBL_MAX : p0->mass;
								double m1 = p1->pinned ? DBL_MAX : p1->mass;
								double m2 = p2->pinned ? DBL_MAX : p2->mass;

								double J = -glm::min(h * repulsionStiffness * d, m0 * (0.1 * d / h - relativeToNormal));
								//cout << "2 --- d: " << d << " H: " << thresholdDistance << " actual distance: " << distance << " alpha: " << alpha << endl;

								p0->v = p0->v + (J / m0) * n;
								p1->v = p1->v - (J / m1) * n * (1 - alpha);
								p2->v = p2->v - (J / m2) * n * (alpha);
							}
						}
					}
					

					
					// - make sure your impulse is in the correct direction!
					// - distribute impulse to the three particles involved in the appropriate manner
				}
			}
		}
		

	}
	*/
	/**
	* Checks all collisions in interval t to t+h
	* @param h
	* @param system
	* @return true if all collisions resolved
	*/
	/*bool checkCollision(double h, std::vector<Voxel*>* particles) {
		// TODO: Most of your assignment code will go here
		// - put your iterative solving code here!
		// - you can use the nextCollision function below to deal with an individual particle-edge pair

		// For each particle-edge pair, find the roots for when the three particles are
		// co-linear, and then pick the first root on (0,h] which corresponds to an actually collision.
		// compute a collision response.  That is, compute an appropriate collision normal, compute the 
		// impulse, and then apply the impulse to the associated particles.  Be sure to deal with pinning 
		// constraints!
		iterations = 0;
		bool hadCollision = true;
		while (hadCollision && iterations < maxIterations)
		{
			hadCollision = false;
			for (int i = 0; i < particles->size(); i++)
			{
				for (int j = 0; j < springs->size(); j++)
				{
					Particle* p0 = particles->at(i);
					Particle* p1 = springs->at(j)->p1;
					Particle* p2 = springs->at(j)->p2;

					if (p1 != p0 && p2 != p0)
					{
						if (!nextCollision(h, restitutionValue, p0, p1, p2))
						{
							hadCollision = true;
						}
					}
				}
			}
			iterations++;
		}


		if (highestIterations < iterations)
		{
			highestIterations = iterations;
		}

		if (iterations >= maxIterations)
		{
			return false;
		}
		return true;
	}
	*/
	/**
	* Processes next collision on (0,h] if it exists, between p0 and segment connecting p1 and p2.
	*
	* Watch out for
	* - normal direction
	* - epsilon tests
	* - segment connecting p1 and p2 at zero length! (unlikely?)
	*
	* @param h				timestep
	* @param restitution	bouncyness parameter
	* @param p0				particle
	* @param p1				line segment end point
	* @param p2				line segment end point
	* @return true if collision processed
	*/
	/*bool CollisionHandler::nextCollision(double h, double restitution, Particle* p0, Particle* p1, Particle* p2) {
		//TODO: 
		// * - finds roots
		double A1 = p0->p.x;
		double A2 = p0->p.y;

		double B1 = p1->p.x;
		double B2 = p1->p.y;

		double C1 = p2->p.x;
		double C2 = p2->p.y;

		double Adot1 = p0->v.x;
		double Adot2 = p0->v.y;

		double Bdot1 = p1->v.x;
		double Bdot2 = p1->v.y;

		double Cdot1 = p2->v.x;
		double Cdot2 = p2->v.y;

		// code generated by MATLAB
		double a = Adot1*Bdot2-Adot2*Bdot1-Adot1*Cdot2+Adot2*Cdot1+Bdot1*Cdot2-Bdot2*Cdot1;
		double b = A1 * Bdot2 - A2 * Bdot1 + Adot1 * B2 - Adot2 * B1 - A1 * Cdot2 + A2 * Cdot1 - Adot1 * C2 + Adot2 * C1 + B1 * Cdot2 - B2 * Cdot1 + Bdot1 * C2 - Bdot2 * C1;
		double c = A1 * B2 - A2 * B1 - A1 * C2 + A2 * C1 + B1 * C2 - B2 * C1;

		double t = -1;
		// if quadratic
		if (-eps > a && a > eps)
		{
			double delta = glm::sqrt(b * b - 4 * a * c);
			double t1 = (-b + delta) / (2 * a);
			double t2 = (-b - delta) / (2 * a);

			if (t1 >= -eps && t2 >= -eps)
			{
				t = glm::min(t1, t2);
			}
			else if (t1 >= -eps)
			{
				t = t1;
			}
			else if (t2 >= -eps)
			{
				t = t2;
			}
		}
		// otherwise
		else
		{
			if (-eps <= b && b <= eps && -eps <= c && c <= eps)
			{
				t = 0;
			}
			else
			{
				t = -c / b;
			}
		}

		

		// Then there is collision
		if (t >= -eps)
		{
			if (t > h + eps)
			{
				return true; // Collision doesn't happen in time interval
			}
			// else collision does occur in time interval
		}
		else
		{
			return true; // Collision doesn't happen
		}


		// * - checks that p0 actually falls on segment
		// Compute alpha
		glm::vec2 p1h = t * p1->v + p1->p;
		glm::vec2 p2h = t * p2->v + p2->p;
		glm::vec2 p0h = t * p0->v + p0->p;

		glm::vec2 repVec = p2h - p1h;
		double denominatorValue = glm::dot(repVec, repVec);

		double alpha = glm::dot((repVec), (p0h - p1h))/(denominatorValue);
		//cout << "alpha " << alpha << endl;

		if (-eps <= alpha && alpha <= 1 + eps)
		{ 
			/*p0->color = glm::vec3(0.0f, 0.0f, 1.0f);
			p1->color = glm::vec3(1.0f, 0.0f, 1.0f);
			p2->color = glm::vec3(1.0f, 0.0f, 1.0f);
			cout << "a: " << a << endl;
			cout << "b: " << b << endl;
			cout << "c: " << c << endl;
			cout << "204, alpha: " << alpha << endl;
			cout << "p0: " << p0->p.x << " " << p0->p.y << endl;
			cout << "p1: " << p1->p.x << " " << p1->p.y << endl;
			cout << "p2: " << p2->p.x << " " << p2->p.y << endl;
			cout << "p0dot: " << p0->v.x << " " << p0->v.y << endl;
			cout << "p1dot: " << p1->v.x << " " << p1->v.y << endl;
			cout << "p2dot: " << p2->v.x << " " << p2->v.y << endl;
			cout << "p0h: " << p0h.x << " " << p0h.y << endl;
			cout << "p1h: " << p1h.x << " " << p1h.y << endl;
			cout << "p2h: " << p2h.x << " " << p2h.y << endl;
			cout << "t: " << t << endl;

			// * - processes collision by computing and applying impulses
			// First we obtain the normal
			vec2 lineToPoint = p1->p - p0->p;
			vec2 u0 = p1->p;
			vec2 u1 = p2->p - p1->p;
			vec2 projectionLineToPoint = glm::dot(lineToPoint - u0, u1)/glm::dot(u1, u1) * u1 + u0;
			vec2 n =  projectionLineToPoint - lineToPoint;
			n = n / glm::sqrt(glm::dot(n, n)); // normalize

			// relative speed 
			vec2 vrel = p0->v - ((1 - alpha) * p1->v + (alpha)*p2->v);

			// obtain impulse maginitude
			// handle infinite mass (for pinned particles)

			double m0 = p0->pinned ? DBL_MAX : p0->mass;
			double m1 = p1->pinned ? DBL_MAX : p1->mass;
			double m2 = p2->pinned ? DBL_MAX : p2->mass;

			// The following is just to handle "infinite mass" without running into overflow
			double denominatorPoint = 1 / m0;
			if (p0->pinned)
				denominatorPoint = 0;

			double denominatorSpring = ((1 - alpha) * m1 + alpha * m2);
			if (denominatorSpring < 0 || denominatorSpring > 1000)
				denominatorSpring = 0;
			else
				denominatorSpring = 1 / denominatorSpring;
			/*
			if (p0->pinned && p1->pinned && !p2->pinned)
			{
				overallPointMass = 0;
				overallSpringMass = ((1 - alpha) * m1 + alpha * m2);
			}
			else if (p0->pinned && !p1->pinned && p2->pinned)
			{
				overallPointMass = 0;
				overallSpringMass = 1 / m1;
			}
			else if (!p0->pinned && p1->pinned && p2->pinned)
			{
				overallPointMass = 1 / m0;
				overallSpringMass = 0;
			}
			else if (!p0->pinned && !p1->pinned && p2->pinned)
			{
				overallPointMass = 1 / m0;
				overallSpringMass = 1 / m1;
			}
			else if (!p0->pinned && p1->pinned && !p2->pinned)
			{
				overallPointMass = 1 / m0;
				overallSpringMass = 1 / m2;
			}
			else if (p0->pinned && !p1->pinned && !p2->pinned)
			{
				overallPointMass = 0;
			}
			else
			{

			}

			double J = -(1 + restitution) * (glm::dot(vrel, n)) / glm::max(denominatorPoint + denominatorSpring, 1.0);
			if (-J > 100000 || J > 100000)
			{
				cout << "J: " << J << " m0: " << m0 << " m1: " << m1 << " m2: " << m2 << " dot product: " << glm::dot(vrel, n) << " Overall Spring Mass: " << denominatorSpring << endl;
			}
			// apply forces to particles
			p0->v = p0->v + (J / m0) * n;
			p1->v = p1->v - (J / m1) * n * (1 - alpha);
			p2->v = p2->v - (J / m2) * n * (alpha);
		}
		else
		{
			return true; //no collision
		}

		// * - returns true if collision processed
		//cout << "detected a collision, alpha: " << alpha << endl;
		return false;
	}
	*/
};

