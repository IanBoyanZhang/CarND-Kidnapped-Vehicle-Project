/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "map.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//num_particles = 1000;
	num_particles = 1;
	// Taking into account Gaussian Sensor Noise around initial heading estimation
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Initialize all particles to first position
  for (auto i = 0; i < num_particles; i+=1) {
		Particle particle;
		particle.id = i;
		particle.x = generateGaussianVariable(x, std_x);
		particle.y = generateGaussianVariable(y, std_y);
		particle.theta = generateGaussianVariable(theta, std_theta);
		weights.push_back(1);
    particle.weight = 1;
		particles.push_back(particle);
	}
	weights.resize(num_particles, 1.0);
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	/**
	 * How to write this in lambda function and callback style,
	 * though it may be detrimental to runtime efficiency?
	 */
	/**
	 * TODO: refactor duplicated code
	 * TODO: Cache random number generator?
	 */
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	// TODO: Add gaussian noise to velocity and yaw_rate
	double threshold = 1e-4;
	for (auto i = 0; i < num_particles; i+=1) {
		if (yaw_rate < threshold) {
			linearMotionParticleProgress(&particles[i], velocity, delta_t);
			continue;
		}
		nonLinearMotionParticleProgress(&particles[i], velocity, delta_t, yaw_rate);
	}

	for (auto &p:particles) {
		p.x = generateGaussianVariable(p.x, std_x);
		p.y = generateGaussianVariable(p.y, std_y);
		p.theta = generateGaussianVariable(p.theta, std_theta);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// associate landmark to this observation (chose the landmark with the shortest)
	// distance to observation
	// Factor below part into data association function?
	// By default associate with first landmark
				// Filtering landmark range
/*        if (distance >= sensor_range) {
					continue;
				}
				if (distance < min) {
					min = distance;
					min_id = l.id_i;
				}
			}
			p.associations.push_back(min_id);*/
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// You would need to normalize the weights for calculations

  weights.clear();
	// TODO: sensor_range for filtering
	// Time complexity M (#_of_particles) * N (#_of_observations) * Q (#_of_landmarks)
	for (auto& p: particles) {
    p.associations.clear();
		p.sense_x.clear();
		p.sense_y.clear();
		for (const auto &o: observations) {
			// Coordinate transform
			//   http://planning.cs.uiuc.edu/node99.html

			// Observation obs in global/map coordinate
      // Factoring below code to data observation
			double o_x_map = o.x * cos(p.theta) - o.y * sin(p.theta) + p.x;
			double o_y_map = o.x * sin(p.theta) + o.y * cos(p.theta) + p.y;

			// Filtering landmark
			double min = numeric_limits<double>::max();
			double distance = 0;
			double distance_p_to_lm = 0;
      p.associations.push_back(-1);
			p.sense_x.push_back(o_x_map);
			p.sense_y.push_back(o_y_map);
			for (const auto &lm: map_landmarks.landmark_list) {
				distance = dist(o_x_map, o_y_map, lm.x_f, lm.y_f);
				distance_p_to_lm = dist(p.x, p.y, lm.x_f, lm.y_f);
				if (distance_p_to_lm < sensor_range && distance < min) {
					min = distance;
					p.associations[o.id] = lm.id_i;
          p.sense_x[o.id] = o_x_map;
					p.sense_y[o.id] = o_y_map;
				}
			}
		}
	}

  // The calculate Multi-variate Gaussian probability (x, y) for each observation of
	// particle

	// calculate particle weight
	const double std_x = std_landmark[0];
	const double std_y = std_landmark[1];

	const double scalar = 1.0/( 2.0 * M_PI * std_x * std_y);
	const double dx_divider = 2.0*pow(std_x,2);
	const double dy_divider = 2.0*pow(std_y,2);

  double inter = 0;
	// 2D Gaussian
	for (auto& p:particles) {
    p.weight = 1;
		for (const auto& o:observations) {
			//const double dx2 = pow(p.sense_x[o.id] - map_landmarks.landmark_list[p.associations[o.id] - 1].x_f ,2);
			const double dx2 = pow(p.sense_x[o.id] - map_landmarks.landmark_list[p.associations[o.id]].x_f ,2);
			//const double dy2 = pow(p.sense_y[o.id] - map_landmarks.landmark_list[p.associations[o.id] - 1].y_f, 2);
			const double dy2 = pow(p.sense_y[o.id] - map_landmarks.landmark_list[p.associations[o.id]].y_f, 2);
			cout << "dx2: " << dx2 << endl;
      cout << "dy2: " << dy2 << endl;

			inter = scalar*exp(-(dx2/dx_divider + dy2/dy_divider));
			cout << "intermediate: " << inter << endl;
			p.weight *= inter;
		}
		weights.push_back(p.weight);
    cout << "P.weight: " << p.weight << endl;
	}

	// Summation
  double sum_of_weights = accumulate(weights.begin(), weights.end(), 0.0);

	cout << "sum_of_weights" << sum_of_weights << endl;

	// Normalization
	for (auto& w:weights) {
		w/=sum_of_weights;
	}
	weights = weights;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  random_device rd;
	mt19937_64 gen(rd());
	discrete_distribution<double> d(weights.begin(), weights.end());
	vector<Particle> new_particles;
  for (auto i = 0; i < num_particles; i+=1) {
    auto index = d(gen);
		new_particles.push_back(std::move(particles[index]));
	}
  particles = std::move(new_particles);
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

void ParticleFilter::linearMotionParticleProgress(Particle *particle, const double v, const double dt) {
	double x_0 = particle->x;
	double y_0 = particle->y;
  double yaw_0 = particle->theta;

	double x_f = x_0 + v * dt * cos(yaw_0);
	double y_f = y_0 + v * dt * sin(yaw_0);

	particle->x = x_f;
  particle->y = y_f;
}

void ParticleFilter::nonLinearMotionParticleProgress(Particle *particle, const double v,
																										 const double dt, const double yaw_rate) {
  double x_0 = particle->x;
	double y_0 = particle->y;
  double yaw_0 = particle->theta;

  double yaw_f = yaw_0 + yaw_rate * dt;
	double x_f = x_0 + v / yaw_rate * (sin(yaw_f) - sin(yaw_0));
	double y_f = y_0 + v / yaw_rate * (cos(yaw_0) - cos(yaw_f));

	particle->x = x_f;
	particle->y = y_f;
	particle->theta = yaw_f;
}

double ParticleFilter::generateGaussianVariable(const double var_, const double std_) {
	default_random_engine gen;
	// Create a normal (gaussian distribution)
  normal_distribution<double> dist_var(var_, std_);
	return dist_var(gen);
}

