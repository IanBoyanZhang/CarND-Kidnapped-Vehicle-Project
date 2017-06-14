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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 1000;
	for (auto i = 0; i < num_particles; i+=1) {
		weights.push_back(1);
	}
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
		particles.push_back(particle);
	}
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

	for (auto i = 0; i < num_particles; i+=1) {
		Particle particle = particles[i];
		particle.id = i;
		particle.x = generateGaussianVariable(particle.x, std_x);
		particle.y = generateGaussianVariable(particle.y, std_y);
		particle.theta = generateGaussianVariable(particle.theta, std_theta);
		particles[i] = particle;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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

	// TODO: sensor_range for filtering
	// Time complexity M (#_of_particles) * N (#_of_observations) * Q (#_of_landmarks)
	for (auto& p: particles) {
    for (const auto& o: observations) {
			// Coordinate transform
			//   http://planning.cs.uiuc.edu/node99.html
			p.sense_x.push_back(o.x * cos(p.theta) - o.y * sin(p.theta) + p.x);
			p.sense_y.push_back(o.x * sin(p.theta) + o.y * sin(p.theta) + p.y);
			// associate landmark to this observation (chose the landmark with the shortest)
			// distance to observation
      // Factor below part into data association function?
			double min = numeric_limits<double>::max();
			// By default associate with first landmark
			double min_id = 0;
			double distance = 0;
			for (const auto& l: map_landmarks.landmark_list) {
				distance = dist(p.sense_x[o.id], p.sense_y[o.id], l.x_f, l.y_f);
				if (distance < min) {
					min = distance;
					min_id = l.id_i;
				}
			}
			p.associations.push_back(min_id);
//      cout << p.associations[o.id] << endl;
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

 	vector<double> weights_;

	// 2D Gaussian
	for (auto& p:particles) {
    p.weight = 1;
		for (const auto& map_obs:observations) {
			const double dx2 = pow(p.sense_x[map_obs.id] - map_landmarks.landmark_list[p.associations[map_obs.id]].x_f ,2);
			const double dy2 = pow(p.sense_y[map_obs.id] - map_landmarks.landmark_list[p.associations[map_obs.id]].y_f, 2);
			p.weight *= scalar*exp(-(dx2/dx_divider + dy2/dy_divider));
		}
		weights_.push_back(p.weight);
	}

	// Summation
  double sum_of_weights = accumulate(weights_.begin(), weights_.end(), 0.0);

	// Normalization
	for (auto& w:weights_) {
		w/=sum_of_weights;
	}
	weights = weights_;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  random_device rd;
	mt19937_64 gen(rd());
	discrete_distribution<double> d(weights.begin(), weights.end());
	map<int, int> m;
	vector<Particle> new_particles;
  for (auto i = 0; i < num_particles; i+=1) {
		new_particles.push_back(particles[d(gen)]);
	}
  particles = new_particles;
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

