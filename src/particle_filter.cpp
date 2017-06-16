/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
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
	num_particles = 10;

	// Taking into account Gaussian Sensor Noise around initial heading estimation
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	// Initialize all particles to first position
	for (auto i = 0; i < num_particles; i+=1) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
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
	//Initialize all the particles to the (first)GPS position as mean
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	std::default_random_engine gen;
	double threshold = 1e-5;

	for (auto &p:particles){
		double x = p.x;
		double y = p.y;
		double theta = p.theta;
		//Add measurements using the physical Model
		if (fabs(yaw_rate) < threshold){
			linearMotionParticleProgress(&p, velocity, delta_t);
		} else{
			nonLinearMotionParticleProgress(&p, velocity, delta_t, yaw_rate);
		}

		// Add Noise - creates a normal (Gaussian) distribution for x,y and theta
		std::normal_distribution<double> dist_x(p.x, std_x);
		std::normal_distribution<double> dist_y(p.y, std_y);
		std::normal_distribution<double> dist_psi(p.theta, std_theta);

		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta  = dist_psi(gen);
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
	vector<double> sense_x_;
	vector<double> sense_y_;

	sense_x_.resize(observations.size());
	sense_y_.resize(observations.size());

	weights.clear();

	for (auto &p : particles) {
		p.associations.clear();
		p.associations.resize(observations.size());

		sense_x_.clear();
		sense_y_.clear();

		for (int i = 0; i < observations.size(); i +=1) {
			double min_length = numeric_limits<double>::max();
      // observation in car to observation in map
			double x = cos(p.theta)*observations[i].x - sin(p.theta)*observations[i].y + p.x;
			double y = sin(p.theta)*observations[i].x + cos(p.theta)*observations[i].y + p.y;

      // Data association
			for (auto &lm: map_landmarks.landmark_list){
				double distance = dist(x, y, lm.x_f, lm.y_f);
				double distance_2 = dist(p.x, p.y, lm.x_f, lm.y_f);
				if (distance < min_length && distance_2 < sensor_range){
					min_length = distance;
					// id_i (1 based)
					p.associations[i] = lm.id_i - 1;
					sense_x_[i] = x;
					sense_y_[i] = y;
				}
			}
		}

    // Calc total probability for each particle
		double total_prob = 1.0;

    const double std_x = std_landmark[0];
    const double std_y = std_landmark[1];

    const double scalar = 1.0/( 2.0 * M_PI * std_x * std_y);
    const double dx_divider = 2.0*pow(std_x,2);
    const double dy_divider = 2.0*pow(std_y,2);

		for (int m =0; m < observations.size();m++){

			const double lm_x = map_landmarks.landmark_list[p.associations[m]].x_f;
			const double lm_y = map_landmarks.landmark_list[p.associations[m]].y_f;

      const double dx2 = pow(lm_x - sense_x_[m], 2);
      const double dy2 = pow(lm_y - sense_y_[m], 2);
      total_prob *= scalar * exp (-(dx2/dx_divider + dy2/dy_divider));
		}
		p.weight = total_prob;
		weights.push_back(total_prob);
	}
}

void ParticleFilter::resample() {
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::discrete_distribution<int> d(weights.begin(), weights.end()); // Define a discrete distribution
	std::vector<Particle> new_particles; // Resampled particles holder
	std::default_random_engine gen3;

	for (int i = 0; i< num_particles; i++){
		auto index = d(gen3);
		new_particles.push_back(std::move(particles[index]));
	}
	//assign the particles from holder to the original
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
