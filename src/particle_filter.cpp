/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <sstream>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to first
   *   position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 200;
  
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle particle;
    particle.weight = 1.;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);

    particles.push_back(particle);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_gen useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_gen/
   */
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);

    for (int i = 0; i < num_particles; i++) {

    	if (fabs(yaw_rate) < 0.00001) {
    		particles[i].x += velocity * delta_t * cos(particles[i].theta);
    		particles[i].y += velocity * delta_t * sin(particles[i].theta);
    	} else {
    		particles[i].x += velocity / yaw_rate * (- sin(particles[i].theta) + sin(particles[i].theta + yaw_rate * delta_t));
    		particles[i].y += velocity / yaw_rate * (  cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
        	particles[i].theta += yaw_rate * delta_t;
        }

        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
     } //  end i

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
  /**
   * Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
  for (unsigned int i = 0; i < observations.size(); i++){

      double min_dist = std::numeric_limits<double>::max();
      int landmark_id;
      for (unsigned j = 0; j < predicted.size(); j++) {
        double euclidean_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

        if (euclidean_dist < min_dist){
          min_dist = euclidean_dist;
          landmark_id = predicted[j].id;
        }
        observations[i].id = landmark_id;
      } //end j
    } //end i
} //end fct


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // for each particle
  for (int i = 0; i < num_particles; i++){
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    // for each landmark
    vector<LandmarkObs> landmark_predictions;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      double landmark_x = map_landmarks.landmark_list[j].x_f;
      double landmark_y = map_landmarks.landmark_list[j].y_f;

      if (fabs(landmark_x-particle_x) <= sensor_range && fabs(landmark_y-particle_y) <= sensor_range) {
        landmark_predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
      }
    } // end j

    //for each observation
    vector<LandmarkObs> observations_map;
    for (unsigned int k = 0; k < observations.size(); k++) {
      int observations_id = observations[k].id;
      double observations_x = particle_x + cos(particle_theta) * observations[k].x - sin(particle_theta) * observations[k].y;
      double observations_y = particle_y + sin(particle_theta) * observations[k].x + cos(particle_theta) * observations[k].y;
      observations_map.push_back(LandmarkObs{observations_id, observations_x, observations_y});
    } // end k
    
    dataAssociation(landmark_predictions, observations_map);

    particles[i].weight = 1.0;
    // for each obs
    for (unsigned int m = 0; m < observations_map.size(); m++) {
      double obs_x;
      double obs_y;
      double predict_x;
      double predict_y;
      obs_x = observations_map[m].x;
      obs_y = observations_map[m].y;

      for (unsigned int k = 0; k < landmark_predictions.size(); k++) {
        if (landmark_predictions[k].id == observations_map[m].id) {
          predict_x = landmark_predictions[k].x;
          predict_y = landmark_predictions[k].y;
        }
     } // end k

      double diff_x2 = pow((obs_x - predict_x), 2);
      double diff_y2 = pow((obs_y - predict_y), 2);
      double sigma_x = std_landmark[0];
      double sigma_y = std_landmark[1];
      double norm = (1. / (2. * M_PI * sigma_x * sigma_y));
      double prob = norm * exp( -( diff_x2/(2*pow(sigma_x, 2)) + (diff_y2/(2*pow(sigma_y, 2))) ) );

      particles[i].weight *= prob;
    } // end m
  } // end i
} // end fct

void ParticleFilter::resample()
{
  /**
   * Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;

  vector<double> weights;
  double weight_max = std::numeric_limits<double>::min();

  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if (particles[i].weight > weight_max) {
      weight_max = particles[i].weight;
    }
  } // end i

  std::uniform_real_distribution<double> dist_double(0.0, weight_max);
  std::uniform_int_distribution<int> dist_int(0, num_particles-1);
  
  int index = dist_int(gen);
  double beta = 0.;
  vector<Particle> resampled_particles;
  
  for (int j = 0; j < num_particles; j++) {
      beta += dist_double(gen) * 2.;
      while (beta > weights[index]) {
         beta -= weights[index];
         index = (index + 1) % num_particles;
      }
      resampled_particles.push_back(particles[index]);
  }  // end j
  particles = resampled_particles;

} // end fct

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X"){
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}