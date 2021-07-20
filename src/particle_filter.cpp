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
#define EPS 0.001
using std::string;
using std::vector;
using namespace std;

using std::normal_distribution;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  if (is_initialized) {
    return;
  }
  num_particles = 10;  // TODO: Set the number of particles
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x,std[0]);
  std::normal_distribution<double> dist_y(y,std[1]);
  std::normal_distribution<double> orient_theta(theta,std[2]);
  
  for(int i=0; i<num_particles; i++) {
    Particle P;
    P.id = i+1;
    P.x = dist_x(gen);
    P.y = dist_y(gen);
    P.theta = orient_theta(gen);
    P.weight = 1.0;
    particles.push_back(P);
    weights.push_back(1.0);
  } 
  is_initialized = true;
  //weights.resize(num_particles);
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */  
  
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0,std_pos[0]);
  std::normal_distribution<double> dist_y(0,std_pos[1]);
  std::normal_distribution<double> orient_theta(0,std_pos[2]);
  
   for (int i=0;i<num_particles;i++) {
   
    //To check if yaw_rate is too small to be insignificant
    if (fabs(yaw_rate) < EPS) { 
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      
    }else{
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate * delta_t;
    }
     //Adding Gaussian noise
     particles[i].x += dist_x(gen);
     particles[i].y += dist_y(gen);
     particles[i].theta += orient_theta(gen);
   }

}
  

void ParticleFilter::dataAssociation(vector<Map::single_landmark_s> predicted, 
                                     vector<LandmarkObs> &observation) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
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
  double stdLandmarkRange = std_landmark[0];
  double stdLandmarkBearing = std_landmark[1];
  
  for (int i=0;i<particles.size;i++) {
    
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
    double weight = 1.0;
    
    vector<Map::single_landmark_s> Landmark_inrange;
    for (unsigned int j=0;j<map_landmarks.landmark_list.size();j++) {
      float x_landmark = map_landmarks.landmark_list[j].x_f;
      float y_landmark = map_landmarks.landmark_list[j].y_f;
      int landmarkId = map_landmarks.landmark_list[j].id_i;
      //double range_distance = dist(particles[i].x, particles[i].y, x_landmark, y_landmark);
      if (fabs(x_landmark - p_x) <= sensor_range && fabs(y_landmark - p_y) <= sensor_range) {
        
      	Landmark_inrange.push_back(Map::single_landmark_s{landmarkId,x_landmark,y_landmark});
      }
    }
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
 
    for (unsigned int j=0;j<observations.size();j++) {
      LandmarkObs obs;
      obs.x = (observations[j].x * cos(p_theta)) - (observations[j].y * sin(p_theta)) + p_x;
      obs.y = (observations[j].y * cos(p_theta)) + (observations[j].x * sin(p_theta)) + p_y;
      obs.id = observations[j].id;
      
      double minDistance = numeric_limits<double>::max();
  	  int mapId = -1;
      double lmarkx, lmarky;
      for (unsigned int k=0;k<Landmark_inrange.size();k++) {
      
        Map::single_landmark_s p = Landmark_inrange[k];
        double distance = dist(obs.x, obs.y, p.x_f, p.y_f);
        if(distance<minDistance) {
          minDistance = distance;
          mapId = p.id_i;
          lmarkx = p.x_f;
          lmarky = p.y_f;
        }
      }
      associations.push_back(mapId);
      sense_x.push_back(lmarkx);
      sense_y.push_back(lmarky);
      
      double dX = obs.x - lmarkx;
      double dY = obs.y - lmarky;

      weight *= ( (1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp( -( dX*dX/(2*stdLandmarkRange*stdLandmarkRange)) - (dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)) ) );
    }
    
    SetAssociations(particles[i], associations, sense_x, sense_y);
    particles[i].weight = weight;
    weights.at(i) = weight;
    
    }
  
  }

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> weights;
  double maxWeight = numeric_limits<double>::min();
  for(int i=0;i<num_particles;i++){
    weights.push_back(particles[i].weight);
    if(maxWeight<particles[i].weight){
      maxWeight = particles[i].weight;
    }
  }
  
  //std::default_random_engine gen;
  std::uniform_real_distribution<double> randWeight(0.0,maxWeight);
  std::uniform_int_distribution<int> randVal(0,(num_particles-1));
 
  int index = randVal(gen);
  double randomWeight = randWeight(gen);
  double x_weight = 0.0;
  vector<Particle> resampledParticles;
  for(int i=0;i<num_particles;i++){
    x_weight += 2*randomWeight;
    while(weights[index] < x_weight){
      x_weight -= weights[index];
      index = (index+1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }
  particles = resampledParticles;
}



void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  //particle.associations.clear();
  //particle.sense_x.clear();
  //particle.sense_y.clear();
  
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
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
