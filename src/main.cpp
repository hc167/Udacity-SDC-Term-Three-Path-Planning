// Author: Hiu Chan
// Bosch Challenge: Path Planning self driving car
//

#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
double mph2ms(double x) { return x/2.23694; }
double ms2mph(double x) { return x*2.23694; }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};
}

void rotatePoint(double x, double y, double angle, double & target_x, double & target_y)
{
  target_x = x*cos(angle) - y*sin(angle);
  target_y = x*sin(angle) + y*cos(angle);
}


tk::spline getSmoothCurvePathPlan(double x1, double x2, double y1, double y2, double angle, double car_s, double lane, 
				  vector<double> & map_waypoints_s, vector<double> & map_waypoints_x, vector<double> & map_waypoints_y)
{

  // create a list of widely spaced (x, y) waypoints.
  // later we will interoplate these waypoint with spline and fill it in with more pointss than control speed
  vector<double> ptsx;
  vector<double> ptsy;

  // Start next point with car+35, that way we can reduce the AccN significantly.
  vector<double> next_wp0 = getXY(car_s+35, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  vector<double> next_wp1 = getXY(car_s+50, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  vector<double> next_wp2 = getXY(car_s+65, 2+4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

  ptsx.push_back(x1);
  ptsx.push_back(x2);
  
  ptsx.push_back(next_wp0[0]);
  ptsx.push_back(next_wp1[0]);
  ptsx.push_back(next_wp2[0]);

  ptsy.push_back(y1);
  ptsy.push_back(y2);

  ptsy.push_back(next_wp0[1]);
  ptsy.push_back(next_wp1[1]);
  ptsy.push_back(next_wp2[1]);

  for(int i=0; i< ptsx.size(); ++i) {
    double shift_x = ptsx[i] - x2;
    double shift_y = ptsy[i] - y2;

    // negative to rotate counter-clockwise                                                                                    
    rotatePoint(shift_x, shift_y, -angle, ptsx[i], ptsy[i]);
  }

  tk::spline s;
  s.set_points(ptsx, ptsy);

  return s;
}

int getLane(double lane)
{
  for (int i=0; i < 3; ++i){
    if (lane >= i*4 && lane < i*4+4)
      return i;
  }
  // return -1 if the lane is unknown. Who care in this case?
  return -1;
}

bool noCarBehindThisLane(int checklane, double mySpeed, double car_behind[], double car_behind_speed[], double scale = 0.16, double scale_other_car = 1.2,
			 double min_dist_with_car = 3.0, double min_dist_with_slower_car = 3.5, const double min_dist = 20)
{
  const double min_maintained_distance = min_dist_with_car + scale*mySpeed;
  return car_behind[checklane] == -1 || (car_behind[checklane] > min_maintained_distance && car_behind_speed[checklane]*scale_other_car <= mySpeed) || 
    car_behind[checklane] > min_dist || (car_behind[checklane] > min_dist_with_slower_car && car_behind_speed[checklane] <= mySpeed);
}

bool safeToProceedLane(int lane, double speed, double car_ahead[], double scale = 1, double min_dist = 35)
{
  return (car_ahead[lane] == -1 || car_ahead[lane] > min_dist);
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  double prev_s_val = -0.1;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
	
	if(s > prev_s_val){
	  prev_s_val = s;
	  map_waypoints_x.push_back(x);
	  map_waypoints_y.push_back(y);
	  map_waypoints_s.push_back(s);
	  map_waypoints_dx.push_back(d_x);
	  map_waypoints_dy.push_back(d_y);  
	}
	else {
	  break;
	}
  }

  // start from lane 1
  int lane = 1;

  // reference velocity to the target speed
  double ref_vel = 0;

  int change_lane_duration = 0;

  // Acceleration is calculated by comparing the rate of change of average speed over .2 second intervals. 
  // In this case total acceleration at one point was as high as 75 m/s^2. Jerk was also very high. 
  // The jerk is calculated as the average acceleration over 1 second intervals. 
  // In order for the passenger to have an enjoyable ride both the jerk and total acceleration should not reach values over 10. 


  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, 
	       &lane, &ref_vel, &change_lane_duration](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

		change_lane_duration++;

		bool too_close = false;  
		double car_speed_ms = mph2ms(car_speed);
		const double topspeed = mph2ms(49.70);
		double smooth_factor = 1.0;

		// If I am driving 20 m/s, I will need 20 meters of distance to stop completely at a rate of 10m/s^2 deceleration
		// Multiple by 1.2 for the lookahead distance to preserve 20% additional distance to avoid being too close to the car in front.
		double lookahead_distance = car_speed_ms * 1.2; 

		int best_lane = lane;
		double speed_reduce = 0;

		if (change_lane_duration > 150){

		  double car_ahead[3]; // Distance of the car ahead of me for each lane. -1 mean there is no car at the front
		  double car_behind[3]; // Distance of the car behind of me for each lane. -1 mean there is no car at the front
		  double car_behind_speed[3]; // Speed of the car behind my car.

		  car_ahead[0] = -1;
		  car_ahead[1] = -1;
		  car_ahead[2] = -1;
		  
		  car_behind[0] = -1;
		  car_behind[1] = -1;
		  car_behind[2] = -1;

		  for(int i=0; i < sensor_fusion.size(); ++i){
		    float d = sensor_fusion[i][6]; // 
		    double vx = sensor_fusion[i][3];
		    double vy = sensor_fusion[i][4];
		    double check_speed = sqrt(vx*vx + vy*vy);
		    double check_car_s = sensor_fusion[i][5];
		    
		    // Cars at the front from the sensor data
		    if ( check_car_s > car_s ){
		      double car_dist = check_car_s - car_s;
		      if (car_ahead[getLane(d)] == -1){
			car_ahead[getLane(d)] = car_dist;
		      }
		      else{
			car_ahead[getLane(d)] = (car_dist < car_ahead[getLane(d)]) ? car_dist : car_ahead[getLane(d)];
		      }
		    }
		    // Cars beside or behind me from the sensor data
		    else{
		      double car_dist = car_s - check_car_s;
		      if (car_behind[getLane(d)] == -1){
			car_behind[getLane(d)] = car_dist;
			car_behind_speed[getLane(d)] = check_speed;
		      }
		      else{
			car_behind[getLane(d)] = (car_dist < car_behind[getLane(d)]) ? car_dist : car_behind[getLane(d)];
			car_behind_speed[getLane(d)] = (car_dist < car_behind[getLane(d)]) ? check_speed : car_behind_speed[getLane(d)];
		      }
		    }
		  }
		  

		  // Calculate the best lane index number
		  if (car_ahead[lane] == -1) // Maintain at current lane
		    best_lane = lane;
		  else{
		    if ((lane == 0 || lane == 2) && car_ahead[1] == -1) // If car is in lane 0, or 2, and lane 1 is -1, it makes sense to switch to lane 1
		      best_lane = 1;
		    else if (lane == 0 && car_ahead[2] == -1)
		      best_lane = 2;
		    else if (lane == 2 && car_ahead[0] == -1)
		      best_lane = 0;
		    else if (lane == 1 && car_ahead[0] == -1)
		      best_lane = 0;
		    else if (lane == 1 && car_ahead[2] == -1)
		      best_lane = 2;
		    else{
		      best_lane = (car_ahead[0] > car_ahead[1] ) ? 0 : 1;
		      best_lane = (car_ahead[best_lane] > car_ahead[2] ) ? best_lane : 2;
		    }			
		  }


		  if(car_ahead[lane] != - 1 && best_lane != lane){ // If there is car in front of me or I am not in the best lane, then I will consider changing lane
		    
		    bool toZero = false;
		    bool toOne = false;

		    switch(lane){

		    case 0:
		      if( best_lane == 1 && noCarBehindThisLane(1, car_speed_ms, car_behind, car_behind_speed)){
			lane = 1;
			change_lane_duration = 0;
		      }
		      // this scenario is if lane 2 may made more sense. for instance, lane 1 may slightly slow my car down. But entering lane 1 allow me to
		      // move to lane two, which is the best possible lane at current time, once lane become one, assuming the condition remain the same in the
		      // next few seconds, it will switch to lane 2
		      
		      else if (best_lane == 2 && (safeToProceedLane(1, car_speed_ms, car_ahead) || car_ahead[1] > car_ahead[0])&& 
			       noCarBehindThisLane(1, car_speed_ms, car_behind, car_behind_speed) &&
			       noCarBehindThisLane(2, car_speed_ms, car_behind, car_behind_speed)){
			lane = 1;
			change_lane_duration = 0;
			}
		      break;

		    case 1:
		      toZero = false;
		      toOne = false;
		      if((car_ahead[0] > car_ahead[1] || car_ahead[0] == -1) && 
			 noCarBehindThisLane(0, car_speed_ms, car_behind, car_behind_speed)
			 ){
			toZero = true;
			change_lane_duration = 0;
		      }
		      if((car_ahead[2] > car_ahead[1] || car_ahead[2] == -1) && 
			 noCarBehindThisLane(2, car_speed_ms, car_behind, car_behind_speed)
			 ){
			toOne = true;
			change_lane_duration = 0;
		      }
		      if (toZero && toOne){
			// if lane zero and lane two are both good, we select the best one.
			lane = best_lane;
		      }
		      else if (toZero)
			lane = 0;
		      else if (toOne)
			lane = 2;
		      break;
		    case 2:
		      if( best_lane == 1 && noCarBehindThisLane(1, car_speed_ms, car_behind, car_behind_speed)){
			lane = 1;
			change_lane_duration = 0;
		      }
		      // this scenario is if lane 0 may made more sense. for instance, lane 1 may slightly slow my car down. But entering lane 1 allow me to
		      // move to lane zero, which is the best possible lane at current time, once lane become one, assuming the condition remain the same in the
		      // next few seconds, it will switch to lane 0
		      
		      else if (best_lane == 0 && (safeToProceedLane(1, car_speed_ms, car_ahead) || car_ahead[1] > car_ahead[2]) && 
			       noCarBehindThisLane(1, car_speed_ms, car_behind, car_behind_speed) &&
			       noCarBehindThisLane(0, car_speed_ms, car_behind, car_behind_speed)){
			lane = 1;
			change_lane_duration = 0;
			}
		      break;
		    }
		  }
		}

		// Looking into the sersor information and check the car in front of me
		// If the car is slow, I better drive slow.

		double additional_dist = 0;

		if (abs(lane - best_lane) == 2)
			additional_dist = 3;

		for(int i=0; i < sensor_fusion.size(); ++i){
		  // car is in my lane
		  float d = sensor_fusion[i][6]; // 
		  double vx = sensor_fusion[i][3];
		  double vy = sensor_fusion[i][4];
		  double check_speed = sqrt(vx*vx + vy*vy);
		  double check_car_s = sensor_fusion[i][5];

		  if(getLane(d) == lane){

		    // There is a car in front of me that is very close, I either maintain the same speed or slow down
		    if((check_car_s > car_s) && ((check_car_s-car_s) < lookahead_distance*1.1)){
		      too_close = true;
		      
		      // If the car is really close and slower than my current speed, then I really need to slow down
		      if (check_speed < car_speed_ms && ((check_car_s-car_s) < lookahead_distance + additional_dist)){
			speed_reduce = 0.9;
			// Smooth is used in case my car is following another car and I start to accelerate again
			// I want to accelerate slowly in case I get too close and need to slow down again, this help made the speed more stable
			smooth_factor = 0.1;
		      }
		    }
		  }
		}

		// I just change lane, I can drive faster
		if (change_lane_duration == 0)
		  smooth_factor = 1;
		  
		if (too_close){
		    ref_vel = car_speed_ms - speed_reduce;
		}
		else if (car_speed_ms < topspeed-1.5){ // Trying to accelerate at 1m/100ms, which is 10m/s
		  ref_vel = car_speed_ms + smooth_factor; // smooth_factor is 1 most of the time. we update once every 100ms. therefore the acceleration is 10m/s^2
		  if (smooth_factor < 1.0)
		    smooth_factor += 0.1; // too smooth out the speed of the car in case it is following another driving car.
		}
		else if(car_speed < topspeed - 0.5) // If I am close to top speed I can allow (50mph), I want to slow down my acceleration.
		  ref_vel = car_speed_ms + 0.3;
		else if (car_speed_ms < topspeed)
		  ref_vel = car_speed_ms + 0.07;


		int prev_size = previous_path_x.size();

		// reference x, y, yaw states
		// either we will reference the starting point as where the car is or at the previous paths end point
		double ref_x = car_x;
		double ref_y = car_y;
		double ref_yaw = deg2rad(car_yaw);

		double prev_car_x;
		double prev_car_y;

		if(prev_size < 2){
		  // use two points that make the path tangent to the car
		  prev_car_x = car_x - cos(ref_yaw);
		  prev_car_y = car_y - sin(ref_yaw);
		}
		else{
		  // redefine reference state as previous path end point

		  // just preserve the first 0.1 seconds for our new path. Yes, we only keep the 0.1 second of the previous planned path.
		  // and since we increase at 1m/s for our acceleration, and update every 0.1 seconds, then we are doing a 10m/s^2 acceleration.
		  if (prev_size > 5 )
		    prev_size = 5;

		  ref_x = previous_path_x[prev_size-1];
		  ref_y = previous_path_y[prev_size-1];
  
		  prev_car_x = previous_path_x[prev_size - 2];
		  prev_car_y= previous_path_y[prev_size - 2];
		  ref_yaw = atan2(ref_y - prev_car_y, ref_x - prev_car_x);
		}

		tk::spline s = getSmoothCurvePathPlan(prev_car_x, ref_x, prev_car_y, ref_y, ref_yaw, car_s, lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

		// define the planned path
		vector<double> next_x_vals;
		vector<double> next_y_vals;

		// start with all the previous path point from up to the first 0.1 seconds. Yes, we only keep the 0.1 second of the previous planned path.
		// and since we increase at 1m/s for our acceleration, and update every 0.1 seconds, then we are doing a 10m/s^2 acceleration.
		for(int i=0; i<prev_size; ++i) {
		  next_x_vals.push_back(previous_path_x[i]);
		  next_y_vals.push_back(previous_path_y[i]);
		}

		double target_x = 25.0; // fastest we can go is 50 mph, which is 22.352m/s. Looking ahead 25 meters is more than enough
		double target_y = s(target_x);
		double target_dist = sqrt(target_x*target_x + target_y*target_y);

		double add_on = 0;
		double N = target_dist/(0.02*ref_vel);

		for(int i=0; i< 50-prev_size; ++i) {

		  double x_point = add_on + target_x/N;
		  double y_point = s(x_point);

		  add_on = x_point;

		  // rotate back to normal after rotating it eariler
		  double tmp_x;
		  double tmp_y;
		  
		  // positive to rotate clockwise back to normal
		  rotatePoint(x_point, y_point, ref_yaw, tmp_x, tmp_y);
		  
		  // shift back to original
		  tmp_x += ref_x;
		  tmp_y += ref_y;

		  next_x_vals.push_back(tmp_x);
		  next_y_vals.push_back(tmp_y);
		}

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
