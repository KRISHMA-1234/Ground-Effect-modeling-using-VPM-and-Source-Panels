#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "kernel.h"
#include "exafmm.h"
#include "surface.hpp"
using namespace exafmm;
#include <fstream>
//lifting line is along x axis, centre at 0,0 and a particular height from the ground (z axis)
class blade {
 public:
  const int number_of_panels=50;
  const double airDensity = 1.225; // kg/m^3
  double wingSpan ;
  double e=1;
  double panel_center[50];
  double panel_edge[51];
  double aoa[50]={0.0};
  double total_lift = 0;
  double chord_length[50]={0.0};
  double lift[50]={0.0};
  double circulation[50] = {0.0};
  Eigen::Vector3d leading_edge[50];
  Eigen::Vector3d trailing_edge[50];
  double C_l[50];
  double height;
  double store_value;
  Eigen::Vector3d freestream;
  double prev_circulation[50] = {0.0};
  Eigen::Vector3d velocityatcenter[50];
  Eigen::Vector3d velocityatedge[51];
  Eigen::Vector3d chordvector[50];
  double drag = 0;
  double set_angle;
  double induced=0;
  double AR = 2.5464790894703255;
  double area = 1/2.5464790894703255;
 int amplitude = 10; // degrees
int time_steps = 25; // number of time steps for half a period
int full_period = time_steps * 2 ; // full period of the sine wave



//double heights[50];
double prevheight=0;


  double get_cl(double aoa){
   
       return 2*M_PI*aoa;
    
    
  }
  void setpanelcenter(){
    for(int i=0;i<number_of_panels;i++){
      panel_center[i] = (i+0.5)*(wingSpan/number_of_panels);
    }
  
  }

  void setpaneledge(){
    for(int i =0; i<number_of_panels+1;i++){
        panel_edge[i] = i*(wingSpan/number_of_panels);
    }
  }
  // initializing the initial aoa(angle of the wing with the free stream, considering no wing twist)
 /* for(int i=0;i<100;i++){
      aoa[i] = 5; //in degrees 
      circulation[i]=0;
      prev_circulation[i]=0;

    }*/
 void calculatechordlength(){
    for (int i=0;i<number_of_panels;i++){
        double c_root=0.5;
        //chord_length[i]=0.5;  //box wing
        //chord_length[i]= c_root * std::sqrt(1 - std::pow((2 * panel_center[i]) / 2*wingSpan, 2)) ; //elliptic wing 
        chord_length[i] = c_root* std::sqrt(1 - std::pow((2 * panel_center[i] - wingSpan) / wingSpan, 2));
    }
 }
 void calculatevelocityatcenter(Bodies bodies){
     for (int i=0;i<number_of_panels;i++){
        velocityatcenter[i]= freestream + P2cell( panel_center[i],0,height , bodies) ;//not added source panels for ground effect yet
        }
 }
 void calculatevelocityatcenterpanel(Bodies bodies, Surface surface){
     for (int i=0;i<number_of_panels;i++){
        Eigen::Vector3d panel_vel= surface.calculateVelocityAtPointcoord(panel_center[i],0,height );
        velocityatcenter[i]= freestream + P2cell( panel_center[i],0,height , bodies) + panel_vel;// added source panels for ground effect yet
        }
 }
 void calculatevelocityatedge(Bodies bodies){
    for (int i=0;i<number_of_panels+1;i++){
        velocityatedge[i]= freestream + P2cell( panel_edge[i],0,height , bodies); //not added source panels for ground effect yet
        }
 }
  void calculatevelocityatedgepanel(Bodies bodies, Surface surface){
    for (int i=0;i<number_of_panels+1;i++){
        Eigen::Vector3d panel_vel= surface.calculateVelocityAtPointcoord(panel_edge[i],0,height );
        velocityatedge[i]= freestream + P2cell( panel_edge[i],0,height , bodies) + panel_vel; // added source panels for ground effect yet
        }
 }

///////////////////////////////////////////////////////////////////////
 void calculatevelocityatcenterheave(Bodies bodies, double addvel){
  Eigen::Vector3d heavvel(0, 0, addvel);
     for (int i=0;i<number_of_panels;i++){
        velocityatcenter[i]= freestream + P2cell( panel_center[i],0,height , bodies) + heavvel;//not added source panels for ground effect yet
        }
 }
 void calculatevelocityatcenterpanelheave(Bodies bodies, Surface surface, double addvel){
    Eigen::Vector3d heavvel(0, 0, addvel);
     for (int i=0;i<number_of_panels;i++){
        Eigen::Vector3d panel_vel= surface.calculateVelocityAtPointcoord(panel_center[i],0,height );
        velocityatcenter[i]= freestream + P2cell( panel_center[i],0,height , bodies) + panel_vel + heavvel;// added source panels for ground effect yet
        }
 }

 void calculatevelocityatedgeheave(Bodies bodies, double addvel){
   Eigen::Vector3d heavvel(0, 0, addvel);
    for (int i=0;i<number_of_panels+1;i++){
        velocityatedge[i]= freestream + P2cell( panel_edge[i],0,height , bodies) + heavvel; //not added source panels for ground effect yet
        }
 }
  void calculatevelocityatedgepanelheave(Bodies bodies, Surface surface, double addvel){
      Eigen::Vector3d heavvel(0, 0, addvel);
    for (int i=0;i<number_of_panels+1;i++){
        Eigen::Vector3d panel_vel= surface.calculateVelocityAtPointcoord(panel_edge[i],0,height );
        velocityatedge[i]= freestream + P2cell( panel_edge[i],0,height , bodies) + panel_vel + heavvel; // added source panels for ground effect yet
        }
 }

////////////////////////////////////////////////////////////////////////////////////
void setleadingtrailingedge(){
    for(int i = 0; i<number_of_panels; i++){
        leading_edge[i][0]=panel_center[i];
        trailing_edge[i][0]=panel_center[i];
        leading_edge[i][1]=  -(0.25)*chord_length[i]*cos(set_angle);
        trailing_edge[i][1]=  (0.75)*chord_length[i]*cos(set_angle);
        leading_edge[i][2]= height + (0.25)*chord_length[i]*sin(set_angle);
        trailing_edge[i][2]= height - (0.75)*chord_length[i]*sin(set_angle);
    }
}
void setchordvector(){
    for(int i = 0; i<number_of_panels; i++){
      chordvector[i]=leading_edge[i]-trailing_edge[i];

    }
}
 //void initializeblade();
 void calculateaoa(){
    for (int i=0;i<number_of_panels;i++){
        //if(j==0){
            //aoa[i]=set_angle;
        //}
        //else {
      Eigen::Vector3d chordNormalized = chordvector[i].normalized();
      //std::cout<<"normalized vector is "<< chordNormalized.transpose()<<std::endl;
      double dotProduct = chordNormalized.dot(velocityatcenter[i]);
      //std::cout<<"dot product is "<< dotProduct<<std::endl;
      Eigen::Vector3d velocityChordProjection = dotProduct * chordNormalized;
      Eigen::Vector3d velocityPerpendicular = velocityatcenter[i] - velocityChordProjection;
      Eigen::Vector3d ProjectedVector ;
      double normchord = std::sqrt(chordvector[i][0]*chordvector[i][0] + chordvector[i][1]*chordvector[i][1] + chordvector[i][2]*chordvector[i][2]);
      //double normproject = ProjectedVector[0]*ProjectedVector[0] + ProjectedVector[1]*ProjectedVector[1] + ProjectedVector[2]*ProjectedVector[2];
      ProjectedVector<<0,velocityatcenter[i][1],velocityatcenter[i][2];
      double normproject = std::sqrt(ProjectedVector[0]*ProjectedVector[0] + ProjectedVector[1]*ProjectedVector[1] + ProjectedVector[2]*ProjectedVector[2]);
      if (set_angle == 0){
        if(velocityatcenter[i][2]<0){
          aoa[i] = -(M_PI - std::acos((chordvector[i].dot(ProjectedVector))/(normchord*normproject)));
        }
        else if (velocityatcenter[i][2]>=0){
            aoa[i] = M_PI - std::acos((chordvector[i].dot(ProjectedVector))/(normchord*normproject));
        }
      }
      else if (set_angle>0){
        aoa[i] = M_PI - std::acos((chordvector[i].dot(ProjectedVector))/(normchord*normproject));
      }
      else{
        aoa[i] = -(M_PI - std::acos((chordvector[i].dot(ProjectedVector))/(normchord*normproject)));
      }
      
      //aoa[i] = std::atan2(velocityPerpendicular.z(), velocityChordProjection.norm());
       // }
   } 
 } //???????????
 void calculateliftinpanel(){
    for (int i=0;i<number_of_panels;i++){
        double cl= get_cl(aoa[i]);
    lift[i] = 0.5*airDensity*((velocityatcenter[i]).norm()*(velocityatcenter[i]).norm())*(wingSpan/number_of_panels)*chord_length[i]*cl ;
    //lift[i] = 0.5*airDensity*((freestream).norm()*(freestream).norm())*(wingSpan/number_of_panels)*chord_length[i]*cl ;//??????????
    }
 }
 void calculatecirculationinpanel(){
 for (int i=0;i<number_of_panels;i++){
   prev_circulation[i]= circulation[i];
   double cl= get_cl(aoa[i]);
    circulation[i] = (velocityatcenter[i].norm())*chord_length[i]*cl/2;
   //circulation[i] = lift[i]/(airDensity*(velocityatcenter[i].norm()));
 }
}

double calculatetotallift(){
    total_lift=0;
    for (int i=0;i<number_of_panels;i++){
        total_lift += lift[i];
    }
    return total_lift;
}
double getdrag(){
    for(int i=0;i<number_of_panels;i++){
      double cl=get_cl(aoa[i]);
      double AR = 2.5464790894703255;
     drag += 0.5*airDensity*(velocityatcenter[i].norm())*(velocityatcenter[i].norm())*chord_length[i]*(wingSpan/number_of_panels)*(cl*cl/M_PI*e*AR);
    }
    return drag;
  }
double getinduced(){
    for(int i=0;i<number_of_panels ;i++){
      induced += -(airDensity*velocityatcenter[i][2]*circulation[i]*wingSpan/number_of_panels);
    }
    return induced;
  }
//std::vector<double> pitchaoa = {0, 2, 5, 10,15,20,25,30,35,40,45,50,55,60,65,60,55,50,45,40,35,30,25,20,15,10, 5, 2, 0, -2, -5, -10,-15,-20,-25,-30,-35,-40,-45,-50,-55,-60,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10, -5, -2, 0};
std::vector<double> pitchaoa = {2,5,10,5,2,0,-2,-5,-10,-5,-2,0};
//std::vector<double> pitchaoa = {5};
double sinosoidangle(int t){


//t = np.arange(0, full_period)  # time steps
double setangle;
if (t == full_period){
   setangle =0;
}
else{
  setangle = amplitude * std::sin(2 * M_PI * t / full_period);
}

return setangle;
}
   
    
int total_points = 50;

  const int heave_total_steps = 52;        // Total number of time steps
const double heave_amplitude = 0.05;      // Amplitude of the heaving motion (max height)
const double heave_period = 0.5;         // Period of one full oscillation (in arbitrary time units)
const double heave_omega = 2 * M_PI / heave_period; // Angular frequency (radians per time unit)
const double heave_time_step = heave_period / heave_total_steps; // Time step size

// Function to calculate the position of the wing at a given time
double calculateWingPosition_heave(double time) {
    // Sinusoidal motion: h(t) = H * sin(omega * t)
    return heave_amplitude * sin(heave_omega * time);
}

// Function to calculate the velocity of the wing at a given time
double calculateWingVelocity_heave(double time) {
    // Velocity: v(t) = H * omega * cos(omega * t)
    return heave_amplitude * heave_omega * cos(heave_omega * time);
}


}

;