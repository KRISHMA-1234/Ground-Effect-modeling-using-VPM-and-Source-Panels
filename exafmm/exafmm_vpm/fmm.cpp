#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
// #if EXAFMM_EAGER
#include "traverse_eager.h"
#include "surface.hpp"
// #include "thirdparty/eigen"  // ???????
#include <Eigen/Core>
#include "rapidobj/rapidobj.hpp" // ???????
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include "blade.h"
#include <vector>
// #elif EXAFMM_LAZY
// #include "traverse_lazy.h"
// #endif
using namespace exafmm;
#include <fstream>

int main(int argc, char **argv)
{
   Surface surface;
   surface.initFromOBJFile("../translatedground1.obj");
   
   blade blade;
   blade.wingSpan = 1;
   blade.height = 0;
   blade.set_angle = 5*M_PI/180; //degrees
   blade.setpanelcenter(); //initilalzing wing geometry 
   blade.setpaneledge();
  
   for (int i =0; i<blade.number_of_panels;i++){
    std::cout<<"the panel center is"<<blade.panel_center[i]<<std::endl;
   }
   for (int i =0; i<blade.number_of_panels +1;i++){
    std::cout<<"the panel edge is"<<blade.panel_edge[i]<<std::endl;
   }
   double number = blade.wingSpan/blade.number_of_panels;
   std::cout<<number<<std::endl;
   blade.calculatechordlength();
   blade.setleadingtrailingedge();
   blade.setchordvector();
   
   for (int i =0; i<blade.number_of_panels;i++){
    std::cout<<"chord length "<<i<<"is "<< blade.chord_length[i]<<std::endl;
    //std::cout<<"chord vector "<<i<<"is "<< blade.chordvector[i].transpose()<<std::endl;
    //std::cout<<"leading_edge vector "<<i<<"is "<< blade.leading_edge[i].transpose()<<std::endl;
    //std::cout<<"trailing_edge vector "<<i<<"is "<< blade.trailing_edge[i].transpose()<<std::endl;
   }
   blade.freestream<<0,1,0;
   blade.calculatechordlength();
   const int                    npanels = surface.getCellCenters().size();
  auto cell_centers= surface.getCellCenters();
  

  const int numBodies = 1800;  // Number of bodies
  const int numTargets = 1'000; // Number of targets for checking answer
  P = 4;                       // Order of expansions
  ncrit = 64;                   // Number of bodies per leaf cell
  theta = 0.4;                  // Multipole acceptance criterion

  printf("--- %-16s ------------\n", "FMM Profiling"); // Start profiling
  //! Initialize bodies

  //Bodies bodies(numBodies); // Initialize bodies
  Bodies bodies;
  Bodies sudobodies;
  Bodies wakebodies;
    auto writeTovtk = [&bodies](int step) {
   std::ofstream file;
file.open("output" + std::to_string(step) + ".vtk");
file << "# vtk DataFile Version 3.0\n";
file << "vtk output\n";
file << "ASCII\n";
file << "DATASET POLYDATA\n";
file << "POINTS " << bodies.size() << " float\n";
for (size_t b = 0; b < bodies.size(); b++)
{
  file << bodies[b].X[0] << " " << bodies[b].X[1] << " " << bodies[b].X[2] << "\n";
}
// Add points as vertices
file << "VERTICES " << bodies.size() << " " << 2 * bodies.size() << "\n";
for (size_t b = 0; b < bodies.size(); b++)
{
  file << "1 " << b << "\n";
}
file << "POINT_DATA " << bodies.size() << "\n";

// Add alpha vectors
file << "VECTORS alpha float\n";
for (const auto &body : bodies)
{
  file << body.alpha[0] << " " << body.alpha[1] << " " << body.alpha[2] << "\n";
}

// Add velocity vectors
file << "VECTORS velocity float\n";
for (const auto &body : bodies)
{
  file << body.velocity[0] << " " << body.velocity[1] << " " << body.velocity[2] << "\n";
}

// Add scalar radius values
file << "SCALARS radius float 1\n";
file << "LOOKUP_TABLE default\n";
for (const auto &body : bodies)
{
  file << body.radius << "\n";
}

file.close();
  };

  auto writeTovtkwake = [&wakebodies](int step) {
    std::ofstream file;
    file.open("output" + std::to_string(step) + ".vtk");
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << wakebodies.size() << " float\n";
    for (size_t b = 0; b < wakebodies.size(); b++) {
      file << wakebodies[b].X[0] << " " << wakebodies[b].X[1] << " " << wakebodies[b].X[2] << "\n";
    }
    // Add points as vertices
    file << "VERTICES " << wakebodies.size() << " " << 2 * wakebodies.size() << "\n";
    for (size_t b = 0; b < wakebodies.size(); b++) {
      file << "1 " << b << "\n";
    }
    file << "POINT_DATA " << wakebodies.size() << "\n";

    // Add alpha vectors
    file << "VECTORS alpha float\n";
    for (const auto &wakebody : wakebodies) {
      file << wakebody.alpha[0] << " " << wakebody.alpha[1] << " " << wakebody.alpha[2] << "\n";
    }

    // Add velocity vectors
    file << "VECTORS velocity float\n";
    for (const auto &wakebody : wakebodies) {
      file << wakebody.velocity[0] << " " << wakebody.velocity[1] << " " << wakebody.velocity[2] << "\n";
    }
    file.close();
  };


  auto advance_euler = [&bodies, &surface , &blade](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      //std::cout<<"The position at beggining of advancing body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
      //std::cout<<"The velocity addition due to particles on body "<<b<<"is"<<bodies[b].velocity[0]<<" "<<bodies[b].velocity[1]<<" "<<bodies[b].velocity[2]<<std::endl;
      //std::cout<<"The alpha at beggining of advancing body "<<b<<"is"<<bodies[b].alpha[0]<<" "<<bodies[b].alpha[1]<<" "<<bodies[b].alpha[2]<<std::endl;
      //std::cout<<"The dadt addition due to particles on body "<<b<<"is"<<bodies[b].dadt[0]<<" "<<bodies[b].dadt[1]<<" "<<bodies[b].dadt[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {
        
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        bodies[b].X[d] += bodies[b].velocity[d] * dt +blade.freestream[d]*dt;
        bodies[b].prev_velocity[d] = bodies[b].velocity[d]  ;
        
        bodies[b].velocity[d] = 0;                    //  Clear velocity
        
        bodies[b].alpha[d] += bodies[b].dadt[d] * dt;
        bodies[b].prev_dadt[d] = bodies[b].dadt[d] ;
        bodies[b].dadt[d] = 0 ;
      }
      //std::cout<<"The position addition due to particles on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };
  auto advance_euler_wake = [&wakebodies, &surface , &blade](real_t dt)
  {
    for (size_t b = 0; b < wakebodies.size(); b++)
    { // Loop over bodies
      //std::cout<<"The position at beggining of advancing body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
      //std::cout<<"The velocity addition due to particles on body "<<b<<"is"<<bodies[b].velocity[0]<<" "<<bodies[b].velocity[1]<<" "<<bodies[b].velocity[2]<<std::endl;
      //std::cout<<"The alpha at beggining of advancing body "<<b<<"is"<<bodies[b].alpha[0]<<" "<<bodies[b].alpha[1]<<" "<<bodies[b].alpha[2]<<std::endl;
      //std::cout<<"The dadt addition due to particles on body "<<b<<"is"<<bodies[b].dadt[0]<<" "<<bodies[b].dadt[1]<<" "<<bodies[b].dadt[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {
        
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        wakebodies[b].X[d] += wakebodies[b].velocity[d] * dt +blade.freestream[d]*dt;
        wakebodies[b].prev_velocity[d] = wakebodies[b].velocity[d]  ;
        
        wakebodies[b].velocity[d] = 0;                    //  Clear velocity
        
        wakebodies[b].alpha[d] += wakebodies[b].dadt[d] * dt;
        wakebodies[b].prev_dadt[d] = wakebodies[b].dadt[d] ;
        wakebodies[b].dadt[d] = 0 ;
      }
      //std::cout<<"The position addition due to particles on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };
auto advance_euler_panel = [&bodies, &surface](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      Eigen::Vector3d extra_vel;
      extra_vel = surface.calculateVelocityAtPointcoord(bodies[b].X[0],bodies[b].X[1], bodies[b].X[2]);
      bodies[b].vel_panel[0] = extra_vel[0];
      bodies[b].vel_panel[1] = extra_vel[1];
      bodies[b].vel_panel[2] = extra_vel[2];
      Eigen::Vector3d stretch_panel; 
      stretch_panel = surface.getstretchingterm(bodies[b].X[0],bodies[b].X[1],bodies[b].X[2],bodies[b].alpha[0],bodies[b].alpha[1],bodies[b].alpha[2]);
      bodies[b].dadt_panel[0] = stretch_panel[0];
      bodies[b].dadt_panel[1] = stretch_panel[1];
      bodies[b].dadt_panel[2] = stretch_panel[2];
      //std::cout<<"The position at beggining of advancing body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
      //std::cout<<"The velocity addition due to panel on body "<<b<<"is"<<bodies[b].vel_panel[0]<<" "<<bodies[b].vel_panel[1]<<" "<<bodies[b].vel_panel[2]<<std::endl;
      //std::cout<<"The alpha at beggining of advancing body "<<b<<"is"<<bodies[b].alpha[0]<<" "<<bodies[b].alpha[1]<<" "<<bodies[b].alpha[2]<<std::endl;
      //std::cout<<"The dadt addition due to panel on body "<<b<<"is"<<bodies[b].dadt_panel[0]<<" "<<bodies[b].dadt_panel[1]<<" "<<bodies[b].dadt_panel[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {
        
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        bodies[b].X[d] += bodies[b].vel_panel[d] * dt;
        bodies[b].prev_velocity_panel[d] =  bodies[b].vel_panel[d];
        
        bodies[b].vel_panel[d] = 0;        
        
        bodies[b].alpha[d] += bodies[b].dadt_panel[d] * dt;
        bodies[b].prev_dadt_panel[d] =  bodies[b].dadt_panel[d];
        bodies[b].dadt_panel[d] = 0;                 //  Clear velocity
        
      }
      //std::cout<<"The position addition due to panel on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };

  auto advance_euler_panel_nostretch = [&bodies, &surface](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      Eigen::Vector3d extra_vel;
      extra_vel = surface.calculateVelocityAtPointcoord(bodies[b].X[0],bodies[b].X[1], bodies[b].X[2]);
      bodies[b].vel_panel[0] = extra_vel[0];
      bodies[b].vel_panel[1] = extra_vel[1];
      bodies[b].vel_panel[2] = extra_vel[2];
      //Eigen::Vector3d stretch_panel; 
      //stretch_panel = surface.getstretchingterm(bodies[b].X[0],bodies[b].X[1],bodies[b].X[2],bodies[b].alpha[0],bodies[b].alpha[1],bodies[b].alpha[2]);
      //bodies[b].dadt_panel[0] = stretch_panel[0];
      //bodies[b].dadt_panel[1] = stretch_panel[1];
      //bodies[b].dadt_panel[2] = stretch_panel[2];
      //std::cout<<"The position at beggining of advancing body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
      //std::cout<<"The velocity addition due to panel on body "<<b<<"is"<<bodies[b].vel_panel[0]<<" "<<bodies[b].vel_panel[1]<<" "<<bodies[b].vel_panel[2]<<std::endl;
      //std::cout<<"The alpha at beggining of advancing body "<<b<<"is"<<bodies[b].alpha[0]<<" "<<bodies[b].alpha[1]<<" "<<bodies[b].alpha[2]<<std::endl;
      //std::cout<<"The dadt addition due to panel on body "<<b<<"is"<<bodies[b].dadt_panel[0]<<" "<<bodies[b].dadt_panel[1]<<" "<<bodies[b].dadt_panel[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {
        
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        bodies[b].X[d] += bodies[b].vel_panel[d] * dt;
        bodies[b].prev_velocity_panel[d] =  bodies[b].vel_panel[d];
        
        bodies[b].vel_panel[d] = 0;        
        
        //bodies[b].alpha[d] += bodies[b].dadt_panel[d] * dt;
        //bodies[b].prev_dadt_panel[d] =  bodies[b].dadt_panel[d];
        //bodies[b].dadt_panel[d] = 0;                 //  Clear velocity
        
      }
      //std::cout<<"The position addition due to panel on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };

  //auto calculatepanelstretching = [&bodies, &surface]{
     // Eigen::Vector3d stretch_panel; 
     // stretch_panel = surface.getstretchingterm(bodies.X[0],bodies.X[1],bodies.X[2],bodies.alpha[0],bodies.alpha[1],bodies.alpha[2]);
      
      
   // }

    auto advance_adambash_panel = [&bodies, &surface](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      Eigen::Vector3d extra_vel;
      extra_vel = surface.calculateVelocityAtPointcoord(bodies[b].X[0],bodies[b].X[1], bodies[b].X[2]);
      Eigen::Vector3d stretch_panel; 
      stretch_panel = surface.getstretchingterm(bodies[b].X[0],bodies[b].X[1],bodies[b].X[2],bodies[b].alpha[0],bodies[b].alpha[1],bodies[b].alpha[2]);
      
      //std::cout<<"velocity due to panels is"<<extra_vel.transpose()<<std::endl;
      bodies[b].vel_panel[0] = extra_vel[0];
      bodies[b].vel_panel[1] = extra_vel[1];
      bodies[b].vel_panel[2] = extra_vel[2];
      bodies[b].dadt_panel[0] = stretch_panel[0];
      bodies[b].dadt_panel[1] = stretch_panel[1];
      bodies[b].dadt_panel[2] = stretch_panel[2];
      //std::cout<<"The velocity addition due to panel on body "<<b<<"is"<<bodies[b].vel_panel[0]<<" "<<bodies[b].vel_panel[1]<<" "<<bodies[b].vel_panel[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        bodies[b].X[d] += 3*dt*bodies[b].vel_panel[d]/2 - dt*bodies[b].prev_velocity_panel[d]/2 ;
        bodies[b].prev_velocity_panel[d] = bodies[b].vel_panel[d];
        //std::cout<<"The velocity addition due to panel on body "<<<b<<"is"<<bodies[b].velocity[0]<<" "<<bodies[b].velocity[1]<<" "bodies[b].velocity[2]<<std::endl;
        bodies[b].vel_panel[d] = 0;                    //  Clear velocity
        
        bodies[b].alpha[d] += 3*dt*bodies[b].dadt_panel[d]/2 - dt*bodies[b].prev_dadt_panel[d]/2 ;
        bodies[b].prev_dadt_panel[d] = bodies[b].dadt_panel[d] ;
        bodies[b].dadt_panel[d] = 0;
      
      }
      //std::cout<<"The position addition due to panel on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };

   auto advance_adambash = [&bodies, &surface](real_t dt)
  {
    for (size_t b = 0; b < bodies.size(); b++)
    { // Loop over bodies
      //std::cout<<"The velocity addition due to particles on body "<<b<<"is"<<bodies[b].velocity[0]<<" "<<bodies[b].velocity[1]<<" "<<bodies[b].velocity[2]<<std::endl;
      for (int d = 0; d < 3; d++)
      {  
        //bodies[b].X[d] += bodies[b].velocity[d] * dt + vel_add_panel[d] * dt; //  Update position
        bodies[b].X[d] += 3*dt*bodies[b].velocity[d]/2 - dt*bodies[b].prev_velocity[d]/2 ;
        bodies[b].prev_velocity[d] = bodies[b].velocity[d] ;
        
        bodies[b].velocity[d] = 0;                    //  Clear velocity
       
        bodies[b].alpha[d] += 3*dt*bodies[b].dadt[d]/2 - dt*bodies[b].prev_dadt[d]/2 ;
        bodies[b].prev_dadt[d] = bodies[b].dadt[d] ;
        bodies[b].dadt[d] = 0;
      }
       //std::cout<<"The position addition due to particles on body "<<b<<"is"<<bodies[b].X[0]<<" "<<bodies[b].X[1]<<" "<<bodies[b].X[2]<<std::endl;
    } // End loop over bodies
  };

    
  // auto store_prev_vel = [&bodies, &surface](){
   // Eigen::Vector3d vel_prev;
    // vel_prev[0] = bodies[]
   // }


  // adding velocity due to source panel. 
  // strength of the source panel does not change due to the vortex particle ??????
  // How to caculate the velocity effect of each vortex particle on the source panel? for no penetration condition

  //! Build tree

  //! FMM evaluation

  auto boundparticle = [&bodies,&blade](int k) {
   for (int i=0;i<blade.number_of_panels;i++){
    Body bodyi;
    bodyi.prev_dadt_panel[0]=0;
    bodyi.prev_dadt_panel[1]=0;
    bodyi.prev_dadt_panel[2]=0;
    bodyi.dadt_panel[0]=0;
    bodyi.dadt_panel[1]=0;
    bodyi.dadt_panel[2]=0;
    bodyi.prev_velocity_panel[0]=0;
    bodyi.prev_velocity_panel[1]=0;
    bodyi.prev_velocity_panel[2]=0;
    bodyi.vel_panel[0]=0;
    bodyi.vel_panel[1]=0;
    bodyi.vel_panel[2]=0;
    bodyi.velocity[0]=0;
    bodyi.velocity[1]=0;
    bodyi.velocity[2]=0;
    bodyi.dadt[0]=0;
    bodyi.dadt[1]=0;
    bodyi.dadt[2]=0;
    bodyi.prev_velocity[0]=0;
    bodyi.prev_velocity[1]=0;
    bodyi.prev_velocity[2]=0;
    bodyi.prev_dadt[0]=0;
    bodyi.prev_dadt[1]=0;
    bodyi.prev_dadt[2]=0;
    bodyi.radius = 7.9060000000e-02;
    bodyi.alpha[0]= blade.circulation[i]*blade.wingSpan/blade.number_of_panels;
    bodyi.alpha[1]= 0;
    bodyi.alpha[2]= 0;
    bodyi.X[0]=  blade.panel_center[i];//(i+0.5)*blade.wingSpan/blade.number_of_panels;
    bodyi.X[1]= 0;//since lifting line is at y = 0, the bound particles are originated from y = 0
    bodyi.X[2]= blade.height;//originate at the height of the lifting line
    if (k==0){
      bodies.push_back(bodyi);
    }
   else {
    bodies[i] = bodyi; //updating the value of the ith body; 
   }
   }
  
  };
   
  auto trailshedparticle = [&bodies,&blade,&wakebodies](double dt){
    for (int i=0;i<blade.number_of_panels +1;i++){ //NOTE THE +1 in this , STARTING FROM I=1
    if(i > 0 && i < blade.number_of_panels){
    Body bodyi;
    bodyi.prev_dadt_panel[0]=0;
    bodyi.prev_dadt_panel[1]=0;
    bodyi.prev_dadt_panel[2]=0;
    bodyi.dadt_panel[0]=0;
    bodyi.dadt_panel[1]=0;
    bodyi.dadt_panel[2]=0;
    bodyi.prev_velocity_panel[0]=0;
    bodyi.prev_velocity_panel[1]=0;
    bodyi.prev_velocity_panel[2]=0;
    bodyi.vel_panel[0]=0;
    bodyi.vel_panel[1]=0;
    bodyi.vel_panel[2]=0;
    bodyi.velocity[0]=0;
    bodyi.velocity[1]=0;
    bodyi.velocity[2]=0;
    bodyi.dadt[0]=0;
    bodyi.dadt[1]=0;
    bodyi.dadt[2]=0;
    bodyi.prev_velocity[0]=0;
    bodyi.prev_velocity[1]=0;
    bodyi.prev_velocity[2]=0;
    bodyi.prev_dadt[0]=0;
    bodyi.prev_dadt[1]=0;
    bodyi.prev_dadt[2]=0;
    bodyi.radius = 7.9060000000e-02;
    /*
    if(i>blade.number_of_panels/2){
       bodyi.alpha[0]= (-(blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][0]*dt)*-1;
    bodyi.alpha[1]= (-(blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][1]*dt)*-1;
    blade.store_value = bodyi.alpha[1];
    bodyi.alpha[2]= (-(blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][2]*dt)*-1;
    bodyi.X[0]= blade.panel_edge[i] + blade.velocityatedge[i][0]*dt*0.5 ;
    bodyi.X[1]= blade.velocityatedge[i][1]*dt*0.5;
    bodyi.X[2]=  blade.velocityatedge[i][2]*dt*0.5 + blade.height ;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
    }
    */
    /*
    else if(i==blade.number_of_panels/2){
    //bodyi.alpha[0]= (-(blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][0]*dt)*-1;
    bodyi.alpha[0]= 0;
    bodyi.alpha[1]= 0;
    bodyi.alpha[2]= 0;
    bodyi.X[0]= blade.panel_edge[i] + blade.velocityatedge[i][0]*dt*0.5 ;
    bodyi.X[1]= blade.velocityatedge[i][1]*dt*0.5;
    bodyi.X[2]=  blade.velocityatedge[i][2]*dt*0.5 + blade.height ;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
    }
    */
    //else{
     bodyi.alpha[0]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][0]*dt)*-1;
    bodyi.alpha[1]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][1]*dt)*-1;
    bodyi.alpha[2]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][2]*dt)*-1;
    //bodyi.alpha[0]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][0]*dt);
    //bodyi.alpha[1]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][1]*dt);
    //bodyi.alpha[2]= ((blade.circulation[i-1]-blade.circulation[i])*blade.velocityatedge[i][2]*dt);
    bodyi.X[0]= blade.panel_edge[i] + blade.velocityatedge[i][0]*dt*0.5 ;
    bodyi.X[1]= blade.velocityatedge[i][1]*dt*0.5;
    bodyi.X[2]=  blade.velocityatedge[i][2]*dt*0.5 + blade.height ;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
    //}
    
      
   }
   else if(i==0){
      Body bodyi;
    bodyi.prev_dadt_panel[0]=0;
    bodyi.prev_dadt_panel[1]=0;
    bodyi.prev_dadt_panel[2]=0;
    bodyi.dadt_panel[0]=0;
    bodyi.dadt_panel[1]=0;
    bodyi.dadt_panel[2]=0;
    bodyi.prev_velocity_panel[0]=0;
    bodyi.prev_velocity_panel[1]=0;
    bodyi.prev_velocity_panel[2]=0;
    bodyi.vel_panel[0]=0;
    bodyi.vel_panel[1]=0;
    bodyi.vel_panel[2]=0;
    bodyi.velocity[0]=0;
    bodyi.velocity[1]=0;
    bodyi.velocity[2]=0;
    bodyi.dadt[0]=0;
    bodyi.dadt[1]=0;
    bodyi.dadt[2]=0;
    bodyi.prev_velocity[0]=0;
    bodyi.prev_velocity[1]=0;
    bodyi.prev_velocity[2]=0;
    bodyi.prev_dadt[0]=0;
    bodyi.prev_dadt[1]=0;
    bodyi.prev_dadt[2]=0;
    bodyi.radius = 7.9060000000e-02;
    bodyi.alpha[0]= ((-blade.circulation[i])*blade.velocityatedge[i][0]*dt)*-1;
    bodyi.alpha[1]= ((-blade.circulation[i])*blade.velocityatedge[i][1]*dt)*-1;
    bodyi.alpha[2]= ((-blade.circulation[i])*blade.velocityatedge[i][2]*dt)*-1;
    //bodyi.alpha[0]= ((-blade.circulation[i])*blade.velocityatedge[i][0]*dt);
    //bodyi.alpha[1]= ((-blade.circulation[i])*blade.velocityatedge[i][1]*dt);
    //bodyi.alpha[2]= ((-blade.circulation[i])*blade.velocityatedge[i][2]*dt);
    bodyi.X[0]= blade.panel_edge[i] + blade.velocityatedge[i][0]*dt*0.5 ;
    bodyi.X[1]= blade.velocityatedge[i][1]*dt*0.5;
    bodyi.X[2]=  blade.velocityatedge[i][2]*dt*0.5 + blade.height ;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
   }
   else if(i==blade.number_of_panels){
    Body bodyi;
    bodyi.prev_dadt_panel[0]=0;
    bodyi.prev_dadt_panel[1]=0;
    bodyi.prev_dadt_panel[2]=0;
    bodyi.dadt_panel[0]=0;
    bodyi.dadt_panel[1]=0;
    bodyi.dadt_panel[2]=0;
    bodyi.prev_velocity_panel[0]=0;
    bodyi.prev_velocity_panel[1]=0;
    bodyi.prev_velocity_panel[2]=0;
    bodyi.vel_panel[0]=0;
    bodyi.vel_panel[1]=0;
    bodyi.vel_panel[2]=0;
    bodyi.velocity[0]=0;
    bodyi.velocity[1]=0;
    bodyi.velocity[2]=0;
    bodyi.dadt[0]=0;
    bodyi.dadt[1]=0;
    bodyi.dadt[2]=0;
    bodyi.prev_velocity[0]=0;
    bodyi.prev_velocity[1]=0;
    bodyi.prev_velocity[2]=0;
    bodyi.prev_dadt[0]=0;
    bodyi.prev_dadt[1]=0;
    bodyi.prev_dadt[2]=0;
    bodyi.radius = 7.9060000000e-02;
    bodyi.alpha[0]= (-(blade.circulation[i-1])*blade.velocityatedge[i][0]*dt);
    bodyi.alpha[1]= (-(blade.circulation[i-1])*blade.velocityatedge[i][1]*dt);
    bodyi.alpha[2]= (-(blade.circulation[i-1])*blade.velocityatedge[i][2]*dt);
    //bodyi.alpha[0]= (-(blade.circulation[i-1])*blade.velocityatedge[i][0]*dt)*-1;
    //bodyi.alpha[1]= (-(blade.circulation[i-1])*blade.velocityatedge[i][1]*dt)*-1;
    //bodyi.alpha[2]= (-(blade.circulation[i-1])*blade.velocityatedge[i][2]*dt)*-1;
    bodyi.X[0]= blade.panel_edge[i] + blade.velocityatedge[i][0]*dt*0.5 ;
    bodyi.X[1]= blade.velocityatedge[i][1]*dt*0.5;
    bodyi.X[2]=  blade.velocityatedge[i][2]*dt*0.5 + blade.height ;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
   }
   
    
   }
    
  };

  auto spanshedparticle = [&bodies,&blade,&wakebodies](double dt){
    for (int i=0;i<blade.number_of_panels;i++){
    Body bodyi;
    bodyi.prev_dadt_panel[0]=0;
    bodyi.prev_dadt_panel[1]=0;
    bodyi.prev_dadt_panel[2]=0;
    bodyi.dadt_panel[0]=0;
    bodyi.dadt_panel[1]=0;
    bodyi.dadt_panel[2]=0;
    bodyi.prev_velocity_panel[0]=0;
    bodyi.prev_velocity_panel[1]=0;
    bodyi.prev_velocity_panel[2]=0;
    bodyi.vel_panel[0]=0;
    bodyi.vel_panel[1]=0;
    bodyi.vel_panel[2]=0;
    bodyi.velocity[0]=0;
    bodyi.velocity[1]=0;
    bodyi.velocity[2]=0;
    bodyi.dadt[0]=0;
    bodyi.dadt[1]=0;
    bodyi.dadt[2]=0;
    bodyi.prev_velocity[0]=0;
    bodyi.prev_velocity[1]=0;
    bodyi.prev_velocity[2]=0;
    bodyi.prev_dadt[0]=0;
    bodyi.prev_dadt[1]=0;
    bodyi.prev_dadt[2]=0;
    bodyi.radius = 7.9060000000e-02;
    bodyi.alpha[0]=((blade.prev_circulation[i]-blade.circulation[i])*blade.wingSpan/blade.number_of_panels)*-1;
    //bodyi.alpha[0]=((blade.prev_circulation[i]-blade.circulation[i])*blade.wingSpan/blade.number_of_panels);
    bodyi.alpha[1]=0;
    bodyi.alpha[2]=0;
    bodyi.X[0]= blade.panel_center[i] + (blade.velocityatedge[i][0] + blade.velocityatedge[i+1][0])*0.5*dt;
    bodyi.X[1]=(blade.velocityatedge[i][1] + blade.velocityatedge[i+1][1])*0.5*dt;
    bodyi.X[2]=(blade.velocityatedge[i][2] + blade.velocityatedge[i+1][2])*0.5*dt + blade.height;
    bodies.push_back(bodyi);
    wakebodies.push_back(bodyi);
   }
  };
  auto BiotSavart = [&bodies]()
  {
    // start("Build tree");             // Start timer
    Cells cells = buildTree(bodies); // Build tree
    // stop("Build tree");              // Stop timer

    // start("P2M & M2M");           // Start timer
    initKernel();                 // Initialize kernel
    upwardPass(cells);            // Upward pass for P2M, M2M
    // stop("P2M & M2M");            // Stop timer
    // start("M2L & P2P");           // Start timer
    horizontalPass(cells, cells); // Horizontal pass for M2L, P2P
    // stop("M2L & P2P");            // Stop timer
    // start("L2L & L2P");           // Start timer
    downwardPass(cells);          // Downward pass for L2L, L2P
    // stop("L2L & L2P");            // Stop timer
  };

  // writeTovtk(0);
    auto advancewakebodies=[&wakebodies,&bodies](){
      for (int i=10;i<bodies.size();i++){
        wakebodies[i-10]=bodies[i];

      }
};
//initializing with the flow field and blade geom
 // /* comment this section for undulated surfaces
    blade.calculatevelocityatcenter(bodies);
    blade.calculatevelocityatedge(bodies);
    blade.calculateaoa();
    std::cout<<"first size of bodies array is "<<bodies.size()<<std::endl;
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    double tot_first=blade.calculatetotallift();
     //std::cout<<"the total lift in iteration 0 "<<"is "<<tot_first<<std::endl;
    for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the first circulation for panel "<< v<< "is "<<blade.circulation[v]<<std::endl;
     //std::cout<<"the first circulation for panel "<< v<< "is "<<blade.prev_circulation[v]<<std::endl;
     //std::cout<<"the first lift for panel "<< v<< "is "<<blade.lift[v]<<std::endl;
    }
   //boundparticle(0);
    for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the first central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     //std::cout<<"chord vector "<<v<<"is "<< blade.chordvector[v].transpose()<<std::endl;
     //std::cout<<"the first aoa at panel centre "<< v<<"is"<<blade.aoa[v]*180/M_PI<<std::endl;
    }
//without source panel 
 double prev_lift = tot_first;
 //  */
////////////////////////////////////////////////////////////////////////////


// /*
for (int k = 1; k < 401; k++)
   { 
    //intialize wing by making array of centres and panel edges 

    //find original clmax by initial aoa set by user
    
    
    //find original circulation 
    //boundparticle();
    blade.total_lift=0;
    trailshedparticle(0.01);
    spanshedparticle(0.01);
     std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
    std::cout<<"the size of the wake body array is "<<wakebodies.size()<<std::endl;
    //checkpoint
    Eigen::Vector3d vel1=P2cell(0.5,0,0,bodies);
    Eigen::Vector3d vel2=P2cell(0.5,0,0,bodies);
    std::cout<<"VEL1"<<vel1.transpose()<<std::endl;
    std::cout<<"VEL2"<<vel2.transpose()<<std::endl;

   
    std::cout<<"writing file"<<k<<std::endl;
    writeTovtk(k);
    
    //blade.calculatevelocityatcenter(wakebodies);
    //blade.calculatevelocityatedge(wakebodies);
    blade.calculatevelocityatcenter(bodies);
    blade.calculatevelocityatedge(bodies);
    blade.calculateaoa();
    
    for (int v=0; v<blade.number_of_panels;v++){
     std::cout<<"the central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     std::cout<<"the aoa at panel centre "<< v<<"is"<<blade.aoa[v]*180/M_PI<<std::endl;
    }
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    for (int v=0; v<blade.number_of_panels;v++){
     std::cout<<"the lift for panel "<< v<< "is "<<blade.lift[v]<<std::endl;
     std::cout<<"the circulation for panel "<< v<< "is "<<blade.circulation[v]<<std::endl;
    }
    //shred initial particles
    //write first vtk file for shred particles     
    
    BiotSavart();
    //advancewakebodies();
    advance_euler(0.01); //add freestream influence on this
    //advance_euler_wake(0.001);
   
    //boundparticle(k); 
    
    //find velocity by particles on the panel centers and find new aoa 
    //advance particles due to fmm + free stream 
    //shred new particles by this new aoa [append]
    //write new vtk file 
    // find new velocity with all particles , advance original with fmm and free stream 
    // find new particles... cont .. 
     double tot=blade.calculatetotallift();
     double diff_lift = tot-prev_lift;
     prev_lift = tot;
     //std::cout<<"the diff is "<< diff_lift<<std::endl;
     //std::cout<<"the relative diff is "<< diff_lift/tot<<std::endl;
     std::cout<<"the total lift in iteration "<<k<<"is "<<tot<<std::endl;
    if (k==400){
    double tot_drag = blade.getdrag();
    double induced = blade.getinduced();
    std::cout<<"The total drag after 400 iterations is "<<tot_drag<<std::endl;
    std::cout<<"The induced drag after 400 iterations is "<<induced<<std::endl;
    std::cout<<"This is the end "<<std::endl;
    }
   }
// */
///////////////////////////////////////////////////////////////////////////////////////////////////


//with source panel 
//very first step if free stream is parallel to ground, no effect on the lift
//intialize wing by making array of centres and panel edges 
    //initialize the source panel AIC matrix
    //find original clmax by initial aoa set by user 
    //find original circulation 
    //shred initial particles
    //write first vtk file for shred particles 
    //find panel strengths with new freestream + particles combination 
    //find net induced velocity on panel centres and shred new particles
    //convect original particles by each other influence + panel + freestream  


    /* uncomment this for undulated surface
  Eigen::Vector3d initial;
  initial << 0.0, 0.0, 0.0;
 std::vector<Eigen::Vector3d> rhs(npanels, initial);
     for (int j = 1; j < cell_centers.size() ; j++){
      rhs[j] = P2cell( cell_centers[j][0], cell_centers[j][1], cell_centers[j][2], bodies) +blade.freestream;  
    }
    
    surface.calculateSigma(rhs);
    Eigen::VectorXd value= surface.getSigma(rhs);
      std::cout<<"sigma"<<value.transpose()<<std::endl;
    blade.calculatevelocityatcenterpanel(bodies,surface);
    blade.calculatevelocityatedgepanel(bodies,surface);
    blade.calculateaoa();
    std::cout<<"first size of bodies array is "<<bodies.size()<<std::endl;
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    double tot_first=blade.calculatetotallift();
     //std::cout<<"the total lift in iteration 0 "<<"is "<<tot_first<<std::endl;
    for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the first circulation for panel "<< v<< "is "<<blade.circulation[v]<<std::endl;
     //std::cout<<"the first circulation for panel "<< v<< "is "<<blade.prev_circulation[v]<<std::endl;
     //std::cout<<"the first lift for panel "<< v<< "is "<<blade.lift[v]<<std::endl;
    }
   //boundparticle(0);
    for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the first central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     //std::cout<<"chord vector "<<v<<"is "<< blade.chordvector[v].transpose()<<std::endl;
     //std::cout<<"the first aoa at panel centre "<< v<<"is"<<blade.aoa[v]*180/M_PI<<std::endl;
    }

 double prev_lift = tot_first;

    for (int k = 1; k < 401; k++)
   { 
    
    blade.total_lift=0;
    trailshedparticle(0.01);
    spanshedparticle(0.01);
    std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
    std::cout<<"writing file"<<k<<std::endl;
    writeTovtk(k);
    Eigen::Vector3d initial;
    initial << 0.0, 0.0, 0.0;
     for (int v=0; v<blade.number_of_panels;v++){
     std::cout<<"the lift for panel "<< v<< "is "<<blade.lift[v]<<std::endl;
     std::cout<<"the circulation for panel "<< v<< "is "<<blade.circulation[v]<<std::endl;
    }
    std::vector<Eigen::Vector3d> rhs(npanels, initial);
     for (int v=0; v<blade.number_of_panels;v++){
     std::cout<<"the central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     std::cout<<"the aoa at panel centre "<< v<<"is"<<blade.aoa[v]*180/M_PI<<std::endl;
    }

   for (int j = 1; j < cell_centers.size() ; j++){
      rhs[j] = P2cell( cell_centers[j][0], cell_centers[j][1], cell_centers[j][2], bodies) +blade.freestream;  
    }
    
    surface.calculateSigma(rhs);
    Eigen::VectorXd value= surface.getSigma(rhs);
    
    blade.calculatevelocityatcenterpanel(bodies,surface);
    blade.calculatevelocityatedgepanel(bodies,surface);
    blade.calculateaoa();
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    BiotSavart();
    advance_euler_panel_nostretch(0.01);
    advance_euler(0.01);//freestream influence added in this. 
    double tot=blade.calculatetotallift();
    std::cout<<"the total lift in iteration "<<k<<"is "<<tot<<std::endl;
    if (k==400){
    double tot_drag = blade.getdrag();
    double induced = blade.getinduced();
    std::cout<<"The total drag after 400 iterations is "<<tot_drag<<std::endl;
    std::cout<<"The induced drag after 400 iterations is "<<induced<<std::endl;
    std::cout<<"This is the end "<<std::endl;
    }
    }
 */
 
 Eigen::MatrixXd clarray(50, 6);

    // Initialize the matrix to zero
    //clarray.setZero();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //for pitching motion 

   /*
     blade.wingSpan = 1;
   blade.height = 0;
   blade.set_angle = 0; //degrees
   blade.setpanelcenter(); //initilalzing wing geometry 
   blade.setpaneledge();
    blade.freestream<<0,1,0;
   blade.calculatechordlength();
   blade.setleadingtrailingedge();
   blade.setchordvector();
    blade.calculatevelocityatcenter(bodies);
    blade.calculatevelocityatedge(bodies);
    blade.calculateaoa();
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
     tot_first=blade.calculatetotallift();
    blade.total_lift=0;
*/
    //without ground effect 
/*
    const int MAX_SIZE = 20000; // Set the maximum size limit for the vector
    const int REMOVE_COUNT = 1010;
for (int k = 1; k < 4; k++){
 std::cout<<"cycle "<<k<<std::endl;
   for (int s = 1; s<blade.full_period +1;s++){
     
        double tot=blade.calculatetotallift();
     std::cout<<"the total lift in cycle "<<k<<"and iteration "<<s<<"is"<<tot<<std::endl;
     double tot_cl = tot/(0.5*blade.airDensity*blade.freestream.norm()*blade.freestream.norm()*blade.area);
     std::cout<<"the total cl in cycle "<<k<<"and iteration "<<s<<"is"<<tot_cl<<std::endl;
     clarray(s-1,k-1) = tot_cl;
    std::cout<<"the current set aoa is "<<blade.set_angle*180/M_PI<<std::endl;
    std::cout<<"the current actual aoa is "<<blade.aoa[25]*180/M_PI<<std::endl;
    std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
    //std::cout<<"the size of the wake body array is "<<wakebodies.size()<<std::endl;
    
    
    
    //blade.calculatevelocityatcenter(wakebodies);
    //blade.calculatevelocityatedge(wakebodies);
    //std::cout<<"outside loop"<<s<<std::endl;
    
    //for (int i =0;i<100;i++){ //convergence for single aoa
    //std::cout<<"inside loop"<<i<<std::endl;
    //std::cout<<"outside loop"<<s<<std::endl;
       trailshedparticle(0.01);
       spanshedparticle(0.01);
       std::cout<<"writing file"<<s<<std::endl;
       int l = 100*k + s;
       writeTovtk(l);
       blade.set_angle= blade.sinosoidangle(s)*M_PI/180;
       blade.setleadingtrailingedge();
       blade.setchordvector(); //degrees
       blade.calculatevelocityatcenter(bodies);
       blade.calculatevelocityatedge(bodies);
       blade.calculateaoa();
       blade.calculateliftinpanel();
       blade.calculatecirculationinpanel();
       
       BiotSavart();    
       advance_euler(0.01);
       
       //std::cout<<"the current aoa is "<<blade.aoa[25]*180/M_PI<<std::endl;
        
   // }
    
     //if (bodies.size() > MAX_SIZE) {
            // Erase the first REMOVE_COUNT elements
            //bodies.erase(bodies.begin(), bodies.begin() + REMOVE_COUNT);
            //std::cout << "Removed first " << REMOVE_COUNT << " elements. New size: " << vec.size() << std::endl;
        //}
        Eigen::Vector3d vel1=P2cell(0.5,0,0,bodies);
    Eigen::Vector3d vel2=P2cell(0.5,0,0,bodies);
    //std::cout<<"VEL1"<<vel1.transpose()<<std::endl;
    //std::cout<<"VEL2"<<vel2.transpose()<<std::endl;
    //blade.set_angle = blade.pitchaoa[s]*M_PI/180;
    
     //for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     //std::cout<<"the aoa at panel centre "<< v<<"is "<<blade.aoa[v]*180/M_PI<<std::endl;
    //}

    
    //blade.calculateliftinpanel();
    //blade.calculatecirculationinpanel();
   
    //shred initial particles
    //write first vtk file for shred particles     
    
    //BiotSavart();
    //advancewakebodies();
    //advance_euler(0.01); //add freestream influence on this
    //advance_euler_wake(0.001);
   
   
     //double diff_lift = tot-prev_lift;
     //prev_lift = tot;
     //std::cout<<"the diff is "<< diff_lift<<std::endl;
     //std::cout<<"the relative diff is "<< diff_lift/tot<<std::endl;
     

   }
   }
   */
  ////////////////////////////////////////////   
   //pitching with ground effect
  /*
   for (int k = 1; k < 4; k++){
 std::cout<<"cylce "<<k<<std::endl;
   for (int s = 1; s<blade.full_period +1;s++){
   
   double tot=blade.calculatetotallift();
     std::cout<<"the total lift in cycle "<<k<<"and iteration "<<s<<"is"<<tot<<std::endl;
     double tot_cl = tot/(0.5*blade.airDensity*blade.freestream.norm()*blade.freestream.norm()*blade.area);
     std::cout<<"the total cl in cycle "<<k<<"and iteration "<<s<<"is"<<tot_cl<<std::endl;
     clarray(s-1,k-1) = tot_cl;

    trailshedparticle(0.01);
    spanshedparticle(0.01);
     std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
    //std::cout<<"the size of the wake body array is "<<wakebodies.size()<<std::endl;
   std::cout<<"the current set aoa is "<<blade.set_angle<<std::endl;
    //std::cout<<"the current aoa is "<<aoa[25]<<std::endl;
   
    std::cout<<"writing file"<<s<<std::endl;
    int i = 100*k + s;
    writeTovtk(i);
    Eigen::Vector3d initial;
    initial << 0.0, 0.0, 0.0;
    std::vector<Eigen::Vector3d> rhs(npanels, initial);
    for (int j = 1; j < cell_centers.size() ; j++){
      rhs[j] = P2cell( cell_centers[j][0], cell_centers[j][1], cell_centers[j][2], bodies) +blade.freestream;  
    }
    
    surface.calculateSigma(rhs);
    //blade.calculatevelocityatcenter(wakebodies);
    //blade.calculatevelocityatedge(wakebodies);
    //blade.calculatevelocityatcenter(bodies);
    //blade.calculatevelocityatedge(bodies);
     blade.set_angle= blade.sinosoidangle(s)*M_PI/180;
     blade.setleadingtrailingedge();
     blade.setchordvector(); //degrees

    blade.calculatevelocityatcenterpanel(bodies,surface);
    blade.calculatevelocityatedgepanel(bodies,surface);
    //blade.set_angle = blade.pitchaoa[s]; //degrees
   // blade.calculateaoa(k);
     blade.calculateaoa();
    
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
   
    //shred initial particles
    //write first vtk file for shred particles     
    
    BiotSavart();
    //advancewakebodies();
    advance_euler_panel_nostretch(0.01);
    advance_euler(0.01);//freestream influence added in this. 
    //advance_euler_wake(0.001);
   
     //double tot=blade.calculatetotallift();
     //double diff_lift = tot-prev_lift;
     //prev_lift = tot;
     //std::cout<<"the diff is "<< diff_lift<<std::endl;
     //std::cout<<"the relative diff is "<< diff_lift/tot<<std::endl;
     //std::cout<<"the total lift in iteration "<<k<<"is "<<tot<<std::endl;

   }
   }
   
*/
  /////////////////////////////////HEAVING MOTION /////////////////////////////////////////////////////
 /*
///without source panel
clarray(blade.heave_total_steps, 6);
std::vector<double> positions(blade.heave_total_steps+1);
std::vector<double> velocities(blade.heave_total_steps+1);
   for (int i = 0; i < blade.heave_total_steps+1; ++i) {
        double time = i * blade.heave_time_step; // Calculate the current time
         
        // Calculate the position and velocity at this time step
        positions[i] = blade.calculateWingPosition_heave(time);
        velocities[i] = -1*blade.calculateWingVelocity_heave(time);

    }

    //initlializing at time = 0
     blade.height = positions[0];
    blade.setleadingtrailingedge();
       blade.setchordvector(); //degrees
blade.calculatevelocityatcenterheave(bodies,velocities[0]);
    blade.calculatevelocityatedgeheave(bodies,velocities[0]);
    blade.calculateaoa();
    std::cout<<"first size of bodies array is "<<bodies.size()<<std::endl;
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    double tot_first=blade.calculatetotallift();
for (int k = 1; k < 4; k++){
 std::cout<<"cycle "<<k<<std::endl;
   for (int s = 1; s<blade.heave_total_steps+1;s++){
      std::cout<<"writing file"<<s<<std::endl;
        double tot=blade.calculatetotallift();
         std::cout<<"the current aoa at centre panel is "<<blade.aoa[25]*180/M_PI<<std::endl;
         std::cout<<"the current lift at centre panel is "<<blade.lift[25]<<std::endl;
        // std::cout<<"the current velocity due to wake and freestream is "<<blade.calculatevelocityatcenter(bodies)<<std::endl;
         // std::cout<<"the current velocity due to wake and freestream and body is "<<blade.calculatevelocityatcenterheave(bodies,velocities[s-1])<<std::endl;
     std::cout<<"the total lift in cycle "<<k<<"and iteration "<<s<<"is"<<tot<<std::endl;
     double tot_cl = tot/(0.5*blade.airDensity*blade.freestream.norm()*blade.freestream.norm()*blade.area);
     std::cout<<"the total cl in cycle "<<k<<"and iteration "<<s<<"is"<<tot_cl<<std::endl;
     clarray(s-1,k-1) = tot_cl;
      Eigen::Vector3d vel1=P2cell(blade.panel_center[25],0,blade.height,bodies);
   
    std::cout<<"VEL1"<<vel1.transpose()<<std::endl;

     std::cout<<"the added velocity is "<<velocities[s-1]<<std::endl;
    std::cout<<"the current set height is "<<blade.height<<std::endl;
    std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
       trailshedparticle(blade.heave_time_step);
       spanshedparticle(blade.heave_time_step);
      
       int l = 100*k + s;
       writeTovtk(l);
      
       blade.prevheight = blade.height;
       blade.height = positions[s];
       double extraheavevel= velocities[s];
   
      blade.setleadingtrailingedge();
       blade.setchordvector(); //degrees
       blade.calculatevelocityatcenterheave(bodies,extraheavevel);
       blade.calculatevelocityatedgeheave(bodies,extraheavevel);
     
       blade.calculateaoa();
       blade.calculateliftinpanel();
       blade.calculatecirculationinpanel();
       
       BiotSavart();    
       advance_euler(blade.heave_time_step);
       
       std::cout<<"End of file"<<s<<std::endl;
        
   // }
    
     //if (bodies.size() > MAX_SIZE) {
            // Erase the first REMOVE_COUNT elements
            //bodies.erase(bodies.begin(), bodies.begin() + REMOVE_COUNT);
            //std::cout << "Removed first " << REMOVE_COUNT << " elements. New size: " << vec.size() << std::endl;
        //}

    //std::cout<<"VEL1"<<vel1.transpose()<<std::endl;
    //std::cout<<"VEL2"<<vel2.transpose()<<std::endl;
    //blade.set_angle = blade.pitchaoa[s]*M_PI/180;
    
     //for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     //std::cout<<"the aoa at panel centre "<< v<<"is "<<blade.aoa[v]*180/M_PI<<std::endl;
    //}

    
    //blade.calculateliftinpanel();
    //blade.calculatecirculationinpanel();
   
    //shred initial particles
    //write first vtk file for shred particles     
    
    //BiotSavart();
    //advancewakebodies();
    //advance_euler(0.01); //add freestream influence on this
    //advance_euler_wake(0.001);
   
   
     //double diff_lift = tot-prev_lift;
     //prev_lift = tot;
     //std::cout<<"the diff is "<< diff_lift<<std::endl;
     //std::cout<<"the relative diff is "<< diff_lift/tot<<std::endl;
     

   }
   }


//with source panel
 */
 /*
clarray(blade.heave_total_steps, 6);
std::vector<double> positions(blade.heave_total_steps+1);
std::vector<double> velocities(blade.heave_total_steps+1);
   for (int i = 0; i < blade.heave_total_steps+1; ++i) {
        double time = i * blade.heave_time_step; // Calculate the current time
         
        // Calculate the position and velocity at this time step
        positions[i] = blade.calculateWingPosition_heave(time);
        velocities[i] = -1*blade.calculateWingVelocity_heave(time);

    }

    //initlializing at time = 0
     blade.height = positions[0];
    blade.setleadingtrailingedge();
       blade.setchordvector(); //degrees
             Eigen::Vector3d initial_panel;
    initial_panel << 0.0, 0.0, 0.0;
    std::vector<Eigen::Vector3d> rhs_panel(npanels, initial_panel);
    for (int j = 1; j < cell_centers.size() ; j++){
      rhs_panel[j] = P2cell( cell_centers[j][0], cell_centers[j][1], cell_centers[j][2], bodies) +blade.freestream;  
    }
    
    surface.calculateSigma(rhs_panel);
blade.calculatevelocityatcenterpanelheave(bodies,surface,velocities[0]);
    blade.calculatevelocityatedgepanelheave(bodies,surface,velocities[0]);
    blade.calculateaoa();
    std::cout<<"first size of bodies array is "<<bodies.size()<<std::endl;
    blade.calculateliftinpanel();
    blade.calculatecirculationinpanel();
    double tot_first=blade.calculatetotallift();
for (int k = 1; k < 4; k++){
 std::cout<<"cycle "<<k<<std::endl;
   for (int s = 1; s<blade.heave_total_steps+1;s++){
      std::cout<<"writing file"<<s<<std::endl;
        double tot=blade.calculatetotallift();
         std::cout<<"the current aoa at centre panel is "<<blade.aoa[25]*180/M_PI<<std::endl;
         std::cout<<"the current lift at centre panel is "<<blade.lift[25]<<std::endl;
        // std::cout<<"the current velocity due to wake and freestream is "<<blade.calculatevelocityatcenter(bodies)<<std::endl;
         // std::cout<<"the current velocity due to wake and freestream and body is "<<blade.calculatevelocityatcenterheave(bodies,velocities[s-1])<<std::endl;
     std::cout<<"the total lift in cycle "<<k<<"and iteration "<<s<<"is"<<tot<<std::endl;
     double tot_cl = tot/(0.5*blade.airDensity*blade.freestream.norm()*blade.freestream.norm()*blade.area);
     std::cout<<"the total cl in cycle "<<k<<"and iteration "<<s<<"is"<<tot_cl<<std::endl;
     clarray(s-1,k-1) = tot_cl;
      Eigen::Vector3d vel1=P2cell(blade.panel_center[25],0,blade.height,bodies);
   
    std::cout<<"VEL1"<<vel1.transpose()<<std::endl;

     std::cout<<"the added velocity is "<<velocities[s-1]<<std::endl;
    std::cout<<"the current set height is "<<blade.height<<std::endl;
    std::cout<<"the size of the body array is "<<bodies.size()<<std::endl;
       trailshedparticle(blade.heave_time_step);
       spanshedparticle(blade.heave_time_step);
      
       int l = 100*k + s;
       writeTovtk(l);
      
       blade.prevheight = blade.height;
       blade.height = positions[s];
       double extraheavevel= velocities[s];
   
      blade.setleadingtrailingedge();
       blade.setchordvector(); //degrees
          Eigen::Vector3d initial;
    initial << 0.0, 0.0, 0.0;
    std::vector<Eigen::Vector3d> rhs(npanels, initial);
    for (int j = 1; j < cell_centers.size() ; j++){
      rhs[j] = P2cell( cell_centers[j][0], cell_centers[j][1], cell_centers[j][2], bodies) +blade.freestream;  
    }
    
    surface.calculateSigma(rhs);
       blade.calculatevelocityatcenterpanelheave(bodies,surface,extraheavevel);
       blade.calculatevelocityatedgepanelheave(bodies,surface,extraheavevel);
     
       blade.calculateaoa();
       blade.calculateliftinpanel();
       blade.calculatecirculationinpanel();
       
       BiotSavart();    
       advance_euler(blade.heave_time_step);
       
       std::cout<<"End of file"<<s<<std::endl;
        
   // }
    
     //if (bodies.size() > MAX_SIZE) {
            // Erase the first REMOVE_COUNT elements
            //bodies.erase(bodies.begin(), bodies.begin() + REMOVE_COUNT);
            //std::cout << "Removed first " << REMOVE_COUNT << " elements. New size: " << vec.size() << std::endl;
        //}

    //std::cout<<"VEL1"<<vel1.transpose()<<std::endl;
    //std::cout<<"VEL2"<<vel2.transpose()<<std::endl;
    //blade.set_angle = blade.pitchaoa[s]*M_PI/180;
    
     //for (int v=0; v<blade.number_of_panels;v++){
     //std::cout<<"the central velocity for panel "<< v<< "is "<<(blade.velocityatcenter[v]).transpose()<<std::endl;
     //std::cout<<"the aoa at panel centre "<< v<<"is "<<blade.aoa[v]*180/M_PI<<std::endl;
    //}

    
    //blade.calculateliftinpanel();
    //blade.calculatecirculationinpanel();
   
    //shred initial particles
    //write first vtk file for shred particles     
    
    //BiotSavart();
    //advancewakebodies();
    //advance_euler(0.01); //add freestream influence on this
    //advance_euler_wake(0.001);
   
   
     //double diff_lift = tot-prev_lift;
     //prev_lift = tot;
     //std::cout<<"the diff is "<< diff_lift<<std::endl;
     //std::cout<<"the relative diff is "<< diff_lift/tot<<std::endl;
     

   }
   }
 */

//////////////////////////////////////////////////////////////////////////////////////////

/*
   std::ofstream file0("outputcl0.txt");
    if (!file0) {
        std::cerr << "Error opening file outputcl0.txt!" << std::endl;
        return 1;
    }
    for (int i = 1; i < blade.heave_total_steps+1; i++) {
        file0 << clarray(i-1, 0) << "\n";
    }
    file0.close();

    // Write to outputcl1.txt
    std::ofstream file1("outputcl1.txt");
    if (!file1) {
        std::cerr << "Error opening file outputcl1.txt!" << std::endl;
        return 1;
    }
    for (int i = 1; i < blade.heave_total_steps+1; i++) {
        file1 << clarray(i-1, 1) << "\n";
    }
    file1.close();

    // Write to outputcl2.txt
    std::ofstream file2("outputcl2.txt");
    if (!file2) {
        std::cerr << "Error opening file outputcl2.txt!" << std::endl;
        return 1;
    }
    for (int i = 1; i < blade.heave_total_steps+1; i++) {
        file2 << clarray(i-1, 2) << "\n";
    }
    file2.close();
//*/
    // Write to outputcl3.txt


    std::ofstream file0("outputcl0.txt");
    if (!file0) {
        std::cerr << "Error opening file outputcl0.txt!" << std::endl;
        return 1;
    }
    for (int i = 0; i < 50; i++) {
        file0 << clarray(i, 0) << "\n";
    }
    file0.close();

    // Write to outputcl4.txt
    std::ofstream file1("outputcl1.txt");
    if (!file1) {
        std::cerr << "Error opening file outputcl4.txt!" << std::endl;
        return 1;
    }
    for (int i = 0; i < 50; i++) {
        file1 << clarray(i, 1) << "\n";
    }
    file1.close();

    // Write to outputcl5.txt
    std::ofstream file2("outputcl2.txt");
    if (!file2) {
        std::cerr << "Error opening file outputcl5.txt!" << std::endl;
        return 1;
    }
    for (int i = 0; i < 50; i++) {
        file2 << clarray(i, 2) << "\n";
    }
    file2.close();

    return 0;

    
}

