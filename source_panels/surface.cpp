//  ███████████    █████████  ███████████ ████                               ████████  ██████████
// ░░███░░░░░███  ███░░░░░███░░███░░░░░░█░░███                              ███░░░░███░░███░░░░███
//  ░███    ░███ ███     ░░░  ░███   █ ░  ░███   ██████  █████ ███ █████   ░░░    ░███ ░███   ░░███
//  ░██████████ ░███          ░███████    ░███  ███░░███░░███ ░███░░███       ██████░  ░███    ░███
//  ░███░░░░░░  ░███    █████ ░███░░░█    ░███ ░███ ░███ ░███ ░███ ░███      ░░░░░░███ ░███    ░███
//  ░███        ░░███  ░░███  ░███  ░     ░███ ░███ ░███ ░░███████████      ███   ░███ ░███    ███
//  █████        ░░█████████  █████       █████░░██████   ░░████░████      ░░████████  ██████████
// ░░░░░          ░░░░░░░░░  ░░░░░       ░░░░░  ░░░░░░     ░░░░ ░░░░        ░░░░░░░░  ░░░░░░░░░░

// =========================================================================
//  Copyright (C) 2024-2024 Technical University of Munich, ENAC - Ecole Nationale de l'Aviation Civile
//  This file is part of PGFlow3D.

//  Permission is hereby granted, free of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files (the "Software"), to deal in the Software without
//  restriction, including without limitation the rights to use,
//  copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following
//  conditions:

//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.

//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
//  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
//  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//  OTHER DEALINGS IN THE SOFTWARE.

//  Authors:
//           M Kuersat Yurt, MSc. (TUM)
// =========================================================================

#include "surface.hpp"
#include <Eigen/Geometry>
#include <chrono>
#include <filesystem>
#include <iostream>
#include "rapidobj/rapidobj.hpp"

Eigen::Vector3d Surface::transformVectorToInertialFrame(Eigen::Vector3d vectorInLocalFrame, size_t panelId)
{
  Eigen::Vector3d vectorInInertialFrame;

  vectorInInertialFrame(0) = (vectorInLocalFrame(0) * _cellOblique[panelId](0) + vectorInLocalFrame(1) * _cellTangent[panelId](0) + vectorInLocalFrame(2) * _cellNormals[panelId](0));
  vectorInInertialFrame(1) = (vectorInLocalFrame(0) * _cellOblique[panelId](1) + vectorInLocalFrame(1) * _cellTangent[panelId](1) + vectorInLocalFrame(2) * _cellNormals[panelId](1));
  vectorInInertialFrame(2) = (vectorInLocalFrame(0) * _cellOblique[panelId](2) + vectorInLocalFrame(1) * _cellTangent[panelId](2) + vectorInLocalFrame(2) * _cellNormals[panelId](2));
  return vectorInInertialFrame;
}

Eigen::Vector3d Surface::transformPointToPanelLocalFrame(const Eigen::Vector3d &pointInInertialFrame, size_t panelId)
{
  Eigen::Vector3d pointInLocalFrame;
  auto            Diff = pointInInertialFrame - _cellCenters[panelId];
  pointInLocalFrame(0) = Diff.dot(_cellOblique[panelId]);
  pointInLocalFrame(1) = Diff.dot(_cellTangent[panelId]);
  pointInLocalFrame(2) = Diff.dot(_cellNormals[panelId]);
  return pointInLocalFrame;
}

Eigen::Vector3d Surface::calculateVelocityDueToSinglePanel(Eigen::Vector3d point, size_t panelId)
{
  Eigen::Vector3d velocityAtPoint{0, 0, 0};
  // Iterate over each edge of the triangular panel
  bool all_inside{true};

  // Get the coordinates of the three nodes and point in panel frame
  Eigen::Vector3d              point_coo = transformPointToPanelLocalFrame(point, panelId);
  std::vector<Eigen::Vector3d> panel_coo(3);
  for (size_t i = 0; i < 3; i++) {
    panel_coo[i] = transformPointToPanelLocalFrame(_vertices[_connectivity[panelId](i)], panelId);
  }

  for (size_t i = 0; i < 3; i++) {
    // Define the indices of the two nodes that make the edge
    size_t j = (i + 1) % 3;
    // Get the coordinates of the two nodes in panel frame
    const double xi = panel_coo[i](0);
    const double yi = panel_coo[i](1);
    const double xj = panel_coo[j](0);
    const double yj = panel_coo[j](1);

    //  Calculate the terms needed for the velocity components
    const double d_ij   = std::sqrt((yj - yi) * (yj - yi) + (xj - xi) * (xj - xi));
    const double S_ij   = (yj - yi) / d_ij;
    const double C_ij   = (xj - xi) / d_ij;
    const double s_ij_i = (xi - point_coo(0)) * C_ij + (yi - point_coo(1)) * S_ij;
    const double s_ij_j = (xj - point_coo(0)) * C_ij + (yj - point_coo(1)) * S_ij;
    const double R_ij   = (point_coo(0) - xi) * S_ij - (point_coo(1) - yi) * C_ij;
    const double r_i    = std::sqrt((point_coo(0) - xi) * (point_coo(0) - xi) + (point_coo(1) - yi) * (point_coo(1) - yi) + (point_coo(2) * point_coo(2)));
    const double r_j    = std::sqrt((point_coo(0) - xj) * (point_coo(0) - xj) + (point_coo(1) - yj) * (point_coo(1) - yj) + (point_coo(2) * point_coo(2)));
    const double Q_ij   = std::log((r_i + r_j + d_ij) / (r_i + r_j - d_ij + 1e-14));
    all_inside          = all_inside && (R_ij > 0.0);
    //const double e_i= point_coo(2)*point_coo(2) + (point_coo(0)-xi)*(point_coo(0)-xi);
    //const double h_i= (point_coo(1)-yi)*(point_coo(0)-xi);
   // const double e_j= point_coo(2)*point_coo(2) + (point_coo(0)-xj)*(point_coo(0)-xj);
  //  const double h_j= (point_coo(1)-yj)*(point_coo(0)-xj);
  //  const double m_ij= (yj-yi)/(xj-xi);


    const double J_ij = std::atan2(R_ij * std::abs(point_coo(2)) * (r_i * s_ij_j - r_j * s_ij_i), r_i * r_j * R_ij * R_ij + point_coo(2) * point_coo(2) * s_ij_j * s_ij_i);
    velocityAtPoint(0) -= S_ij * Q_ij;
    velocityAtPoint(1) += C_ij * Q_ij;
    velocityAtPoint(2) -= std::copysign(1.0, point_coo(2)) * J_ij;
    //velocityAtPoint(2) = std::atan2((m_ij*e_i-h_i),(point_coo(2)*r_i)) - std::atan2((m_ij*e_j-h_j),(point_coo(2)*r_j));
  }
  if (all_inside && std::abs(point_coo(2)) < 1e-12) { // Very close to zero
    velocityAtPoint(2) += 2 * M_PI;
  } else if (all_inside) { // Inside panel but not close to zero
    velocityAtPoint(2) += std::copysign(1.0, point_coo(2)) * 2 * M_PI;
  }

  velocityAtPoint /= 4 * M_PI;

  // Transform the velocity to the inertial frame
  velocityAtPoint = transformVectorToInertialFrame(velocityAtPoint, panelId);
  return velocityAtPoint;
}

void Surface::initFromOBJFile(const std::string &filename)
{
  // Check if the file exist or not
  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("File does not exist " + filename);
  }

  // Get the extension of the file
  std::filesystem::path pathOfFile(filename);
  std::string           extension = pathOfFile.extension().string();
  // Read the file based on the extension
  if (extension != ".obj") {
    throw std::runtime_error("File format not supported");
  }
  auto objFile = rapidobj::ParseFile(pathOfFile);

  if (objFile.error) {
    std::cout << objFile.error.code.message() << "\n";
    if (!objFile.error.line.empty()) {
      std::cout << "On line " << objFile.error.line_num << ": \"" << objFile.error.line << "\"\n";
    }
  }

  // Get the number of cells and points
  size_t nPoints = objFile.attributes.positions.size() / 3;

  auto nCells = size_t(0);
  for (const auto &shape : objFile.shapes) {
    nCells += shape.mesh.num_face_vertices.size();
  }

  // Resize the vectors
  _vertices.resize(nPoints);
  _connectivity.resize(nCells);
  _cellCenters.resize(nCells);
  _cellNormals.resize(nCells);
  _cellTangent.resize(nCells);
  _cellOblique.resize(nCells);
  _AICMatrix.resize(nCells, nCells);
  _AICMatrixInv.resize(nCells, nCells);
  _sigma.resize(nCells);

  // Get the points
  for (size_t i = 0; i < nPoints; i++) {
    _vertices[i] = Eigen::Vector3d(objFile.attributes.positions[3 * i], objFile.attributes.positions[3 * i + 1], objFile.attributes.positions[3 * i + 2]);
  }

  // Get the cells
  for (const auto &shape : objFile.shapes) {
    for (unsigned char num_face_vertice : shape.mesh.num_face_vertices) {
      if (num_face_vertice != 3) {
        throw std::runtime_error("Only triangular cells are supported");
      }
    }
    auto &indices = shape.mesh.indices;
    for (size_t i = 0; i < shape.mesh.num_face_vertices.size(); i++) {
      _connectivity[i] = Eigen::Vector3i(indices[3 * i].position_index, indices[3 * i + 1].position_index, indices[3 * i + 2].position_index);
    }
  }

  // Calculate the cell centers
  for (size_t i = 0; i < nCells; i++) {
    _cellCenters[i] = (_vertices[_connectivity[i](0)] + _vertices[_connectivity[i](1)] + _vertices[_connectivity[i](2)]) / 3;
  }

  // Calculate the cell tangents,normals and obliques
  for (size_t i = 0; i < nCells; i++) {
    Eigen::Vector3d v1 = _vertices[_connectivity[i](1)] - _vertices[_connectivity[i](2)];
    Eigen::Vector3d v2 = _vertices[_connectivity[i](1)] - _vertices[_connectivity[i](0)];
    _cellTangent[i]    = v1.normalized();
    _cellNormals[i]    = v1.cross(v2).normalized();
    _cellOblique[i]    = _cellNormals[i].cross(_cellTangent[i]);
  }

  calculateAICMatrix();
  calculateAICMatrixInv();
}

void Surface::calculateAICMatrix()
{
  std::cout << "Calculating AIC matrix" << std::endl;

  const auto start{std::chrono::steady_clock::now()};

  size_t nCells = _connectivity.size();
#pragma omp parallel for collapse(2)
  for (size_t i = 0; i < nCells; i++) {
    for (size_t j = 0; j < nCells; j++) {
      const Eigen::Vector3d AIC = calculateVelocityDueToSinglePanel(_cellCenters[i], j);
      _AICMatrix(i, j)          = AIC.dot(_cellNormals[i]);
    }
  }
  std::cout <<nCells << std::endl;
  //for (size_t i = 0; i < nCells; i++) {
   // for (size_t j = 0; j < nCells; j++) {
    //  if (j==i+1){
    //  std::cout <<"The value at" <<i<<"and"<<j<<"is:"<< _AICMatrix(i,j) << std::endl;
     // }
  // }
  //}
  const auto                          end{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{end - start};
  std::cout << elapsed_seconds.count() << "s\n";
  
         
      
    
}

void Surface::calculateAICMatrixInv()
{
  std::cout << "Calculating AIC matrix inverse" << std::endl;
  const auto start{std::chrono::steady_clock::now()};
  _AICMatrixInv = _AICMatrix.inverse();
  const auto                          end{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{end - start};
  std::cout << elapsed_seconds.count() << "s\n";
}

void Surface::calculateSigma(const std::vector<Eigen::Vector3d> &rhs )
{
  // Multiply RHS Velocity by normal vectors to construct actual RHS
  Eigen::VectorXd b(_cellCenters.size());
  for (size_t i = 0; i < _cellCenters.size(); i++) {
    b(i) = -(_cellNormals[i].dot(rhs[i]));
  }
  _sigma = _AICMatrixInv * b;
  //return _sigma;
}
Eigen::VectorXd Surface::getSigma(const std::vector<Eigen::Vector3d> &rhs )
{
  // Multiply RHS Velocity by normal vectors to construct actual RHS
  Eigen::VectorXd b(_cellCenters.size());
  for (size_t i = 0; i < _cellCenters.size(); i++) {
    b(i) = -(_cellNormals[i].dot(rhs[i]));
  }
  sigmaget = _AICMatrixInv * b;
  return sigmaget;
}
Eigen::Vector3d Surface::calculateVelocityAtPoint(const Eigen::Vector3d &point)
{
  Eigen::Vector3d velocity{0, 0, 0};
  for (size_t i = 0; i < _cellCenters.size(); i++) {
    velocity += calculateVelocityDueToSinglePanel(point, i) * _sigma(i);
  }
  return velocity;
}

Eigen::Vector3d Surface::calculateVelocityAtPointcoord(const double x_coord,const  double y_coord, const double z_coord)
{
  Eigen::Vector3d velocity{0, 0, 0};
  Eigen::Vector3d point;
      point[0]=x_coord;
      point[1]=y_coord;
      point[2]=z_coord;
  for (size_t i = 0; i < _cellCenters.size(); i++) {
    velocity += calculateVelocityDueToSinglePanel(point, i) * _sigma(i);
  }
  return velocity;
}

std::vector<Eigen::Vector3d> &Surface::getCellCenters()
{
  return _cellCenters;
}
/*
Eigen::Vector3d paneleffectonparticle(double x_coord, double y_coord, double z_coord){
      
      Eigen::Vector3d point;
      point[0]=x_coord;
      point[1]=y_coord;
      point[2]=z_coord;

      Eigen::Vector3d Vel_due_panel;
      for( int i=0;i<_cellCenters.size();i++){

        Vel_due_panel += calculateVelocityDueToSinglePanel((point, i))*_sigma(i); // iterate over j panel

      }
      
      return Vel_due_panel;
      
}
*/
Eigen::Matrix3d Surface::Calculatederivatesduetosinglepanel(Eigen::Vector3d point, size_t panelId){
   
   Eigen::Matrix3d matrix;
   Eigen::Matrix3d T;
   Eigen::Vector3d              point_coo = transformPointToPanelLocalFrame(point, panelId);
  std::vector<Eigen::Vector3d> panel_coo(3);
  for (size_t i = 0; i < 3; i++) {
    panel_coo[i] = transformPointToPanelLocalFrame(_vertices[_connectivity[panelId](i)], panelId);
  }
   double du_dx=0;
   double du_dy=0;
   double du_dz=0;
   double dv_dx=0;
   double dv_dy=0;
   double dv_dz=0;
   double dw_dx=0;
   double dw_dy=0;
   double dw_dz=0;

  for (size_t i = 0; i < 3; i++) {
    // Define the indices of the two nodes that make the edge
    size_t j = (i + 1) % 3;
    // Get the coordinates of the two nodes in panel frame
    
    const double xi = panel_coo[i](0);
    const double yi = panel_coo[i](1);
    const double xj = panel_coo[j](0);
    const double yj = panel_coo[j](1);
    
    const double d_ij   = std::sqrt((yj - yi) * (yj - yi) + (xj - xi) * (xj - xi));
    const double S_ij   = (yj - yi) / d_ij;
    const double C_ij   = (xj - xi) / d_ij;
    const double s_ij_i = (xi - point_coo(0)) * C_ij + (yi - point_coo(1)) * S_ij;
    const double s_ij_j = (xj - point_coo(0)) * C_ij + (yj - point_coo(1)) * S_ij;
    const double R_ij   = (point_coo(0) - xi) * S_ij - (point_coo(1) - yi) * C_ij;
    const double r_i    = std::sqrt((point_coo(0) - xi) * (point_coo(0) - xi) + (point_coo(1) - yi) * (point_coo(1) - yi) + (point_coo(2) * point_coo(2)));
    const double r_j    = std::sqrt((point_coo(0) - xj) * (point_coo(0) - xj) + (point_coo(1) - yj) * (point_coo(1) - yj) + (point_coo(2) * point_coo(2)));
    const double Q_ij   = std::log((r_i + r_j + d_ij) / (r_i + r_j - d_ij + 1e-14));
    //Calculating Derivatives
    double epsilon = 1e-10;
    //std::cout << "r_i: " << r_i << ", r_j: " << r_j << ", Q_ij: " << Q_ij << std::endl;
    /*
    double term1=(yj-yi)/d_ij;
    double term2 = (point_coo(0) - xi)/(r_i);
    double term3 = (point_coo(0) - xj)/(r_j);
    double term4 = r_i + r_j - d_ij;
    double term5 = r_i + r_j + d_ij;
    double term6 = (point_coo(1) - yi)/(r_i);
    double term7 = (point_coo(1) - yj)/(r_j);
    double term8 = (xi-xj)/d_ij;
    double term9 = point_coo(2)/(r_i) ;
    double term10 = point_coo(2)/(r_j) ;
    double term11 = (yj-yi)/(xj-xi);
    double term12 = (point_coo(1) - yi)*(point_coo(0) - xi);
    double term13 = (point_coo(1) - yj)*(point_coo(0) - xj);
    double term14 = 2*term11*(point_coo(0)-xi) + yi - point_coo(1) ;
    double term15 = term11*((point_coo(0)-xi)*((point_coo(0)-xi)) + point_coo(2)*point_coo(2));
    double term18 = 2*term11*(point_coo(0)-xj) + yj - point_coo(1) ;
    double term19 = term11*((point_coo(0)-xj)*((point_coo(0)-xj)) + point_coo(2)*point_coo(2));
    double term16 = (term14/(point_coo(2)*(r_i)) - ((point_coo(0)- xi)*(term15 - term12))/(point_coo(2)*r_i*r_i*r_i))/(((term15-term12)*(term15-term12)/(point_coo(2)*point_coo(2)*r_i*r_i)) +1);
    double term17 = (term18/(point_coo(2)*(r_j)) - ((point_coo(0)- xj)*(term19 - term13))/(point_coo(2)*r_j*r_j*r_j))/(((term19-term13)*(term19-term13)/(point_coo(2)*point_coo(2)*r_j*r_j)) +1);
    //double term16 = (term14/(point_coo(2)*(r_i)) - ((point_coo(0)- xi)*(term15 - term12))/(point_coo(2)*r_i*r_i*r_i));
    //double term17 = (term18/(point_coo(2)*(r_j)) - ((point_coo(0)- xj)*(term19 - term13))/(point_coo(2)*r_j*r_j*r_j));
     du_dx += term1*(term3 + term2)/term4 - term1*(term3 + term2)/term5  ;
     du_dy += term1*(term6 + term7)/term4 - term1*(term6 + term7)/term5  ;
     du_dz += term1*(term9 + term10)/(term4) - term1*(term9 + term10)/term5;
     dv_dx += term8*(term3 + term2)/term4 - term8*(term3 + term2)/term5  ;
     dv_dy += term8*(term6 + term7)/term4 - term8*(term6 + term7)/term5  ;
     dv_dz += term8*(term9 + term10)/(term4) - term8*(term9 + term10)/term5;
     dw_dx += term16 - term17;
    double term20 = (point_coo(1)- yi)*(term15 - term12)/(point_coo(1)*r_i*r_i*r_i);
    double term21 = (point_coo(0) - xi)/(point_coo(1)*r_i);
    double term22 = (point_coo(1)- yj)*(term19 - term13)/(point_coo(1)*r_j*r_j*r_j);
    double term23 = (point_coo(0) - xj)/(point_coo(1)*r_j);
    double term24 = (-term20-term21)/(((term15 - term12)*(term15 - term12))/((point_coo(2)*point_coo(2)*r_i*r_i)) +1);
    double term25 = (-term22-term23)/(((term19 - term13)*(term19 - term13))/((point_coo(2)*point_coo(2)*r_j*r_j)) +1);
     dw_dy += term24 - term25 ;
    double term26 = 2*(term11)/(r_i);
    double term27 = (term15 - term12)/(point_coo(2)*point_coo(2)*r_i);
    double term28 = (term15 - term12)/(r_i*r_i*r_i);
    double term29 = (term15-term12)*(term15-term12)/(point_coo(2)*point_coo(2)*r_i*r_i) +1;
    double term30 = 2*(term11)/(r_j);
    double term31 = (term19 - term13)/(point_coo(2)*point_coo(2)*r_j);
    double term32 = (term19 - term13)/(r_j*r_j*r_j);
    double term33 = (term19-term13)*(term19-term13)/(point_coo(2)*point_coo(2)*r_j*r_j) +1;
     dw_dz +=  ((term26 - term27 - term28)/(term29)) - ((term30- term31 - term32)/(term33));
*/
        double term1 = (yj - yi) / d_ij;
        double term2 = (point_coo(0) - xi) / (r_i + epsilon);
        double term3 = (point_coo(0) - xj) / (r_j + epsilon);
        double term4 = r_i + r_j - d_ij + epsilon;
        double term5 = r_i + r_j + d_ij + epsilon;
        double term6 = (point_coo(1) - yi) / (r_i + epsilon);
        double term7 = (point_coo(1) - yj) / (r_j + epsilon);
        double term8 = (xi - xj) / d_ij;
        double term9 = point_coo(2) / (r_i + epsilon);
        double term10 = point_coo(2) / (r_j + epsilon);
        double term11 = (yj - yi) / (xj - xi + epsilon);
        double term12 = (point_coo(1) - yi) * (point_coo(0) - xi);
        double term13 = (point_coo(1) - yj) * (point_coo(0) - xj);
        double term14 = 2 * term11 * (point_coo(0) - xi) + yi - point_coo(1);
        double term15 = term11 * ((point_coo(0) - xi) * ((point_coo(0) - xi)) + point_coo(2) * point_coo(2));
        double term18 = 2 * term11 * (point_coo(0) - xj) + yj - point_coo(1);
        double term19 = term11 * ((point_coo(0) - xj) * ((point_coo(0) - xj)) + point_coo(2) * point_coo(2));

        long double term16 = (term14 / (point_coo(2) * (r_i + epsilon)) - ((point_coo(0) - xi) * (term15 - term12)) / (point_coo(2) * r_i * r_i * r_i + epsilon)) / (((term15 - term12) * (term15 - term12) / (point_coo(2) * point_coo(2) * r_i * r_i + epsilon)) + 1);
        long double term17 = (term18 / (point_coo(2) * (r_j + epsilon)) - ((point_coo(0) - xj) * (term19 - term13)) / (point_coo(2) * r_j * r_j * r_j + epsilon)) / (((term19 - term13) * (term19 - term13) / (point_coo(2) * point_coo(2) * r_j * r_j + epsilon)) + 1);
        //long double term16 = (term14 - ((point_coo(0) - xi) * (term15 - term12))) / (((term15 - term12) * (term15 - term12) / (point_coo(2) * point_coo(2) * r_i * r_i + epsilon)) + 1);
        //long double term17 = (term18  - ((point_coo(0) - xj) * (term19 - term13))) / (((term19 - term13) * (term19 - term13) / (point_coo(2) * point_coo(2) * r_j * r_j + epsilon)) + 1);
        du_dx += term1 * (term3 + term2) / term4 - term1 * (term3 + term2) / term5;
        du_dy += term1 * (term6 + term7) / term4 - term1 * (term6 + term7) / term5;
        du_dz += term1 * (term9 + term10) / (term4) - term1 * (term9 + term10) / term5;
        dv_dx += term8 * (term3 + term2) / term4 - term8 * (term3 + term2) / term5;
        dv_dy += term8 * (term6 + term7) / term4 - term8 * (term6 + term7) / term5;
        dv_dz += term8 * (term9 + term10) / (term4) - term8 * (term9 + term10) / term5;
        //dw_dx += term17 ;
        dw_dx += term16 - term17;
        //dw_dx += 0;

        double term20 = (point_coo(1) - yi) * (term15 - term12) / (point_coo(1) * r_i * r_i * r_i + epsilon);
        double term21 = (point_coo(0) - xi) / (point_coo(1) * r_i + epsilon);
        double term22 = (point_coo(1) - yj) * (term19 - term13) / (point_coo(1) * r_j * r_j * r_j + epsilon);
        double term23 = (point_coo(0) - xj) / (point_coo(1) * r_j + epsilon);
        double term24 = (-term20 - term21) / (((term15 - term12) * (term15 - term12)) / ((point_coo(2) * point_coo(2) * r_i * r_i + epsilon)) + 1);
        double term25 = (-term22 - term23) / (((term19 - term13) * (term19 - term13)) / ((point_coo(2) * point_coo(2) * r_j * r_j + epsilon)) + 1);

        dw_dy += term24 - term25;

        double term26 = 2 * (term11) / (r_i + epsilon);
        double term27 = (term15 - term12) / (point_coo(2) * point_coo(2) * r_i + epsilon);
        double term28 = (term15 - term12) / (r_i * r_i * r_i + epsilon);
        double term29 = (term15 - term12) * (term15 - term12) / (point_coo(2) * point_coo(2) * r_i * r_i + epsilon) + 1;
        double term30 = 2 * (term11) / (r_j + epsilon);
        double term31 = (term19 - term13) / (point_coo(2) * point_coo(2) * r_j + epsilon);
        double term32 = (term19 - term13) / (r_j * r_j * r_j + epsilon);
        double term33 = (term19 - term13) * (term19 - term13) / (point_coo(2) * point_coo(2) * r_j * r_j + epsilon) + 1;

        dw_dz += ((term26 - term27 - term28) / (term29)) - ((term30 - term31 - term32) / (term33));
    }


   
   matrix<< du_dx, du_dy, du_dz,
            dv_dx, dv_dy, dv_dz,
            dw_dx, dw_dy, dw_dz;
  //std::cout<<"I am being called:"<<panelId<<std::endl;
    //std::cout<<"The point under consideration is:"<<point_coo.transpose()<<std::endl;
    //std::cout<<"The w derivative values for panel"<<panelId<<"are"<<dw_dx<<" "<<dw_dy<<" "<<dw_dz<<std::endl;
    T<< _cellOblique[0] , _cellTangent[0] , _cellNormals[0] ,
        _cellOblique[1] , _cellTangent[1] , _cellNormals[1] ,
        _cellOblique[2] , _cellTangent[2] , _cellNormals[2] ;
    Eigen::Matrix3d matrixInInertialFrame = T * matrix * T.transpose();

    return matrixInInertialFrame;

}
Eigen::Matrix3d Surface::calculatederivativematrixatpoint(const Eigen::Vector3d &point){
  Eigen::Matrix3d derive_matrix;
   derive_matrix.setZero();
  for (size_t i = 0; i < _cellCenters.size(); i++) {
    derive_matrix +=  _sigma(i)*Calculatederivatesduetosinglepanel(point, i)/(4*M_PI) ;
  }
  return derive_matrix;
}
Eigen::Vector3d Surface::getstretchingterm(const double x_coord,const double y_coord,const double z_coord,const double alpha_x, const double alpha_y, const double alpha_z){
  Eigen::Vector3d changed_vector;
  changed_vector.setZero();

  Eigen::Vector3d vector;
  vector.setZero();
  vector[0]=alpha_x;
  vector[1]=alpha_y;
  vector[2]=alpha_z;
  Eigen::Vector3d point_vector;
  point_vector[0]=x_coord;
  point_vector[1]=y_coord;
  point_vector[2]=z_coord;
  Eigen::Matrix3d matrix = calculatederivativematrixatpoint(point_vector);
  changed_vector[0] = vector[0]*matrix(0,0) + vector[1]*matrix(0,1) + vector[2]*matrix(0,2);
  changed_vector[1] = vector[0]*matrix(0,0) + vector[1]*matrix(0,1) + vector[2]*matrix(0,2);
  changed_vector[2] = vector[0]*matrix(0,0) + vector[1]*matrix(0,1) + vector[2]*matrix(0,2);

return changed_vector;

} 

