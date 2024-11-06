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

#ifndef SURFACE_HPP
#define SURFACE_HPP
#include <Eigen/Core>
#include <string>
#include <vector>


class Surface {
public:
  std::vector<Eigen::Vector3d> _vertices;
  std::vector<Eigen::Vector3i> _connectivity;
  std::vector<Eigen::Vector3d> _cellCenters;
  std::vector<Eigen::Vector3d> _cellNormals;
  std::vector<Eigen::Vector3d> _cellTangent;
  std::vector<Eigen::Vector3d> _cellOblique;
  Eigen::MatrixXd              _AICMatrix;
  Eigen::MatrixXd              _AICMatrixInv;
  Eigen::VectorXd              _sigma;
  Eigen::VectorXd               sigmaget;
  Eigen::Vector3d transformVectorToInertialFrame(Eigen::Vector3d vectorInLocalFrame, size_t panelId);
  Eigen::Vector3d transformPointToPanelLocalFrame(const Eigen::Vector3d& pointInInertialFrame, size_t panelId);
  Eigen::Vector3d calculateVelocityDueToSinglePanel(Eigen::Vector3d point, size_t panelId);
  void            calculateAICMatrix();
  void            calculateAICMatrixInv();
  Eigen::VectorXd getSigma(const std::vector<Eigen::Vector3d> &rhs );
public:
  Eigen::Vector3d calculateVelocityAtPointcoord(const double x_coord,const  double y_coord, const double z_coord);
  Eigen::Matrix3d calculatederivativematrixatpoint(const Eigen::Vector3d &point);
  //Eigen::VectorXd calculateSigma(const std::vector<Eigen::Vector3d> &rhs );
  Eigen::Matrix3d Calculatederivatesduetosinglepanel(Eigen::Vector3d point, size_t panelId);
  void initFromOBJFile(const std::string &filename);
  void translateBuilding(const Eigen::Vector3d &translation);
  void calculateSigma(const std::vector<Eigen::Vector3d> &rhs);
  //Eigen::VectorXd getSigma(const std::vector<Eigen::Vector3d> &rhs);
  std::vector<Eigen::Vector3d> &getCellCenters();
  Eigen::Vector3d               calculateVelocityAtPoint(const Eigen::Vector3d &point);
  Eigen::Vector3d getstretchingterm(const double x_coord,const double y_coord,const double z_coord,const double alpha_x, const double alpha_y, const double alpha_z);
};

#endif // SURFACE_HPP
