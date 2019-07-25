#ifndef _GET_RAW_DATA_H_
#define _GET_RAW_DATA_H_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>
#include<iostream>
#include "Eigen/Dense"
#include "get_raw_data.h"

using namespace Eigen;
using namespace std;

typedef void    (*func_ptr)(const MatrixXd &input, const VectorXd &params, VectorXd &output);
void            get_axis(const MatrixXd &input, const VectorXd &params, VectorXd &output);
void            get_pos(const MatrixXd &input, const VectorXd &params, VectorXd &output);
float           get_angle_acc(Vector3f j1, Vector3f j2, RowVector3f a1, RowVector3f a2, RowVector3f g1, RowVector3f g2, RowVector3f g_dot1, RowVector3f g_dot2, RowVector3f o1, RowVector3f o2);
void            get_jacobian(func_ptr func, const MatrixXd &input, const VectorXd &params, MatrixXd &output );
void            gauss_newton(func_ptr func, const MatrixXd &input, const VectorXd &output, VectorXd &params );
int             imu_joint_pos_data_fit();
int             imu_joint_axis_data_fit();
void            get_raw_data();
void            test_angle();
float           **getData( char filename[] );

#endif