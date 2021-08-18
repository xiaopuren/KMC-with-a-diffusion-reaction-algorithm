//
//  main.cpp
//  TNF-KMC-3D-CPT
//
//  Created by Zhaoqian Su on 11/8/18.
//  Copyright (c) 2018 @AE. All rights reserved.
//

//
//  main.cpp
//  TNF-KMC-3DFinal
//
//  Created by Zhaoqian Su on 4/22/18.
//  Copyright (c) 2018 @AE. All rights reserved.
//

//
//  main.cpp
//  TNF-KMC-04
//
//  Created by Zhaoqian Su on 4/4/18.
//  Copyright (c) 2018 @AE. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <math.h>
#include <iomanip> //needed to use manipulators with parameters (precision, width)
#include <queue>    // storage index of bonded protein index
#include <set>      // storage complex number
#include <vector>
#include <algorithm>


using namespace std;

int simu_step = 20000000 ;//1000000000;
double time_step = 10; //nano second
double distance_step = 10; //Angstrom
double distance_amp;
double cell_range_x = 5773;
double cell_range_y = 5773;
double cell_range_z = 1000;

#define RB_A_tot_num 600
#define protein_A_tot_num 150
#define RB_A_num_per_protein 4
#define RB_A_res_num 4

#define protein_A_tot_num_matrix 151
#define RB_A_num_per_protein_matrix 5
#define RB_A_res_num_matrix 5

#define RB_B_tot_num 200
#define protein_B_tot_num 50
#define RB_B_num_per_protein 4       // 3+1 add a virtual point in the center of protein
#define RB_B_res_num 2

#define protein_B_tot_num_matrix 51
#define RB_B_num_per_protein_matrix 5//  3+1+1 because c++ start from 0, and becuase adding a virtual center of protein RB_B[i][1]
#define RB_B_res_num_matrix 3

#define protein_tot_num 200
#define max_bond_num 900  // 2* protein_A_tot_num

#define protein_tot_num_matrix 201
#define max_bond_num_matrix 901  // 2* protein_A_tot_num +1

double pai = 3.1415926;
double RB_A_radius = 20;     // receptor
double RB_A_D = 1;           // parameter from mebrane protein of plosOne paper
double RB_A_rot_D = 0.0174;

double RB_B_radius = 30;         // Ligands
double RB_B_D = 7.2614;          // Unit: A^2/ns
double RB_B_rot_D = 0.0061209;   // Unit: radian/ns   parameters are calculated in OneNote @AE/TNF

double mono_cis_Ass_Rate = 0.000047;
double mono_cis_Diss_Rate = 0.000000000000112;    // This 0.1 means the monomer cis will dissociate automatically

double cis_D = 0.5;      // two different cis interaction: (1) two single receptors (2)
double cis_rot_D = 0.005;
double cis_Ass_Rate = 0.00096;
double cis_Diss_Rate = 0.000000000000112;

double bond_D = 0.5;
double bond_rot_D = 0.005;
double Ass_Rate = 0.04; //units: per nanosecond
double Diss_Rate = 0.000000000000348;

double bond_dist_cutoff = 18;
double bond_thetapd = 90;
double bond_thetapd_cutoff = 45;
double bond_thetaot = 180;
double bond_thetaot_cutoff = 90;
double cis_thetaot_cutoff = 10;
double cis_dist_cutoff = 15;

double bond_D_cal, bond_rot_D_cal;
double R_x [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];  //matrix size must be constant, so use #define in the above
double R_y [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_z [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_x_0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_y_0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_z_0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_x_new [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];  //matrix in fortran and C++ are different starting from 0 or 1
double R_y_new [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_z_new [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_x_new0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_y_new0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];
double R_z_new0 [protein_tot_num_matrix][RB_A_num_per_protein_matrix][RB_A_res_num_matrix];

int protein_status[protein_tot_num_matrix][RB_A_num_per_protein_matrix];
int protein_status_new[protein_tot_num_matrix][RB_A_num_per_protein_matrix];
int res_nei[protein_tot_num_matrix][7];
int res_nei_new[protein_tot_num_matrix][7];
int seed[protein_tot_num_matrix];
int seed_new[protein_tot_num_matrix];
int visited[protein_tot_num_matrix];
int moved[protein_tot_num_matrix];

int results[protein_tot_num_matrix][protein_tot_num_matrix];  // use BFS to find out all RB index of each complex in the whole system.
int results_prev[protein_tot_num_matrix][protein_tot_num_matrix]; // use results_prev to storage the results of previous step
int results_new[protein_tot_num_matrix][protein_tot_num_matrix]; // use reults_new to find out newly attached proteins on seed
int results_temp[protein_tot_num_matrix][protein_tot_num_matrix]; // use results_prev to storage the results of previous step
int results_temp_size, results_prev_size;
int temp_int;

int complex_index; //complex_index corresponds to row number of results[][]
int complex_sub_index; // subunits index of a complex, from 1 to protein_tot_num
int protein_num_in_Max_Complex;

int bond_num, bond_num_new, bond_num_rl, bond_num_rl_new; //number of bond formed in the simulation
int bond_num_cis, bond_num_cis_new, bond_num_mono_cis, bond_num_mono_cis_new; //number of bond formed in the simulation

double rand2(); //declare random number generator
double gettheta (double point_x[3], double point_y[3], double point_z[3], double test_theta); //declear the function for angle calculation
int i,j,k,n_t, m, n, b, it, jt, kt;
double temp_i,temp_j,temp_k;
double theta,phi,psai,phai;
double dist, dist1, dist2, dist3, dist4;
double cm1_a_x;
double cm1_a_y, cm0_a_y;
double cm1_a_z;
double t[3][3];
double PB_x, PB_y, PB_z;
int initial_simu_time, current_simu_time;
int mc_time_step;
int iteration_mole_step;
int selecting_mole_index; //select one molecule
int RB_A_index, RB_B_index, protein_A_index, protein_B_index, res_B_index1, res_B_index2,res_B_index3,protein_A_index1, protein_A_index2, protein_A_index3, q_index;
double Prob_diff;
int collision_flag;
double point_x[3], point_y[3], point_z[3];
double theta_pd, theta_ot;
double theta_pd2, theta_ot2;
double Prob_Ass, Prob_Diss;
double angle, angle1, angle2;
int selected_bond, selected_pro_A, selected_pro_A2, selected_pro_B, selected_res_B;
double prob;
double temp, angle_x1, angle_x2, angle_y1, angle_y2, dot, det;
int part1,part2,part3,part4,part5,part6,part7,part8,part9;
int complex_size, tot_cluster_num, tot_proteins_in_cluster;    // cluster_size = total proteins in clusters / total cluster_num
double cluster_size;
int protein_A_num_inComplex, protein_B_num_inComplex;
queue<int> q;
std::ofstream parameter ("parameter.log");

#define EPSILON    (1.0E-8)
bool AreSame(double a, double b);



int main (){
    
    //////////////write the input parameter to log file
    std::ofstream parameter ("parameter.log", std::ofstream::app);
    parameter <<setw(25)<< "box size: x y z" <<setw(15)<< cell_range_x << setw(7)<< cell_range_y << setw(7) << cell_range_z<<'\n'<<'\n';
    parameter <<setw(25)<< "protein_A_tot_num" <<setw(15)<< protein_A_tot_num <<'\n';
    parameter <<setw(25)<< "RB_A_tot_num" <<setw(15)<< RB_A_tot_num <<'\n';
    parameter <<setw(25)<< "protein_B_tot_num" <<setw(15)<< protein_B_tot_num <<'\n';
    parameter <<setw(25)<< "RB_B_tot_num" <<setw(15)<< RB_B_tot_num <<'\n'<<'\n';
    
    parameter <<setw(25)<< "RB_A_D" <<setw(15)<< RB_A_D <<'\n';
    parameter <<setw(25)<< "RB_A_rot_D" <<setw(15)<< RB_A_rot_D <<'\n';
    parameter <<setw(25)<< "RB_B_D" <<setw(15)<< RB_B_D <<'\n';
    parameter <<setw(25)<< "RB_B_rot_D" <<setw(15)<< RB_B_rot_D <<'\n'<<'\n';
    
    parameter <<setw(25)<< "R-L interaction:" <<'\n';
    parameter <<setw(25)<< "bond_D" <<setw(15)<< bond_D<<'\n';
    parameter <<setw(25)<< "bond_rot_D" <<setw(15)<< bond_rot_D<<'\n';
    parameter <<setw(25)<< "Ass_Rate" <<setw(15)<< Ass_Rate<<'\n';
    parameter <<setw(25)<< "Diss_Rate" <<setw(15)<< Diss_Rate<<'\n'<<'\n';
    
    parameter <<setw(25)<< "Cis interaction:" <<'\n';
    parameter <<setw(25)<< "cis_D" <<setw(15)<< cis_D<<'\n';
    parameter <<setw(25)<< "cis_rot_D" <<setw(15)<< cis_rot_D<<'\n';
    parameter <<setw(25)<< "mono_cis_Ass_Rate" <<setw(15)<< mono_cis_Ass_Rate<<'\n';
    parameter <<setw(25)<< "mono_cis_Diss_Rate" <<setw(15)<< mono_cis_Diss_Rate<<'\n'<<'\n';
    parameter <<setw(25)<< "cis_Ass_Rate" <<setw(15)<< cis_Ass_Rate<<'\n';
    parameter <<setw(25)<< "cis_Diss_Rate" <<setw(15)<< cis_Diss_Rate<<'\n'<<'\n';
    
    parameter.close();
    
    
    /////////////// initialize bond information
    for (i = 1; i <= protein_A_tot_num; i++){
        protein_status[i][2] = 0;
        protein_status[i][3] = 0;
        res_nei[i][2] = 0;
        res_nei[i][4] = 0;
        res_nei[i][3] = 0;
    }
    for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
        for (j = 1; j <= RB_B_num_per_protein; j++){
            protein_status[i][j] = 0;
            res_nei[i][j] = 0;
        }
    }
    
    
    
    
    std::ifstream check_point ("position.cpt");
    
    if (check_point.is_open()){
        cout << "CPT file is exist\n";
        
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = 1; j <= RB_A_num_per_protein; j++){
                for (k = 1; k <= RB_A_res_num; k++){
                    check_point >> R_x[i][j][k];
                    check_point >> R_y[i][j][k];
                    check_point >> R_z[i][j][k];

                }
            }
            check_point >> protein_status[i][2];
            check_point >> protein_status[i][3];
            check_point >> res_nei[i][2];
            check_point >> res_nei[i][4];
            check_point >> res_nei[i][3];
            
        }
        
        
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            for (j = 1; j <= RB_B_num_per_protein; j++){
                for (k = 1; k <= RB_B_res_num; k++){
                    check_point >> R_x[i][j][k];
                    check_point >> R_y[i][j][k];
                    check_point >> R_z[i][j][k];
                }
                check_point >> protein_status[i][j];
                check_point >> res_nei[i][j];
            }
        }
        
        check_point >> bond_num;
        check_point >> bond_num_rl;
        check_point >> bond_num_cis;
        check_point >> bond_num_mono_cis;
        check_point >> protein_num_in_Max_Complex;
        check_point >> initial_simu_time;
        initial_simu_time = initial_simu_time + 1;    // avoid output same cpt frame twice! previous frame = current frame
    
        
    } /////// end if cpt file exist
      //////  else randomly generate inital configuration
    
    else{
        cout << "CPT file not exist\n";
        std::ofstream ofs ("test.gro");  // open a new test.gro file, delete the previous
        std::ofstream bond ("bond.dat");  // open a new test.gro file, delete the previous
        std::ofstream cluster ("cluster.log");
        std::ofstream check_point ("position.cpt");

        //randomly insert RB_A,  A is membrane molecule in x-y plane
        for (i=1; i <= protein_A_tot_num; i++){
            
        lable1:
            temp_i=rand2()*cell_range_x - cell_range_x/2;
            temp_j=rand2()*cell_range_y - cell_range_y/2;
            temp_k= 0;
            
            //avoid overlap with RB_A
            for (j=1; j <= i-1; j++){
                dist=sqrt( (temp_i-R_x[j][1][1])*(temp_i-R_x[j][1][1])
                          +(temp_j-R_y[j][1][1])*(temp_j-R_y[j][1][1])
                          );
                if (dist <= RB_A_radius + RB_A_radius) {
                    goto lable1;
                }
            }
            
            for (j = 1; j <= RB_A_num_per_protein; j++){
                R_x[i][j][1] = temp_i;
                R_y[i][j][1] = temp_j;
                R_z[i][j][1] = temp_k + (j*2-2)*RB_A_radius; // change radius
                
                R_x_0[i][j][1] = temp_i;
                R_y_0[i][j][1] = temp_j;
                R_z_0[i][j][1] = temp_k + (j*2-2)*RB_A_radius;
                
                R_x_0[i][j][2] = temp_i + RB_A_radius;
                R_y_0[i][j][2] = temp_j ;
                R_z_0[i][j][2] = temp_k + (j*2-2)*RB_A_radius;
                R_x_0[i][j][3] = temp_i - RB_A_radius;
                R_y_0[i][j][3] = temp_j ;
                R_z_0[i][j][3] = temp_k + (j*2-2)*RB_A_radius;
                R_x_0[i][j][4] = temp_i ;
                R_y_0[i][j][4] = temp_j ;
                R_z_0[i][j][4] = temp_k + (j*2-1)*RB_A_radius;
            }
            protein_status[i][2] = 0;
            protein_status[i][3] = 0;
            res_nei[i][2] = 0;  // index of protein_B binding with protein_A
            res_nei[i][4] = 0;  // which rb of protein_B binding with protein_A[i]
            
            res_nei[i][3] = 0;   // cis interaction, index of protein_A
            
            
            
            
            //rotation...
            theta = 0;
            phi = 0;
            psai = (2*rand2() - 1)*pai;
            
            t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
            t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
            t[0][2] = sin(theta)*sin(phi);
            
            t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
            t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
            t[1][2] = -sin(theta)*cos(phi);
            
            t[2][0] = sin(psai)*sin(theta);
            t[2][1] = cos(psai)*sin(theta);
            t[2][2] = cos(theta);
            
            for (j = 1; j <= RB_A_num_per_protein; j++ ){
                for (k = 2; k <= RB_A_res_num; k++){
                    R_x[i][j][k] = t[0][0]*(R_x_0[i][j][k] - R_x[i][j][1]) + t[0][1]*(R_y_0[i][j][k] - R_y[i][j][1]) + t[0][2]*(R_z_0[i][j][k] - R_z[i][j][1]) + R_x[i][j][1];
                    R_y[i][j][k] = t[1][0]*(R_x_0[i][j][k] - R_x[i][j][1]) + t[1][1]*(R_y_0[i][j][k] - R_y[i][j][1]) + t[1][2]*(R_z_0[i][j][k] - R_z[i][j][1]) + R_y[i][j][1];
                    R_z[i][j][k] = t[2][0]*(R_x_0[i][j][k] - R_x[i][j][1]) + t[2][1]*(R_y_0[i][j][k] - R_y[i][j][1]) + t[2][2]*(R_z_0[i][j][k] - R_z[i][j][1]) + R_z[i][j][1];
                }
            }
        }
        
        //randomly insert RB_B
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            
        lable2:
            temp_i=rand2()*cell_range_x - cell_range_x/2;
            temp_j=rand2()*cell_range_y - cell_range_x/2;
            temp_k=rand2()*cell_range_z;
            
            //avoid overlap with RB_A,  use the virtue RB[j][1][1]center
            for (j = 1; j <= protein_A_tot_num; j++){
                for (k = 1; k <= RB_A_num_per_protein; k++){
                    dist=sqrt( (temp_i-R_x[j][k][1])*(temp_i-R_x[j][k][1])
                              +(temp_j-R_y[j][k][1])*(temp_j-R_y[j][k][1])
                              +(temp_k-R_z[j][k][1])*(temp_k-R_z[j][k][1])
                              );
                    if (dist <= RB_A_radius + RB_B_radius*2/sqrt(3) + RB_B_radius) {
                        goto lable2;
                    }
                }
            }
            
            //avoid overlap with RB_B
            for (j = protein_A_tot_num + 1; j <= i-1; j++){
                dist=sqrt( (temp_i-R_x[j][1][1])*(temp_i-R_x[j][1][1])
                          +(temp_j-R_y[j][1][1])*(temp_j-R_y[j][1][1])
                          +(temp_k-R_z[j][1][1])*(temp_k-R_z[j][1][1])
                          );
                if (dist <= RB_B_radius*2/sqrt(3) + RB_B_radius*2/sqrt(3) + 2*RB_B_radius) {
                    goto lable2;
                }
            }
            
            
            R_x[i][1][1] = temp_i;   // virtual point, for calculation only, will not output
            R_y[i][1][1] = temp_j;   // virtual point, for calculation only, will not output
            R_z[i][1][1] = temp_k;   // virtual point, for calculation only, will not output
            
            R_x_0[i][1][2] = temp_i;                 // virtual point, for bonding angle calculation
            R_y_0[i][1][2] = temp_j;                 // virtual point, for bonding angle calculation
            R_z_0[i][1][2] = temp_k + RB_B_radius;   // virtual point, for bonding angle calculation
            
            R_x_0[i][2][1] = temp_i;
            R_y_0[i][2][1] = temp_j + RB_B_radius*2/sqrt(3);
            R_z_0[i][2][1] = temp_k;
            R_x_0[i][3][1] = temp_i - RB_B_radius;
            R_y_0[i][3][1] = temp_j - RB_B_radius/sqrt(3);
            R_z_0[i][3][1] = temp_k;
            R_x_0[i][4][1] = temp_i + RB_B_radius;
            R_y_0[i][4][1] = temp_j - RB_B_radius/sqrt(3);
            R_z_0[i][4][1] = temp_k;
            
            R_x_0[i][2][2] = temp_i;
            R_y_0[i][2][2] = temp_j + RB_B_radius*(2/sqrt(3) + 1);
            R_z_0[i][2][2] = temp_k;
            R_x_0[i][3][2] = temp_i - RB_B_radius*(sqrt(3)/2 + 1);
            R_y_0[i][3][2] = temp_j - RB_B_radius/sqrt(3) - RB_B_radius/2;
            R_z_0[i][3][2] = temp_k;
            R_x_0[i][4][2] = temp_i + RB_B_radius*(sqrt(3)/2 + 1);
            R_y_0[i][4][2] = temp_j - RB_B_radius/sqrt(3) - RB_B_radius/2;
            R_z_0[i][4][2] = temp_k;
            
            for (j = 1; j <= RB_B_num_per_protein; j++ ){
                protein_status[i][j] = 0;
                res_nei[i][j] = 0;   //index of protein_A binding with protein_B[i][j]
            }
            
            
            
            //rotation...
            theta = (2*rand2() - 1)*pai;
            phi = (2*rand2() - 1)*pai;
            psai = (2*rand2() - 1)*pai;
            
            t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
            t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
            t[0][2] = sin(theta)*sin(phi);
            
            t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
            t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
            t[1][2] = -sin(theta)*cos(phi);
            
            t[2][0] = sin(psai)*sin(theta);
            t[2][1] = cos(psai)*sin(theta);
            t[2][2] = cos(theta);
            
            for (j = 1; j <= RB_B_num_per_protein; j++ ){    //j = 2,3,4
                for (k = 1; k <= RB_B_res_num; k++){         //k = 1,2
                    if (j!=1 || k!=1) {                      // skip the central virtual point [i][1][1]  Becarefully! it must be ||, not &&
                        R_x[i][j][k] = t[0][0]*(R_x_0[i][j][k] - R_x[i][1][1]) + t[0][1]*(R_y_0[i][j][k] - R_y[i][1][1]) + t[0][2]*(R_z_0[i][j][k] - R_z[i][1][1]) + R_x[i][1][1];
                        R_y[i][j][k] = t[1][0]*(R_x_0[i][j][k] - R_x[i][1][1]) + t[1][1]*(R_y_0[i][j][k] - R_y[i][1][1]) + t[1][2]*(R_z_0[i][j][k] - R_z[i][1][1]) + R_y[i][1][1];
                        R_z[i][j][k] = t[2][0]*(R_x_0[i][j][k] - R_x[i][1][1]) + t[2][1]*(R_y_0[i][j][k] - R_y[i][1][1]) + t[2][2]*(R_z_0[i][j][k] - R_z[i][1][1]) + R_z[i][1][1];
                    }
                }
            }
        }
        
        
        bond_num = 0;
        bond_num_rl = 0;
        bond_num_cis = 0;
        bond_num_mono_cis = 0;
        protein_num_in_Max_Complex = 0;
        initial_simu_time = 1;
    }
    
    

    /////////////////Begin main loop of Diffusion(part1)-Reaction(part2) simulation
    for (mc_time_step = initial_simu_time; mc_time_step <= simu_step; mc_time_step++){
        
        
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = 1; j <= RB_A_num_per_protein; j++){
                for (k = 1; k <= RB_A_res_num; k++){
                    R_x_new[i][j][k] = R_x[i][j][k];
                    R_y_new[i][j][k] = R_y[i][j][k];
                    R_z_new[i][j][k] = R_z[i][j][k];
                }
            }
            protein_status_new[i][2] = protein_status[i][2];
            protein_status_new[i][3] = protein_status[i][3];
            res_nei_new[i][2] = res_nei[i][2];
            res_nei_new[i][4] = res_nei[i][4];
            
            res_nei_new[i][3] = res_nei[i][3];
            
        }
        
        
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            for (j = 1; j <= RB_B_num_per_protein; j++){
                for (k = 1; k <= RB_B_res_num; k++){
                    R_x_new[i][j][k] = R_x[i][j][k];
                    R_y_new[i][j][k] = R_y[i][j][k];
                    R_z_new[i][j][k] = R_z[i][j][k];
                }
                protein_status_new[i][j] = protein_status[i][j];
                res_nei_new[i][j] = res_nei[i][j];
            }
        }
        
        bond_num_new = bond_num;
        bond_num_rl_new = bond_num_rl;
        bond_num_cis_new = bond_num_cis;
        bond_num_mono_cis_new = bond_num_mono_cis;
        
        
        tot_cluster_num = 0;
        tot_proteins_in_cluster = 0;  // to culculate cluster size
        cluster_size = 0.0;
        
        
        ///////////////////////////////////////////////////////////
        /////////BFS find out all the clusters containing protein_B in the system!
        
        //        cout <<  "///////////////////////////////////////// " <<'\n';
        //        cout <<  "come into matrix " <<'\n';
        
        ///////use matrix instead of vector in BFS
        
        /// initialize the matrix results to 0
        for (i = 0; i <= protein_tot_num; i++){
            for (j = 0; j <= protein_tot_num; j++){
                results[i][j] = 0;
            }
        }
        
        for (i = 1; i <= protein_tot_num; i++){
            visited[i] = 0;
            moved[i] = 0;
        }
        
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            if (visited[i] == 0){
                
                queue<int> qc;   //queue for complex index storage...
                
                j = 0;
                visited[i] = 1;
                qc.push(i);
                
                while(!qc.empty()){
                    
                    j = j+1;
                    results[i][j] = qc.front();    // put front element of queue into matrix results
                    q_index = qc.front();  // find out the index of front element of queue
                    qc.pop();              // pop out this element
                    
                    ////////find neighbors of this element
                    vector<int> neighbors;
                    if (q_index <= protein_A_tot_num){
                        if (res_nei[q_index][2] > 0) {neighbors.push_back(res_nei[q_index][2]);}
                        if (res_nei[q_index][3] > 0) {neighbors.push_back(res_nei[q_index][3]);}
                    }
                    if (q_index > protein_A_tot_num){
                        if (res_nei[q_index][2] > 0) {neighbors.push_back(res_nei[q_index][2]);}
                        if (res_nei[q_index][3] > 0) {neighbors.push_back(res_nei[q_index][3]);}
                        if (res_nei[q_index][4] > 0) {neighbors.push_back(res_nei[q_index][4]);}
                    }
                    
                    ////////iterate neighbors
                    for (std::vector<int>::iterator it = neighbors.begin(); it!=neighbors.end(); ++it){
                        if (visited[*it] == 0){        // if not visited, put in queue
                            visited[*it] = 1;          // mark as visited
                            qc.push(*it);              // put it in queue qc
                        }
                    }
                }
            }
        }
        
        
        /*
         for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
         for (j = 1; j <= protein_tot_num; j++){
         if (results[i][j] != 0){
         cout  << results[i][j]<< "  ";
         }
         }
         cout <<'\n';
         }
         cout <<  "finished results " <<'\n';
         */
        ///////////////// Randomly select one molecule, randomly move every molecule      START part 1 from here!
        for (iteration_mole_step = 1; iteration_mole_step <= protein_A_tot_num + protein_B_tot_num; iteration_mole_step++ ){
            
            selecting_mole_index = iteration_mole_step;
            
            if (selecting_mole_index <= protein_A_tot_num){              //////// if the selected molecule is an enzyme
                protein_A_index = selecting_mole_index;       //make random diffusion for selected molecule  RB_A_index = 0,1,2,3... RB_A_tot_num - 1
                
                if (protein_status_new[protein_A_index][2] == 0 && protein_status_new[protein_A_index][3] == 0){   //if selected molecule is unbound
                    distance_amp = 2*sqrt(RB_A_D*time_step/6)*rand2();  // move this unbound molecule
                    //    theta = rand2()*pai;
                    phai = rand2()*2*pai;
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        for (k = 1; k <= RB_A_res_num; k++){
                            R_x_new0[protein_A_index][j][k] = R_x[protein_A_index][j][k] + distance_amp*cos(phai);
                            R_y_new0[protein_A_index][j][k] = R_y[protein_A_index][j][k] + distance_amp*sin(phai);
                            R_z_new0[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                        }
                    }
                    
                    PB_x = cell_range_x*round(R_x_new0[protein_A_index][1][1]/cell_range_x);
                    PB_y = cell_range_y*round(R_y_new0[protein_A_index][1][1]/cell_range_y);
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        for (k = 1; k <= RB_A_res_num; k++){
                            R_x_new0[protein_A_index][j][k] = R_x_new0[protein_A_index][j][k] - PB_x;
                            R_y_new0[protein_A_index][j][k] = R_y_new0[protein_A_index][j][k] - PB_y;
                        }
                    }
                    
                    
                    ////// rotation
                    theta = 0;
                    phi = 0;
                    psai = (2*rand2() - 1)*sqrt(RB_A_rot_D*time_step);
                    
                    t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
                    t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
                    t[0][2] = sin(theta)*sin(phi);
                    
                    t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
                    t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
                    t[1][2] = -sin(theta)*cos(phi);
                    
                    t[2][0] = sin(psai)*sin(theta);
                    t[2][1] = cos(psai)*sin(theta);
                    t[2][2] = cos(theta);
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        R_x_new[protein_A_index][j][1] = R_x_new0[protein_A_index][j][1];
                        R_y_new[protein_A_index][j][1] = R_y_new0[protein_A_index][j][1];
                        R_z_new[protein_A_index][j][1] = R_z_new0[protein_A_index][j][1];
                        
                        for (k = 2; k <= RB_A_res_num; k++){
                            R_x_new[protein_A_index][j][k] = t[0][0]*(R_x_new0[protein_A_index][j][k] - R_x_new[protein_A_index][j][1]) + t[0][1]*(R_y_new0[protein_A_index][j][k] - R_y_new[protein_A_index][j][1]) + t[0][2]*(R_z_new0[protein_A_index][j][k] - R_z_new[protein_A_index][j][1]) + R_x_new[protein_A_index][j][1];
                            R_y_new[protein_A_index][j][k] = t[1][0]*(R_x_new0[protein_A_index][j][k] - R_x_new[protein_A_index][j][1]) + t[1][1]*(R_y_new0[protein_A_index][j][k] - R_y_new[protein_A_index][j][1]) + t[1][2]*(R_z_new0[protein_A_index][j][k] - R_z_new[protein_A_index][j][1]) + R_y_new[protein_A_index][j][1];
                            R_z_new[protein_A_index][j][k] = t[2][0]*(R_x_new0[protein_A_index][j][k] - R_x_new[protein_A_index][j][1]) + t[2][1]*(R_y_new0[protein_A_index][j][k] - R_y_new[protein_A_index][j][1]) + t[2][2]*(R_z_new0[protein_A_index][j][k] - R_z_new[protein_A_index][j][1]) + R_z_new[protein_A_index][j][1];
                        }
                    }
                    
                    
                    ///////check collision
                    collision_flag = 0;
                    for (i = 1; i <= protein_A_tot_num; i++){      ///check overlap between RB_A and all the other RB_A
                        if (protein_A_index != i){
                            dist=sqrt( (R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])*(R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])
                                      +(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])*(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])
                                      +(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])*(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])
                                      );
                            if (dist < RB_A_radius + RB_A_radius) {
                                collision_flag = 1;
                            }
                        }
                    }
                    
                    for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){        //////check overlap between RB_A with RB_B
                        for (j = 2; j <= RB_B_num_per_protein; j++){
                            for (k = 1; k <= RB_A_num_per_protein; k++){
                                dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])*(R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])
                                          +(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])*(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])
                                          +(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])*(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])
                                          );
                                if (dist < RB_A_radius + RB_B_radius) {
                                    collision_flag = 1;
                                }
                            }
                        }
                    }
                    
                    if (collision_flag == 1){
                        for (j = 1; j <= RB_A_num_per_protein; j++){
                            for (k = 1; k <= RB_A_res_num; k++){
                                R_x_new[protein_A_index][j][k] = R_x[protein_A_index][j][k];
                                R_y_new[protein_A_index][j][k] = R_y[protein_A_index][j][k];
                                R_z_new[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                            }
                        }
                    }
                    
                    
                } //end ---- selected molecule RB_A is unbound
                
                
                //protein_A -- protein_A cis interaction
                
                if (visited[protein_A_index] == 0
                    && res_nei[protein_A_index][2] == 0
                    && protein_A_index == res_nei[res_nei[protein_A_index][3]][3]
                    && res_nei[res_nei[protein_A_index][3]][2] == 0){
                    
                    protein_A_index2 = res_nei[protein_A_index][3];
                    visited[protein_A_index2] = 1;
                    
                    //cout << "proteinA1 " << protein_A_index << "  proteinA2 " << protein_A_index2<<'\n';
                    
                    
                    distance_amp = 2*sqrt(cis_D*time_step/6)*rand2();  // move this unbound molecule
                    //    theta = rand2()*pai;
                    phai = rand2()*2*pai;
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        for (k = 1; k <= RB_A_res_num; k++){
                            R_x_new0[protein_A_index][j][k] = R_x[protein_A_index][j][k] + distance_amp*cos(phai);
                            R_y_new0[protein_A_index][j][k] = R_y[protein_A_index][j][k] + distance_amp*sin(phai);
                            R_z_new0[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                            R_x_new0[protein_A_index2][j][k] = R_x[protein_A_index2][j][k] + distance_amp*cos(phai);
                            R_y_new0[protein_A_index2][j][k] = R_y[protein_A_index2][j][k] + distance_amp*sin(phai);
                            R_z_new0[protein_A_index2][j][k] = R_z[protein_A_index2][j][k];
                        }
                    }
                    
                    
                    PB_x = cell_range_x*round((R_x_new0[protein_A_index][1][1] + R_x_new0[protein_A_index2][1][1])/2/cell_range_x);
                    PB_y = cell_range_y*round((R_y_new0[protein_A_index][1][1] + R_y_new0[protein_A_index2][1][1])/2/cell_range_y);
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        for (k = 1; k <= RB_A_res_num; k++){
                            R_x_new0[protein_A_index][j][k] = R_x_new0[protein_A_index][j][k] - PB_x;
                            R_y_new0[protein_A_index][j][k] = R_y_new0[protein_A_index][j][k] - PB_y;
                            R_x_new0[protein_A_index2][j][k] = R_x_new0[protein_A_index2][j][k] - PB_x;
                            R_y_new0[protein_A_index2][j][k] = R_y_new0[protein_A_index2][j][k] - PB_y;
                        }
                    }
                    
                    
                    
                    ////// rotate protein_A bound in bond
                    theta = 0;
                    phi = 0;
                    psai = (2*rand2() - 1)*sqrt(cis_rot_D*time_step);
                    
                    t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
                    t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
                    t[0][2] = sin(theta)*sin(phi);
                    
                    t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
                    t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
                    t[1][2] = -sin(theta)*cos(phi);
                    
                    t[2][0] = sin(psai)*sin(theta);
                    t[2][1] = cos(psai)*sin(theta);
                    t[2][2] = cos(theta);
                    
                    cm1_a_x = 0;
                    cm1_a_y = 0;
                    cm1_a_z = 0;
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        cm1_a_x = cm1_a_x + R_x_new[protein_A_index][j][1] + R_x_new[protein_A_index2][j][1];
                        cm1_a_y = cm1_a_y + R_y_new[protein_A_index][j][1] + R_y_new[protein_A_index2][j][1];
                        cm1_a_z = cm1_a_z + R_z_new[protein_A_index][j][1] + R_z_new[protein_A_index2][j][1];
                    }
                    
                    // 2 protein_A
                    cm1_a_x = cm1_a_x / (RB_A_num_per_protein*2);
                    cm1_a_y = cm1_a_y / (RB_A_num_per_protein*2);
                    cm1_a_z = cm1_a_z / (RB_A_num_per_protein*2);
                    
                    
                    for (j = 1; j <= RB_A_num_per_protein; j++){
                        for (k = 1; k <= RB_A_res_num; k++){
                            R_x_new[protein_A_index][j][k] = t[0][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[0][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[0][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_x;
                            R_y_new[protein_A_index][j][k] = t[1][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[1][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[1][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_y;
                            R_z_new[protein_A_index][j][k] = t[2][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[2][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[2][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_z;
                            
                            R_x_new[protein_A_index2][j][k] = t[0][0]*(R_x_new0[protein_A_index2][j][k] - cm1_a_x) + t[0][1]*(R_y_new0[protein_A_index2][j][k] - cm1_a_y) + t[0][2]*(R_z_new0[protein_A_index2][j][k] - cm1_a_z) + cm1_a_x;
                            R_y_new[protein_A_index2][j][k] = t[1][0]*(R_x_new0[protein_A_index2][j][k] - cm1_a_x) + t[1][1]*(R_y_new0[protein_A_index2][j][k] - cm1_a_y) + t[1][2]*(R_z_new0[protein_A_index2][j][k] - cm1_a_z) + cm1_a_y;
                            R_z_new[protein_A_index2][j][k] = t[2][0]*(R_x_new0[protein_A_index2][j][k] - cm1_a_x) + t[2][1]*(R_y_new0[protein_A_index2][j][k] - cm1_a_y) + t[2][2]*(R_z_new0[protein_A_index2][j][k] - cm1_a_z) + cm1_a_z;
                        }
                    }
                    
                    
                    /////////////////////////////// relax the complex
                    dist2=sqrt( (R_x_new[protein_A_index][3][3]-R_x_new[protein_A_index2][3][3])
                               *(R_x_new[protein_A_index][3][3]-R_x_new[protein_A_index2][3][3])
                               +(R_y_new[protein_A_index][3][3]-R_y_new[protein_A_index2][3][3])
                               *(R_y_new[protein_A_index][3][3]-R_y_new[protein_A_index2][3][3])
                               );
                    dist1=sqrt( (R_x_new[protein_A_index][3][1]-R_x_new[protein_A_index2][3][1])
                               *(R_x_new[protein_A_index][3][1]-R_x_new[protein_A_index2][3][1])
                               +(R_y_new[protein_A_index][3][1]-R_y_new[protein_A_index2][3][1])
                               *(R_y_new[protein_A_index][3][1]-R_y_new[protein_A_index2][3][1])
                               );
                    dist3 = cis_dist_cutoff/2 + RB_A_radius + RB_A_radius;
                    dist4 = cis_dist_cutoff/2;
                    if ( !AreSame(dist1, dist3) || !AreSame(dist2, dist4)){
                        //      cout << "dist1 =  " <<  dist1 << "   dist2 =  " <<  dist2 <<'\n';
                        //      cout << "dist3 =  " <<  dist3 << "   dist4 =  " <<  dist4 <<'\n';
                        
                        for (j = 1; j <= RB_A_num_per_protein; j++){
                            R_x_new[protein_A_index2][j][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index][3][3] - R_x_new[protein_A_index][3][1]) + R_x_new[protein_A_index][3][3];
                            R_y_new[protein_A_index2][j][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index][3][3] - R_y_new[protein_A_index][3][1]) + R_y_new[protein_A_index][3][3];
                            
                            R_x_new[protein_A_index2][j][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index][3][3] - R_x_new[protein_A_index][3][1]) + R_x_new[protein_A_index][3][3];
                            R_y_new[protein_A_index2][j][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index][3][3] - R_y_new[protein_A_index][3][1]) + R_y_new[protein_A_index][3][3];
                            
                            R_x_new[protein_A_index2][j][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_x_new[protein_A_index][3][3] - R_x_new[protein_A_index][3][1]) + R_x_new[protein_A_index][3][3];
                            R_y_new[protein_A_index2][j][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_y_new[protein_A_index][3][3] - R_y_new[protein_A_index][3][1]) + R_y_new[protein_A_index][3][3];
                            
                            R_x_new[protein_A_index2][j][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index][3][3] - R_x_new[protein_A_index][3][1]) + R_x_new[protein_A_index][3][3];
                            R_y_new[protein_A_index2][j][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index][3][3] - R_y_new[protein_A_index][3][1]) + R_y_new[protein_A_index][3][3];
                        }
                    }
                    
                    
                    
                    ///////check collision
                    collision_flag = 0;
                    
                    for (i = 1; i <= protein_A_tot_num; i++){      ///check overlap between RB_A and all the other RB_A
                        if (protein_A_index != i){
                            dist=sqrt( (R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])*(R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])
                                      +(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])*(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])
                                      +(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])*(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])
                                      );
                            if (dist < RB_A_radius + RB_A_radius) {
                                collision_flag = 1;
                            }
                        }
                        
                        if (protein_A_index2 != i){
                            dist=sqrt( (R_x_new[i][1][1]-R_x_new[protein_A_index2][1][1])*(R_x_new[i][1][1]-R_x_new[protein_A_index2][1][1])
                                      +(R_y_new[i][1][1]-R_y_new[protein_A_index2][1][1])*(R_y_new[i][1][1]-R_y_new[protein_A_index2][1][1])
                                      +(R_z_new[i][1][1]-R_z_new[protein_A_index2][1][1])*(R_z_new[i][1][1]-R_z_new[protein_A_index2][1][1])
                                      );
                            if (dist < RB_A_radius + RB_A_radius) {
                                collision_flag = 1;
                            }
                        }
                    }
                    
                    
                    for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){        //////check overlap between RB_A with RB_B
                        for (j = 2; j <= RB_B_num_per_protein; j++){
                            for (k = 1; k <= RB_A_num_per_protein; k++){
                                dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])*(R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])
                                          +(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])*(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])
                                          +(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])*(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])
                                          );
                                if (dist < RB_A_radius + RB_B_radius) {
                                    collision_flag = 1;
                                }
                                
                                dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_A_index2][k][1])*(R_x_new[i][j][1]-R_x_new[protein_A_index2][k][1])
                                          +(R_y_new[i][j][1]-R_y_new[protein_A_index2][k][1])*(R_y_new[i][j][1]-R_y_new[protein_A_index2][k][1])
                                          +(R_z_new[i][j][1]-R_z_new[protein_A_index2][k][1])*(R_z_new[i][j][1]-R_z_new[protein_A_index2][k][1])
                                          );
                                if (dist < RB_A_radius + RB_B_radius) {
                                    collision_flag = 1;
                                }
                            }
                        }
                    }
                    
                    if (collision_flag == 1){
                        for (j = 1; j <= RB_A_num_per_protein; j++){
                            for (k = 1; k <= RB_A_res_num; k++){
                                R_x_new[protein_A_index][j][k] = R_x[protein_A_index][j][k];
                                R_y_new[protein_A_index][j][k] = R_y[protein_A_index][j][k];
                                R_z_new[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                                
                                R_x_new[protein_A_index2][j][k] = R_x[protein_A_index2][j][k];
                                R_y_new[protein_A_index2][j][k] = R_y[protein_A_index2][j][k];
                                R_z_new[protein_A_index2][j][k] = R_z[protein_A_index2][j][k];
                            }
                        }
                    }
                    
                } //end ---- protein_A -- protein_A cis interaction
                
                
                
            } //end -- selected molecule is an enzyme, all RB_A_index
            
            ///////////////////////////////////////////////////
            ///////////////////////////////////////////////////start protein_B
            
            
            
            
            //// start to traversal each complex
            
            if ((selecting_mole_index > protein_A_tot_num) && (selecting_mole_index <= (protein_A_tot_num + protein_B_tot_num) )){
                
                complex_index = selecting_mole_index;
                complex_size = 0;
                protein_A_num_inComplex = 0;    // prepare to find out how many protein_A in complex
                protein_B_num_inComplex = 0;    // prepare to find out how many protein_B in complex
                
                for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                    if (results[complex_index][complex_sub_index] == 0){break;}
                    if (results[complex_index][complex_sub_index] >  protein_A_tot_num) {protein_B_num_inComplex++;}
                    if (results[complex_index][complex_sub_index] <= protein_A_tot_num) {protein_A_num_inComplex++;}
                    complex_size++;
                }
                //  cout << " complex size := " << complex_size << endl;
                //  cout << " protein_B_num_inComplex := " << protein_B_num_inComplex << endl;
                //  cout << " protein_A_num_inComplex := " << protein_A_num_inComplex << endl;
                
                if (complex_size > protein_num_in_Max_Complex){
                    protein_num_in_Max_Complex = complex_size;
                }
                
                
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // move and rotate one single protein_B
                
                if (complex_size == 1){
                    
                    protein_B_index = results[complex_index][1];
                    
                    distance_amp = 2*sqrt(RB_B_D*time_step/6)*rand2();  // move this unbound molecule
                    theta = rand2()*pai;
                    phai = rand2()*2*pai;
                    
                    for (j = 1; j <= RB_B_num_per_protein; j++){
                        for (k = 1; k <= RB_B_res_num; k++){
                            R_x_new0[protein_B_index][j][k] = R_x[protein_B_index][j][k] + distance_amp*sin(theta)*cos(phai);
                            R_y_new0[protein_B_index][j][k] = R_y[protein_B_index][j][k] + distance_amp*sin(theta)*sin(phai);
                            R_z_new0[protein_B_index][j][k] = R_z[protein_B_index][j][k] + distance_amp*cos(theta);
                        }
                    }
                    
                    PB_x = cell_range_x*round(R_x_new0[protein_B_index][1][1]/cell_range_x);
                    PB_y = cell_range_y*round(R_y_new0[protein_B_index][1][1]/cell_range_y);
                    PB_z = cell_range_z*round(R_z_new0[protein_B_index][1][1]/cell_range_z);
                    
                    if (R_z_new0[protein_B_index][1][1] > cell_range_z || R_z_new0[protein_B_index][1][1] < 0){
                        for (j = 1; j <= RB_B_num_per_protein; j++){
                            for (k = 1; k <= RB_B_res_num; k++){
                                R_z_new0[protein_B_index][j][k] = -R_z_new0[protein_B_index][j][k] + 2*PB_z;
                            }
                        }
                    }
                    
                    for (j = 1; j <= RB_B_num_per_protein; j++){
                        for (k = 1; k <= RB_B_res_num; k++){
                            R_x_new0[protein_B_index][j][k] = R_x_new0[protein_B_index][j][k] - PB_x;
                            R_y_new0[protein_B_index][j][k] = R_y_new0[protein_B_index][j][k] - PB_y;
                        }
                    }
                    
                    
                    ////  rotate unbound RA_B
                    theta = (2*rand2() - 1)*sqrt(RB_B_rot_D*time_step);
                    phi =   (2*rand2() - 1)*sqrt(RB_B_rot_D*time_step);
                    psai =  (2*rand2() - 1)*sqrt(RB_B_rot_D*time_step);
                    
                    t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
                    t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
                    t[0][2] = sin(theta)*sin(phi);
                    
                    t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
                    t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
                    t[1][2] = -sin(theta)*cos(phi);
                    
                    t[2][0] = sin(psai)*sin(theta);
                    t[2][1] = cos(psai)*sin(theta);
                    t[2][2] = cos(theta);
                    
                    R_x_new[protein_B_index][1][1] = R_x_new0[protein_B_index][1][1];
                    R_y_new[protein_B_index][1][1] = R_y_new0[protein_B_index][1][1];
                    R_z_new[protein_B_index][1][1] = R_z_new0[protein_B_index][1][1];
                    
                    for (j = 1; j <= RB_B_num_per_protein; j++){
                        for (k = 1; k <= RB_B_res_num; k++){
                            R_x_new[protein_B_index][j][k] = t[0][0]*(R_x_new0[protein_B_index][j][k] - R_x_new[protein_B_index][1][1]) + t[0][1]*(R_y_new0[protein_B_index][j][k] - R_y_new[protein_B_index][1][1]) + t[0][2]*(R_z_new0[protein_B_index][j][k] - R_z_new[protein_B_index][1][1]) + R_x_new[protein_B_index][1][1];
                            R_y_new[protein_B_index][j][k] = t[1][0]*(R_x_new0[protein_B_index][j][k] - R_x_new[protein_B_index][1][1]) + t[1][1]*(R_y_new0[protein_B_index][j][k] - R_y_new[protein_B_index][1][1]) + t[1][2]*(R_z_new0[protein_B_index][j][k] - R_z_new[protein_B_index][1][1]) + R_y_new[protein_B_index][1][1];
                            R_z_new[protein_B_index][j][k] = t[2][0]*(R_x_new0[protein_B_index][j][k] - R_x_new[protein_B_index][1][1]) + t[2][1]*(R_y_new0[protein_B_index][j][k] - R_y_new[protein_B_index][1][1]) + t[2][2]*(R_z_new0[protein_B_index][j][k] - R_z_new[protein_B_index][1][1]) + R_z_new[protein_B_index][1][1];
                        }
                    }
                } //// end single protein_B
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // rotate one single protein_B complex
                
                if (complex_size > 1){
                    
                    tot_cluster_num = tot_cluster_num + 1;
                    tot_proteins_in_cluster = tot_proteins_in_cluster + complex_size;
                    
                    PB_x = 0;
                    PB_y = 0;
                    PB_z = 0;
                    
                    //bond_D = bond_D - complex_size*0.2;  // the bigger the complex_size, the smaller the diffusion coefficient
                    if (protein_B_num_inComplex == 1){ bond_D_cal = bond_D;}
                    if (protein_B_num_inComplex >= 2){ bond_D_cal = 0;}
                    //      if (protein_B_num_inComplex == 3){ bond_D = 1.0;}
                    //      if (protein_B_num_inComplex >  3){ bond_D = 0;}
                    
                    
                    distance_amp = 2*sqrt(bond_D_cal*time_step/6)*rand2();  // move this bond together
                    phai = rand2()*2*pai;
                    
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {  ///////move protein_A
                            
                            protein_A_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                for (k = 1; k <= RB_A_res_num; k++){
                                    R_x_new0[protein_A_index][j][k] = R_x[protein_A_index][j][k] + distance_amp*cos(phai);
                                    R_y_new0[protein_A_index][j][k] = R_y[protein_A_index][j][k] + distance_amp*sin(phai);
                                    R_z_new0[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                                }
                            }
                            PB_x = PB_x + R_x_new0[protein_A_index][1][1];
                            PB_y = PB_y + R_y_new0[protein_A_index][1][1];
                        }
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num){   ////move protein_B
                            
                            protein_B_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_B_num_per_protein; j++){
                                for (k = 1; k <= RB_B_res_num; k++){
                                    R_x_new0[protein_B_index][j][k] = R_x[protein_B_index][j][k] + distance_amp*cos(phai);
                                    R_y_new0[protein_B_index][j][k] = R_y[protein_B_index][j][k] + distance_amp*sin(phai);
                                    R_z_new0[protein_B_index][j][k] = R_z[protein_B_index][j][k];
                                }
                            }
                            PB_x = PB_x + R_x_new0[protein_B_index][1][1];
                            PB_y = PB_y + R_y_new0[protein_B_index][1][1];
                        }
                        
                    }// end traversal all A and B, and move them together
                    
                    PB_x = cell_range_x*round(PB_x/(protein_A_num_inComplex + protein_B_num_inComplex)/cell_range_x);
                    PB_y = cell_range_y*round(PB_y/(protein_A_num_inComplex + protein_B_num_inComplex)/cell_range_y);
                    
                    cm1_a_x = 0;
                    cm1_a_y = 0;
                    cm1_a_z = 0;
                    
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {
                            protein_A_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                for (k = 1; k <= RB_A_res_num; k++){
                                    R_x_new0[protein_A_index][j][k] = R_x_new0[protein_A_index][j][k] - PB_x;
                                    R_y_new0[protein_A_index][j][k] = R_y_new0[protein_A_index][j][k] - PB_y;
                                }
                            }
                            
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                cm1_a_x = cm1_a_x + R_x_new0[protein_A_index][j][1];
                                cm1_a_y = cm1_a_y + R_y_new0[protein_A_index][j][1];
                                cm1_a_z = cm1_a_z + R_z_new0[protein_A_index][j][1];
                            }
                        }
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num){
                            protein_B_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_B_num_per_protein; j++){
                                for (k = 1; k <= RB_B_res_num; k++){
                                    R_x_new0[protein_B_index][j][k] = R_x_new0[protein_B_index][j][k] - PB_x;
                                    R_y_new0[protein_B_index][j][k] = R_y_new0[protein_B_index][j][k] - PB_y;
                                }
                            }
                            for (j = 1; j <= RB_B_num_per_protein; j++){
                                cm1_a_x = cm1_a_x + R_x_new0[protein_B_index][j][1];
                                cm1_a_y = cm1_a_y + R_y_new0[protein_B_index][j][1];
                                cm1_a_z = cm1_a_z + R_z_new0[protein_B_index][j][1];
                            }
                        }
                    }  // end PBC checking
                    
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    
                    //find rotation center//center of mass of complex
                    cm1_a_x = cm1_a_x / (RB_A_num_per_protein*protein_A_num_inComplex + RB_B_num_per_protein*protein_B_num_inComplex);
                    cm1_a_y = cm1_a_y / (RB_A_num_per_protein*protein_A_num_inComplex + RB_B_num_per_protein*protein_B_num_inComplex);
                    cm1_a_z = cm1_a_z / (RB_A_num_per_protein*protein_A_num_inComplex + RB_B_num_per_protein*protein_B_num_inComplex);
                    
                    ////// prepare for rotational move
                    //bond_rot_D = bond_rot_D/complex_size;   // the bigger the complex_size, the smaller the diffusion coefficient
                    if (protein_B_num_inComplex == 1){ bond_rot_D_cal = bond_rot_D;}
                    if (protein_B_num_inComplex >= 2){ bond_rot_D_cal = 0;}
                    //        if (protein_B_num_inComplex == 3){ bond_rot_D = 0.001;}
                    //        if (protein_B_num_inComplex >  3){ bond_rot_D = 0;}
                    
                    theta = 0;
                    phi = 0;
                    psai = (2*rand2() - 1)*sqrt(bond_rot_D_cal*time_step);
                    
                    t[0][0] = cos(psai)*cos(phi) - cos(theta)*sin(phi)*sin(psai);
                    t[0][1] = -sin(psai)*cos(phi) - cos(theta)*sin(phi)*cos(psai);
                    t[0][2] = sin(theta)*sin(phi);
                    
                    t[1][0] = cos(psai)*sin(phi) + cos(theta)*cos(phi)*sin(psai);
                    t[1][1] = -sin(psai)*sin(phi) + cos(theta)*cos(phi)*cos(psai);
                    t[1][2] = -sin(theta)*cos(phi);
                    
                    t[2][0] = sin(psai)*sin(theta);
                    t[2][1] = cos(psai)*sin(theta);
                    t[2][2] = cos(theta);
                    
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {
                            protein_A_index = results[complex_index][complex_sub_index];
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                for (k = 1; k <= RB_A_res_num; k++){
                                    R_x_new[protein_A_index][j][k] = t[0][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[0][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[0][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_x;
                                    R_y_new[protein_A_index][j][k] = t[1][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[1][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[1][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_y;
                                    R_z_new[protein_A_index][j][k] = t[2][0]*(R_x_new0[protein_A_index][j][k] - cm1_a_x) + t[2][1]*(R_y_new0[protein_A_index][j][k] - cm1_a_y) + t[2][2]*(R_z_new0[protein_A_index][j][k] - cm1_a_z) + cm1_a_z;
                                }
                            }
                        }
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num){
                            protein_B_index = results[complex_index][complex_sub_index];
                            for (j = 1; j <= RB_B_num_per_protein; j++){
                                for (k = 1; k <= RB_B_res_num; k++){
                                    R_x_new[protein_B_index][j][k] = t[0][0]*(R_x_new0[protein_B_index][j][k] - cm1_a_x) + t[0][1]*(R_y_new0[protein_B_index][j][k] - cm1_a_y) + t[0][2]*(R_z_new0[protein_B_index][j][k] - cm1_a_z) + cm1_a_x;
                                    R_y_new[protein_B_index][j][k] = t[1][0]*(R_x_new0[protein_B_index][j][k] - cm1_a_x) + t[1][1]*(R_y_new0[protein_B_index][j][k] - cm1_a_y) + t[1][2]*(R_z_new0[protein_B_index][j][k] - cm1_a_z) + cm1_a_y;
                                    R_z_new[protein_B_index][j][k] = t[2][0]*(R_x_new0[protein_B_index][j][k] - cm1_a_x) + t[2][1]*(R_y_new0[protein_B_index][j][k] - cm1_a_y) + t[2][2]*(R_z_new0[protein_B_index][j][k] - cm1_a_z) + cm1_a_z;
                                }
                            }
                        }
                        
                    }  // end rotation
                    ///////////////////////////////////////////////////////////////////////////////////
                    /////////////finished translational and rotational move
                }  // end if complex_size > 1
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// put protein_B down and align protein_A with cis interaction
                
                if (complex_size > 1){
                    
                    if (protein_B_num_inComplex == 1
                        && (R_z_new[protein_B_index][1][2]!= (R_z_new[protein_B_index][1][1] + RB_B_radius))){
                        // if protein_B is not laid down.  Just come in this loop in the 1st time !!! afterwards it's already laid down
                        //// protein_B_index already got its value in previous translational and rotational move
                        
                        for (j = 1; j <= RB_B_num_per_protein; j++){
                            for (k = 1; k <= RB_B_res_num; k++){
                                R_z_new[protein_B_index][j][k] = R_z_new[protein_A_index][3][1];
                            }
                        }
                        R_z_new[protein_B_index][1][2] = R_z_new[protein_A_index][3][1] + RB_B_radius;
                        
                        angle = atan2((R_x_new[protein_B_index][2][1] - R_x_new[protein_B_index][1][1]), (R_y_new[protein_B_index][2][1] - R_y_new[protein_B_index][1][1])) + pai;
                        //         cout <<  angle <<'\n';
                        
                        
                        // ghost protein_B
                        R_x_0[protein_B_index][1][1] = 0;
                        R_y_0[protein_B_index][1][1] = 0;
                        
                        R_x_0[protein_B_index][1][2] = 0;
                        R_y_0[protein_B_index][1][2] = 0;
                        
                        R_x_0[protein_B_index][2][1] = 0;
                        R_y_0[protein_B_index][2][1] = RB_B_radius*2/sqrt(3);
                        
                        R_x_0[protein_B_index][2][2] = 0;
                        R_y_0[protein_B_index][2][2] = RB_B_radius*(2/sqrt(3) + 1);  // done, not move
                        
                        R_x_0[protein_B_index][3][1] = - RB_B_radius;
                        R_y_0[protein_B_index][3][1] = - RB_B_radius/sqrt(3);
                        
                        R_x_0[protein_B_index][3][2] = - RB_B_radius*(sqrt(3)/2 + 1);
                        R_y_0[protein_B_index][3][2] = - RB_B_radius/sqrt(3) - RB_B_radius/2;
                        
                        R_x_0[protein_B_index][4][1] =   RB_B_radius;
                        R_y_0[protein_B_index][4][1] = - RB_B_radius/sqrt(3);
                        
                        R_x_0[protein_B_index][4][2] =   RB_B_radius*(sqrt(3)/2 + 1);
                        R_y_0[protein_B_index][4][2] = - RB_B_radius/sqrt(3) - RB_B_radius/2;
                        
                        cm1_a_x = R_x_new[protein_B_index][1][1];
                        cm1_a_y = R_y_new[protein_B_index][1][1];
                        
                        for (j = 1; j <= RB_B_num_per_protein; j++){
                            for (k = 1; k <= RB_B_res_num; k++){
                                R_x_new[protein_B_index][j][k] = R_x_0[protein_B_index][j][k]*cos(angle) - R_y_0[protein_B_index][j][k]*sin(angle) + cm1_a_x;
                                R_y_new[protein_B_index][j][k] = R_x_0[protein_B_index][j][k]*sin(angle) + R_y_0[protein_B_index][j][k]*cos(angle) + cm1_a_y;
                            }
                        }
                        
                        
                    } // end if protein_B_num_inComplex == 1 and attached with protein_A
                } // end if complex_size > 1, lay down one single protein_B process
                /////////////////// finished lay down protein_B only, next step : align attached protein_A
                
                if (complex_size > 1 && protein_B_num_inComplex == 1) {
                    
                    for (j = 2; j <= RB_B_num_per_protein; j++){           // becareful!! just for j = 2 3 4
                        if (res_nei_new[protein_B_index][j] != 0){        ///////////// find out attached protein_A index
                            protein_A_index1 = res_nei_new[protein_B_index][j];
                            
                            //cout << "proteinA1 " << protein_A_index1 <<'\n';
                            
                            
                            dist2=sqrt( (R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                       *(R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                       +(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                       *(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                       );
                            dist1=sqrt( (R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                       *(R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                       +(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                       *(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                       );
                            if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                for (k = 1; k <= RB_A_num_per_protein; k++){
                                    R_x_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                    R_y_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                    
                                    R_x_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                    R_y_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                    
                                    R_x_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                    R_y_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                    
                                    R_x_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                    R_y_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                }
                            } // end move protein_A_index1
                        }   // end if (res_nei_new[protein_B_index][j] != 0)
                        
                    } // end traveral each residue of protein_B [2] [3] [4]
                }// end     if (complex_size > 1 && protein_B_num_inComplex == 1)
                ////////////////// finished align attached protein_A, next step to consider protein_A's cis interaction
                
                //////////////////// start to align cis protein_A (2nd attached protein_A)
                if (complex_size > 1 && protein_B_num_inComplex == 1) {
                    
                    for (j = 2; j <= RB_B_num_per_protein; j++){           // becareful!! just for j = 2 3 4
                        if (res_nei_new[protein_B_index][j] != 0 && res_nei_new[res_nei_new[protein_B_index][j]][3] != 0){        ///////////// find out attached protein_A index
                            
                            protein_A_index1 = res_nei_new[protein_B_index][j];
                            protein_A_index2 = res_nei_new[protein_A_index1][3];
                            
                            dist2=sqrt( (R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                       *(R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                       +(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                       *(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                       );
                            dist1=sqrt( (R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                       *(R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                       +(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                       *(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                       );
                            if ( !AreSame(dist1,cis_dist_cutoff/2 + RB_A_radius + RB_A_radius) || !AreSame(dist2, cis_dist_cutoff/2)){
                                for (k = 1; k <= RB_A_num_per_protein; k++){
                                    R_x_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                    R_y_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                    
                                    R_x_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                    R_y_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                    
                                    R_x_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                    R_y_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                    
                                    R_x_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                    R_y_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                }
                            } // finished move cis protein_A
                            
                            
                        }
                    }
                }
                /////////// finished align cis protein_A
                
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                //step 0 : put protein_A attached with it's protein_B
                if (protein_B_num_inComplex > 1) {
                    random_shuffle(&results[complex_index][1],&results[complex_index][complex_size]);  /// to avoid collision during alignment, randomly iterate this array
                    
                    for ( complex_sub_index= 1; complex_sub_index <= complex_size; complex_sub_index++){
                        //if (results[complex_index][complex_sub_index] == 0){continue;}
                        //cout << "randomly iterate ?   " << results[complex_index][complex_sub_index] <<'\n';
                        
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {  ///////move protein_A
                            
                            protein_A_index1 = results[complex_index][complex_sub_index];
                            
                            if (res_nei_new[protein_A_index1][2]!=0){
                                protein_B_index = res_nei_new[protein_A_index1][2];
                                j = res_nei_new[protein_A_index1][4]; // corresopnding protein_B residue index
                                
                                dist2=sqrt( (R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                           *(R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                           +(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                           *(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                           );
                                dist1=sqrt( (R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                           *(R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                           +(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                           *(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                           );
                                if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                    moved[protein_A_index1] = 1;  //////////////////////// mark moved!
                                    
                                    for (k = 1; k <= RB_A_num_per_protein; k++){
                                        R_x_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                    }
                                    
                                    
                                }
                            }
                        }
                    }
                }
                
                
                
                
                
                
                
                
            lable3:
                /////////// step 1 : align two proteinA with cis interaction
                if (protein_B_num_inComplex > 1) {
                    
                    random_shuffle(&results[complex_index][1],&results[complex_index][complex_size]);  /// to avoid collision during alignment, randomly iterate this array
                    
                    for ( complex_sub_index= 1; complex_sub_index <= complex_size; complex_sub_index++){
                        //if (results[complex_index][complex_sub_index] == 0){continue;}
                        //cout << "randomly iterate ?   " << results[complex_index][complex_sub_index] <<'\n';
                        
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {  ///////move protein_A
                            
                            protein_A_index = results[complex_index][complex_sub_index];
                            
                            if (    res_nei_new[protein_A_index][2]!=0
                                && (res_nei_new[protein_A_index][3]!=0)
                                && (res_nei_new[res_nei_new[protein_A_index][3]][2]!=0)
                                &&  moved[protein_A_index] == 0
                                ){
                                
                                protein_A_index1 = results[complex_index][complex_sub_index];
                                protein_A_index2 = res_nei_new[protein_A_index1][3];
                                
                                moved[protein_A_index1] = 1;
                                moved[protein_A_index2] = 1;
                                
                                //            for (j = 2; j <= RB_B_num_per_protein; j++){
                                //                moved[res_nei_new[res_nei_new[protein_A_index][2]][j]] = 1;
                                //            }
                                
                                dist1=sqrt( (R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                           *(R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                           +(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                           *(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                           );
                                
                                dist2=sqrt( (R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                           *(R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                           +(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                           *(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                           );
                                
                                if (!AreSame(dist1, cis_dist_cutoff/2 + RB_A_radius + RB_A_radius) || !AreSame(dist2, cis_dist_cutoff/2)){
                                    //cout << "inside lable3   " << protein_B_num_inComplex <<'\n';
                                    /////////////////// mark mistake !!! should be : dist1, cis_dist_cutoff/2 + RB_A_radius + RB_A_radius
                                    
                                    
                                    for (k = 1; k <= RB_A_num_per_protein; k++){
                                        R_x_new[protein_A_index1][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index2][3][3] - R_x_new[protein_A_index2][3][1]) + R_x_new[protein_A_index2][3][3];
                                        R_y_new[protein_A_index1][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index2][3][3] - R_y_new[protein_A_index2][3][1]) + R_y_new[protein_A_index2][3][3];
                                        
                                        R_x_new[protein_A_index1][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index2][3][3] - R_x_new[protein_A_index2][3][1]) + R_x_new[protein_A_index2][3][3];
                                        R_y_new[protein_A_index1][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index2][3][3] - R_y_new[protein_A_index2][3][1]) + R_y_new[protein_A_index2][3][3];
                                        
                                        R_x_new[protein_A_index1][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_x_new[protein_A_index2][3][3] - R_x_new[protein_A_index2][3][1]) + R_x_new[protein_A_index2][3][3];
                                        R_y_new[protein_A_index1][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_y_new[protein_A_index2][3][3] - R_y_new[protein_A_index2][3][1]) + R_y_new[protein_A_index2][3][3];
                                        
                                        R_x_new[protein_A_index1][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index2][3][3] - R_x_new[protein_A_index2][3][1]) + R_x_new[protein_A_index2][3][3];
                                        R_y_new[protein_A_index1][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index2][3][3] - R_y_new[protein_A_index2][3][1]) + R_y_new[protein_A_index2][3][3];
                                    }
                                } // finished move cis protein_A
                            }
                        }  // end traversal and find each protein_A_index
                    } /// end traversal this complex, traversal and find each protein_A_index
                } //////////  end if protein_B_num_inComplex > 1
                //////////////// finished align protein_A_index with cis interaction !!!
                
                //   lable4:
                ///////////// step 2 : put protein_B_index attached with protein_A having cis interaction
                if (protein_B_num_inComplex > 1) {
                    
                    random_shuffle(&results[complex_index][1],&results[complex_index][complex_size]);  /// to avoid collision during alignment, randomly iterate this array
                    for ( complex_sub_index= 1; complex_sub_index <= complex_size; complex_sub_index++){
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num) {  ///////move protein_A
                            
                            protein_B_index = results[complex_index][complex_sub_index];
                            for (j = 2; j <= RB_B_num_per_protein; j++){           // becareful!! just for j = 2 3 4
                                if (   res_nei_new[protein_B_index][j] != 0         // protein_B has a neighbor protein_A1
                                    && res_nei_new[res_nei_new[protein_B_index][j]][3] != 0     // protein_A1 has a cis neighbor protein_A2  // protein_A2 has the other neighbor protein_B
                                    && res_nei_new[res_nei_new[res_nei_new[protein_B_index][j]][3]][2] != 0
                                    && moved[protein_B_index] == 0
                                    ){        ///////////// find out attached protein_A index
                                    
                                    protein_A_index1 = res_nei_new[protein_B_index][j];
                                    
                                    dist2=sqrt( (R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                               *(R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                               +(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                               *(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                               );
                                    dist1=sqrt( (R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                               *(R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                               +(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                               *(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                               );
                                lable4:
                                    if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                        
                                        moved[protein_B_index] = 1;
                                        
                                        for (k = 1; k <= RB_B_res_num; k++){
                                            R_z_new[protein_B_index][1][k] = R_z_new[protein_A_index1][3][1];
                                            R_z_new[protein_B_index][2][k] = R_z_new[protein_A_index1][3][1];
                                            R_z_new[protein_B_index][3][k] = R_z_new[protein_A_index1][3][1];
                                            R_z_new[protein_B_index][4][k] = R_z_new[protein_A_index1][3][1];
                                            
                                        }
                                        R_z_new[protein_B_index][1][2] = R_z_new[protein_A_index1][3][1] + RB_B_radius;
                                        
                                        // ghost protein_B
                                        R_x_0[protein_B_index][1][1] = 0;
                                        R_y_0[protein_B_index][1][1] = 0;
                                        
                                        R_x_0[protein_B_index][1][2] = 0;
                                        R_y_0[protein_B_index][1][2] = 0;
                                        
                                        R_x_0[protein_B_index][2][1] = 0;
                                        R_y_0[protein_B_index][2][1] = RB_B_radius*2/sqrt(3);
                                        
                                        R_x_0[protein_B_index][2][2] = 0;
                                        R_y_0[protein_B_index][2][2] = RB_B_radius*(2/sqrt(3) + 1);  // done, not move
                                        
                                        R_x_0[protein_B_index][3][1] = - RB_B_radius;
                                        R_y_0[protein_B_index][3][1] = - RB_B_radius/sqrt(3);
                                        
                                        R_x_0[protein_B_index][3][2] = - RB_B_radius*(sqrt(3)/2 + 1);
                                        R_y_0[protein_B_index][3][2] = - RB_B_radius/sqrt(3) - RB_B_radius/2;
                                        
                                        R_x_0[protein_B_index][4][1] =   RB_B_radius;
                                        R_y_0[protein_B_index][4][1] = - RB_B_radius/sqrt(3);
                                        
                                        R_x_0[protein_B_index][4][2] =   RB_B_radius*(sqrt(3)/2 + 1);
                                        R_y_0[protein_B_index][4][2] = - RB_B_radius/sqrt(3) - RB_B_radius/2;
                                        
                                        /////////////modify the following:
                                        // angle between  point 1 counterclokewise to point 2 !!
                                        angle_x1 = R_x_0[protein_B_index][j][1];
                                        angle_y1 = R_y_0[protein_B_index][j][1];
                                        angle_x2 = R_x_new[protein_A_index1][3][1] - R_x_new[protein_A_index1][3][2];
                                        angle_y2 = R_y_new[protein_A_index1][3][1] - R_y_new[protein_A_index1][3][2];
                                        dot = angle_x1*angle_x2 + angle_y1*angle_y2 ;
                                        det = angle_x1*angle_y2 - angle_y1*angle_x2 ;
                                        
                                        angle = atan2(-det, -dot) + pai;
                                        
                                        
                                        
                                        
                                        cm1_a_x = (bond_dist_cutoff/2 + RB_B_radius*2/sqrt(3) + RB_B_radius)/RB_A_radius*(R_x_new[protein_A_index1][3][2]-R_x_new[protein_A_index1][3][1])
                                        + R_x_new[protein_A_index1][3][2];
                                        cm1_a_y = (bond_dist_cutoff/2 + RB_B_radius*2/sqrt(3) + RB_B_radius)/RB_A_radius*(R_y_new[protein_A_index1][3][2]-R_y_new[protein_A_index1][3][1])
                                        + R_y_new[protein_A_index1][3][2];;
                                        
                                        for (m = 1; m <= RB_B_num_per_protein; m++){
                                            for (n = 1; n <= RB_B_res_num; n++){
                                                R_x_new[protein_B_index][m][n] = R_x_0[protein_B_index][m][n]*cos(angle) - R_y_0[protein_B_index][m][n]*sin(angle) + cm1_a_x;
                                                R_y_new[protein_B_index][m][n] = R_x_0[protein_B_index][m][n]*sin(angle) + R_y_0[protein_B_index][m][n]*cos(angle) + cm1_a_y;
                                            }
                                        }
                                        
                                        ///////////// finished moving protein_B
                                        ////////////////////////////////////////////////////////////////////////////////////
                                        
                                        ////////////////////////////////////////////////////////
                                        /////////// start to align protein_A (res_nei_new[protein_B_index][j]) attached with protein_B_j
                                        for (m = 2; m <= RB_B_num_per_protein; m++){
                                            
                                            
                                            protein_A_index1 = res_nei_new[protein_B_index][m];
                                            
                                            
                                            if (res_nei_new[protein_A_index1][2]!=0){
                                                protein_B_index = res_nei_new[protein_A_index1][2];
                                                n = res_nei_new[protein_A_index1][4]; // corresopnding protein_B residue index
                                                
                                                dist2=sqrt( (R_x_new[protein_B_index][n][2]-R_x_new[protein_A_index1][3][2])
                                                           *(R_x_new[protein_B_index][n][2]-R_x_new[protein_A_index1][3][2])
                                                           +(R_y_new[protein_B_index][n][2]-R_y_new[protein_A_index1][3][2])
                                                           *(R_y_new[protein_B_index][n][2]-R_y_new[protein_A_index1][3][2])
                                                           );
                                                dist1=sqrt( (R_x_new[protein_B_index][n][1]-R_x_new[protein_A_index1][3][1])
                                                           *(R_x_new[protein_B_index][n][1]-R_x_new[protein_A_index1][3][1])
                                                           +(R_y_new[protein_B_index][n][1]-R_y_new[protein_A_index1][3][1])
                                                           *(R_y_new[protein_B_index][n][1]-R_y_new[protein_A_index1][3][1])
                                                           );
                                                if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                                    moved[protein_A_index1] = 1;  //////////////////////// mark moved!
                                                    
                                                    
                                                    for (k = 1; k <= RB_A_num_per_protein; k++){
                                                        R_x_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][n][2] - R_x_new[protein_B_index][n][1]) + R_x_new[protein_B_index][n][2];
                                                        R_y_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][n][2] - R_y_new[protein_B_index][n][1]) + R_y_new[protein_B_index][n][2];
                                                        
                                                        R_x_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][n][2] - R_x_new[protein_B_index][n][1]) + R_x_new[protein_B_index][n][2];
                                                        R_y_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][n][2] - R_y_new[protein_B_index][n][1]) + R_y_new[protein_B_index][n][2];
                                                        
                                                        R_x_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][n][2] - R_x_new[protein_B_index][n][1]) + R_x_new[protein_B_index][n][2];
                                                        R_y_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][n][2] - R_y_new[protein_B_index][n][1]) + R_y_new[protein_B_index][n][2];
                                                        
                                                        R_x_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_x_new[protein_B_index][n][2] - R_x_new[protein_B_index][n][1]) + R_x_new[protein_B_index][n][2];
                                                        R_y_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_y_new[protein_B_index][n][2] - R_y_new[protein_B_index][n][1]) + R_y_new[protein_B_index][n][2];
                                                    }
                                                }
                                                
                                                //// inside protein_A_index1 loop check(align) 2nd protein_A (cis protein_A)
                                                if (res_nei_new[protein_A_index1][3]!=0){
                                                    protein_A_index2 = res_nei_new[protein_A_index1][3];
                                                    
                                                    dist2=sqrt( (R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                                               *(R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                                               +(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                                               *(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                                               );
                                                    dist1=sqrt( (R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                                               *(R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                                               +(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                                               *(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                                               );
                                                    if (!AreSame(dist1, cis_dist_cutoff/2 + RB_A_radius + RB_A_radius) || !AreSame(dist2, cis_dist_cutoff/2)){
                                                        moved[protein_A_index2] = 1;  //////////////////////// mark moved!
                                                        
                                                        for (k = 1; k <= RB_A_num_per_protein; k++){
                                                            R_x_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                                            R_y_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                                            
                                                            R_x_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                                            R_y_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                                            
                                                            R_x_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                                            R_y_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                                            
                                                            R_x_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                                            R_y_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                                        }
                                                        
                                                    }
                                                    
                                                }//finished align second protein_A within 1st protein_A loop
                                            } // finished align 1st protein_A
                                            
                                        }/////////// finished align protein_A (res_nei_new[protein_B_index][n]) attached with protein_B_n
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                ///////////////////finished: put protein_B attached (match) with protein_A having cis interaction
                
                
                //////////////////////repeat lable 4
                if (protein_B_num_inComplex > 1) {
                    
                    random_shuffle(&results[complex_index][1],&results[complex_index][complex_size]);  /// to avoid collision during alignment, randomly iterate this array
                    for ( complex_sub_index= 1; complex_sub_index <= complex_size; complex_sub_index++){
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num) {
                            
                            protein_B_index = results[complex_index][complex_sub_index];
                            for (j = 2; j <= RB_B_num_per_protein; j++){           // becareful!! just for j = 2 3 4
                                if (   res_nei_new[protein_B_index][j] != 0         // protein_B has a neighbor protein_A1
                                    && res_nei_new[res_nei_new[protein_B_index][j]][3] != 0     // protein_A1 has a cis neighbor protein_A2  // protein_A2 has the other neighbor protein_B
                                    && res_nei_new[res_nei_new[res_nei_new[protein_B_index][j]][3]][2] != 0
                                    && moved[protein_B_index] == 0
                                    ){        ///////////// find out attached protein_A index
                                    
                                    protein_A_index1 = res_nei_new[protein_B_index][j];
                                    
                                    dist2=sqrt( (R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                               *(R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                               +(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                               *(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                               );
                                    dist1=sqrt( (R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                               *(R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                               +(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                               *(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                               );
                                    if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                        //              cout << "prepare to repeat   " << protein_B_num_inComplex <<'\n';
                                        //              cout << "dist1 =  " <<  dist1 << "   dist2 =  " <<  dist2 <<'\n';
                                        //              cout << "dist =  " <<  cis_dist_cutoff/2 + RB_A_radius + RB_A_radius << "   dist =  " <<  cis_dist_cutoff/2 <<'\n';
                                        //              cout << "Bindex " <<  protein_B_index << "   Aindex1  " <<  protein_A_index1 <<'\n';
                                        
                                        goto lable4;
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                
                //////////////////////////////////////////////////////////////////////////////////////////////
                // already finished the alignment !! But, cometime the complex freeze due to collision problem
                // to solve this proble, we change two protein_A, to make the alignment the other way !
                
                ////////////////// finished: put 1st protein_A(without cis interaction) back to protein_B
                ////////////////// finished: put 2nd protein_A (with cis interaction, but not attached to protein_B) back to its cis protein_A
                
                
                if (protein_B_num_inComplex > 1) {
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {  ///////move protein_A
                            
                            protein_A_index1 = results[complex_index][complex_sub_index];
                            
                            if (res_nei_new[protein_A_index1][2]!=0){
                                protein_B_index = res_nei_new[protein_A_index1][2];
                                j = res_nei_new[protein_A_index1][4]; // corresopnding protein_B residue index
                                
                                dist2=sqrt( (R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                           *(R_x_new[protein_B_index][j][2]-R_x_new[protein_A_index1][3][2])
                                           +(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                           *(R_y_new[protein_B_index][j][2]-R_y_new[protein_A_index1][3][2])
                                           );
                                dist1=sqrt( (R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                           *(R_x_new[protein_B_index][j][1]-R_x_new[protein_A_index1][3][1])
                                           +(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                           *(R_y_new[protein_B_index][j][1]-R_y_new[protein_A_index1][3][1])
                                           );
                                if (!AreSame(dist1, bond_dist_cutoff/2 + RB_A_radius + RB_B_radius) || !AreSame(dist2, bond_dist_cutoff/2)){
                                    moved[protein_A_index1] = 1;  //////////////////////// mark moved!
                                    
                                    for (k = 1; k <= RB_A_num_per_protein; k++){
                                        R_x_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][1] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][4] = (bond_dist_cutoff/2 + RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][3] = (bond_dist_cutoff/2 + 2*RB_A_radius)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                        
                                        R_x_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_x_new[protein_B_index][j][2] - R_x_new[protein_B_index][j][1]) + R_x_new[protein_B_index][j][2];
                                        R_y_new[protein_A_index1][k][2] = (bond_dist_cutoff/2)/RB_B_radius * (R_y_new[protein_B_index][j][2] - R_y_new[protein_B_index][j][1]) + R_y_new[protein_B_index][j][2];
                                    }
                                }
                            }
                        }
                    }
                }
                ////////////////// finished: put 1st protein_A(without cis interaction) back to protein_B
                ////////////////// next step: put 2nd protein_A (with cis interaction, but not attached to protein_B) back to its cis protein_A
                
                if (protein_B_num_inComplex > 1) {
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {  ///////move protein_A
                            
                            protein_A_index1 = results[complex_index][complex_sub_index];
                            
                            if (res_nei_new[protein_A_index1][2]!=0
                                && res_nei_new[protein_A_index1][3]!=0    // make sure this cis protein_A not attached with protein_B
                                && res_nei_new[res_nei_new[protein_A_index1][3]][2] == 0){
                                protein_A_index2 = res_nei_new[protein_A_index1][3];
                                
                                dist2=sqrt( (R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                           *(R_x_new[protein_A_index1][3][3]-R_x_new[protein_A_index2][3][3])
                                           +(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                           *(R_y_new[protein_A_index1][3][3]-R_y_new[protein_A_index2][3][3])
                                           );
                                dist1=sqrt( (R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                           *(R_x_new[protein_A_index1][3][1]-R_x_new[protein_A_index2][3][1])
                                           +(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                           *(R_y_new[protein_A_index1][3][1]-R_y_new[protein_A_index2][3][1])
                                           );
                                if (!AreSame(dist1, cis_dist_cutoff/2 + RB_A_radius + RB_A_radius) || !AreSame(dist2, cis_dist_cutoff/2)){
                                    for (k = 1; k <= RB_A_num_per_protein; k++){
                                        R_x_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                        R_y_new[protein_A_index2][k][1] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                        
                                        R_x_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                        R_y_new[protein_A_index2][k][4] = (cis_dist_cutoff/2 + RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                        
                                        R_x_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                        R_y_new[protein_A_index2][k][3] = (cis_dist_cutoff/2              )/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                        
                                        R_x_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_x_new[protein_A_index1][3][3] - R_x_new[protein_A_index1][3][1]) + R_x_new[protein_A_index1][3][3];
                                        R_y_new[protein_A_index2][k][2] = (cis_dist_cutoff/2 + 2*RB_A_radius)/RB_A_radius * (R_y_new[protein_A_index1][3][3] - R_y_new[protein_A_index1][3][1]) + R_y_new[protein_A_index1][3][3];
                                    }
                                } // finished move cis protein_A
                            }
                        }
                    }
                }
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////check collision
                collision_flag = 0;
                
                for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                    if (results[complex_index][complex_sub_index] == 0){break;}
                    
                    if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {
                        protein_A_index = results[complex_index][complex_sub_index];
                        
                        for (i = 1; i <= protein_A_tot_num; i++){      ///check overlap between RB_A and all the other RB_A
                            if (protein_A_index != i){
                                dist=sqrt( (R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])*(R_x_new[i][1][1]-R_x_new[protein_A_index][1][1])
                                          +(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])*(R_y_new[i][1][1]-R_y_new[protein_A_index][1][1])
                                          +(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])*(R_z_new[i][1][1]-R_z_new[protein_A_index][1][1])
                                          );
                                if (dist < RB_A_radius + RB_A_radius) {
                                    collision_flag = 1;
                                }
                            }
                        }
                        
                        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){        //////check overlap between RB_A with RB_B
                            for (j = 2; j <= RB_B_num_per_protein; j++){
                                for (k = 1; k <= RB_A_num_per_protein; k++){
                                    dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])*(R_x_new[i][j][1]-R_x_new[protein_A_index][k][1])
                                              +(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])*(R_y_new[i][j][1]-R_y_new[protein_A_index][k][1])
                                              +(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])*(R_z_new[i][j][1]-R_z_new[protein_A_index][k][1])
                                              );
                                    if (dist < RB_A_radius + RB_B_radius) {
                                        collision_flag = 1;
                                    }
                                }
                            }
                        }
                    }
                    
                    if ( results[complex_index][complex_sub_index] > protein_A_tot_num){
                        protein_B_index = results[complex_index][complex_sub_index];
                        
                        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){      ///check overlap between RB_B and all the other RB_B
                            for (j = 2; j <= RB_B_num_per_protein; j++){
                                for (k = 2; k <= RB_B_num_per_protein; k++){
                                    if (protein_B_index != i){
                                        dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_B_index][k][1])*(R_x_new[i][j][1]-R_x_new[protein_B_index][k][1])
                                                  +(R_y_new[i][j][1]-R_y_new[protein_B_index][k][1])*(R_y_new[i][j][1]-R_y_new[protein_B_index][k][1])
                                                  +(R_z_new[i][j][1]-R_z_new[protein_B_index][k][1])*(R_z_new[i][j][1]-R_z_new[protein_B_index][k][1])
                                                  );
                                        if (dist < RB_B_radius + RB_B_radius) {
                                            collision_flag = 1;
                                        }
                                    }
                                }
                            }
                        }
                        
                        for (i = 1; i <= protein_A_tot_num; i++){        //////check overlap between RB_A with RB_B
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                for (k = 2; k <= RB_B_num_per_protein; k++){        //////check overlap between RB_A with RB_B
                                    dist=sqrt( (R_x_new[i][j][1]-R_x_new[protein_B_index][k][1])*(R_x_new[i][j][1]-R_x_new[protein_B_index][k][1])
                                              +(R_y_new[i][j][1]-R_y_new[protein_B_index][k][1])*(R_y_new[i][j][1]-R_y_new[protein_B_index][k][1])
                                              +(R_z_new[i][j][1]-R_z_new[protein_B_index][k][1])*(R_z_new[i][j][1]-R_z_new[protein_B_index][k][1])
                                              );
                                    if (dist < RB_A_radius + RB_B_radius) {
                                        collision_flag = 1;
                                    }
                                }
                            }
                        }
                    }
                }  ////////end  check collision, and next step : put them back when flag = 1
                
                
                if (collision_flag == 1){
                    for ( complex_sub_index= 1; complex_sub_index <= protein_tot_num; complex_sub_index++){
                        if (results[complex_index][complex_sub_index] == 0){break;}
                        
                        
                        if ( results[complex_index][complex_sub_index] <= protein_A_tot_num) {
                            protein_A_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_A_num_per_protein; j++){
                                for (k = 1; k <= RB_A_res_num; k++){
                                    R_x_new[protein_A_index][j][k] = R_x[protein_A_index][j][k];
                                    R_y_new[protein_A_index][j][k] = R_y[protein_A_index][j][k];
                                    R_z_new[protein_A_index][j][k] = R_z[protein_A_index][j][k];
                                }
                            }
                        }
                        
                        if ( results[complex_index][complex_sub_index] > protein_A_tot_num){
                            protein_B_index = results[complex_index][complex_sub_index];
                            
                            for (j = 1; j <= RB_B_num_per_protein; j++){
                                for (k = 1; k <= RB_B_res_num; k++){
                                    R_x_new[protein_B_index][j][k] = R_x[protein_B_index][j][k];
                                    R_y_new[protein_B_index][j][k] = R_y[protein_B_index][j][k];
                                    R_z_new[protein_B_index][j][k] = R_z[protein_B_index][j][k];
                                }
                            }
                        }
                    }
                }/// put the collision molecules back when flag = 1
                
            }  // end traversal each protein_B, selecting_mole_index > protein_A_tot_num) && (selecting_mole_index <= (protein_A_tot_num + protein_B_tot_num)
            
            
            
            
            
            
            
            /////////////////////////////////////////////////////////////// End protein_B
            
        } // end iteration_mole_step, end diffusion! end part1 !
        
        /////////////// START Reaction, start part2
        
        // (1) RB_A and RB_B associate into a bond
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = protein_A_tot_num + 1; j <= protein_tot_num; j++){    /// j = protein_B_index
                for (k = 2; k <= RB_B_num_per_protein; k++){      // go over each binding site of RB_B 2,3,4
                    if (protein_status_new[i][2] == 0 ){           // for protein A binding site [i][3][2]
                        if (protein_status_new[j][k] == 0){     // for protein B binding site [j][3][2]
                            dist=sqrt( (R_x_new[j][k][2]-R_x_new[i][3][2])*(R_x_new[j][k][2]-R_x_new[i][3][2])
                                      +(R_y_new[j][k][2]-R_y_new[i][3][2])*(R_y_new[j][k][2]-R_y_new[i][3][2])
                                      +(R_z_new[j][k][2]-R_z_new[i][3][2])*(R_z_new[j][k][2]-R_z_new[i][3][2])
                                      );
                            if (dist < bond_dist_cutoff){
                                
                                theta_ot = 0;
                                point_x[0] = R_x_new[i][3][1] - R_x_new[i][3][2];
                                point_y[0] = R_y_new[i][3][1] - R_y_new[i][3][2];
                                point_z[0] = R_z_new[i][3][1] - R_z_new[i][3][2];
                                point_x[1] = 0;
                                point_y[1] = 0;
                                point_z[1] = 0;
                                point_x[2] = R_x_new[j][k][1] - R_x_new[j][k][2];
                                point_y[2] = R_y_new[j][k][1] - R_y_new[j][k][2];
                                point_z[2] = R_z_new[j][k][1] - R_z_new[j][k][2];
                                theta_ot2 =  gettheta ( point_x,  point_y,  point_z,  theta_ot);
                                //         cout << " theta_pd2 = " << theta_ot2 <<'\n';
                                
                                theta_pd = 0;
                                point_x[0] = R_x_new[i][3][1] - R_x_new[i][3][4];
                                point_y[0] = R_y_new[i][3][1] - R_y_new[i][3][4];
                                point_z[0] = R_z_new[i][3][1] - R_z_new[i][3][4];
                                point_x[1] = 0;
                                point_y[1] = 0;
                                point_z[1] = 0;
                                point_x[2] = R_x_new[j][1][1] - R_x_new[j][1][2];  //// adding RB_B [j][1][2] to protein B...
                                point_y[2] = R_y_new[j][1][1] - R_y_new[j][1][2];
                                point_z[2] = R_z_new[j][1][1] - R_z_new[j][1][2];
                                theta_pd2 =  gettheta ( point_x,  point_y,  point_z,  theta_pd);
                                //        cout << "theta_ot2 = " << theta_pd2 <<'\n';
                                
                                
                                if ((abs(theta_pd2) < bond_thetapd_cutoff) && (abs(theta_ot2 - 180) < bond_thetaot_cutoff)){
                                    ////////////////////mark mistake!!  angle - 180 !!
                                    
                                    Prob_Ass = Ass_Rate*time_step;
                                    prob = rand2();
                                    
                                    if ((prob < Prob_Ass)) {
                                        
                                        protein_status_new[i][2] = 1;       // protein_A binding site [i][3][2]  with protein_B
                                        protein_status_new[j][k] = 1;
                                        
                                        res_nei_new[j][k] = i;   // j = protein_B_index here, j = protein_A_tot_num + 1; j <= protein_tot_num;
                                        res_nei_new[i][2] = j;   // connected protein_B index is j;
                                        res_nei_new[i][4] = k;    // Becareful, matrix[i][3] is taken by cis interaction!
                                        // so we use [i][4]  connected protein_B_res index is k; here i <= protein_A_tot_num;
                                        
                                        bond_num_new = bond_num_new + 1;
                                        bond_num_rl_new = bond_num_rl_new + 1;
                                        
                                        // corner case: if protein A1 and protein A2 are monomer, after R-L association  ==> monomer cis (-1) became into complex cis(+1)
                                        selected_pro_A2 = res_nei_new[i][3];
                                        if ( selected_pro_A2!= 0 && protein_status_new[selected_pro_A2][2] == 0){
                                            bond_num_mono_cis_new = bond_num_mono_cis_new - 1;
                                            bond_num_cis_new = bond_num_cis_new + 1;
                                        }
                                        
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }/// end RB_A and RB_B associate into a bond
        
        /////  type (1) single R -- single R  cis interaction occurs !!
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = 1; j <= protein_A_tot_num; j++){
                if ( i != j
                    && protein_status_new[i][3] == 0
                    && protein_status_new[j][3] == 0
                    && protein_status_new[i][2] == 0
                    && protein_status_new[j][2] == 0
                    ){           // for protein A binding site [i][3][2]  mistake! should be [i][3][3]
                    dist=sqrt( (R_x_new[j][3][3]-R_x_new[i][3][3])*(R_x_new[j][3][3]-R_x_new[i][3][3])
                              +(R_y_new[j][3][3]-R_y_new[i][3][3])*(R_y_new[j][3][3]-R_y_new[i][3][3])
                              +(R_z_new[j][3][3]-R_z_new[i][3][3])*(R_z_new[j][3][3]-R_z_new[i][3][3])
                              );
                    if (dist < cis_dist_cutoff){
                        
                        theta_ot = 0;
                        point_x[0] = R_x_new[i][3][1] - R_x_new[i][3][3];
                        point_y[0] = R_y_new[i][3][1] - R_y_new[i][3][3];
                        point_z[0] = R_z_new[i][3][1] - R_z_new[i][3][3];
                        point_x[1] = 0;
                        point_y[1] = 0;
                        point_z[1] = 0;
                        point_x[2] = R_x_new[j][3][1] - R_x_new[j][3][3];
                        point_y[2] = R_y_new[j][3][1] - R_y_new[j][3][3];
                        point_z[2] = R_z_new[j][3][1] - R_z_new[j][3][3];
                        theta_ot2 =  gettheta ( point_x,  point_y,  point_z,  theta_ot);
                        //     cout << "cis theta_ot2 = " << theta_ot2<<'\n';
                        
                        
                        
                        if (abs(theta_ot2 - 180) < cis_thetaot_cutoff){
                            ///// mark mistake!!  angle - 180!!
                            
                            Prob_Ass = mono_cis_Ass_Rate*time_step;
                            prob = rand2();
                            
                            if ((prob < Prob_Ass)) {
                                
                                protein_status_new[i][3] = 1;       // protein_A binding site [i][3][2]  with protein_B
                                protein_status_new[j][3] = 1;
                                bond_num_new = bond_num_new + 1;
                                bond_num_mono_cis_new = bond_num_mono_cis_new + 1;
                                
                                res_nei_new[j][3] = i;
                                res_nei_new[i][3] = j;
                                
                                //            cout << "cis i = " << bond_nb_ctg_pro_idx_new[bond_num_new][3]<< "  cis j = " << bond_nb_ctg_pro_idx_new[bond_num_new][4]<<'\n';
                            }
                        }
                    }
                }
            }
        }
        
        
        /////  type (2) NOT (single R -- single R) cis interaction occurs !!
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = 1; j <= protein_A_tot_num; j++){
                if ( i != j
                    && protein_status_new[i][3] == 0
                    && protein_status_new[j][3] == 0
                    && (protein_status_new[j][2] == 1 || protein_status_new[i][2] == 1)
                    ){           // for protein A binding site [i][3][2]
                    dist=sqrt( (R_x_new[j][3][3]-R_x_new[i][3][3])*(R_x_new[j][3][3]-R_x_new[i][3][3])
                              +(R_y_new[j][3][3]-R_y_new[i][3][3])*(R_y_new[j][3][3]-R_y_new[i][3][3])
                              +(R_z_new[j][3][3]-R_z_new[i][3][3])*(R_z_new[j][3][3]-R_z_new[i][3][3])
                              );
                    if (dist < cis_dist_cutoff){
                        
                        theta_ot = 0;
                        point_x[0] = R_x_new[i][3][1] - R_x_new[i][3][3];
                        point_y[0] = R_y_new[i][3][1] - R_y_new[i][3][3];
                        point_z[0] = R_z_new[i][3][1] - R_z_new[i][3][3];
                        point_x[1] = 0;
                        point_y[1] = 0;
                        point_z[1] = 0;
                        point_x[2] = R_x_new[j][3][1] - R_x_new[j][3][3];
                        point_y[2] = R_y_new[j][3][1] - R_y_new[j][3][3];
                        point_z[2] = R_z_new[j][3][1] - R_z_new[j][3][3];
                        theta_ot2 =  gettheta ( point_x,  point_y,  point_z,  theta_ot);
                        //     cout << "cis theta_ot2 = " << theta_ot2<<'\n';
                        
                        
                        
                        if (abs(theta_ot2 - 180) < cis_thetaot_cutoff){
                            ///// mark mistake!!  angle - 180!!
                            
                            Prob_Ass = cis_Ass_Rate*time_step;
                            prob = rand2();
                            
                            if ((prob < Prob_Ass)) {
                                
                                protein_status_new[i][3] = 1;       // protein_A binding site [i][3][2]  with protein_B
                                protein_status_new[j][3] = 1;
                                bond_num_new = bond_num_new + 1;
                                bond_num_cis_new = bond_num_cis_new + 1;
                                
                                
                                res_nei_new[j][3] = i;
                                res_nei_new[i][3] = j;
                                
                                //            cout << "cis i = " << bond_nb_ctg_pro_idx_new[bond_num_new][3]<< "  cis j = " << bond_nb_ctg_pro_idx_new[bond_num_new][4]<<'\n';
                            }
                        }
                    }
                }
            }
        }
        
        
        
        // (2) bond dissociate into a RB_A and RB_B
        for (i = 1; i <= protein_A_tot_num; i++){
            if (protein_status_new[i][2] == 1){
                selected_pro_A = i;
                selected_pro_B = res_nei_new[i][2];
                selected_res_B = res_nei_new[i][4];
                
                Prob_Diss = Diss_Rate * time_step;
                prob = rand2();
                if (prob < Prob_Diss){
                    protein_status_new[selected_pro_A][2] = 0;
                    protein_status_new[selected_pro_B][selected_res_B] = 0;
                    
                    res_nei_new[selected_pro_A][2] = 0;
                    res_nei_new[selected_pro_A][4] = 0;
                    
                    res_nei_new[selected_pro_B][selected_res_B] = 0;
                    
                    bond_num_new = bond_num_new - 1;
                    bond_num_rl_new = bond_num_rl_new - 1;
                    
                    // corner case:  after R-L dissociation  ==> complex cis(-1) may become monomer cis
                    selected_pro_A2 = res_nei_new[i][3];
                    if (selected_pro_A2!= 0 && protein_status_new[selected_pro_A2][2] == 0){
                        bond_num_mono_cis_new = bond_num_mono_cis_new + 1;
                        bond_num_cis_new = bond_num_cis_new - 1;
                    }
                    
                }
            }
        }
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////
        //// type (1)cis dissociation  momoner cis
        for (i = 1; i <= protein_A_tot_num; i++){
            if (protein_status_new[i][3] == 1){
                selected_pro_A = i;
                selected_pro_A2 = res_nei_new[i][3];
                
                if (protein_status_new[selected_pro_A][2] == 0 && protein_status_new[selected_pro_A2][2] == 0){
                    
                    Prob_Diss = mono_cis_Diss_Rate * time_step;
                    prob = rand2();
                    if (prob < Prob_Diss){
                        protein_status_new[selected_pro_A][3] = 0;
                        protein_status_new[selected_pro_A2][3] = 0;
                        res_nei_new[selected_pro_A][3] = 0;
                        res_nei_new[selected_pro_A2][3] = 0;
                        
                        bond_num_new = bond_num_new - 1;
                        bond_num_mono_cis_new = bond_num_mono_cis_new - 1;
                    }
                }
            }
        }
        
        //// type (2)cis dissociation    complex cis
        for (i = 1; i <= protein_A_tot_num; i++){
            if (protein_status_new[i][3] == 1){
                selected_pro_A = i;
                selected_pro_A2 = res_nei_new[i][3];
                
                if (protein_status_new[selected_pro_A][2] == 1 || protein_status_new[selected_pro_A2][2] == 1){
                    
                    Prob_Diss = cis_Diss_Rate * time_step;
                    prob = rand2();
                    if (prob < Prob_Diss){
                        protein_status_new[selected_pro_A][3] = 0;
                        protein_status_new[selected_pro_A2][3] = 0;
                        res_nei_new[selected_pro_A][3] = 0;
                        res_nei_new[selected_pro_A2][3] = 0;
                        
                        bond_num_new = bond_num_new - 1;
                        bond_num_cis_new = bond_num_cis_new - 1;
                        
                    }
                }
            }
        }
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////
        
/*
        // how many complexes in the system?
        set<int> complex_set;   // use c++ data structure set
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            for (j = 1; j <= RB_B_num_per_protein; j++){
                if (protein_status_new[i][j] == 1){
                    complex_set.insert(i);              // "set" data sturcture in c++ only storage unique elements!
                }
            }
        }
        
        //cout<<"Set Size = "<<complex_set.size()<<endl;
*/
        
        
        ///////update coordinates for all molecules
        
        for (i = 1; i <= protein_A_tot_num; i++){
            for (j = 1; j <= RB_A_num_per_protein; j++){
                for (k = 1; k <= RB_A_res_num; k++){
                    R_x[i][j][k] = R_x_new[i][j][k];
                    R_y[i][j][k] = R_y_new[i][j][k];
                    R_z[i][j][k] = R_z_new[i][j][k];
                }
            }
            protein_status[i][2] = protein_status_new[i][2];
            protein_status[i][3] = protein_status_new[i][3];
            res_nei[i][2] = res_nei_new[i][2];
            res_nei[i][4] = res_nei_new[i][4];
            res_nei[i][3] = res_nei_new[i][3];
            
        }
        
        
        for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
            for (j = 1; j <= RB_B_num_per_protein; j++){
                for (k = 1; k <= RB_B_res_num; k++){
                    R_x[i][j][k] = R_x_new[i][j][k];
                    R_y[i][j][k] = R_y_new[i][j][k];
                    R_z[i][j][k] = R_z_new[i][j][k];
                }
                protein_status[i][j] = protein_status_new[i][j];
                res_nei[i][j] = res_nei_new[i][j];
            }
        }
        
        
        
        bond_num = bond_num_new;
        bond_num_rl = bond_num_rl_new;
        bond_num_cis = bond_num_cis_new;
        bond_num_mono_cis = bond_num_mono_cis_new;
        
        if (tot_cluster_num != 0){
            cluster_size = (double)tot_proteins_in_cluster/tot_cluster_num;
        }
        
        
///////////Check Point file:
        if (mc_time_step % 5000 == 0){
            std::ofstream check_point ("position.cpt", std::ofstream::trunc);
            check_point.setf(ios::fixed, ios::floatfield);  // make output numbers look neat
            check_point << setprecision(3);
            
            for (i = 1; i <= protein_A_tot_num; i++){
                for (j = 1; j <= RB_A_num_per_protein; j++){
                    for (k = 1; k <= RB_A_res_num; k++){
                        check_point <<setw(10)<< R_x[i][j][k] <<setw(10)<< R_y[i][j][k] <<setw(10)<< R_z[i][j][k] <<'\n';
                    }
                }
                check_point <<setw(8)<< protein_status_new[i][2];
                check_point <<setw(8)<< protein_status_new[i][3];
                check_point <<setw(8)<< res_nei_new[i][2];
                check_point <<setw(8)<< res_nei_new[i][4];
                check_point <<setw(8)<< res_nei_new[i][3]<<'\n';
                
            }
            
            
            for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
                for (j = 1; j <= RB_B_num_per_protein; j++){
                    for (k = 1; k <= RB_B_res_num; k++){
                        check_point <<setw(10)<< R_x[i][j][k] <<setw(10)<< R_y[i][j][k] <<setw(10)<< R_z[i][j][k] <<'\n';
                    }
                    check_point <<setw(8)<< protein_status_new[i][j];
                    check_point <<setw(8)<< res_nei_new[i][j]<<'\n';
                }
            }
            
            check_point << bond_num_new<<'\n';
            check_point << bond_num_rl_new<<'\n';
            check_point << bond_num_cis_new<<'\n';
            check_point << bond_num_mono_cis_new<<'\n';
            check_point << protein_num_in_Max_Complex<<'\n';
            check_point << mc_time_step<<'\n';

            check_point.close();
        }
        
        
        if (mc_time_step % 5000 == 0){
            std::ofstream bond ("bond.dat", std::ofstream::app);
            bond.setf(ios::fixed, ios::floatfield);  // make output numbers look neat
            bond << setprecision(3);
            bond <<setw(15)<< mc_time_step*time_step<<setw(5)<< bond_num_rl<<setw(5)<< bond_num_mono_cis<<setw(10)<< bond_num_cis<<setw(10)<<bond_num<<setw(10)<< cluster_size <<setw(10)<< protein_num_in_Max_Complex<<'\n';
            bond.close();
        }
        
        
        
        
        if (mc_time_step % 5000 == 0){
            std::ofstream ofs ("test.gro", std::ofstream::app); //open file and append to write
            ofs.setf(ios::fixed, ios::floatfield);  // make output numbers look neat
            ofs << setprecision(3);
            ofs << "Hello Gro!" << ", t="<< mc_time_step*time_step <<'\n';
            ofs << protein_A_tot_num*RB_A_num_per_protein + protein_B_tot_num*(RB_B_num_per_protein - 1) <<'\n';
            //  ofs << protein_A_tot_num*RB_A_num_per_protein*(RB_A_res_num) + protein_B_tot_num*(RB_B_num_per_protein)*RB_B_res_num <<'\n';
            for (i=1; i<= protein_A_tot_num; i++){
                for (j = 1; j <= RB_A_num_per_protein; j++){
                    ofs <<setw(5)<<i<<"ALA" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][1]/10 <<setw(8)<< R_y[i][j][1]/10 <<setw(8)<< R_z[i][j][1]/10 <<'\n';
                    //   ofs <<setw(5)<<i<<"VAL" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][2]/10 <<setw(8)<< R_y[i][j][2]/10 <<setw(8)<< R_z[i][j][2]/10 <<'\n';
                    //   ofs <<setw(5)<<i<<"VAL" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][3]/10 <<setw(8)<< R_y[i][j][3]/10 <<setw(8)<< R_z[i][j][3]/10 <<'\n';
                    //   ofs <<setw(5)<<i<<"VAL" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][4]/10 <<setw(8)<< R_y[i][j][4]/10 <<setw(8)<< R_z[i][j][4]/10 <<'\n';
                    //          cout << " Ax = " << R_x[i][j][4] << " Ay = " << R_y[i][j][4] <<" Az = " << R_z[i][j][4] <<endl;
                    
                }
            }
            
            for (i = protein_A_tot_num+1; i <= protein_tot_num; i++){
                for (j=2; j<= RB_B_num_per_protein; j++){   //j = 2,3,4
                    ofs <<setw(5)<<i<<"LEU" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][1]/10 <<setw(8)<< R_y[i][j][1]/10 <<setw(8)<< R_z[i][j][1]/10 <<'\n';
                    //   ofs <<setw(5)<<i<<"GLY" <<setw(7)<< "CA" << setw(5)<< i <<setw(8)<< R_x[i][j][2]/10 <<setw(8)<< R_y[i][j][2]/10 <<setw(8)<< R_z[i][j][2]/10 <<'\n';
                    //         cout << " Bx = " << R_x[i][j][1] << " By = " << R_y[i][j][1] <<" Bz = " << R_z[i][j][1] <<endl;
                }
            }
            
            ofs << setw(8)<< cell_range_x/10 << setw(12)<< cell_range_y/10 << setw(12)<< cell_range_z/10<<'\n';
            ofs.close();
            
        } // end writing coordinates every 10 time_step
        
        
        
        if (mc_time_step % 5000 == 0){
            std::ofstream cluster ("cluster.log", std::ofstream::app); //open file and append to write
            cluster << "Hello Cluster!" << ", t="<< mc_time_step*time_step <<'\n';
            for (i = protein_A_tot_num + 1; i <= protein_tot_num; i++){
                for (j = 1; j <= protein_tot_num; j++){
                    if (results[i][j] != 0){
                        cluster  << results[i][j]<< "  ";
                    }
                }
                cluster <<'\n';
            }
            
            cluster.close();
            
        } // end writing coordinates every 10 time_step
        
        
    } //end mc_time_step, end diffusion-reaction loop
} // end main loop



double rand2(){
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    // ready to generate random numbers
    
    double currentRandomNumber = unif(rng);
    return currentRandomNumber;
    //    std::cout << currentRandomNumber << std::endl;
}


double gettheta (double point_x[3], double point_y[3], double point_z[3], double test_theta){
    double lx[2],ly[2],lz[2],lr[2];
    double conv, doth1, doth2;
    int in;
    
    for ( in=0; in < 2; in++){
        lx[in] = 0;
        ly[in] = 0;
        lz[in] = 0;
        lr[in] = 0;
    }
    
    lx[0]= point_x[1] - point_x[0];
    ly[0]= point_y[1] - point_y[0];
    lz[0]= point_z[1] - point_z[0];
    lr[0]= sqrt(lx[0]*lx[0]+ly[0]*ly[0]+lz[0]*lz[0]);
    
    lx[1]= point_x[2] - point_x[1];
    ly[1]= point_y[2] - point_y[1];
    lz[1]= point_z[2] - point_z[1];
    lr[1]= sqrt(lx[1]*lx[1]+ly[1]*ly[1]+lz[1]*lz[1]);
    
    test_theta = 0;
    
    conv = 180/3.14159;
    
    doth1 = -(lx[1]*lx[0] + ly[0]*ly[1] + lz[0]*lz[1]);
    doth2 = doth1/(lr[1]*lr[0]);
    if (doth2 > 1){
        doth2 = 1;
    }
    if (doth2 < -1){
        doth2 = -1;
    }
    
    test_theta = acos(doth2)* conv;
    return test_theta;
}

bool AreSame(double a, double b)
{
    return fabs(a - b) < EPSILON;
}


///////// put z coordinates into right value

//angle between  RB_B [j][1][1] -------> RB_A[i][3][1] and [0][0] ---> [0][1] in x-y plane
//   x2 = R_x_new[i][2][1] - R_x_new[j][1][1];
//   y2 = R_y_new[i][2][1] - R_y_new[j][1][1];
//     x1 = 0;
//     y1 = 1;
//   dot = x1*x2 + y1*y2;    // dot product between [x1, y1] and [x2, y2]
//   det = x1*y2 - y1*x2;    // determinant
//   angle = atan2(det, dot)   // this is wrong! because the range is -180 to 180, but we need 0 to 360
// atan2 usually is in the range [-180°,180°]. To get [0°,360°] without a case distinction, one can replace atan2(y,x) with atan2(-y,-x) + 180°
// from point 1 counterclokewise to point 2 !!

//angle = atan2((R_x_new[protein_B_index][2][1] - R_x_new[protein_B_index][1][1]), (R_y_new[protein_B_index][2][1] - R_y_new[protein_B_index][1][1])) + pai;

// for 2nd case, angle between







