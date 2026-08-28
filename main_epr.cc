/*------------------------------------------------------------------------------------------
 
 Simulation of an active Brownian particle following a given navigation policy.
 The goal is to compute the Entropy Production Rate (EPR) of the particle along its trajectory.
 
 Compile the code with the associated Makefile.
 
 Date created: 24-06-2023
 Last update: 24-06-2023
 
 Author: Lorenzo Piro
 
 -------------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cfloat>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string.h>
#include <random>

using namespace std;

enum{DIR_X,DIR_Y,DIMENSION};
enum{DIR_XX,DIR_XY,DIR_YX,DIR_YY,RANK2_DIMENSION};
enum{P_A,P_kappa,P_Pe_r,P_Pe_m,P_mu,P_dt,P_flow,P_policy,P_epsilon,P_theta0,P_Rcutoff,NParameters}; //Parameters
enum{Const}; //Flows
enum{SP,AAP}; //Policies

#define vect vector<double>

//Particle structure
struct particle{
  vect pos;
  vect vel;
  double theta;
  double thetaMeas;
};

/* ------------------------- Global Variables ------------------------------ */

//Particle variables
double kappa, Pe_r, Pe_m, mu, sigma_m, theta0, R_cutoff; //Torque strength, Peclet numbers, Measurement strength, Initial orientation, Cut-off radius
vect init_pos; //Initial position
particle P; //Particle

//Simulation variables
int N, Nsave_trajs, t_save; //# of simulated trajectories, # of trajectories to save and at which freq
double tmax, t_steady, dt; //Max run time, transient time and time step

//Flow variables
double A; //Flow field amplitude and flow wavelength
string name_flow; //Name of the flow field

//Policy variables
double eps; //Sharpness of the aligning potential
string policy; // Name of the policy used

//Parameters array
double Parameters[NParameters];

//Task of the particle
bool target_point;

//Vectors containing quantities to be saved
vect stoch_traj, epr, moments; 

//File names and saving booleans
string name_params, name_file_traj, name_save;
bool save_trajs,save_epr,save_moments; //=1 if trajectories/EPR/moments are saved
int run_index; //Index such that parallel runs don't overwrite the data

//Profiling
clock_t start_time, end_time;

#include "../include_local/utility.h"
#include "../include_local/flow.h"
#include "../include_local/dynamics.h"
#include "../include_local/saving.h"

/* ---------------------------------- Main --------------------------------- */

int main(int arg, char *args[])
{
    //Local variables
    int seed;
    double x_n, y_n, deltaTheta1_n, deltaTheta2_n, deltaTheta1_np1, deltaTheta2_np1, theta_star_n, theta_star_np1;
    double noise_theta_n, noise_thetaMeas_n;
    vect Flow_real_n, Flow_np1;
     
    //Start time for profiling
    start_time = clock();
    
    //Set the variables to their default values and import from the command line
    Set_default_values();
    Command_line(arg,args);
    
    //Fix the parameters of the simulation and print configuration                                              
    Set_parameters(arg,args);
    Print_configuration();

    //Seed rng
    seed = 1716217304 + run_index*17; //time(NULL); 
    srand48(seed);
    
    //Set output folders
    init_output_folders();
    
    //Loop over different trajectories (or episodes if you prefer)
    for(int k=0;k<N;k++){
    
      //Initialize particle position and orientation
      Initialize_particle(&P,init_pos,theta0,Parameters);
      
      //Initialize values to zero before starting time
      long long int t = 0;
      double entr_deter = 0., entr_stoch_ito = 0., entr_stoch_strat = 0., entr_act = 0., entr_sens = 0.;
      double msd_n_cum = 0., div_u_star_cum = 0., grad_thetaStar2_cum = 0.;
      double cosdTheta_n_cum = 0., cosdThetaStar_n_cum = 0., sindThetaStar2_n_cum = 0., sindTheta2_n_cum = 0.;
      double sindThetam_n_cum = 0., sindThetam2_n_cum = 0., cosdThetam_n_cum = 0., cosdThetam2_n_cum = 0.;
      bool active;
      double sensing;
      
      //Loop over time
      while(t*dt <= tmax){	  
	
    	//Compute current value of theta_star from policy
    	Flow_real_n = Compute_flow_real(&P,Parameters);

        theta_star_n = Compute_theta_star(&P,Flow_real_n,Parameters);
	
    	deltaTheta1_n = P.theta - P.thetaMeas; 
        deltaTheta2_n = P.thetaMeas - theta_star_n;
    	noise_theta_n = rand_normal(0.,1.);
    	noise_thetaMeas_n = rand_normal(0.,1.);
    	x_n = P.pos[DIR_X];
    	y_n = P.pos[DIR_Y];
	
    	//Introduce sensing cutoff close to the point target
    	if(target_point)
    	  active = (x_n*x_n+y_n*y_n > R_cutoff*R_cutoff);
    	else
    	  active = 1;
    	
    	sensing = active*sin(deltaTheta2_n);
    	
    	//Cumulate stuff we want to save
    	if(t*dt >= t_steady){

    	  //sin^2(theta-theta_m) and cos(theta-theta_m) to check general expression of EPR 
    	  sindTheta2_n_cum += dt*sin(deltaTheta1_n)*sin(deltaTheta1_n);
    	  cosdTheta_n_cum += dt*cos(deltaTheta1_n);
    	  //sin^2(theta-theta^*) and cos(theta-theta*) to check general expression of EPR
    	  sindThetaStar2_n_cum += dt*sin(deltaTheta1_n+deltaTheta2_n)*sin(deltaTheta1_n+deltaTheta2_n);
    	  cosdThetaStar_n_cum += dt*cos(deltaTheta1_n+deltaTheta2_n);
    	  //sin(theta_m-theta^*) and cos(theta_m-theta^*) to check moment expansion approximation
    	  sindThetam_n_cum += active*dt*sin(deltaTheta2_n);
    	  sindThetam2_n_cum += active*dt*sin(deltaTheta2_n)*sin(deltaTheta2_n);
    	  cosdThetam_n_cum += active*dt*cos(deltaTheta2_n);
    	  cosdThetam2_n_cum += active*dt*cos(deltaTheta2_n)*cos(deltaTheta2_n);

    	  //Mean squared distance from target or straight line depending on the task
    	  if(target_point){
    	    msd_n_cum += dt*(y_n*y_n+x_n*x_n);
    	    if(active){
    	      div_u_star_cum += dt/sqrt(y_n*y_n+x_n*x_n);
    	      grad_thetaStar2_cum += dt/(y_n*y_n+x_n*x_n);
    	    }
    	  }
          else{
    	    msd_n_cum += dt*y_n*y_n;
    	    div_u_star_cum += Compute_policy_divergence(&P,Flow_real_n,theta_star_n,Parameters);
    	    grad_thetaStar2_cum -= Compute_policy_divergence(&P,Flow_real_n,theta_star_n,Parameters)/Parameters[P_epsilon];
    	  }
    	}
	
    	//Update dynamics
    	P.vel[DIR_X] = cos(P.theta);
    	P.vel[DIR_Y] = sin(P.theta);
    
    	P.pos[DIR_X] += dt*P.vel[DIR_X] + Flow_real_n[DIR_X] + Parameters[P_Pe_r]*rand_normal(0.,1.);
    	P.pos[DIR_Y] += dt*P.vel[DIR_Y] + Flow_real_n[DIR_Y] + Parameters[P_Pe_r]*rand_normal(0.,1.); 
    	P.theta += -Parameters[P_kappa]*sin(deltaTheta1_n) + sqrt(2.*dt)*noise_theta_n;
    	P.thetaMeas += -Parameters[P_mu]*sensing + Parameters[P_Pe_m]*noise_thetaMeas_n;
    	
    	//Compute new value of theta_star from policy
    	Flow_np1 = Compute_flow_real(&P,Parameters);
    	theta_star_np1 = Compute_theta_star(&P,Flow_np1,Parameters);
    	deltaTheta1_np1 = P.theta - P.thetaMeas;
    	deltaTheta2_np1 = P.thetaMeas - theta_star_np1;

    	//Compute entropy production rate
    	if(t*dt >= t_steady){
    	  entr_deter += dt*(pow(kappa,2)*pow(sin(deltaTheta1_n),2) + Pe_m*pow(mu,2)*active*pow(sin(deltaTheta2_n),2)); //Should be OK
    	  
    	  entr_stoch_ito += - sqrt(2.*pow(kappa,2))*(sqrt(dt)*sin(deltaTheta1_n)*noise_theta_n + 0.5*dt*sqrt(2.)*cos(deltaTheta1_n)) \
    	    - sqrt(2.*Pe_m*pow(mu,2))*(sqrt(dt)*active*sin(deltaTheta2_n)*noise_thetaMeas_n + 0.5*dt*sqrt(2./Pe_m)*active*cos(deltaTheta2_n));  //Should be OK
	  
    	  entr_act += dt*(pow(kappa,2)*pow(sin(deltaTheta1_n),2)) - sqrt(2.*pow(kappa,2))*(sqrt(dt)*sin(deltaTheta1_n)*noise_theta_n + 0.5*dt*sqrt(2.)*cos(deltaTheta1_n));
    	  entr_sens += dt*Pe_m*pow(mu,2)*active*pow(sin(deltaTheta2_n),2)-sqrt(2.*Pe_m*pow(mu,2))*(sqrt(dt)*active*sin(deltaTheta2_n)*noise_thetaMeas_n \
		       + active*0.5*dt*sqrt(2./Pe_m)*cos(deltaTheta2_n));
    	}
    	
    	t++;
    	
    	//Store traj if it needs to be saved
    	if(save_trajs && k < Nsave_trajs && t%t_save==0 && t*dt>=t_steady){
    	  stoch_traj.push_back(entr_sens);
    	  stoch_traj.push_back(entr_deter + entr_stoch_ito);
    	  stoch_traj.push_back(grad_thetaStar2_cum);
    	  stoch_traj.push_back(msd_n_cum);
    	  stoch_traj.push_back(div_u_star_cum);
    	  stoch_traj.push_back(P.pos[DIR_X]);
    	  stoch_traj.push_back(P.pos[DIR_Y]);
    	  stoch_traj.push_back(P.theta);
    	  stoch_traj.push_back(P.thetaMeas);
    	  stoch_traj.push_back(theta_star_n);
    	}
	
      }
      //End loop over time 

      cout << "EPR det: " << entr_deter/(tmax-t_steady) << endl;
      cout << "EPR ito: " << entr_stoch_ito/(tmax-t_steady) << endl;
      cout << "MSD    : " << msd_n_cum/(tmax-t_steady) << endl;
      cout << "dUstar : " << div_u_star_cum/(tmax-t_steady) << endl;
	
      //Store in a vector EPR and other useful quantities to test theory prediction
      if(save_epr){
    	epr.push_back(entr_deter/(tmax-t_steady));
    	epr.push_back(entr_stoch_ito/(tmax-t_steady));
    	epr.push_back(entr_stoch_strat/(tmax-t_steady)); 
    	epr.push_back(msd_n_cum/(tmax-t_steady));
    	epr.push_back(div_u_star_cum/(tmax-t_steady));
    	epr.push_back(cosdTheta_n_cum/(tmax-t_steady));
    	epr.push_back(sindTheta2_n_cum/(tmax-t_steady));
    	epr.push_back(cosdThetaStar_n_cum/(tmax-t_steady));
    	epr.push_back(sindThetaStar2_n_cum/(tmax-t_steady));
    	epr.push_back(grad_thetaStar2_cum/(tmax-t_steady));
    	epr.push_back(entr_act/(tmax-t_steady));
    	epr.push_back(entr_sens/(tmax-t_steady));
      }

      //Store in a vector sine and cosine moments to check accuracy of moments expansion
      if(save_moments){
        moments.push_back(cosdThetam_n_cum/(tmax-t_steady));
        moments.push_back(cosdThetam2_n_cum/(tmax-t_steady));
        moments.push_back(sindThetam_n_cum/(tmax-t_steady));
        moments.push_back(sindThetam2_n_cum/(tmax-t_steady));
      }
      
      //Save traj on file
      if(save_trajs && k < Nsave_trajs){	
    	name_save = "../output/" + name_flow + "/Trajectories/" + name_params + "_tsave" + double_to_string(t_save*dt,4) + "_n" + int_to_string(run_index) + ".bin";
    	Save_binary_vect(name_save,stoch_traj);
    	stoch_traj.clear();
      }
    }
    
    //Save EPR and related stuff on file
    if(save_epr){
      name_save = "../output/"+ name_flow +"/EPR/" + name_params + "_n" + int_to_string(run_index) + ".dat";
      Save_vector(name_save,epr,12);
      epr.clear();
    }

    //Save moments of sine and cosine on file to check moments expansion
    if(save_moments){
      name_save = "../output/"+ name_flow +"/Moments/" + name_params + "_n" + int_to_string(run_index) + ".dat";
      Save_vector(name_save,moments,4);
      moments.clear();
    }
  
    //Print total CPU time of the simulation
    end_time = clock();
    print_time(start_time,end_time);
    
    return 0;
}















