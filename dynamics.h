//Initial condition
void Initialize_particle(particle *P, const vect init_pos, const double theta0, const double *Parameters)
{
    P->pos = init_pos;
    
    vect e0;
    e0.push_back(cos(theta0));
    e0.push_back(sin(theta0));
    P->vel = e0;
    P->theta = theta0;
    P->thetaMeas = theta0;
}

//Update theta^* according to the chosen policy
double Compute_theta_star(particle *P, const vect Flow, const double *Parameters)
{
  double theta_star;
  
  if(target_point)
      theta_star = M_PI + atan2(P->pos[DIR_Y],P->pos[DIR_X]+DBL_EPSILON);
  else{
      if(Parameters[P_policy] == SP)
        theta_star = 0.;
      else{
        double Dx = P->pos[DIR_Y];
        double epsilon = Parameters[P_epsilon];
        double alpha,alpha0;

    	alpha0 = - tanh(Dx/epsilon);
    
        if(alpha0 < -1.)
          alpha = -1.;
        else if(alpha0 > 1.)
          alpha = 1.;
        else
          alpha = alpha0;
        
        theta_star = asin(alpha);
      }
  }
  return theta_star;
}

//Compute the divergence of theta^* depending on the chosen policy
double Compute_policy_divergence(particle *P, const vect Flow, const double theta_star, const double *Parameters){
  
  double der_u_star[DIMENSION]={0.,0.};
  double div_u_star;
  double check, Dx = P->pos[DIR_Y], eps = Parameters[P_epsilon];

  if(Parameters[P_policy]==AAP){ //AAP
      der_u_star[DIR_X] = 0.;
      der_u_star[DIR_Y] =  - Parameters[P_dt]*sech(Dx/eps)*sech(Dx/eps)/eps;
  }
  else if(Parameters[P_policy]==SP){ //SP
    der_u_star[DIR_X] = 0.;
    der_u_star[DIR_Y] = 0.; 
  }
  
  div_u_star = der_u_star[DIR_X] + der_u_star[DIR_Y];

  return div_u_star;
} 


