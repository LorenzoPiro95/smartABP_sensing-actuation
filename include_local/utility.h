//Some potentially useful routines
static inline string int_to_string(const int a)
{
  ostringstream str;
  str << a;
  return str.str();
}

static inline string double_to_string(const double a, const int precision)
{
  ostringstream str;
  str << fixed << setprecision(precision) << a;
  return str.str();
}

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
	  x = 2.0*drand48() - 1;
	  y = 2.0*drand48() - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double sech(double x) {
    return 1.0/cosh(x);
}

/* ------------------------------------------------------------------------ */

//Print config of the simulation on standard output
static inline void Print_configuration()
{
  cout << "#-----------------------------------------------------------#" << endl;
  cout << "                      Simulation parameters                  " << endl;
  cout << "#-----------------------------------------------------------#" << endl;
  cout << "Flow intensity:           A = " << A << endl;
  cout << "Flow name:                    " << name_flow << endl;
  cout << "Peclet number (pos):   Pe_r = " << Pe_r << endl;
  cout << "Peclet number (meas):  Pe_m = " << Pe_m << endl;
  cout << "Torque strength:          k = " << kappa << endl;
  cout << "Measurement strength:    mu = " << mu << endl;
  if(target_point){
    cout << "Target point" << endl;
    cout << "Cut-off radius:           R = " << R_cutoff << endl;
  }
  else
    cout << "Target path" << endl;
  cout << "# of simulated trajs:     N = " << N << endl;
  cout << "Policy:                       " << policy << endl;
  if(Parameters[P_policy]==AAP)
      cout << "Sensitivity:            eps = " << eps << endl;
  cout << "Initial point:      (xi,yi) = (" << init_pos[DIR_X] << "," << init_pos[DIR_Y] << ")" << endl;
  cout << "Initial orientation:   chi0 = " << theta0 << endl;
  cout << "Saving EPR?                   " << save_epr << endl;
  cout << "Saving moments?               " << save_moments << endl;
  cout << "Saving trajectories?          " << save_trajs << endl;
}

//Read some variables from the command line
void Command_line(int arg, char *args[]){
  int Argument_cmd=1;
  while (Argument_cmd<arg)
    {
      if(strcmp(args[Argument_cmd],"-N")==0)
        {
    	  Argument_cmd += 1;
    	  N = atoi(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-A")==0)
        {
    	  Argument_cmd += 1;
    	  A = atof(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-Pe_r")==0)
        {
          Argument_cmd += 1;
          Pe_r = atof(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-kappa")==0)
        {
    	  Argument_cmd += 1;
    	  kappa = atof(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-mu")==0)
        {
    	  Argument_cmd += 1;
    	  mu = atof(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-sigma_m")==0)
        {
          Argument_cmd += 1;
          sigma_m = atof(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-theta0")==0)
        {
          Argument_cmd += 1;
          theta0 = atof(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-dt")==0)
        {
    	  Argument_cmd += 1;
    	  dt = atof(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-eps")==0)
    	{
    	  Argument_cmd += 1;
    	  eps = atof(args[Argument_cmd]);
    	  Argument_cmd += 1;
    	}
      else if(strcmp(args[Argument_cmd],"-tmax")==0)
        {
    	  Argument_cmd += 1;
    	  tmax = atof(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-n")==0)
    	{
    	  Argument_cmd += 1;
    	  run_index = atoi(args[Argument_cmd]);
    	  Argument_cmd += 1;
    	}
      else if(strcmp(args[Argument_cmd],"-save_trajs")==0)
    	{
    	  Argument_cmd += 1;
    	  save_trajs = 1;
    	}
      else if(strcmp(args[Argument_cmd],"-save_epr")==0)
        {
          Argument_cmd += 1;
          save_epr = 1;
        }
      else if(strcmp(args[Argument_cmd],"-save_moments")==0)
        {
          Argument_cmd += 1;
          save_moments = 1;
        }
      else if(strcmp(args[Argument_cmd],"-Nsave_trajs")==0)
    	{
    	  Argument_cmd += 1;
    	  Nsave_trajs = atoi(args[Argument_cmd]);
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-tsave")==0)
        {
          Argument_cmd += 1;
          t_save = atoi(args[Argument_cmd]);
          Argument_cmd += 1;
        }	
      else if(strcmp(args[Argument_cmd],"-t_steady")==0)
        {
          Argument_cmd += 1;
          t_steady = atoi(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-R_cutoff")==0)
        {
          Argument_cmd += 1;
          R_cutoff = atof(args[Argument_cmd]);
          Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-name_flow")==0)
        {
    	  Argument_cmd += 1;
    	  name_flow = args[Argument_cmd];
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-policy")==0)
        {
    	  Argument_cmd += 1;
    	  policy = args[Argument_cmd];
    	  Argument_cmd += 1;
        }
      else if(strcmp(args[Argument_cmd],"-target_point")==0)
        {
          Argument_cmd += 1;
          target_point = 1;
        }
      else
        {
    	  Argument_cmd += 1;
        }
    }
}

//Default values of the sim parameters if not initialized from command line
static inline void Set_default_values(){
    A = 0.;
    kappa = 1.;
    Pe_r = 10.;
    mu = 1.;
    sigma_m = 0.1;
    Pe_m = 1./(sigma_m*mu);
    N = 1;
    R_cutoff = 0.;

    name_flow = "Const";
    policy = "SP";
    
    tmax = 100.;
    dt = 1.e-4;
    t_save = 100;
    t_steady = 0;

    save_trajs = 0;
    Nsave_trajs = 1;
    save_epr = 0;
    save_moments = 0;

    theta0 = 0.;
    init_pos.push_back(0.1);
    init_pos.push_back(0.1);
    eps = 1.;

    target_point = 0;
    
    run_index = 0;
}


//Routine to set all sim parameters as well as the name structure of the output files
static inline void Set_parameters(int arg, char *args[])
{
  Pe_m = 1./(sigma_m*mu);
  Parameters[P_A] = A*dt;
  Parameters[P_dt] = dt;
  Parameters[P_epsilon] = eps;
  Parameters[P_theta0] = theta0;
  Parameters[P_kappa] = kappa*dt;
  Parameters[P_mu] = mu*dt;
  Parameters[P_Pe_r] = sqrt(2.*dt/Pe_r);
  Parameters[P_Pe_m] = sqrt(2.*dt/Pe_m);
  Parameters[P_Rcutoff] = R_cutoff;
  
  if(name_flow == "Const")
    Parameters[P_flow] = Const;
  else
    cout << "Wrong flow name" << endl;
  
  if(policy == "SP")
    Parameters[P_policy] = SP;
  else if(policy == "AAP")
    Parameters[P_policy] = AAP;
  else
    cout << "Selected policy does not exist" << endl;

  if(target_point)
    name_params = "Target_point/";
  else
    name_params = "Target_path/";
  
  if(Parameters[P_policy]==SP){
    name_params += policy + "_Pe_r" + double_to_string(Pe_r,4) + "_kappa" + double_to_string(kappa,4) + "_sigma_m" + double_to_string(sigma_m,4) + "_mu" \
      + double_to_string(mu,4) + "_A" + double_to_string(A,2) + "_dt" + double_to_string(dt,6) + "_theta" + double_to_string(theta0,3) + "_tmax" + int_to_string(tmax);
  }
  else if(Parameters[P_policy]==AAP){
    name_params += policy + "_Pe_r" + double_to_string(Pe_r,4) + "_kappa" + double_to_string(kappa,4) + "_sigma_m" + double_to_string(sigma_m,4) + "_mu" \
      + double_to_string(mu,4) + "_A" + double_to_string(A,2) + "_dt" + double_to_string(dt,6) + "_theta" + double_to_string(theta0,3) + "_tmax" + int_to_string(tmax) + "_eps" + double_to_string(eps,2);
  }

  if(target_point)
    name_params += "_Rcutoff" + double_to_string(R_cutoff,3);
  
  if(t_steady>0)
    name_params += "_tsteady" + int_to_string(t_steady);
}

//Routine to print CPU time at the end of the simulation
void print_time(const clock_t start_time, const clock_t end_time)
{
    double length = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    
    printf("Total running time = ");
    
    int Nhours = floor(length)/3600;
    
    if(Nhours > 0)
    {
        int Nmin = floor(length - 3600*Nhours)/60;
        printf("%dh%dm%ds\n",Nhours,Nmin,(int)(length - 3600*Nhours - 60*Nmin));
    }
    
    else
    {
        int Nmin = floor(length)/60;
        if(Nmin > 0)
        printf("%dm%ds\n",Nmin,(int)(length - 60*Nmin));
        else
        {
            if(length >= 1)
            printf("%1.2lfs\n",length);
            else
           printf("%1.2es\n",length);
        }
    }
}
