/* ----------------------------- Output folders ---------------------------- */

void init_output_folders(){

    string fname; 
    
    fname="../output/" + name_flow + "/Target_path/Trajectories";
    system(("mkdir -p "+fname).c_str());
    fname="../output/" + name_flow + "/Trajectories/Target_path";
    system(("mkdir -p "+fname).c_str());
    fname="../output/" + name_flow + "/EPR/Target_path";
    system(("mkdir -p "+fname).c_str());
    fname="../output/" + name_flow + "/EPR/Target_point";
    system(("mkdir -p "+fname).c_str());
    fname="../output/" + name_flow + "/Moments/Target_path";
    system(("mkdir -p "+fname).c_str());
    fname="../output/" + name_flow + "/Moments/Target_point";
    system(("mkdir -p "+fname).c_str());
}

/* --------------------------- Saving routines ----------------------------- */

void Save_vector(const string name, const vect traj, const int col){
  ofstream out;
  out.open(name.c_str(),ios::out); 

  if(out.is_open())
    {
      int Npoints = traj.size()/col;
      for(int i=0;i<Npoints;i++){
	for(int j=0;j<col;j++)
	  out << traj[col*i+j] << "   ";
	out << endl;
      }
      out.close();
    }
  else
    {
      cout << "ERROR: Unable to open the file: " << name << endl;
    }    
}


// Function to write a vector of doubles to a binary file
void Save_binary_vect(const string filename, const vect& data) {
  
  ofstream out;

  out.open(filename.c_str(),ios::out|ios::binary);

  // Write the vector data to the file
  out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));

  // Close the file
  out.close();
}
