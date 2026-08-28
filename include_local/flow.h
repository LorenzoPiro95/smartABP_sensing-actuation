//Constant flow
static inline vect Const_flow(const vect pos, const double *Parameters)
{
    vect F(DIMENSION);
    
    F[DIR_X] = Parameters[P_A];
    F[DIR_Y] = 0.;
    
    return F;
}


//Flow selection (real value)
static inline vect Compute_flow_real(const particle *P, const double *Parameters){
  vect F;
  
  if(Parameters[P_flow] == Const)
    F = Const_flow(P->pos,Parameters);

  return F;
}
