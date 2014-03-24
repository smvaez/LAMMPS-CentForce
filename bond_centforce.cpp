/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_centforce.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondCentforce::BondCentforce(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondCentforce::~BondCentforce()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(Rr0);
    memory->destroy(cf_eq_sign);
  }
}

/* ---------------------------------------------------------------------- */

void BondCentforce::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx0,dely0,delz0,fx,fy,fz;
  double delx,dely,delz,r0x,r0y,r0z,dir0x,dir0y,dir0z,dirlen,ebond,fbond;
  double r,rk,cdot;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **x0 = atom->x0;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);
    
    delx0 = x0[i1][0] - x0[i2][0];
    dely0 = x0[i1][1] - x0[i2][1];
    delz0 = x0[i1][2] - x0[i2][2];
    domain->minimum_image(delx0,dely0,delz0);
    
    r0x = Rr0[type][0];
    r0y = Rr0[type][1];
    r0z = Rr0[type][2];
    domain->minimum_image(r0x,r0y,r0z);
    
    cdot = delx0*r0x + dely0*r0y + delz0*r0z;
    
    if(cf_eq_sign[type] == 1){
      if (cdot < 0.0 ) {
        r0x *= -1;
        r0y *= -1;
        r0z *= -1;
      }
    }
    else if (cdot >= 0.0 ) {
        r0x *= -1;
        r0y *= -1;
        r0z *= -1;
    }
      
    dirlen = sqrt(Rr0[type][3]*Rr0[type][3]+Rr0[type][4]*Rr0[type][4]+Rr0[type][5]*Rr0[type][5]);
    dir0x = Rr0[type][3]/dirlen;
    dir0y = Rr0[type][4]/dirlen;
    dir0z = Rr0[type][5]/dirlen;

    r = (delx-r0x)*dir0x + (dely-r0y)*dir0y + (delz-r0z)*dir0z;
    
    // force & energy

    rk = k[type] * r;

    fbond = -2.0*rk;

    if (eflag) ebond = rk*r;

    // apply force to each of 2 atoms
    
    fx = dir0x*fbond;
    fy = dir0y*fbond;
    fz = dir0z*fbond;

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fx;
      f[i1][1] += fy;
      f[i1][2] += fz;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fx;
      f[i2][1] -= fy;
      f[i2][2] -= fz;
    }

    if (evflag) ev_tally_bondcf_xyz(i1,i2,nlocal,newton_bond,ebond,fx,fy,fz,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondCentforce::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(Rr0,n+1,6,"bond:Rr0");
  memory->create(cf_eq_sign,n+1,"bond:cf_eq_sign");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondCentforce::coeff(int narg, char **arg)
{
  double *arr0;
  int  cf_eq_sign_one;
  
  if (!atom->angxo_flag && !atom->bondxo_flag) error->all(FLERR,"Bond style centforce requires atom style bondxo or angxo");
  if (narg != 9) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();
  
  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  
  memory->create(arr0,6,"bond:arr0");

  double k_one = force->numeric(arg[1]);
  arr0[0] = force->numeric(arg[2]);
  arr0[1] = force->numeric(arg[3]);
  arr0[2] = force->numeric(arg[4]);
  arr0[3] = force->numeric(arg[5]);
  arr0[4] = force->numeric(arg[6]);
  arr0[5] = force->numeric(arg[7]);
  cf_eq_sign_one = force->inumeric(arg[8]);
 
  if (cf_eq_sign_one != 1 && cf_eq_sign_one != -1) 
    error->all(FLERR,"Incorrect args for bond coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    cf_eq_sign[i] = cf_eq_sign_one;
    for (int j = 0; j <= 5; j++) Rr0[i][j] = arr0[j];
    setflag[i] = 1;
    count++;
  }
  
  memory->destroy(arr0);
  
  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondCentforce::equilibrium_distance(int i)
{
  double dx0 , dy0 , dz0;
  
  dx0 = Rr0[i][0];
  dy0 = Rr0[i][1];
  dz0 = Rr0[i][2];
  domain->minimum_image(dx0,dy0,dz0);

  return (sqrt(dx0*dx0+dy0*dy0+dz0*dz0));
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondCentforce::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  for (int i = 1; i <= atom->nbondtypes; i++) 
    fwrite(Rr0[i],sizeof(double),6,fp);
  fwrite(&cf_eq_sign[1],sizeof(int),atom->nbondtypes,fp);
 
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondCentforce::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    for (int i = 1; i <= atom->nbondtypes; i++) 
      fread(Rr0[i],sizeof(double),6,fp);
    fread(&cf_eq_sign[1],sizeof(int),atom->nbondtypes,fp);
    
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  for (int i = 1; i <= atom->nbondtypes; i++) 
    MPI_Bcast(Rr0[i],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&cf_eq_sign[1],atom->nbondtypes,MPI_INT,0,world);
  
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondCentforce::single(int type, double rsq, int i, int j)
{
  double delx ,dely ,delz ,r0x ,r0y ,r0z ,dir0x ,dir0y ,dir0z ,dirlen ,dr ,rk ;
  double delx0 ,dely0 ,delz0 ,cdot ;
  
  double **x = atom->x;
  double **x0 = atom->x0;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  domain->minimum_image(delx,dely,delz);
  
  delx0 = x0[i][0] - x0[j][0];
  dely0 = x0[i][1] - x0[j][1];
  delz0 = x0[i][2] - x0[j][2];
  domain->minimum_image(delx0,dely0,delz0);
  
  r0x = Rr0[type][0];
  r0y = Rr0[type][1];
  r0z = Rr0[type][2];
  domain->minimum_image(r0x,r0y,r0z);
  
  cdot = delx0*r0x + dely0*r0y + delz0*r0z;
    
  if(cf_eq_sign[type] != 1){
      if (cdot < 0.0 ) {
        r0x *= -1;
        r0y *= -1;
        r0z *= -1;
      }
  }
  else if (cdot >= 0.0 ) {
        r0x *= -1;
        r0y *= -1;
        r0z *= -1;
  }
    
  dirlen = sqrt(Rr0[type][3]*Rr0[type][3]+Rr0[type][4]*Rr0[type][4]+Rr0[type][5]*Rr0[type][5]);
  dir0x = Rr0[type][3]/dirlen;
  dir0y = Rr0[type][4]/dirlen;
  dir0z = Rr0[type][5]/dirlen;

  dr = (delx-r0x)*dir0x + (dely-r0y)*dir0y + (delz-r0z)*dir0z;

  rk = k[type] * dr;

  return rk*dr;
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void BondCentforce::ev_tally_bondcf_xyz(int i, int j, int nlocal, int newton_bond,
                    double ebond, double fx, double fy, double fz,
                        double delx, double dely, double delz)
{
  double ebondhalf,vv[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy += ebond;
      else {
        ebondhalf = 0.5*ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5*ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    vv[0] = delx*fx;
    vv[1] = dely*fy;
    vv[2] = delz*fz;
    vv[3] = delx*fy;
    vv[4] = delx*fz;
    vv[5] = dely*fz;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += vv[0];
        virial[1] += vv[1];
        virial[2] += vv[2];
        virial[3] += vv[3];
        virial[4] += vv[4];
        virial[5] += vv[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*vv[0];
          virial[1] += 0.5*vv[1];
          virial[2] += 0.5*vv[2];
          virial[3] += 0.5*vv[3];
          virial[4] += 0.5*vv[4];
          virial[5] += 0.5*vv[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*vv[0];
          virial[1] += 0.5*vv[1];
          virial[2] += 0.5*vv[2];
          virial[3] += 0.5*vv[3];
          virial[4] += 0.5*vv[4];
          virial[5] += 0.5*vv[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5*vv[0];
        vatom[i][1] += 0.5*vv[1];
        vatom[i][2] += 0.5*vv[2];
        vatom[i][3] += 0.5*vv[3];
        vatom[i][4] += 0.5*vv[4];
        vatom[i][5] += 0.5*vv[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5*vv[0];
        vatom[j][1] += 0.5*vv[1];
        vatom[j][2] += 0.5*vv[2];
        vatom[j][3] += 0.5*vv[3];
        vatom[j][4] += 0.5*vv[4];
        vatom[j][5] += 0.5*vv[5];
      }
    }
  }
}
