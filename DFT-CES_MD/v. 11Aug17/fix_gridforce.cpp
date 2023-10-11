/* ----------------------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_gridforce.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "grid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGridForce::FixGridForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix gridforce command");
  scalar_flag = 1;
  extscalar = 1;
  
  weight = atof(arg[3]);
  sfactor = atoi(arg[4]);

  if(comm->me == 0) printf("DFT-CES: weight factor    = %f\n", weight  );
  if(comm->me == 0) printf("DFT-CES: supercell factor = %d\n", sfactor );

  force_flag = 0;
  Egrid = 0.0;
  tmp3[0]=tmp3[1]=tmp3[2]=0;
  tmp4[0]=tmp4[1]=tmp4[2]=tmp4[3]=0;
}

/* ---------------------------------------------------------------------- */

FixGridForce::~FixGridForce()
{
  return;
}

/* ---------------------------------------------------------------------- */

int FixGridForce::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGridForce::init()
{
  // check variables
  if (domain->grid->natoms == 0) error->all(FLERR,"Grid data is unavilable");
}

/* ---------------------------------------------------------------------- */

void FixGridForce::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGridForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double fx, fy, fz;

  Egrid = 0.0;
  force_flag = 0;
  tmp3[0]=tmp3[1]=tmp3[2]=0;
  tmp4[0]=tmp4[1]=tmp4[2]=tmp4[3]=0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
        triInter(q[i], x[i][0], x[i][1], x[i][2], tmp4);
        Egrid += weight*23.06092*q[i]*tmp4[0];
        fx = weight*23.06092*q[i]*tmp4[1];
        fy = weight*23.06092*q[i]*tmp4[2];
        fz = weight*23.06092*q[i]*tmp4[3];
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
    }
  }
}

double FixGridForce::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(&Egrid,&Egrid_all,1,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return Egrid_all;
}

void FixGridForce::triInter(double q, double x, double y, double z, double* tmp4){
  int i;
  double* gvin = domain->grid->gvin;
  double* gvout = domain->grid->gvout;
  double px, py, pz, xd, yd, zd;
  double p000, p100, p010, p001, p110, p101, p011, p111;
  double e000[3], e100[3], e010[3], e001[3], e110[3], e101[3], e011[3], e111[3];
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gsx = domain->grid->gx[0];
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];
  double gvol;
  if(x<0){
    px =(fmod(x,gnx*gsx)+gnx*gsx)/gsx;
  }else{
    px = fmod(x,gnx*gsx)/gsx;
  }
  if(y<0){
    py =(fmod(y,gny*gsy)+gny*gsy)/gsy;
  }else{
    py = fmod(y,gny*gsy)/gsy;
  }
  if(z<0){
    pz =(fmod(z,gnz*gsz)+gnz*gsz)/gsz;
  }else{
    pz = fmod(z,gnz*gsz)/gsz;
  }
  xd=(double)(px-(int)px);
  yd=(double)(py-(int)py);
  zd=(double)(pz-(int)pz);
  
  tmp3[0]=tmp3[1]=tmp3[2]=0;
  tmp4[0]=tmp4[1]=tmp4[2]=tmp4[3]=0;    

  p000 = gridValue(gvin, (int)px, (int)py, (int)pz);
  p100 = gridValue(gvin, 1+(int)px, (int)py, (int)pz);
  p010 = gridValue(gvin, (int)px, 1+(int)py, (int)pz);
  p001 = gridValue(gvin, (int)px, (int)py, 1+(int)pz);
  p110 = gridValue(gvin, 1+(int)px, 1+(int)py, (int)pz);
  p101 = gridValue(gvin, 1+(int)px, (int)py, 1+(int)pz);
  p011 = gridValue(gvin, (int)px, 1+(int)py, 1+(int)pz);
  p111 = gridValue(gvin, 1+(int)px, 1+(int)py, 1+(int)pz);
  tmp4[0] = triInterValue(p000,p100,p010,p001,p110,p101,p011,p111,xd,yd,zd);
  gridGrad(gvin, ((int)px)%gnx, ((int)py)%gny, ((int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e000[i] = tmp3[i];
  gridGrad(gvin, (1+(int)px)%gnx, ((int)py)%gny, ((int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e100[i] = tmp3[i];
  gridGrad(gvin, ((int)px)%gnx, (1+(int)py)%gny, ((int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e010[i] = tmp3[i];
  gridGrad(gvin, ((int)px)%gnx, ((int)py)%gny, (1+(int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e001[i] = tmp3[i];
  gridGrad(gvin, (1+(int)px)%gnx, (1+(int)py)%gny, ((int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e110[i] = tmp3[i];
  gridGrad(gvin, (1+(int)px)%gnx, ((int)py)%gny, (1+(int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e101[i] = tmp3[i];
  gridGrad(gvin, ((int)px)%gnx, (1+(int)py)%gny, (1+(int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e011[i] = tmp3[i];
  gridGrad(gvin, (1+(int)px)%gnx, (1+(int)py)%gny, (1+(int)pz)%gnz, tmp3);
  for(i=0;i<3;i++) e111[i] = tmp3[i];
  tmp4[1] = triInterValue(e000[0],e100[0],e010[0],e001[0],e110[0],e101[0],e011[0],e111[0],xd,yd,zd);
  tmp4[2] = triInterValue(e000[1],e100[1],e010[1],e001[1],e110[1],e101[1],e011[1],e111[1],xd,yd,zd);
  tmp4[3] = triInterValue(e000[2],e100[2],e010[2],e001[2],e110[2],e101[2],e011[2],e111[2],xd,yd,zd);

  // inverse trilinear interpolation for saving MD rho
  p000=(1-xd)*(1-yd)*(1-zd)*abs(weight)*q;
  p100=xd*(1-yd)*(1-zd)*abs(weight)*q;
  p010=(1-xd)*yd*(1-zd)*abs(weight)*q;
  p110=xd*yd*(1-zd)*abs(weight)*q;
  p001=(1-xd)*(1-yd)*zd*abs(weight)*q;
  p101=xd*(1-yd)*zd*abs(weight)*q;
  p011=(1-xd)*yd*zd*abs(weight)*q;
  p111=xd*yd*zd*abs(weight)*q;

  gvol = gsx*gsy*gsz;

  gvout[((int)pz)%gnz+(((int)py)%gny)*gnz+(((int)px)%gnx)*gnz*gny] += p000/gvol/sfactor;
  gvout[((int)pz)%gnz+(((int)py)%gny)*gnz+(((int)px+1)%gnx)*gnz*gny] += p100/gvol/sfactor;
  gvout[((int)pz)%gnz+(((int)py+1)%gny)*gnz+(((int)px)%gnx)*gnz*gny] += p010/gvol/sfactor;
  gvout[((int)pz+1)%gnz+(((int)py)%gny)*gnz+(((int)px)%gnx)*gnz*gny] += p001/gvol/sfactor;
  gvout[((int)pz)%gnz+(((int)py+1)%gny)*gnz+(((int)px+1)%gnx)*gnz*gny] += p110/gvol/sfactor;
  gvout[((int)pz+1)%gnz+(((int)py)%gny)*gnz+(((int)px+1)%gnx)*gnz*gny] += p101/gvol/sfactor;
  gvout[((int)pz+1)%gnz+(((int)py+1)%gny)*gnz+(((int)px)%gnx)*gnz*gny] += p011/gvol/sfactor;
  gvout[((int)pz+1)%gnz+(((int)py+1)%gny)*gnz+(((int)px+1)%gnx)*gnz*gny] += p111/gvol/sfactor;

  return;
}

double FixGridForce::triInterValue(double c000,double c100, double c010, double c001, double c110, double c101, double c011, double c111, double xd,double yd,double zd){
  double c00, c10, c01, c11, c0, c1, c;
  c00 = c000*(1-xd) + c100*xd;
  c01 = c001*(1-xd) + c101*xd;
  c10 = c010*(1-xd) + c110*xd;
  c11 = c011*(1-xd) + c111*xd;
  c0 = c00*(1-yd) + c10*yd;
  c1 = c01*(1-yd) + c11*yd;
  c = c0*(1-zd) + c1*zd;
  return c;
}

double FixGridForce::gridValue(double* gvin, int x, int y, int z){
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gValue=0;

  gValue=gvin[(z%gnz)+(y%gny)*gnz+(x%gnx)*gnz*gny];

 return gValue;
}

void FixGridForce::gridGrad(double* gvin, int x, int y, int z, double* tmp3){
  int gnx = domain->grid->gnx;
  int gny = domain->grid->gny;
  int gnz = domain->grid->gnz;
  double gsx = domain->grid->gx[0];
  double gsy = domain->grid->gy[1];
  double gsz = domain->grid->gz[2];

  tmp3[0]=tmp3[1]=tmp3[2]=0;

  tmp3[0] = -1*(gvin[z+y*gnz+((x+1)%gnx)*gnz*gny]-gvin[z+y*gnz+(x?(x-1):(gnx-1))*gnz*gny])/(2*gsx);
  tmp3[1] = -1*(gvin[z+((y+1)%gny)*gnz+x*gnz*gny]-gvin[z+(y?(y-1):(gny-1))*gnz+x*gnz*gny])/(2*gsy);
  tmp3[2] = -1*(gvin[((z+1)%gnz)+y*gnz+x*gnz*gny]-gvin[(z?(z-1):(gnz-1))+y*gnz+x*gnz*gny])/(2*gsz);

  return;
}
