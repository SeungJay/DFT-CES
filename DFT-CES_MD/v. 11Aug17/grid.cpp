/* ----------------------------------------------------------------------
   DFT-CES core subroutines. Written by H.-K. Lim
   Copyright (C) 2016 M-design group @ KAIST
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "fix_gridforce.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

Grid::Grid(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  int me = comm -> me;
  int i;

  natoms = 0;
  atomns = NULL;
  basis  = NULL;
  gvin = NULL;
  gvout = NULL;
  gvout_all = NULL;

  if (narg != 1) error->all(FLERR,"Illegal grid command");

  if (strcmp(arg[0],"empty") != 0 ){ 
    // read header until atom information
    read_header(arg[0]);

    if(gx[1]+gx[2]+gy[0]+gy[2]+gz[0]+gz[1] != 0) error->all(FLERR,"Grid file is not orthorhombic");
    if(me == 0){
      printf("Grid file has been parsed: %s\n",arg[0]);
      printf("DFT-CES: # of atoms in grid: %d\n", natoms);
      printf("DFT-CES: # of grid points: %d %d %d\n", gnx, gny, gnz);
      printf("DFT-CES: grid spacing: %f %f %f Angs\n", gx[0], gy[1], gz[2]);
    }
    memory->grow(atomns, natoms, "grid:atomns");
    memory->grow(basis, natoms, 3, "grid:basis");
    memory->grow(gvin, gnx*gny*gnz, "grid:gvin");
    memory->grow(gvout, gnx*gny*gnz, "grid:gvout");
    if(me==0)memory->grow(gvout_all, gnx*gny*gnz, "grid:gvout");

    // read atom and grid value information
    read_content(arg[0]);
  }
  else{
    return;
  }
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->destroy(atomns);
  memory->destroy(basis);
  memory->destroy(gvin);
  memory->destroy(gvout);
  memory->destroy(gvout_all);
}

/* ---------------------------------------------------------------------- */

void Grid::read_header(char *filename)
{
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];
  double temp[4];

  if (me == 0) {
    fptr = fopen(filename,"r");
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open grid file %s",filename);
      error->one(FLERR,str);
    }

    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr); // skip first two rows
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lf %lf %lf", &natoms, &temp[0], &temp[1], &temp[2]);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lf %lf %lf", &gnx, &gx[0], &gx[1], &gx[2]);
    gx[0] *= 0.52917721;
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lf %lf %lf", &gny, &gy[0], &gy[1], &gy[2]);
    gy[1] *= 0.52917721;
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lf %lf %lf", &gnz, &gz[0], &gz[1], &gz[2]);
    gz[2] *= 0.52917721;
    
    fclose(fptr);
  }

  MPI_Bcast(&natoms,1,MPI_INT,0,world);
  MPI_Bcast(&gnx,1,MPI_INT,0,world);
  MPI_Bcast(&gny,1,MPI_INT,0,world);
  MPI_Bcast(&gnz,1,MPI_INT,0,world);
  MPI_Bcast(&gx[0],3,MPI_DOUBLE,0,world);
  MPI_Bcast(&gy[0],3,MPI_DOUBLE,0,world);
  MPI_Bcast(&gz[0],3,MPI_DOUBLE,0,world);

}

void Grid::read_content(char *filename)
{
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE], *str_ptr;
  double temp;
  int i, cnt;

  if (me == 0) {
    fptr = fopen(filename,"r");
    for(i=0; i<6; i++) fgets(line,MAXLINE,fptr);
    for(i=0; i<natoms; i++){
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lf %lf %lf %lf", &atomns[i], &temp, &basis[i][0], &basis[i][1], &basis[i][2]);
      basis[i][0] *= 0.52917721;
      basis[i][1] *= 0.52917721;
      basis[i][2] *= 0.52917721;
    }
    cnt = 0;
    while(fgets(line,MAXLINE,fptr)!=NULL){
      str_ptr = strtok(line," ");
      for(; str_ptr!=NULL; cnt++) {
        sscanf(str_ptr,"%lf",&gvin[cnt]);
        gvin[cnt] *= 13.60569253;              // Ry to eV
        str_ptr = strtok(NULL," ");
      }
    }    
    fclose(fptr);
  }

  MPI_Bcast(&atomns[0],natoms,MPI_INT,0,world);
  MPI_Bcast(&basis[0][0],3*natoms,MPI_DOUBLE,0,world);
  MPI_Bcast(&gvin[0],gnx*gny*gnz,MPI_DOUBLE,0,world);

}

void Grid::save_grid(char* filename, int nsteps)
{
  int me = comm->me;
  FILE *fptr;
  int i, j, k, cnt;
  
  MPI_Reduce(gvout,gvout_all,gnx*gny*gnz,MPI_DOUBLE,MPI_SUM,0,world);

  if (me == 0 ){
    fptr = fopen(filename, "w");
    fprintf(fptr,"MDrho from DFT-CES\n");
    fprintf(fptr,"ZYX, Bohr unit\n");
    fprintf(fptr,"% 5d    0.000000    0.000000    0.000000\n",natoms);
    fprintf(fptr,"% 5d% 12.6lf    0.000000    0.000000\n",gnx,gx[0]/0.52917721);
    fprintf(fptr,"% 5d    0.000000% 12.6lf    0.000000\n",gny,gy[1]/0.52917721);
    fprintf(fptr,"% 5d    0.000000    0.000000% 12.6lf\n",gnz,gz[2]/0.52917721);
    for(i=0; i<natoms; i++) {
      fprintf(fptr,"% 5d% 12.6lf% 12.6lf% 12.6lf% 12.6lf\n",atomns[i],(double)atomns[i],basis[i][0]/0.52917721,basis[i][1]/0.52917721,basis[i][2]/0.52917721);
    }
    for(i=0;i<gnx;i++){
      for(j=0;j<gny;j++){
        cnt = 0;
        for(k=0;k<gnz;k++){
          fprintf(fptr,"% 13.5lE",gvout_all[k+j*gnz+i*gnz*gny]*pow(0.52917721,3)/(nsteps+1));
          if(cnt%6 == 5 && k < gnz-1) fprintf(fptr,"\n");
          cnt++;
        }
        fprintf(fptr,"\n");
      }
    }
    fclose(fptr);
    printf("DFT-CES: MDrho.cube has been successfully saved\n");
  }
}

