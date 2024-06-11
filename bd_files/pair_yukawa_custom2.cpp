// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_yukawa_custom2.h"

#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairYukawaCustom2::PairYukawaCustom2(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairYukawaCustom2::~PairYukawaCustom2()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(rad);
    memory->destroy(cut);
    memory->destroy(a);
    memory->destroy(b);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairYukawaCustom2::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,radi,radj;
  double rsq,r,rinv,screening1,screening2,forceyukawa,hx,hy,uvdw,factor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        screening1 = exp(-kappa1*(r-(radi+radj)));
        screening2 = exp(-kappa2*(r-(radi+radj)));
        forceyukawa = (a[itype][jtype] * screening1) + (b[itype][jtype] * screening2);

        hx = (r-(radi+radj))/(2*radi);
        hy = radi/radj;
        
        uvdw = -(hamaker/12) * ( ( hy / ( (hx*hx) + (hx*hy) + hx )) + ( hy / ( (hx*hx) + (hx*hy) + hx + hy )) + ( 2 * log( ( (hx*hx) + (hx*hy) + hx ) / ( (hx*hx) + (hx*hy) + hx + hy ) ) ));

        fpair = (factor*forceyukawa * rinv)+(factor*uvdw * rinv);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = (a[itype][jtype]/kappa1 * screening1 - offset[itype][jtype]) + (b[itype][jtype]/kappa2 * screening2 - offset[itype][jtype]) + (uvdw-offset[itype][jtype]);
          evdwl *= factor;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
    init specific to this pair style
------------------------------------------------------------------------- */

void PairYukawaCustom2::init_style()
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair yukawacustom requires atom style sphere");

  neighbor->request(this,instance_me);

  // require that atom radii are identical within each type, I have removed

}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairYukawaCustom2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(rad,n+1,"pair:rad");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(b,n+1,n+1,"pair:b");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairYukawaCustom2::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  kappa1 = utils::numeric(FLERR,arg[0],false,lmp);
  kappa2 = utils::numeric(FLERR,arg[1],false,lmp);
  hamaker = utils::numeric(FLERR,arg[2],false,lmp);
  cut_global = utils::numeric(FLERR,arg[3],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairYukawaCustom2::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a_one = utils::numeric(FLERR,arg[2],false,lmp);
  double b_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      b[i][j] = b_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairYukawaCustom2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a[i][j] = mix_energy(a[i][i],a[j][j],1.0,1.0);
    b[i][j] = mix_energy(b[i][i],b[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag && (cut[i][j] > 0.0)) {
    double screening1 = exp(-kappa1 * (cut[i][j] - (rad[i]+rad[j])));
    double screening2 = exp(-kappa2 * (cut[i][j] - (rad[i]+rad[j])));
    double hx = (cut[i][j]-(rad[i]+rad[j]))/(2*rad[i]);
    double hy = rad[i]/rad[j];
    //double hx = (cut[i][j] - (rad[i]+rad[j]))/4.82;
    //double hy = 1;
    double uvdw = -(hamaker/12) * ( ( hy / ( (hx*hx) + (hx*hy) + hx )) + ( hy / ( (hx*hx) + (hx*hy) + hx + hy )) + ( 2 * log( ( (hx*hx) + (hx*hy) + hx ) / ( (hx*hx) + (hx*hy) + hx + hy ) ) ));

    offset[i][j] = (a[i][j]/kappa1 * screening1) + (b[i][j]/kappa2 * screening2) + uvdw;
  } else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  b[j][i] = b[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairYukawaCustom2::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&b[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairYukawaCustom2::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&a[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&b[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairYukawaCustom2::write_restart_settings(FILE *fp)
{
  fwrite(&kappa1,sizeof(double),1,fp);
  fwrite(&kappa2,sizeof(double),1,fp);
  fwrite(&hamaker,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairYukawaCustom2::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&kappa1,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&kappa2,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&hamaker,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&kappa1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&hamaker,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairYukawaCustom2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,a[i][i],b[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairYukawaCustom2::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,a[i][j],b[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairYukawaCustom2::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                          double /*factor_coul*/, double factor_lj,
                          double &fforce)
{
  double r,rinv,screening1,screening2,forceyukawa,hx,hy,uvdw,phi;

  r = sqrt(rsq);
  rinv = 1.0/r;
  screening1 = exp(-kappa1*(r-(rad[itype]+rad[jtype])));
  screening2 = exp(-kappa2*(r-(rad[itype]+rad[jtype])));;
  forceyukawa = (a[itype][jtype] * screening1) + (b[itype][jtype] * screening2);
  //uvdw = -(2*((rad[itype]*rad[jtype])/(rad[itype]+rad[jtype])))/(6.0*(r-(rad[itype]+rad[jtype])));
  hx = (r-(rad[itype]+rad[jtype]))/4.82; 
  hy = 1;
  uvdw = -(hamaker/12) * ( ( hy / ( (hx*hx) + (hx*hy) + hx )) + ( hy / ( (hx*hx) + (hx*hy) + hx + hy )) + ( 2 * log( ( (hx*hx) + (hx*hy) + hx ) / ( (hx*hx) + (hx*hy) + hx + hy ) ) )); 

  fforce = (factor_lj*forceyukawa * rinv)+(factor_lj*uvdw * rinv);

  phi = (a[itype][jtype]/kappa1 * screening1) + (b[itype][jtype]/kappa2 * screening2) + uvdw - offset[itype][jtype];
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairYukawaCustom2::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "a") == 0) return (void *) a;
  if (strcmp(str, "b") == 0) return (void *) b;
  return nullptr;
}
