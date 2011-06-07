/*************************************************************
 *
 *  findvel.c - subroutines for velocity finding in part-nd
 *
 *  Copyright (C) 2000-2003,2005  Mark J. Stock, mstock@umich.edu
 *
 *  This file is part of part-nd.
 *
 *  part-nd is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  part-nd is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with part-nd; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 *********************************************************** */


#include "structs.h"

int find_new_vels(sim_ptr,cell_ptr,cell_ptr,int);
int find_acc_on_this_cells_parts(sim_ptr,cell_ptr,cell_ptr,int);
int find_grav_acc_on_this_part(sim_ptr,cell_ptr,particle_ptr);
int find_direct_grav_acc_on_this_part(sim_ptr,cell_ptr,particle_ptr);
int find_cont_acc_on_this_part(sim_ptr,cell_ptr,particle_ptr);
int find_acc_on_this_cell(sim_ptr,cell_ptr,cell_ptr,FLOAT*);
FLOAT find_vmax(cell_ptr,cell_ptr,FLOAT);


/* 
 *  Find the new velocity of each particle in the simulation
 */
int find_new_vels(sim_ptr sim,cell_ptr top,cell_ptr curr_cell,int uga) {

   int i;
   particle_ptr this;

   // if (sim->step > 472) {
   // fprintf(stdout,"this cell is at level %d, has %d particles\n",curr_cell->level,curr_cell->num);
   // fflush;
   // }

   /* first, set all particle accelerations to zero */
   // this is done in advect_nodes, in main.c


   if (curr_cell->has_subcells) {

      /* if this cell has subcells, run this subroutine on those */
      for (i=0;i<NCHILD;i++)
         find_new_vels(sim,top,curr_cell->s[i],uga);

   } else {

      /* otherwise, find the accelerations on each particle in this cell
       * as induced by all other particles in the system */

      find_acc_on_this_cells_parts(sim,top,curr_cell,uga);

   }

   return(0);
}


/*
 *  Find the new velocity of a specific particle in the simulation
 */
int find_acc_on_this_cells_parts(sim_ptr sim,cell_ptr top,cell_ptr curr_cell,int update_grav) {

   int i,d,xdir,ydir;
   FLOAT true_dist,rel_vel,distsq,sacc;
   FLOAT kacc,cacc,save[DIM];
   particle_ptr curr = curr_cell->first;

   #pragma omp parallel private(i,d,xdir,ydir,true_dist,rel_vel,distsq,sacc,kacc,cacc,save)
   {
   #pragma omp single nowait
   {
   while (curr) {

      #pragma omp task firstprivate(curr)
      {

      /* new way: find contribution of accelerations from each other
       * cell, one at a time */

#ifndef GRAV_ONLY
      /* apply uniform gravity */
      if (sim->use_uniform_grav) for (i=0;i<DIM;i++) curr->a[i] += sim->g[i];

      /* apply wall contact forces */
      if (sim->use_contact) {
         for (i=0;i<DIM;i++) {
            if (sim->bdry[i][0] == WALL) {
               // true_dist = curr->x[0][i]-top->min[i]-curr->rad;
               /* true_dist is now penetration */
               true_dist = top->min[i]+curr->rad-curr->x[0][i];
               // if (curr->x[0][i]-curr->rad < top->min[i]) {
               if (true_dist > 0.0) {
                  /* fprintf(stdout,"ball hits floor\n"); */
                  // true_dist = curr->x[0][i]-top->min[i]-curr->rad;
                  rel_vel = -curr->u[0][i];
                  /* now using k*x, not k*x*x */
                  // distsq = true_dist*true_dist;
                  /* sacc = ( sim->rk * distsq * distsq + sim->rc * rel_vel ) / curr->mass; */
                  /* now, the damping is merely a damping on acceleration, so we don't have
                   * to change the coefficient when mass changes! */
                  /* sacc = ( sim->rk * distsq + sim->rc * rel_vel ) / curr->mass; */
                  // sacc = sim->rk * distsq / curr->mass + sim->rc * rel_vel;
                  if (curr->parent)
                     kacc = curr->parent->k * true_dist / curr->mass;
                  else
                     kacc = sim->rk * true_dist / curr->mass;
                  // sacc = sim->rk * true_dist / curr->mass + sim->rc * rel_vel;
                  sacc = kacc + sim->rc*rel_vel;
                  /* fprintf(stdout,"  forces are %g from spring, %g from damping\n",sim->rk*distsq*distsq,sim->rc*rel_vel); */
                  /* fprintf(stdout,"  making acceleration %g\n",sacc); */
                  /* fflush(stdout); */
                  /* exit(0); */
                  curr->a[i] += sacc;

                  // now, for sliding friction force
                  // two valid directions
/*
                  if (DIM == 3) {
                     xdir = (i+1)%3;
                     ydir = (i+2)%3;
                     // normalize the velocities (to allow direction finding in 2D)
                     distsq = sqrt(curr->u[0][xdir]*curr->u[0][xdir] + curr->u[0][ydir]*curr->u[0][ydir]);
                     if (distsq > 1.0e-4) {
                        // acceleration is a fraction of the normal acceleration (so only occurs when penetration occurs)
                        // start mu at 0.1
                        curr->a[xdir] -= 0.4*sacc*curr->u[0][xdir]/distsq;
                        curr->a[ydir] -= 0.4*sacc*curr->u[0][ydir]/distsq;
                     }
                  } else if (DIM == 2) {
                     xdir = (i+1)%2;
                     // normalize the velocities (to allow direction finding in 2D)
                     distsq = curr->u[0][xdir];
                     if (distsq > 1.0e-4) {
                        // acceleration is a fraction of the normal acceleration (so only occurs when penetration occurs)
                        // start mu at 0.1
                        curr->a[xdir] -= 0.4*sacc;
                     }
                  } else {
*/
                     // fprintf(stderr,"ERROR (find_acc..): not programmed\n");
                     // exit(0);
                     // normalize the velocities (to allow direction finding in 2D)
                     distsq = 0.;
                     for (d=0;d<DIM;d++) {
                        if (d != i) {
                           distsq += curr->u[0][d]*curr->u[0][d];
                        }
                     }
                     distsq = sqrt(distsq);
                     if (distsq > 1.0e-4) {
                        for (d=0;d<DIM;d++) {
                           if (d != i) {
                              curr->a[d] -= 0.4*sacc*curr->u[0][d]/distsq;
                           }
                        }
                     }
//                }
               }
            }
            // and the opposite wall
            if (sim->bdry[i][1] == WALL) {
               true_dist = curr->x[0][i] + curr->rad - top->max[i];
               // if (curr->x[0][i]+curr->rad > top->max[i]) {
               if (true_dist > 0.0) {
                  // true_dist = top->max[i]-curr->x[0][i]-curr->rad;
                  rel_vel = curr->u[0][i];
                  // distsq = true_dist*true_dist;
                  /* sacc = ( sim->rk * distsq * distsq + sim->rc * rel_vel ) / curr->mass; */
                  /* sacc = ( sim->rk * distsq + sim->rc * rel_vel ) / curr->mass; */
                  // sacc = sim->rk * distsq /curr->mass + sim->rc * rel_vel;
                  if (curr->parent)
                     kacc = curr->parent->k * true_dist / curr->mass;
                  else
                     kacc = sim->rk * true_dist / curr->mass;
                  // sacc = sim->rk * true_dist / curr->mass + sim->rc * rel_vel;
                  sacc = kacc + sim->rc * rel_vel;
                  curr->a[i] -= sacc;

                  // now, for sliding friction force
                  // two valid directions
                  distsq = 0.;
                  for (d=0;d<DIM;d++) {
                     if (d != i) {
                        distsq += curr->u[0][d]*curr->u[0][d];
                     }
                  }
                  distsq = sqrt(distsq);
                  if (distsq > 1.0e-4) {
                     for (d=0;d<DIM;d++) {
                        if (d != i) {
                           curr->a[d] -= 0.4*sacc*curr->u[0][d]/distsq;
                        }
                     }
                  }
/*
                  fprintf(stderr,"ERROR (find_acc..): not programmed\n");
                  exit(0);
                  xdir = (i+1)%3;
                  ydir = (i+2)%3;
                  // normalize the velocities (to allow direction finding in 2D)
                  distsq = sqrt(curr->u[0][xdir]*curr->u[0][xdir] + curr->u[0][ydir]*curr->u[0][ydir]);
                  if (distsq > 1.0e-4) {
                     // acceleration is a fraction of the normal acceleration (so only occurs when penetration occurs)
                     // start mu at 0.1
                     curr->a[xdir] -= 0.4*sacc*curr->u[0][xdir]/distsq;
                     curr->a[ydir] -= 0.4*sacc*curr->u[0][ydir]/distsq;
                  }
*/
               }
            }
         }
      }
#endif

      /* apply inter-particle gravitational and contact forces */
      if (sim->use_self_grav) {
         // compute new ga[3]
         if (update_grav) {
            for (i=0;i<DIM;i++) curr->ga[i] = 0.0;
            find_grav_acc_on_this_part(sim,top,curr);
            if (sim->compare_to_direct) {
               for (i=0;i<DIM;i++) save[i] = curr->ga[i];
               for (i=0;i<DIM;i++) curr->ga[i] = 0.0;
               find_direct_grav_acc_on_this_part(sim,top,curr);
               // error by count
               for (i=0;i<DIM;i++) sim->sserr += pow(curr->ga[i]-save[i],2);
               for (i=0;i<DIM;i++) sim->ssval += save[i]*save[i];
               // error by mass
               for (i=0;i<DIM;i++) sim->smerr += curr->mass*pow(curr->ga[i]-save[i],2);
               for (i=0;i<DIM;i++) sim->smval += curr->mass*save[i]*save[i];
               sim->scnt += DIM;
            }
         }
#ifndef GRAV_ONLY
         // if grav_only, ga is the only acceleration
         for (i=0;i<DIM;i++) curr->a[i] += curr->ga[i];
#endif
         // fprintf(stdout,"acc after grav is %g %g %g\n",curr->a[0],curr->a[1],curr->a[2]); fflush(stdout);
      }

#ifndef GRAV_ONLY
      if (sim->use_contact) {
         find_cont_acc_on_this_part(sim,top,curr);
         // fprintf(stdout,"acc after contact is %g %g %g\n",curr->a[0],curr->a[1],curr->a[2]); fflush(stdout);
      }

      /* apply strand forces (neighbor stretching and bending) */
      if (sim->num_strands > 0)
         find_strand_acc_on_this_part(sim,curr);
#endif

      /* if (sim->step > 473) {
         fprintf(stdout,"  acc after find_acc_on_this_part is %g %g %g\n",curr->a[0],curr->a[1],curr->a[2]);
         fflush;
      } */

      // fprintf(stdout,"  vel at pt %g %g %g is %g %g %g\n",curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->u[1][0],curr->u[1][1],curr->u[1][2]);

      } // end of pragma omp task

      curr = curr->next;
   }
   } // end of pragma omp single nowait
   } // end of pragma omp parallel

   // free(acc);

   return(0);
}


#ifndef GRAV_ONLY
/*
 *  Find the acceleration on point from this cell due to strand forces only
 */
find_strand_acc_on_this_part(sim_ptr sim,particle_ptr this){

   int i;
   // FLOAT bending_stiffness = 1.657;
   FLOAT bending_stiffness = this->parent->b;
   FLOAT distsq,dist,sacc,width,rel_vel[DIM];
   FLOAT fd[DIM],bd[DIM],fdist,bdist,dotp;
   FLOAT theta,moment,m[DIM],mlen,f1[DIM],f2[DIM],f0[DIM];
   particle_ptr curr;

   // here's where we connect the strand together

   // the backwards segment
   if (this->s_last) {
      curr = this->s_last;
      for (i=0;i<DIM;i++) bd[i] = curr->x[0][i] - this->x[0][i];
      distsq = 0.;
      for (i=0;i<DIM;i++) distsq += bd[i]*bd[i];
      dist = curr->s_rad+this->s_rad;
      bdist = sqrt(distsq);	// distance of the centers
      /* fprintf(stdout,"true dist is %g rad are %g %g\n",bdist,curr->rad,this->rad); */
      /* fprintf(stdout,"contact, %d %d\n",curr->index,this->index); */
      /* square the penetration, this is not really the distance squared, tho */
      /* screw it, use k*x now, not k*x*x */
      // distsq = pow(bdist-dist,2);
      distsq = dist - bdist;
      /* normalize the direction vector */
      for (i=0;i<DIM;i++) bd[i] /= bdist;
      /* find the relative normal velocity */
      for (i=0;i<DIM;i++) rel_vel[i] = (curr->u[0][i] - this->u[0][i]) * bd[i];
      /* compute the spring force */
      // sacc = -1.0 * sim->rk * distsq;
      sacc = -1.0 * this->parent->k * distsq;
      /* compute the total change in acceleration due to spring and damping forces */
      /* for (i=0;i<DIM;i++) this->a[i] += (sacc + sim->rc*rel_vel[i]) * bd[i] / this->mass; */
      /* now, the damping is merely a damping on acceleration, so we don't have
       * to change the coefficient when mass changes! */
      for (i=0;i<DIM;i++) this->a[i] += (sacc/this->mass + sim->rc*rel_vel[i]) * bd[i];
      /* fprintf(stdout,"contact dist %g, acc is %g %g %g\n",bdist,this->a[0],this->a[1],this->a[2]); */
   }
   // fprintf(stdout,"acc is %g %g %g\n",this->a[0],this->a[1],this->a[2]);

   // and the forward segment
   if (this->s_next) {
      curr = this->s_next;
      for (i=0;i<DIM;i++) fd[i] = curr->x[0][i] - this->x[0][i];
      distsq = 0.;
      for (i=0;i<DIM;i++) distsq += fd[i]*fd[i];
      dist = curr->s_rad+this->s_rad;
      fdist = sqrt(distsq);	// distance of the centers
      distsq = dist - fdist;
      for (i=0;i<DIM;i++) fd[i] /= fdist;
      for (i=0;i<DIM;i++) rel_vel[i] = (curr->u[0][i] - this->u[0][i]) * fd[i];
      // sacc = -1.0 * sim->rk * distsq;
      sacc = -1.0 * this->parent->k * distsq;
      for (i=0;i<DIM;i++) this->a[i] += (sacc/this->mass + sim->rc*rel_vel[i]) * fd[i];
   }
   // fprintf(stdout,"  acc is %g %g %g\n",this->a[0],this->a[1],this->a[2]);

   // now, do the restoring moment - but only for the middle particles
   if (this->s_next && this->s_last) {
      for (i=0;i<DIM;i++) bd[i] = this->s_last->x[0][i] - this->x[0][i];
      for (i=0;i<DIM;i++) fd[i] = this->s_next->x[0][i] - this->x[0][i];
      // fprintf(stdout,"bdist,fdist %g %g\n",bdist,fdist);
      // theta is in radians, as it should be
      dotp = 0.;
      for (i=0;i<DIM;i++) dotp += bd[i]*fd[i];
      // theta = acos(-1.0e-12 -(bd[0]*fd[0]+bd[1]*fd[1]+bd[2]*fd[2])/(bdist*fdist));
      theta = acos(-1.0e-12 -(dotp)/(bdist*fdist));
      if (isnan(theta)) {
         // fprintf(stdout,"ERROR: theta is nan, fixing\n");
         theta = 0.0;
      }
      moment = bending_stiffness * (theta/2.0) / (bdist+fdist);
      // moment = bending_stiffness * sin(theta/2) / (bdist+fdist);
      // moment = 1.0e+2 * sin(theta/2);
      // fprintf(stdout,"theta is %g\n",theta);
      // fprintf(stdout,"moment is %g\n",moment);
      // fprintf(stdout,"force is %g\n",moment*bdist);

      // this constant needs to be 0.01 of the constant above
      // if (moment > 0.01*bending_stiffness) {
      if (theta > 0.001) {

         // fprintf(stdout,"theta %g, moment %g, force %g\n",theta,moment,moment*bdist);

         if (DIM == 3) {
            // find vector moment
            m[0] = fd[1]*bd[2] - fd[2]*bd[1];
            m[1] = fd[2]*bd[0] - fd[0]*bd[2];
            m[2] = fd[0]*bd[1] - fd[1]*bd[0];
            mlen = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
            for (i=0;i<3;i++) m[i] *= moment/mlen;
            // fprintf(stdout,"  m is %g %g %g\n",m[0],m[1],m[2]);

            // now, find the forces
            f1[0] = m[1]*bd[2] - m[2]*bd[1];
            f1[1] = m[2]*bd[0] - m[0]*bd[2];
            f1[2] = m[0]*bd[1] - m[1]*bd[0];
            f2[0] = -m[1]*fd[2] + m[2]*fd[1];
            f2[1] = -m[2]*fd[0] + m[0]*fd[2];
            f2[2] = -m[0]*fd[1] + m[1]*fd[0];
            for (i=0;i<DIM;i++) f0[i] = -f1[i]-f2[i];
            // fprintf(stdout,"  f0 is %g %g %g\n",f0[0],f0[1],f0[2]);

            // lastly, apply them as accelerations
            // acceleration of "this" particle
            for (i=0;i<DIM;i++) this->a[i] += f0[i]/this->mass;
            // acceleration of neighbor particles
            for (i=0;i<DIM;i++) this->s_next->a[i] += f2[i]/this->s_next->mass;
            for (i=0;i<DIM;i++) this->s_last->a[i] += f1[i]/this->s_last->mass;

            // fprintf(stdout,"  acc is %g %g %g\n",this->a[0],this->a[1],this->a[2]);
            if (isnan(this->a[0])) {
               fprintf(stderr,"\nDied because m is %g %g %g\n",m[0],m[1],m[2]);
               fprintf(stderr,"  fd is %g %g %g  bd is %g %g %g\n",fd[0],fd[1],fd[2],bd[0],bd[1],bd[2]);
               fprintf(stdout,"  theta %g (acos %g), moment %g, force %g\n",theta,(bd[0]*fd[0]+bd[1]*fd[1]+bd[2]*fd[2])/(bdist*fdist),moment,moment*bdist);
               exit(0);
            }
         } else if (DIM == 2) {
            // find vector moment
            mlen = fd[0]*bd[1] - fd[1]*bd[0];
            mlen = moment/mlen;
            // fprintf(stdout,"  m is %g\n",mlen);

            // now, find the forces
            f1[0] = -mlen*bd[1];
            f1[1] = mlen*bd[0];
            f2[0] = mlen*fd[1];
            f2[1] = -mlen*fd[0];
            for (i=0;i<DIM;i++) f0[i] = -f1[i]-f2[i];
            // fprintf(stdout,"  f0 is %g %g\n",f0[0],f0[1]);

            // lastly, apply them as accelerations
            // acceleration of "this" particle
            for (i=0;i<DIM;i++) this->a[i] += f0[i]/this->mass;
            // acceleration of neighbor particles
            for (i=0;i<DIM;i++) this->s_next->a[i] += f2[i]/this->s_next->mass;
            for (i=0;i<DIM;i++) this->s_last->a[i] += f1[i]/this->s_last->mass;

            // fprintf(stdout,"  acc is %g %g\n",this->a[0],this->a[1]);
            if (isnan(this->a[0])) {
               fprintf(stderr,"\nDied because m is %g\n",mlen);
               fprintf(stderr,"  fd is %g %g  bd is %g %g\n",fd[0],fd[1],bd[0],bd[1]);
               fprintf(stdout,"  theta %g (acos %g), moment %g, force %g\n",theta,(bd[0]*fd[0]+bd[1]*fd[1])/(bdist*fdist),moment,moment*bdist);
               exit(0);
            }
         } else {
            // dimensions 1,4+ not supported
            fprintf(stderr,"ERROR (find_strand_acc..): strands of dim 1 or >3 unsupported\n");
            exit(0);
         }
      }
   }

   // ah, heck, it had to have gone well
   return(0);
}
#endif


/*
 *  Find the acceleration on point from this cell; if the point is a
 *  particle, include the pointer.
 *
 *  HACK!!! This is the 3-dimensional form of Biot-Savart! Must make
 *  correctly multi-dimensional!
 */
int find_grav_acc_on_this_part(sim_ptr sim,cell_ptr curr_cell,particle_ptr this){

   int i;
   FLOAT distsq,dist,sacc,d[DIM],width;
   particle_ptr curr;

   /* under what conditions can we lump particles, or do a whole cell at a time? */

   /* first, and here is the crucial block that should speed this code up
    * immensely, if the cell is far enough away and small enough, count
    * the cell as ONE particle for purposes of computing the acceleration */
   for (i=0;i<DIM;i++) d[i] = curr_cell->cm[i] - this->x[0][i];
   distsq = 0.0;
   for (i=0;i<DIM;i++) distsq += d[i]*d[i];
   // the old way: pure gravitation by box
   // dist = sqrt(distsq);
   // the new way: add delta to box distance - MUCH better
   dist = sqrt(distsq) + sim->delta;
   distsq = dist*dist;

   // this is the spatial width
   width = curr_cell->max[0] - curr_cell->min[0];

   // this is the acceleration check right here:
   if (dist > width*sim->theta) {
   // if (dist > width*sim->theta && dist > 0.866*width+10.0*sim->delta) {
      sacc = sim->G * curr_cell->mass / distsq / dist;
      for (i=0;i<DIM;i++) this->ga[i] += sacc*d[i];
      /* and DON'T check the subcells, just return here */
      // fprintf(stdout,"  cell at level %d has dist %g, width %g\n",curr_cell->level,dist,width);
      return(0);
   }

   /* if this cell has subcells, run the subroutine on those */
   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         find_grav_acc_on_this_part(sim,curr_cell->s[i],this);

   /* otherwise, find the acceleration on "this" particle from all of the particles
    * in the current cell */
   } else {
      curr = curr_cell->first;
      while (curr) {
         // modify the "if" for spaghetti sim
#ifdef GRAV_ONLY
         if (curr != this) {
#else
         if (curr != this && curr != this->s_next && curr != this->s_last) {
#endif
            for (i=0;i<DIM;i++) d[i] = curr->x[0][i] - this->x[0][i];
            distsq = 0.0;
            for (i=0;i<DIM;i++) distsq += d[i]*d[i];
            dist = sqrt(distsq) + sim->delta;
            distsq = dist*dist;
            sacc = sim->G * curr->mass / distsq / dist;
            for (i=0;i<DIM;i++) this->ga[i] += sacc*d[i];
            // fprintf(stdout,"dist %g, grav acc is %g %g %g\n",dist,this->ga[0],this->ga[1],this->ga[2]);
         }
         curr = curr->next;
      }
   }

   // ah, heck, it had to have gone well
   return(0);
}


/*
 *  Find the acceleration on point from this cell using direct summation
 */
int find_direct_grav_acc_on_this_part(sim_ptr sim,cell_ptr curr_cell,particle_ptr this){

   int i;
   FLOAT distsq,dist,sacc,d[DIM];
   particle_ptr curr;

   /* if this cell has subcells, run the subroutine on those */
   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         find_direct_grav_acc_on_this_part(sim,curr_cell->s[i],this);

   /* otherwise, find the acceleration on "this" particle from all of the particles
    * in the current cell */
   } else {
      curr = curr_cell->first;
      while (curr) {
         // modify the "if" for spaghetti sim
#ifdef GRAV_ONLY
         if (curr != this) {
#else
         if (curr != this && curr != this->s_next && curr != this->s_last) {
#endif
            for (i=0;i<DIM;i++) d[i] = curr->x[0][i] - this->x[0][i];
            distsq = 0.0;
            for (i=0;i<DIM;i++) distsq += d[i]*d[i];
            dist = sqrt(distsq) + sim->delta;
            distsq = dist*dist;
            sacc = sim->G * curr->mass / distsq / dist;
            for (i=0;i<DIM;i++) this->ga[i] += sacc*d[i];
            // fprintf(stdout,"dist %g, grav acc is %g %g %g\n",dist,this->ga[0],this->ga[1],this->ga[2]);
         }
         curr = curr->next;
      }
   }

   // ah, heck, it had to have gone well
   return(0);
}


#ifndef GRAV_ONLY
/*
 *  Find the acceleration on point from this cell; if the point is a
 *  particle, include the pointer.
 */
int find_cont_acc_on_this_part(sim_ptr sim,cell_ptr curr_cell,particle_ptr this){

   int i;
   FLOAT distsq,dist,sacc,d[DIM],width,true_dist,rel_vel[DIM];
   particle_ptr curr;

   // if (curr_cell->level==0) fprintf(stderr,"particle %d at %g %g %g\n",this->index,this->x[0][0],this->x[0][1],this->x[0][2]); fflush(stderr);
   // fprintf(stderr,"testing cell level %d (%g %g %g)\n",curr_cell->level,curr_cell->min[0],curr_cell->min[1],curr_cell->min[2]); fflush(stderr);

   /* if the subcell's boundaries are such that no particles *can*
    * contact the current particle, then skip it completely */
   dist = 2.0*this->rad;
   for (i=0;i<DIM;i++) {
      if (this->x[0][i] < curr_cell->min[i] - dist) return(0);
      if (this->x[0][i] > curr_cell->max[i] + dist) return(0);
   }
   // fprintf(stderr,"  it could be here\n"); fflush(stderr);

   /* if this cell has subcells, run the subroutine on those */
   if (curr_cell->has_subcells) {
      // fprintf(stderr,"  we have subcells\n"); fflush(stderr);
      for (i=0;i<NCHILD;i++)
         find_cont_acc_on_this_part(sim,curr_cell->s[i],this);

   /* otherwise, find the acceleration on "this" particle from all of the particles
    * in the current cell */
   } else {
      // fprintf(stderr,"  we have no subcells, trying %d particles\n",curr_cell->num); fflush(stderr);
      curr = curr_cell->first;
      while (curr) {
         // fprintf(stderr,"    compare to particle %d\n",curr->index); fflush(stderr);
         // modify the "if" for spaghetti sim
         if (curr != this && curr != this->s_next && curr != this->s_last) {
            for (i=0;i<DIM;i++) d[i] = curr->x[0][i] - this->x[0][i];
            dist = curr->rad+this->rad;
            if (fabs(d[0]) < dist) {
             // if (fabs(d[1]) < dist) {
              // if (fabs(d[2]) < dist) {
               distsq = 0.0;
               for (i=0;i<DIM;i++) distsq += d[i]*d[i];
               true_dist = sqrt(distsq);
               // fprintf(stdout,"true dist is %g rad are %g %g\n",true_dist,curr->rad,this->rad);
               if (true_dist < dist) {
                  /* fprintf(stdout,"contact, %d %d\n",curr->index,this->index); */
                  /* square the penetration, this is not really the distance squared, tho */
                  /* screw it, use k*x now, not k*x*x */
                  // distsq = pow(true_dist-dist,2);
                  // distsq now represents the penetration distance
                  distsq = dist - true_dist;
                  /* normalize the direction vector */
                  for (i=0;i<DIM;i++) d[i] /= true_dist;
                  /* find the relative normal velocity */
                  for (i=0;i<DIM;i++) rel_vel[i] = (curr->u[0][i] - this->u[0][i]) * d[i];
                  /* compute the spring force */
                  sacc = -1.0 * sim->rk * distsq;
                  /* compute the total change in acceleration due to spring and damping forces */
                  /* for (i=0;i<DIM;i++) this->a[i] += (sacc + sim->rc*rel_vel[i]) * d[i] / this->mass; */
                  /* now, the damping is merely a damping on acceleration, so we don't have
                   * to change the coefficient when mass changes! */
                  for (i=0;i<DIM;i++) this->a[i] += (sacc/this->mass + sim->rc*rel_vel[i]) * d[i];
                  // fprintf(stdout,"contact dist %g, acc is %g %g %g\n",true_dist,this->a[0],this->a[1],this->a[2]);
               }
              // }
             // }
            }
         }
         curr = curr->next;
      }
   }

   // ah, heck, it had to have gone well
   return(0);
}
#endif


/*
 *  Find the acceleration on this_cell from all other cells
 */
int find_acc_on_this_cell(sim_ptr sim,cell_ptr curr_cell,cell_ptr this_cell,FLOAT* acc){

   int i;
   FLOAT distsq,dist,sacc,d[3];

   /* loop thru cells, ignoring curr_cell */
   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         find_acc_on_this_cell(sim,curr_cell->s[i],this_cell,acc);
   } else {
      if (curr_cell != this_cell) {
         for (i=0;i<DIM;i++) d[i] = curr_cell->cm[i] - this_cell->cm[i];
         distsq = 0.0;
         for (i=0;i<DIM;i++) distsq += d[i]*d[i];
         dist = sqrt(distsq);
         if (sim->use_self_grav) {
            sacc = sim->G * curr_cell->mass / dist / dist / dist;
            for (i=0;i<DIM;i++) acc[i] += sacc*d[i];
         }
      }
   }

   return(0);
}


/*
 * find maximum velocity of particles in cell
 */
FLOAT find_vmax(cell_ptr top,cell_ptr curr_cell,FLOAT max_vel){

   int i;
   FLOAT this_vel;
   particle_ptr curr = curr_cell->first;

   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         max_vel = find_vmax(top,curr_cell->s[i],max_vel);
   } else {
      while (curr) {

         this_vel = 0.0;
         for (i=0;i<DIM;i++) this_vel += curr->u[0][i]*curr->u[0][i];
         this_vel = sqrt(this_vel);

         if (this_vel > max_vel) max_vel = this_vel;

         /* jump to the next one */
         curr = curr->next;

      }  /* end while (curr->next) */
   }
   return(max_vel);
}

