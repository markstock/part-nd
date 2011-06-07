/*************************************************************
 *
 *  part-nd.c - Three-dimensional particle simulator
 *
 *  Copyright (C) 2000-2005  Mark J. Stock, mstock@umich.edu
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
 * Revision information
 *
 * 2000-12-27 MJS  v0.0  First compile, no octree
 * 2000-12-31 MJS  v0.1  With octrees
 * 2001-01-02 MJS  v0.2  With working octrees
 * 2001-01-09 MJS  v0.3  First bug-free working version
 * 2001-02-17 MJS  v0.4  Now called part3d, self-contact forces
 * 2001-08-29 MJS  v0.6  Damping fixed, density field added
 * 2002-05-09 MJS  v0.7  Now called spag3d, add strand forces
 * 2002-08-22 MJS  v0.8  Back to part3d, merge of both codes
 * 2002-09-03 MJS  v0.9  Finally GPL'd
 * 2003-05-31 MJS  v1.0  self-grav fix, accel cap, better output
 * 2005-02-01 MJS  v1.1  Now called grav-nd, arbitrary-dimensional,
 *                       gravity-only
 * 2005-02-09 MJS  v1.2  Now called part-nd, arbitrary-dimensional,
 *                       supports all other forces in arbitrary dims
 * 2005-02-23 MJS  v1.3  Now writes png natively, removed GIF support
 *
 *********************************************************** */


#include "structs.h"

int run_sim(fileprop_ptr,sim_ptr,cell_ptr);
int set_dt(sim_ptr,cell_ptr);
int advect_nodes(sim_ptr,cell_ptr,cell_ptr);
FLOAT calculate_total_energy(FLOAT,sim_ptr,cell_ptr,cell_ptr);
FLOAT calculate_momentum(sim_ptr,cell_ptr,cell_ptr,FLOAT*,FLOAT);


int main(int argc,char **argv) {

   FILEP file_props;				/* initialize struct for file props */
   fileprop_ptr file = &file_props;
   SIMP sim_props;				/* initialize struct for simulation props */
   sim_ptr sim = &sim_props;
   CELL top_cell;				/* initialize struct for octree */
   cell_ptr top = &top_cell;

   // set some important pointers
   file->out = &file->out_field;
   sim->ff2 = &sim->density_field2;
   sim->ff3 = &sim->density_field3;

   /* Parse command-line args */
   parse_args(argc,argv,file);

   /* Initialize the top-level cell and participating particles */
   // initialize_system(file,sim,top,ff);
   initialize_system(file,sim,top);

   /* Do the time integration */
   // run_sim(file,sim,top,ff);
   run_sim(file,sim,top);

   fprintf(stderr,"\nDone.\n");
   exit(0);
}



/*
 *  Run the time integration of the flowfield simulation
 */
// int run_sim(fileprop_ptr file,sim_ptr sim,cell_ptr top,field2_ptr ff) {
int run_sim(fileprop_ptr file,sim_ptr sim,cell_ptr top) {

   int i,retval,update_grav;
   FLOAT momentum[DIM],mass,energy;
   double dtime;
   unsigned long int stop_tics, start_tics;

   /* until we get output working better, temp filename */
   char outfile[80],outfile2[80],command[160];

   // write particle count and cell utilization
   if (0) {
      /* write_particle_count(top); */
      fprintf(stdout,"  %d particles\n",top->num);
      for (i=0; i<sim->max_levels; i++) sim->cell_count[i] = 0;
      count_cells(sim,top);
      fprintf(stdout,"  ");
      for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%d ",sim->cell_count[i]);
      fprintf(stdout,"cells\n  utilization ");
      for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%g ",sim->cell_count[i]/pow(8,i));
      fprintf(stdout,"\n");
   }

   // make a good tree
   top = build_new_tree(file,sim,top);

   /* iterate over time */
   for (sim->step=0;sim->step<sim->maxstep;sim->step++) {

      // This was after the two ifs
      find_all_cells_cm(sim,top);
      // was 1000 and 100
      fprintf(stdout,"Step %d at time %g, next print %g, with %d particles\n",sim->step,sim->time,sim->next_output_time,top->num);
      // if (sim->step%1000 == 0) fprintf(stdout,"\nStep %d at time %g, next print %g, with %d particles",sim->step,sim->time,sim->next_output_time,top->num);
      // if (sim->step%100 == 0) { fprintf(stdout,"."); fflush(stdout); }

      /* write particle count and cell utilization */
      if (0) {
         /* write_particle_count(top); */
         fprintf(stdout,"  %d particles\n",top->num);
         for (i=0; i<sim->max_levels; i++) sim->cell_count[i] = 0;
         count_cells(sim,top);
         fprintf(stdout,"  ");
         for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%d ",sim->cell_count[i]);
         fprintf(stdout,"cells\n  utilization ");
         for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%g ",sim->cell_count[i]/pow(8,i));
         fprintf(stdout,"\n");
      }

      /* check current time vs. next_output_time, write output if necessary */
      // if (fabs(sim->time - sim->next_output_time) < 1e-5) {
      // if (fabs(sim->time - sim->next_output_time) < sim->dt/2.0) {
      if ((sim->time - sim->next_output_time) > sim->dt/2.0) {
         fprintf(stdout,"  writing output step %d\n",sim->next_output_index);

         /* actually write the output here */
         write_output(file,sim,top);

         /* compute some other stats */
         if (0) {
            for (i=0;i<DIM;i++) momentum[i] = 0.0;
            mass = calculate_momentum(sim,top,top,momentum,0.0);
            fprintf(stdout,"    avg. momentum %g %g %g\n",momentum[0]/mass,momentum[1]/mass,momentum[2]/mass);
         }

         sim->next_output_time += sim->output_dt;
         sim->next_output_index++;
      }

      /* check current time vs. end_time, return if neccessary */
      /* return code 0 - successfully concluded */
      // if (sim->time > sim->end_time) return(0);
      if (sim->time+0.1*sim->dt > sim->end_time) {
         fprintf(stdout,"Out of time.\n");
         return(0);
      }
      if (top->num < 1) {
         fprintf(stdout,"Out of particles.\n");
         return(0);
      }

      /* Okay: contact and strand forces get updated at every time
       * step; gravity acts much more slowly, so when gravity is 
       * included in a simulation that has strands or contact, update
       * the gravitational forces only once per 100 steps
       */
#ifdef GRAV_ONLY
      update_grav = TRUE;
#else
      if (sim->num_strands > 0 || sim->use_contact) {
         if (sim->step%sim->contact_per_grav == 0) update_grav = TRUE;
         else update_grav = FALSE;
      } else {
         // if there's no contact or strands, update gravity every time
         update_grav = TRUE;
      }
#endif

      // one full forward advection step
      {

      /* find the velocities of the particles */
      if (sim->compare_to_direct) {
         sim->ssval = 0.0;
         sim->sserr = 0.0;
         sim->smval = 0.0;
         sim->smerr = 0.0;
         sim->scnt = 0;
      }
      start_tics = clock();
      find_new_vels(sim,top,top,update_grav);
      stop_tics = clock();
      dtime = ((double)stop_tics-(double)start_tics)/CLOCKS_PER_SEC;
         fprintf(stdout,"    took %g sec\n",dtime);
         fflush(stdout);
      if (sim->compare_to_direct) {
         sim->rmserror = sqrt(sim->sserr/sim->ssval);
         sim->rmsmerror = sqrt(sim->smerr/sim->smval);
         fprintf(stdout,"    at %g error %g by mass %g\n",sim->time,sim->rmserror,sim->rmsmerror);
         fflush(stdout);
      }

      /* find an appropriate dt */
      /* set_dt(sim,top); */

      /* advect the nodes, this doesn't move particles to their appropriate cells */
      retval = advect_nodes(sim,top,top);
      // fprintf(stdout,"    advect nodes threw %d\n",retval);

      /* grow strands, if neccessary */
#ifndef GRAV_ONLY
      if (sim->num_strands > 0) grow_strands(sim);
#endif

      /* replace the particles that are not in the correct cell */
      /* it doesn't recalculate the proper number of particles per cell */
      if (retval < 0)
         top = build_new_tree(file,sim,top);
      else
         replace_all_particles(sim,top,top);

      /* remove any cells with no particles, and combine cells with
       * too few particles, it also recalculates the proper number of
       * particles per cell */
      clean_up_all_cells(sim,top);

      }

      /* calculate total energy of the system */
      if (0) {
         energy = calculate_total_energy(0.0,sim,top,top);
         fprintf(stdout,"  total energy is %g\n",energy);
      }

      sim->time += sim->dt;

   /* jump back to head */
   }

   /* return code 1 - ran out of steps */
   fprintf(stdout,"Simulation completed.\n");
   return(1);

}


/*
 *  Calculate and print total momentum of the system
 */
FLOAT calculate_momentum(sim_ptr sim,cell_ptr top,cell_ptr curr_cell,FLOAT *mom,FLOAT mass){

   int i,j,k;
   particle_ptr curr;

   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         mass = calculate_momentum(sim,top,curr_cell->s[i],mom,mass);
   } else {
      curr = curr_cell->first;
      while (curr) {
         for (i=0;i<DIM;i++) mom[i] += curr->mass*curr->u[0][i];
         mass += curr->mass;
         curr = curr->next;
      }
   }

   return(mass);
}


/*
 *  Find the total energy of the system
 */
FLOAT calculate_total_energy(FLOAT energy,sim_ptr sim,cell_ptr top,cell_ptr curr_cell){

   int i,j,k;
   FLOAT velsq;
   particle_ptr curr;

   /* for each particle, see if it needs to be moved, if so,
    * move it back up to the top cell where it will reallocate
    * downward to an appropriate cell */

   if (curr_cell->has_subcells) {

      for (i=0;i<NCHILD;i++)
         energy = calculate_total_energy(energy,sim,top,curr_cell->s[i]);

   } else {

      curr = curr_cell->first;
      while (curr) {

         fprintf(stderr,"ERROR (calculate_total_e..): not coded\n");
         exit(0);
         /* add this particle's energy to the total */
         velsq = curr->u[0][0]*curr->u[0][0] + curr->u[0][1]*curr->u[0][1] + curr->u[0][2]*curr->u[0][2];
         energy += curr->mass*(sim->g[2]*(top->min[2]-curr->x[0][2]) + 0.5*velsq);
         // fprintf(stdout,"energy %g mass %g g %g floor %g elev %g vel %g\n",energy,curr->mass,sim->g[2],top->min[2],curr->x[0][2],sqrt(velsq));

         /* jump to the next one */
         curr = curr->next;
      }
   }

   return(energy);
}


/*
 *  Find an appropriate dt
 */
int set_dt(sim_ptr sim,cell_ptr top) {

   FLOAT dt1;

   /* find maximum dt to allow nodes to move 1/2 grid cell */
   sim->vmax = find_vmax(top,top,0.0);
   /* change this 0.001 number and see what affect it has on both max vel
      and average velocities! */
   dt1 = 0.001*top->max[0]/(1.1*sim->vmax);
   fprintf(stdout,"  vmax is %g gives dt1 of %g ",sim->vmax,dt1);

   /* set dt to either max_dt, or dt1, whichever is shorter */
   if (dt1 < sim->max_dt) {
      sim->dt = dt1;
   } else {
      sim->dt = sim->max_dt;
   }
   /* fprintf(stdout,"  sim->dt is %g\n",sim->dt); */

   /* Now, let the dt adjust itself so that time always
    * lands directly on an output time. Say, if time=0.91 and dt=0.1,
    * but printtime=1.0, reset dt to be 0.9. Note: give this block the
    * ability to stretch dt just a little. If time=0.899, then set dt 
    * to be 0.101, even though it is greater than dtmax. */

   if (sim->time + 1.1*sim->dt > sim->next_output_time) {
      sim->dt = sim->next_output_time - sim->time;
   }

   fprintf(stdout,"and dt used is %g\n",sim->dt);

   /* return code 0 - everything worked */
   return(0);
}


/*
 * advect all of the nodes
 */
int advect_nodes(sim_ptr sim,cell_ptr top,cell_ptr curr_cell) {

   int i;
   int retval = 0;
   int delete_it = FALSE;
   int cap_accel = FALSE;
   static int cap_cnt;
   FLOAT newu[DIM];
   FLOAT accel;
   static FLOAT accelcap,accelcapsq,max_accelsq;
   particle_ptr curr = curr_cell->first;
   particle_ptr last = curr;

   if (curr_cell == top) {
      cap_cnt = 0;
      max_accelsq = 0.0;
      if (cap_accel) {
         accelcap = 0.001/sim->dt;
         accelcapsq = accelcap*accelcap;
      }
   }

   if (curr_cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         retval = retval + advect_nodes(sim,top,curr_cell->s[i]);

   } else {

      while (curr) {

         // impose an acceleration cap
         if (cap_accel) {
            fprintf(stderr,"ERROR (advect_nodes): not coded\n");
            exit(0);
            accel = 0.;
#ifdef GRAV_ONLY
            for (i=0;i<DIM;i++) accel += curr->ga[i]*curr->ga[i];
#else
            for (i=0;i<DIM;i++) accel += curr->a[i]*curr->a[i];
#endif
            if (accel > max_accelsq) max_accelsq = accel;
            if (accel > accelcapsq) {
               cap_cnt++;
               accel = sqrt(accel);
#ifdef GRAV_ONLY
               for (i=0;i<DIM;i++) curr->ga[i] = accelcap*curr->ga[i]/accel;
#else
               for (i=0;i<DIM;i++) curr->a[i] = accelcap*curr->a[i]/accel;
#endif
            }
         }

         // old way, but I guess we're keeping it...
         if (!curr->stationary) {
#ifdef GRAV_ONLY
            for (i=0;i<DIM;i++) newu[i] = curr->u[0][i] + curr->ga[i] * sim->dt;
#else
            for (i=0;i<DIM;i++) newu[i] = curr->u[0][i] + curr->a[i] * sim->dt;
#endif
            for (i=0;i<DIM;i++) curr->x[0][i] += 0.5 * (curr->u[0][i] + newu[i]) * sim->dt;
            for (i=0;i<DIM;i++) curr->u[0][i] = newu[i];
         }

         // and reset the acceleration to zero
#ifndef GRAV_ONLY
         for (i=0;i<DIM;i++) curr->a[i] = 0.0;
#endif

         /* ===================== nah, do that later ========================================== */
         /* if the particle moves outside of the bounds of it's parent cell,
          * move it to this cell's parent, and check there. */
         /* ===================== nah, do that later ========================================== */
         /* if it moved out of bounds, remove it */
         delete_it = FALSE;
         // Do not delete any particles anymore! Just flag for tree rebuilding
         for (i=0;i<DIM;i++)
            if (curr->x[0][i] < top->min[i] && sim->bdry[i][0] == OPEN) {
               delete_it = TRUE;
               // fprintf(stdout,"  i, curr->x[0][i], top->min are %d %g %g\n",i,curr->x[0][i],top->min[i]);
               i+=3;
            }
         if (!delete_it) for (i=0;i<DIM;i++)
            if (curr->x[0][i] > top->max[i] && sim->bdry[i][1] == OPEN) {
               delete_it = TRUE;
               // fprintf(stdout,"  i, curr->x[0][i], top->max are %d %g %g\n",i,curr->x[0][i],top->max[i]);
               i+=3;
            }

         if (delete_it) {
            retval = -1;
         }

         /* jump to the next one */
         last = curr;
         curr = curr->next;

      }  /* end while (curr->next) */

      /*
      if (max_accel > accelcap) {
         fprintf(stdout,"max accel %g\n",max_accel);
         fflush(stdout);
         getchar();
      }
      */

   }

   /*
   if (curr_cell == top && cap_accel) {
      fprintf(stderr,"  %d of %d particles hit accel cap of %g, max %g\n\n",cap_cnt,top->num,accelcap,sqrt(max_accelsq));
   }
   */

   /* if all went well */
   return(retval);
}


#ifndef GRAV_ONLY
/*
 * grow all of the strands: force second particle to move, place
 * new particle if neccessary
 */
int grow_strands(sim_ptr sim) {

   int i,j;
   FLOAT dist[DIM],absdist;
   particle_ptr end;
   particle_ptr second;
   particle_ptr third;
   particle_ptr newpart;

   for (i=0; i<sim->num_strands; i++) {
      if (sim->strand[i].grows) {

         // first, save pointers to the first three particles
         end = sim->strand[i].start;
         second = end->s_next;                // must exist
         third = second->s_next;      // may be NULL

         // move the second particle, regardless of force required
         for (j=0;j<DIM;j++) second->x[0][j] += sim->strand[i].growdir[j] * sim->dt;

         // check to see if there's enough room to spawn a new particle
         for (j=0;j<DIM;j++) dist[j] = second->x[0][j] - end->x[0][j];

         fprintf(stderr,"ERROR (grow_strands): not coded\n");
         exit(0);
         absdist = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);

         if (absdist > 1.0) {

         }
      }
   }

   return(0);
}
#endif

