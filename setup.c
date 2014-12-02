/*************************************************************
 *
 *  setup.c - subroutines for setting up simulations
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
 *********************************************************** */


#include "structs.h"
#include <time.h>

int initialize_system(fileprop_ptr,sim_ptr,cell_ptr);
int add_particles_from_file(sim_ptr,cell_ptr,int,char);
int add_box_of_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT*);
int add_box_of_ni_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT*);
int add_sphere_of_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT);
int add_sphere_of_ni_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT);
int add_strand_of_particles(sim_ptr,cell_ptr,int,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT);
int is_this_space_open(cell_ptr,FLOAT*,FLOAT);
int is_this_space_open2(cell_ptr,FLOAT*,FLOAT);
//int add_galaxy(sim_ptr,cell_ptr,int,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT);
int add_galaxy (sim_ptr,cell_ptr,int,FLOAT*,FLOAT*,FLOAT,FLOAT*);
void pick_point_on_sphere (FLOAT*,FLOAT);
particle_ptr new_particle(int,FLOAT,FLOAT,FLOAT*,FLOAT*);
particle_ptr new_stationary_particle(int,FLOAT,FLOAT*);
int read_input_file(fileprop_ptr,sim_ptr,cell_ptr);
int parse_args(int,char**,fileprop_ptr);
int Usage(char[80],int);
FLOAT*** allocate_3d_array_F(int,int,int);
FLOAT** allocate_2d_array_F(int,int);
png_byte** allocate_2d_array_pb(int,int,int);


/*
 *  Initialize the system - make static and free element structures, etc.
 */
int initialize_system(fileprop_ptr file,sim_ptr sim,cell_ptr top) {

   int i,j,d,ret_val,num_parts;
   FLOAT loc[DIM],vel[DIM],xpos,ypos,zpos,xvel,yvel,zvel,distsq,mpp;
   FLOAT r,rcubed,orbvel,M;
   particle_ptr newpart;

   // -------------------------------------------------------------------------

   /* set default values for variables */
   // strcpy(file->out_fn_root,"output_");
   strcpy(file->out_fn_root,"");
   file->write_pgm = FALSE;
   // file->write_gif = FALSE;
   file->write_png = FALSE;
   file->write_dot = FALSE;
   file->write_rad = FALSE;
   file->write_part = FALSE;
   file->write_seg = FALSE;
   file->image_depth = 8;
   file->out_img_size = 256;

   sim->step = 0;
   sim->maxstep = 999999999;			/* make this sufficiently high */
   sim->start_time = 0.0;
   sim->end_time = 999999.0;			/* this too */
   sim->dt = 0.001;
   sim->max_dt = 0.001;
   sim->output_dt = 0.01;
   sim->max_levels = 10;
   sim->max_parts_in_cell = 40;

   sim->time = sim->start_time;
   sim->next_output_time = sim->start_time;
   sim->next_output_index = 0;

   sim->use_self_grav = TRUE;
   sim->G = 0.0001;
   sim->delta = 0.001;
   sim->theta = 2.5;
   sim->compare_to_direct = FALSE;
   sim->contact_per_grav = 1;

   sim->use_uniform_grav = FALSE;
   sim->g[0] = 0.0;
   sim->g[1] = 0.0;
   sim->g[2] = -9.8;

   sim->use_contact = FALSE;
   sim->rk = 1.0e+9;
   sim->rn = 4;
   sim->rc = 1.0e-6;
   sim->ak = 1.0;
   sim->an = 4;

   for (d=0;d<DIM;d++) {
      sim->bdry[d][0] = WALL;
      sim->bdry[d][1] = WALL;
   }

   sim->particle_cnt = 1000;
   sim->new_part_rad = 0.005;
   sim->new_part_mass = 1.0;

   sim->num_blocks = 0;
   sim->num_spheres = 0;
   sim->num_galaxy = 0;
   sim->num_indiv_parts = 0;
   sim->num_read_part = 0;
   sim->num_read_stat = 0;
   sim->num_strands = 0;
   sim->use_density_field = FALSE;
   sim->construct_dens_surface = FALSE;


   /* assign default values for the top-level cell */
   top->level = 0;
   top->has_subcells = FALSE;
   top->first = NULL;
   for (i=0;i<NCHILD;i++)
      top->s[i] = NULL;
   for (d=0;d<DIM;d++) {
      top->min[d] = 0.0;
      top->max[d] = 1.0;
   }
   top->num = 0;
   top->mass = 0.0;

   // -------------------------------------------------------------------------

   /* now, if the file listed on the command line exists, read it */
   /* ret_val is 0 if file is OK and read, 1 if not */
   ret_val = read_input_file(file,sim,top);

   if (ret_val == 1) {
      fprintf(stdout,"ERROR (setup): upon reading input file, please check syntax\n");
      exit(1);
   }

   // -------------------------------------------------------------------------

   // I think we only need these for the self-gravitation thing
   for (i=0;i<DIM+1;i++) top->mid[i] = (top->min[i]+top->max[i])/2.0;
   for (i=0;i<DIM;i++) top->cm[i] = (top->min[i]+top->max[i])/2.0;

   // -------------------------------------------------------------------------
   // now, place particles

   for (i=0; i<sim->num_read_part; i++)
      add_particles_from_file(sim,top,i,FALSE);

   for (i=0; i<sim->num_read_stat; i++)
      add_particles_from_file(sim,top,i,TRUE);

#ifdef ALLDONE
#ifndef GRAV_ONLY
   for (i=0; i<sim->num_strands; i++) {
      // fprintf(stdout,"sending %d %lf %lf %lf\n",i,sim->block[i][0],sim->block[i][1],sim->block[i][2]);
      add_strand_of_particles(sim,top,i,sim->block[i][0],
         sim->block[i][1],sim->block[i][2],sim->block[i][3],
         sim->block[i][4],sim->block[i][5],sim->block[i][6],
         sim->block[i][7],sim->block[i][8],sim->block[i][9],
         sim->block[i][10],sim->block[i][11]);
   }
#endif
#endif

   // then, place specific particles
   for (i=0; i<sim->num_indiv_parts; i++) {
      for (d=0;d<DIM;d++) loc[d] = sim->indiv_part[i][d+1];
      for (d=0;d<DIM;d++) vel[d] = sim->indiv_part[i][d+DIM+1];
      newpart = new_particle(i,sim->indiv_part[i][0],sim->new_part_rad,loc,vel);
      add_particle_to_cell(sim,newpart,top);
   }

   // lastly, place particles to fill cubical or spherical volumes
   for (i=0; i<sim->num_blocks; i++)
      add_box_of_particles(sim,top,(int)(sim->block[i][0]),
         &sim->block[i][1],&sim->block[i][1+DIM]);

   // add a box of non-intersecting particles!
   for (i=0; i<sim->num_niblocks; i++)
      add_box_of_ni_particles(sim,top,(int)(sim->niblock[i][0]),
         &sim->niblock[i][1],&sim->niblock[i][1+DIM]);

   // add a sphere of randomly-placed particles
   for (i=0; i<sim->num_spheres; i++)
      add_sphere_of_particles(sim,top,(int)(sim->sphere[i][0]),
         &sim->sphere[i][1],sim->sphere[i][1+DIM]);

   // add a sphere of non-intersecting particles
   for (i=0; i<sim->num_nispheres; i++)
      add_sphere_of_ni_particles(sim,top,(int)(sim->nisphere[i][0]),
         &sim->nisphere[i][1],sim->nisphere[i][1+DIM]);

   // add a galaxy of randomly-placed particles
   for (i=0; i<sim->num_galaxy; i++)
      add_galaxy(sim,top,(int)(sim->galaxy[i][0]),
         &sim->galaxy[i][1],&sim->galaxy[i][1],1.0,
         &sim->galaxy[i][1+DIM]);

   // -------------------------------------------------------------------------

   if (sim->use_density_field) {
      // set the cell sizes
      for (i=0; i<2; i++) sim->ff2->n[i] = file->out_img_size;
      for (i=0; i<2; i++) sim->ff2->d[i] = (top->max[i]-top->min[i])/sim->ff2->n[i];
      for (i=0; i<2; i++) sim->ff2->o[i] = top->min[i];

      // allocate memory for the array
      // ff->rho = allocate_3d_array_F(ff->n[0],ff->n[1],ff->n[2]);
      // sim->ff2->rho = allocate_2d_array_F(sim->ff2->n[0],sim->ff2->n[2]);
      sim->ff2->rho = allocate_2d_array_F(sim->ff2->n[0],sim->ff2->n[1]);

      // zero the array
      for (i=0;i<sim->ff2->n[0];i++)
         for (j=0;j<sim->ff2->n[1];j++)
            sim->ff2->rho[i][j] = 0.0;
            // for (k=0;k<sim->ff2->n[2];k++)
               // ff->rho[i][j][k] = 0.0;
   }

   // fprintf(stdout,"\nStarting with %d particles",top->num);

   // and the output file (just dots)
   for (i=0; i<2; i++) file->out->d[i] = (top->max[i]-top->min[i])/file->out->n[i];
   if (file->write_dot || sim->use_density_field)
      file->out->rho = allocate_2d_array_F(file->out_img_size,file->out_img_size);

   // malloc the space for png writing, if necessary
   if (file->write_png)
      file->image = allocate_2d_array_pb(file->out_img_size,file->out_img_size,file->image_depth);

   /* return 0, everything OK */
   return(0);
}


/*
 *  Read particle locations from a file, include them in the simulation
 *  Do this for both static and non-static particles
 */
int add_particles_from_file(sim_ptr sim,cell_ptr top,int i,char is_stat) {

   int d;
   int cnt = 0;
   FLOAT loc[DIM],rad,vel[DIM];
   char token[7][16];
   char twochar[2];
   char sbuf[512];
   particle_ptr newpart;
   FILE *infile;

   /* open file for reading */
   if (is_stat) {
      infile = fopen(sim->read_stat[i],"r");
      if (infile==NULL) {
         fprintf(stderr,"Could not open particle file %s\n",sim->read_stat[i]);
         fflush(stderr);
         return(1);
      }
      fprintf(stdout,"Opening file %s\n",sim->read_stat[i]);
      fflush(stdout);
   } else {
      infile = fopen(sim->read_part[i],"r");
      if (infile==NULL) {
         fprintf(stderr,"Could not open particle file %s\n",sim->read_part[i]);
         fflush(stderr);
         return(1);
      }
      fprintf(stdout,"Opening file %s\n",sim->read_part[i]);
      fflush(stdout);
   }

   /* read a line from the input file */
   while (fscanf(infile,"%[^\n]",sbuf) != EOF) {

      /* fprintf(stdout,"%s\n",sbuf); */

      /* grab the line */
      sscanf(sbuf,"%s %s %s %s %s %s %s",token[0],token[1],token[2],token[3],token[4],token[5],token[6]);
      /* fprintf(stdout,"   first token is %s\n",token[0]); */

      if (token[0][0] == '#') {
         /* read a comment line, or some other descriptor */
         fscanf(infile,"%[\n]",twochar);    /* read up to newline */
         /* fprintf(stdout,"%s\n",sbuf);        /* write comment */

      } else {
         /* read the particle parameters */
         /* fprintf(stdout,"   second token is %s\n",token[1]); */

         for (d=0;d<DIM;d++) loc[d] = atof(token[d]);
         rad = atof(token[DIM]);
         for (d=0;d<DIM;d++) vel[d] = atof(token[d+DIM+1]);

         if (is_stat) {
            for (d=0;d<DIM;d++) loc[d] += sim->stat_xform[i][d];
         } else {
            for (d=0;d<DIM;d++) loc[d] += sim->part_xform[i][d];
            for (d=0;d<DIM;d++) vel[d] += sim->part_xform[i][d+DIM];
         }

         newpart = new_particle(cnt,sim->new_part_mass,rad,loc,vel);
         add_particle_to_cell(sim,newpart,top);

         cnt++;
         if (cnt/10000 != (cnt+1)/10000) {
            fprintf(stdout,".");
            fflush(stdout);
         }

         fscanf(infile,"%[\n]",twochar);    /* read newline */
      }
   }

   fprintf(stdout,"\n");
   fprintf(stdout,"Placed %d particles\n",cnt);
   fflush(stdout);

   fclose(infile);

   /* ret_val is 0 is file is OK and read, 1 if not */
   return(0);
}


/*
 *  Create a cube/rectangle of particles (can be intersecting)
 */
int add_box_of_particles(sim_ptr sim,cell_ptr top,int cnt,FLOAT* start,FLOAT* finish) {

   int d,i,keep_trying,ntried,isitfree,quit_altogether;
   FLOAT loc[DIM],size[DIM];
   FLOAT vel[DIM];
   particle_ptr newpart;
   unsigned long int tics,last_tics;
   FLOAT temp;

   quit_altogether = FALSE;
   for (d=0;d<DIM;d++) size[d] = finish[d] - start[d];
   for (d=0;d<DIM;d++) vel[d] = 0.;

   /* create the particles */
   last_tics = clock();
   for (i=0;i<cnt;i++) {

      for (d=0;d<DIM;d++) loc[d] = start[d] + size[d]*(rand()/(RAND_MAX+1.0));
      /* fprintf(stdout,"  testing pt %g %g %g\n",xpos,ypos,zpos); */

      /* this location is good, place the particle */
      newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,loc,vel);
      add_particle_to_cell(sim,newpart,top);
      /* fprintf(stdout,"  placed pt there\n"); */

      if (i/1000 != (i+1)/1000) {
         fprintf(stdout,".");
         fflush(stdout);
      }
   }
   fprintf(stdout,"\n");
   tics = clock();
   temp = ((FLOAT)tics-(FLOAT)last_tics)/CLOCKS_PER_SEC;

   if (i>cnt) fprintf(stdout,"Placed %d of %d particles\n",i-cnt-1,cnt);
   else fprintf(stdout,"Placed %d of %d particles\n",i,cnt);

   fprintf(stdout,"Took %g seconds of cpu time for it\n",temp);

   /* return 0 is all went well */
   return(0);
}


/*
 *  Create a cube/rectangle of non-intersecting particles
 */
int add_box_of_ni_particles(sim_ptr sim,cell_ptr top,int cnt,FLOAT* start,FLOAT* finish) {

   int i,d,keep_trying,ntried,isitfree,quit_altogether;
   FLOAT loc[DIM],size[DIM];
   FLOAT vel[DIM];
   particle_ptr newpart;
   unsigned long int tics,last_tics;
   FLOAT temp;

   quit_altogether = FALSE;
   for (d=0;d<DIM;d++) size[d] = finish[d] - start[d];
   for (d=0;d<DIM;d++) vel[d] = 0.;

   /* create the particles */
   for (i=0;i<cnt;i++) {

      keep_trying = TRUE;
      ntried = 0;

      /* fprintf(stdout,"placing particle %d\n",i+1); */

      /* try to find an open place for the particle */
      while (keep_trying) {

         ntried++;
         isitfree = TRUE;

         for (d=0;d<DIM;d++) loc[d] = start[d] + size[d]*(rand()/(RAND_MAX+1.0));
         //fprintf(stdout,"  testing pt %g %g %g\n",loc[0],loc[1],loc[2]);

         /* is this location open? check walls, then other particles */
         //fprintf(stdout,"pt %g %g %g is free\n",loc[0],loc[1],loc[2]);
         for (d=0;d<DIM;d++) {
            if (loc[d]-start[d] < sim->new_part_rad) isitfree = FALSE;
            //fprintf(stdout," still free %d %g %g %g\n",d,loc[d],start[d],sim->new_part_rad);
            if (start[d]+size[d]-loc[d] < sim->new_part_rad) isitfree = FALSE;
            //fprintf(stdout," still free %d %g %g %g %g\n",d,start[d],size[d],loc[d],sim->new_part_rad);
         }
         if (isitfree) isitfree = is_this_space_open2(top,loc,sim->new_part_rad);
         //if (isitfree) fprintf(stdout," still free\n");
         if (isitfree) keep_trying = FALSE;

         /* have we just tried too many times already? */
         if (keep_trying && ntried > 2000) {
            keep_trying = FALSE;
            quit_altogether = TRUE;
            fprintf(stdout,"  Couldn't fit any more particles\n");
         }
      }

      if (quit_altogether) {
         /* get me out of this shitbag loop */
         i += cnt;
      } else {
         /* this location is good, place the particle */
         newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,loc,vel);
         add_particle_to_cell(sim,newpart,top);
         /* fprintf(stdout,"  placed pt there\n"); */
      }

      if (i/1000 != (i+1)/1000) {
         fprintf(stdout,".");
         fflush(stdout);
      }
   }
   fprintf(stdout,"\n");

   if (i>cnt) fprintf(stdout,"Placed %d of %d particles\n",i-cnt-1,cnt);
   else fprintf(stdout,"Placed %d of %d particles\n",i,cnt);

   /* return 0 is all went well */
   return(0);
}


/*
 *  Create a sphere of randomly-placed particles
 */
int add_sphere_of_particles(sim_ptr sim,cell_ptr top,int cnt,FLOAT* c,FLOAT rad) {

   int i,d,keep_trying,ntried,isitfree,quit_altogether;
   int pancake = TRUE;
   FLOAT pos[DIM],s[DIM],dist;
   FLOAT loc[DIM],size[DIM];
   FLOAT vel[DIM];
   FLOAT tangvel;
   particle_ptr newpart;

   quit_altogether = FALSE;
   for (d=0;d<DIM;d++) s[d] = c[d]-rad;
   for (d=0;d<DIM;d++) vel[d] = 0.0;

   /* create the particles */
   for (i=0;i<cnt;i++) {

      keep_trying = TRUE;
      ntried = 0;

      // fprintf(stdout,"placing particle %d\n",i+1);

      /* try to find an open place for the particle */
      while (keep_trying) {

         ntried++;
         isitfree = TRUE;

         //for (d=0;d<DIM;d++) pos[d] = s[d] + 2.0*rad*(rand()/(RAND_MAX+1.0));
         for (d=0;d<DIM;d++) pos[d] = 2.0*rand()/(RAND_MAX+1.0) - 1.0;
         //if (pancake)
         //   for (d=2;d<DIM;d++) pos[d] = c[d] + 0.1*rad*(2.0*rand()/(RAND_MAX+1.0)-1.0);
         // fprintf(stdout,"  testing pt %g %g %g\n",xpos,ypos,zpos);

         // is this space within the sphere?
         dist = 0.;
         //for (d=0;d<DIM;d++) dist += pow(pos[d]-c[d],2);
         for (d=0;d<DIM;d++) dist += pow(pos[d],2);
         //dist = sqrt(dist);
         //if (dist+sim->new_part_rad > rad) isitfree = FALSE;
         if (dist > 1.0) isitfree = FALSE;
         // fprintf(stdout,"    dist (%g) + newrad (%g) > rad (%g)\n",dist,sim->new_part_rad,rad);

         // if (isitfree) isitfree = is_this_space_open2(top,pos,sim->new_part_rad);
         if (isitfree) keep_trying = FALSE;

         /* have we just tried too many times already? */
         if (keep_trying && ntried > 2000) {
            keep_trying = FALSE;
            quit_altogether = TRUE;
         }
      }

      // flatten it like a pancake
      if (pancake) for (d=2;d<DIM;d++) pos[d] *= 0.1;

      // position to center and radius
      for (d=0;d<DIM;d++) pos[d] = c[d] + rad*pos[d];

      if (quit_altogether) {
         /* get me out here */
         i += cnt;
      } else {
         /* this location is good, place the particle */
         //if (pancake) {
            // find distance from center of local cluster
            dist = 0.0;
            for (d=0;d<DIM;d++) dist += pow(pos[d]-c[d],2);
            dist = sqrt(dist);
            // find mass of all particles within that radius
            tangvel = (FLOAT)cnt * sim->new_part_mass * pow(dist/rad,DIM);
            // find tangential velocity to offset potential
            #if DIM==1
               tangvel = sqrt(2.*sim->G*tangvel);
            #elif DIM==2
               // apparently we still use the 1/r^2 influence even though that's wrong
               //tangvel = sqrt(2.*sim->G*tangvel/log(dist));
               tangvel = sqrt(2.*sim->G*tangvel/dist);
            #elif DIM==3
               tangvel = sqrt(2.*sim->G*tangvel/dist);
            #else
               //tangvel = sqrt(2.*sim->G*tangvel/pow(dist,DIM-2));
               tangvel = sqrt(2.*sim->G*tangvel/dist);
            #endif
            #if DIM==1
              vel[0] = 0.0;
            #else
              vel[0] = -(pos[1]-c[1])*tangvel/dist;
              vel[1] = (pos[0]-c[0])*tangvel/dist;
              for (d=2;d<DIM;d++) vel[d] = 0.;
            #endif
         //}
         newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,pos,vel);
         add_particle_to_cell(sim,newpart,top);
         /* fprintf(stdout,"  placed pt there\n"); */
      }

      if (i/10000 != (i+1)/10000) {
         fprintf(stdout,".");
         fflush(stdout);
      }
   }
   fprintf(stdout,"\n");

   if (i>cnt) fprintf(stdout,"Placed %d of %d particles\n",i-cnt-1,cnt);
   else fprintf(stdout,"Placed %d of %d particles\n",i,cnt);

   /* return 0 is all went well */
   return(0);
}


/*
 *  Create a sphere of non-intersecting particles
 */
int add_sphere_of_ni_particles(sim_ptr sim,cell_ptr top,int cnt,FLOAT* c,FLOAT rad) {

   int i,d,keep_trying,ntried,isitfree,quit_altogether;
   FLOAT pos[DIM],s[DIM],dist;
   FLOAT loc[DIM],size[DIM];
   FLOAT vel[DIM];
   particle_ptr newpart;

   quit_altogether = FALSE;
   for (d=0;d<DIM;d++) s[d] = c[d]-rad;
   for (d=0;d<DIM;d++) vel[d] = 0.0;

   /* create the particles */
   for (i=0;i<cnt;i++) {

      keep_trying = TRUE;
      ntried = 0;

      // fprintf(stdout,"placing particle %d\n",i+1);

      /* try to find an open place for the particle */
      while (keep_trying) {

         ntried++;
         isitfree = TRUE;

         for (d=0;d<DIM;d++) pos[d] = s[d] + 2.0*rad*(rand()/(RAND_MAX+1.0));
         // fprintf(stdout,"  testing pt %g %g\n",pos[0],pos[1]);

         // is this space within the sphere?
         dist = 0.;
         for (d=0;d<DIM;d++) dist += pow(pos[d]-c[d],2);
         dist = sqrt(dist);
         if (dist+sim->new_part_rad > rad) isitfree = FALSE;
         // fprintf(stdout,"    dist (%g) + newrad (%g) > rad (%g)\n",dist,sim->new_part_rad,rad);

         if (isitfree) isitfree = is_this_space_open2(top,pos,sim->new_part_rad);
         if (isitfree) keep_trying = FALSE;

         /* have we just tried too many times already? */
         if (keep_trying && ntried > 2000) {
            keep_trying = FALSE;
            quit_altogether = TRUE;
         }
      }

      if (quit_altogether) {
         /* get me out of this shitbag loop */
         i += cnt;
      } else {
         /* this location is good, place the particle */
         newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,pos,vel);
         add_particle_to_cell(sim,newpart,top);
         /* fprintf(stdout,"  placed pt there\n"); */
      }

      if (i/10000 != (i+1)/10000) {
         fprintf(stdout,".");
         fflush(stdout);
      }
   }
   fprintf(stdout,"\n");

   if (i>cnt) fprintf(stdout,"Placed %d of %d particles\n",i-cnt-1,cnt);
   else fprintf(stdout,"Placed %d of %d particles\n",i,cnt);

   /* return 0 is all went well */
   return(0);
}


/*
 *  Check to see if a point in space can support a non-intersecting particle
 */
int is_this_space_open(cell_ptr cell,FLOAT* x,FLOAT r){

   int i,isitfree;
   int passed_it_on = FALSE;
   FLOAT tr,dist;
   particle_ptr curr = cell->first;

   tr = 2.0*r;

   if (cell->level == 0) isitfree = TRUE;

   if (cell->has_subcells) {
      /* if the possible contact volume can be completely enclosed by one of the
       * subcells, then test that one instead */
//    for (i=0;i<2;i++)
//       if (x-tr > cell->s[i][0][0]->min[0] && x+tr < cell->s[i][0][0]->max[0])
//          for (j=0;j<2;j++)
//             if (y-tr > cell->s[0][j][0]->min[1] && y+tr < cell->s[0][j][0]->max[1])
//                for (k=0;k<2;k++)
//                   if (z-tr > cell->s[0][0][k]->min[2] && z+tr < cell->s[0][0][k]->max[2]) {
      for (i=0;i<NCHILD;i++) {
         isitfree = is_this_space_open(cell->s[i],x,r);
         passed_it_on = TRUE;
      }
   }

   if (!passed_it_on) {
      /* this means that we have to check against all particles in this cell
       * and its subcells */

      if (cell->has_subcells) {
         /* must check all of the particles in all of the subcells */
         for (i=0;i<NCHILD;i++) {
            isitfree = is_this_space_open(cell->s[i],x,r);
            if (!isitfree) return(FALSE);
         }
      } else {
         /* check the individual particles here */
         while (curr) {
            dist = 0.0;
            for (i=0;i<DIM;i++) dist += pow(x[i]-curr->x[0][i],2);
            dist = sqrt(dist);
            if (dist < tr) {
               /* the position is not free, stop right here */
               return(FALSE);
            } else {
               /* the verdict is still out, so keep looking */
               curr = curr->next;
            }
         }
      }
   }

   return(isitfree);
}


/*
 *  Check to see if a point in space can support a non-intersecting particle
 *
 *  This second method is clearly superior to the earlier method,
 *  there is no more banding of particles---the distribution is
 *  even. Now, to see if over time, the banding re-appears!
 */
int is_this_space_open2(cell_ptr cell,FLOAT* x,FLOAT r){

   int i,isitfree;
   FLOAT tr,dist;
   particle_ptr curr;

   tr = 2.0*r;
   //printf("  tr is %g\n",tr);
   // if (cell->level == 0) isitfree = TRUE;

   // first, if it's outside of the possible contact bounds of
   //   this cube, then just return TRUE (meaning there's nothing
   //   in this cell to prevent the placement of the particle
   for (i=0;i<DIM;i++) {
      //printf("    min/max %g / %g\n",cell->min[i]-tr,cell->max[i]+tr);
      if (x[i] < cell->min[i]-tr) return(TRUE);
      if (x[i] > cell->max[i]+tr) return(TRUE);
   }

   // then, if we have subcells, go check them all (yes, inefficient)

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++) {
         //printf("    checking child %d\n",i);
         isitfree = is_this_space_open2(cell->s[i],x,r);
         // if any of these cells finds an obstruction, then we
         //   can abort the rest of the comparisons
         if (!isitfree) return(FALSE);
      }
   } else {
      // check vs. all of the particles in this cell
      curr = cell->first;
      while (curr) {
         dist = 0.0;
         for (i=0;i<DIM;i++) dist += pow(x[i]-curr->x[0][i],2);
         dist = sqrt(dist);
         //printf("    checking distance %g\n",dist);
         if (dist < tr) {
            // the position is not free, stop right here
            return(FALSE);
         } else {
            // the verdict is still out, so keep looking
            curr = curr->next;
         }
      }
   }

   return(TRUE);
}


#ifdef ALLDONE
#ifndef GRAV_ONLY
/*
 *  Create a cube/rectangle of non-intersecting particles
 *  instead of recieving dl from subroutine, recieve radius and other stuff
 */
int add_strand_of_particles(sim_ptr sim,cell_ptr top,int strand_num,FLOAT rad,FLOAT dens,FLOAT E,FLOAT xs,FLOAT ys,FLOAT zs,FLOAT xf,FLOAT yf,FLOAT zf,FLOAT gx,FLOAT gy,FLOAT gz) {

   int i,num;
   FLOAT dl;
   FLOAT xpos,ypos,zpos;
   FLOAT sx,sy,sz;
   FLOAT loc[DIM],size[DIM];
   FLOAT vel[DIM] = {0.0,0.0,0.0};
   FLOAT tot_len,new_len;
   particle_ptr newpart,last_particle;
   strand_ptr this = &sim->strand[strand_num];

   // create and define the new strand
   this->index = strand_num;
   this->dens = dens;
   this->rad = rad;
   this->E = E;

   // fprintf(stdout,"recieved %lf %lf %lf  %lf %lf %lf\n",xs,ys,zs,xf,yf,zf);
   sx = xf-xs;
   sy = yf-ys;
   sz = zf-zs;
   tot_len = sqrt(sx*sx + sy*sy + sz*sz);
   // num is the number of segments, so num+1 is the number of particles to use
   // num = (int)(tot_len/dl);
   num = (int)(tot_len/(1.2*rad));
   this->nseg = num;
   dl = tot_len/num;
   this->dl = dl;
   fprintf(stdout,"Adding strand with length %g, %g and %d segments\n",tot_len,dl,num);

   // set longitudinal and bending stiffness
   this->k = this->E * rad*rad*M_PI / (num/tot_len);
   this->b = this->E * rad*rad*rad*rad*8.0*M_PI / 27.0;
   fprintf(stdout,"nseg %d  dl %g  k %g  b %g  mass %g\n",this->nseg,this->dl,this->k,this->b,dens*dl*rad*rad*M_PI);

   // now, this re-assignment of length is unneccessary
   // new_len = num*dl;
   // redefine sx,sy,sz to be the segment component lengths
   sx /= (FLOAT)num;
   sy /= (FLOAT)num;
   sz /= (FLOAT)num;

   // set the initial position
   loc[0] = xs;
   loc[1] = ys;
   loc[2] = zs;

   // now, the growth part
   if (isinf(gx)) {
      // it's not growing
      this->grows = FALSE;
   } else {
      this->grows = TRUE;
      fprintf(stdout,"   and it grows at %g %g %g\n",gx,gy,gz);
   }


   /* create the particles */

   // the first particle
   // newpart = new_particle(0,sim->new_part_mass,sim->new_part_rad,xpos,ypos,zpos,0.0,0.0,0.0);
   newpart = new_particle(0,dens*dl*rad*rad*M_PI,rad,loc,vel);
   add_particle_to_cell(sim,newpart,top);
   newpart->s_rad = dl/2.0;
   newpart->s_last = NULL;
   newpart->parent = this;
   last_particle = newpart;
   // sim->strand_start[strand_num] = newpart;
   // sim->strand[strand_num]->start = newpart;
   this->start = newpart;

   // add the middle particles
   for (i=1; i<num; i++) {
      // newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,
      //                        xpos+i*sx,ypos+i*sy,zpos+i*sz,0.0,0.0,0.0);
      loc[0] += sx;
      loc[1] += sy;
      loc[2] += sz;
      newpart = new_particle(i,dens*dl*rad*rad*M_PI,rad,loc,vel);
      add_particle_to_cell(sim,newpart,top);
      newpart->s_rad = dl/2.0;
      newpart->s_last = last_particle;
      newpart->parent = this;
      last_particle->s_next = newpart;
      last_particle = newpart;
   }

   // add the last particle
   // newpart = new_particle(num,sim->new_part_mass,sim->new_part_rad,
   //                        xpos+num*sx,ypos+num*sy,zpos+num*sz,0.0,0.0,0.0);
   loc[0] += sx;
   loc[1] += sy;
   loc[2] += sz;
   newpart = new_particle(num,dens*dl*rad*rad*M_PI,rad,loc,vel);
   add_particle_to_cell(sim,newpart,top);
   newpart->s_rad = dl/2.0;
   newpart->s_last = last_particle;
   last_particle->s_next = newpart;
   newpart->s_next = NULL;
   newpart->parent = this;

   // make appropriate modifications if this is a growing strand
   if (this->grows) {
      // first particle is stationary
      this->start->stationary = TRUE;
      // next part is, too
      this->start->s_next->stationary = TRUE;
   }

   // if all went well, return 0
   return (0);
}
#endif


/*
 *  Create a galaxy with cnt particles
 *  call with
 *  add_galaxy(sim,top,cnt,x,y,z,nx,ny,nz,rad,vx,vy,vz);
 */
//int add_galaxy(sim_ptr sim,cell_ptr top,int cnt,FLOAT x,FLOAT y,FLOAT z,FLOAT nx,FLOAT ny,FLOAT nz,FLOAT rad,FLOAT vx,FLOAT vy,FLOAT vz){
int add_galaxy_old(sim_ptr sim,cell_ptr top,int cnt,FLOAT x,FLOAT y,FLOAT z,FLOAT nx,FLOAT ny,FLOAT nz,FLOAT rad,FLOAT vx,FLOAT vy,FLOAT vz){

   int i;
   /* v[2][0:2] is the basis vector for the rotational axis of the galaxy */
   FLOAT v[3][3];
   FLOAT len;
   FLOAT mpp = 1.0/cnt;
   FLOAT theta,r,n;
   FLOAT xpos,ypos,zpos,xvel,yvel,zvel;
   particle_ptr newpart;

   /* normalize the vector normal to the plane of the galaxy */
   len = sqrt(nx*nx+ny*ny+nz*nz);
   v[2][0] = nx/len;
   v[2][1] = ny/len;
   v[2][2] = nz/len;

   /* create the basis vectors v[0:2][] */

   /* create the particles */
   for (i=0;i<cnt;i++) {

      /* theta is nicely random */
      theta = 2.0*M_PI*(rand()/(RAND_MAX+1.0));

      /* determine r using exponential power law */
      /* but not yet! */
      r = rad*(rand()/(RAND_MAX+1.0));

      /* and determine n, the z component, by a power law, too */
      /* but not yet! */
      n = 0.1*rad*(rand()/(RAND_MAX+1.0));

      /* calculate particle initial location */
      xpos = x + r*cos(theta)*v[0][0] + r*sin(theta)*v[1][0] + n*v[2][0];
      ypos = y + r*cos(theta)*v[0][1] + r*sin(theta)*v[1][1] + n*v[2][1];
      zpos = z + r*cos(theta)*v[0][2] + r*sin(theta)*v[1][2] + n*v[2][2];

      /* determine particle's orbital speed */

      /* calculate paricle's initial velocity */
      xvel = vx;
      yvel = vy;
      zvel = vz;

      newpart = new_particle(i,mpp,sim->new_part_rad,pos,vel);
      add_particle_to_cell(sim,newpart,top);
   }

   /* return 0 is all went well */
   return(0);
}
#endif


/*
 * Create a galaxy according to the Plummer model, scaled to units such
 * that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 * 
 * (Is it obvious that I grabbed this code from Josh Barnes?)
 * Thanks, dude!
 *
 * cnt is the number of particles in this galaxy
 * x is the position of the center of mass of the galaxy
 * nx is the vector normal to the plane of the galaxy
 * vx is a speed to add to all particles in the galaxy
 * rad is the smoothing radius of the particle
 */
int add_galaxy (sim_ptr sim, cell_ptr top, int cnt,
   FLOAT *px, FLOAT *nx, FLOAT rad, FLOAT *vx) {

   int i,d;
   FLOAT rsc, vsc, r, v, x, y;
   FLOAT pos[DIM],vel[DIM];
   particle_ptr newpart;

   // choose scaling factors
   rsc = rad * (3 * PI) / 16;
   vsc = sqrt(1.0 / rsc);

   // for each of the new particles...
   for (i=0; i<cnt; i++) {

      // pick mass shell radius
      x = 0.999*(rand()/(RAND_MAX+1.0));
      // convert it to a real radius
      r = 1.0 / sqrt(pow(x, -2.0/3.0) - 1);
      // and choose that position on an N-dimensional sphere
      pick_point_on_sphere(pos,r*rsc);

      // now, choose the velocity in much the same way
      do {
         x = 1.0*(rand()/(RAND_MAX+1.0));
         y = 0.1*(rand()/(RAND_MAX+1.0));
      } while (y > x*x * pow(1 - x*x, 3.5));
      v = x * sqrt(2.0 / sqrt(1 + r*r));
      pick_point_on_sphere(vel,v*vsc);

      // rotate to match the normal vector

      // scale the location and velocity
      for (d=0; d<DIM; d++) pos[d] += px[d];
      for (d=0; d<DIM; d++) vel[d] += vx[d];

      // finally, create the particle
      newpart = new_particle(i,sim->new_part_mass,sim->new_part_rad,pos,vel);
      add_particle_to_cell(sim,newpart,top);
   }

   return(0);
}


/*
 * find a random point on a given sphere shell centered at the origin
 *
 * inputs
 *   rad	radius of shell
 *
 * outputs
 *   loc[DIM]  new point
 */
void pick_point_on_sphere (FLOAT* loc,FLOAT rad) {

   int d;
   FLOAT temp[DIM];
   FLOAT len;

   // then move it due to random motion
#if DIM==1
   temp[0] = rand()/(RAND_MAX+1.0);
   if (temp[0] > 0.5) {
      loc[0] = rad;
   } else {
      loc[0] = -rad;
   }
#elif DIM==2
   // second method uses no iteration at all, because sin/cos is accelerated
   temp[0] = 6.2831853071795864 * rand()/(RAND_MAX+1.0);
   loc[0] = cos(temp[0]);
   loc[1] = sin(temp[0]);
   for (d=0;d<DIM;d++) loc[d] = rad*loc[d];
#elif DIM==3
   // new method needs no sqrt
   len = 2.;
   while (len > 1.) {
      for (d=0;d<DIM;d++) temp[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      len = temp[0]*temp[0]+temp[1]*temp[1];
   }
   loc[2] = 2.*len -1.;
   len = 2.*sqrt(1.-len);
   loc[0] = temp[0]*len;
   loc[1] = temp[1]*len;
   for (d=0;d<DIM;d++) loc[d] = rad*loc[d];
#else
   // else use the standard method
   len = 2.;
   while (len > 1.) {
      for (d=0;d<DIM;d++) loc[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      len = 0.0;
      for (d=0;d<DIM;d++) len += pow(loc[d],2);
      //len = vec_length_sq(loc);
   }
   len = sqrt(len);
   for (d=0;d<DIM;d++) loc[d] = rad*loc[d]/len;
#endif

   return;
}


/*
 *  Create a new particle
 */
particle_ptr new_particle(int i,FLOAT m,FLOAT r,FLOAT* x,FLOAT* u){

   int d;
   particle_ptr newpart;

   newpart = (PARTICLE*)malloc(sizeof(PARTICLE));
   newpart->index = i;
   newpart->mass = m;
   newpart->logmass = log(abs(m));
   newpart->rad = r;
#ifndef GRAV_ONLY
   newpart->s_rad = r;
#endif
   for (d=0;d<DIM;d++) newpart->x[0][d] = x[d];
   for (d=0;d<DIM;d++) newpart->u[0][d] = u[d];
#ifndef GRAV_ONLY
   for (d=0;d<DIM;d++) newpart->a[d] = 0.0;
#endif
   newpart->stationary = FALSE;
   newpart->flag = FALSE;
   newpart->next = NULL;
#ifndef GRAV_ONLY
   newpart->s_last = NULL;
   newpart->s_next = NULL;
   newpart->parent = NULL;
#endif

   // fprintf(stderr,"Added particle %d with mass %g and radius %g\n",newpart->index,newpart->mass,newpart->rad);

   return(newpart);
}


/*
 *  Create a new stationary particle
 */
particle_ptr new_stationary_particle(int i,FLOAT r,FLOAT* x){

   int d;
   particle_ptr newpart;

   newpart = (PARTICLE*)malloc(sizeof(PARTICLE));
   newpart->index = i;
   newpart->mass = 0.0;
   newpart->rad = r;
#ifndef GRAV_ONLY
   newpart->s_rad = r;
#endif
   for (d=0;d<DIM;d++) newpart->x[0][d] = x[d];
   for (d=0;d<DIM;d++) newpart->u[0][d] = 0.0;
#ifndef GRAV_ONLY
   for (d=0;d<DIM;d++) newpart->a[d] = 0.0;
#endif
   newpart->stationary = TRUE;
   newpart->flag = FALSE;
   newpart->next = NULL;
#ifndef GRAV_ONLY
   newpart->s_last = NULL;
   newpart->s_next = NULL;
   newpart->parent = NULL;
#endif

   return(newpart);
}


/*
 * Read the input file and set the listed parameters
 */
int read_input_file(fileprop_ptr file,sim_ptr sim,cell_ptr top) {

   int d;
   char token[13][32];
   char twochar[2];
   char sbuf[512];
   FILE *infile;

   /* open file for reading */
   infile = fopen(file->input_fn,"r");
   if (infile==NULL) {
      fprintf(stderr,"Could not open input file %s\n",file->input_fn);
      fflush(stderr);
      return(1);
   }
   fprintf(stdout,"Opening file %s\n",file->input_fn);
   fflush(stdout);

   /* read a line from the input file */
   while (fscanf(infile,"%[^\n]",sbuf) != EOF) {

      // fprintf(stdout,"%s\n",sbuf);

      /* grab the first word */
      sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s",token[0],token[1],token[2],token[3],token[4],token[5],token[6],token[7],token[8],token[9],token[10],token[11],token[12]);
      //fprintf(stdout,"   first token is %s\n",token[0]);

      if (token[0][0] == '#') {
         /* read a comment line, or some other descriptor */
         fscanf(infile,"%[\n]",twochar);    /* read up to newline */
         /* fprintf(stdout,"%s\n",sbuf);        /* write comment */

      } else {
         /* first word is a parameter name */
         /* read the parameter value */
         /* fprintf(stdout,"   second token is %s\n",token[1]); */

         if (strncmp(token[0],"outfile_root",12) == 0) {
            strcpy(file->out_fn_root,token[1]);
         } else if (strncmp(token[0],"particles",9) == 0) {
            sim->particle_cnt = atoi(token[1]);
         } else if (strncmp(token[0],"start_frame",11) == 0) {
            sim->next_output_index = atoi(token[1]);
         } else if (strncmp(token[0],"write_pgm",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_pgm = TRUE;
            else
               file->write_pgm = FALSE;
         /*} else if (strncmp(token[0],"write_gif",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_gif = TRUE;
            else
               file->write_gif = FALSE; */
         } else if (strncmp(token[0],"write_png",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_png = TRUE;
            else
               file->write_png = FALSE;
         } else if (strncmp(token[0],"write_dot",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_dot = TRUE;
            else
               file->write_dot = FALSE;
         } else if (strncmp(token[0],"write_rad",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_rad = TRUE;
            else
               file->write_rad = FALSE;
         } else if (strncmp(token[0],"write_part",10) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_part = TRUE;
            else
               file->write_part = FALSE;
         } else if (strncmp(token[0],"write_seg",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               file->write_seg = TRUE;
            else
               file->write_seg = FALSE;
         } else if (strncmp(token[0],"end_time",8) == 0) {
            sim->end_time = atof(token[1]);
         } else if (strncmp(token[0],"dt",2) == 0) {
            sim->dt = atof(token[1]);
         } else if (strncmp(token[0],"output_dt",9) == 0) {
            sim->output_dt = atof(token[1]);
         } else if (strncmp(token[0],"max_levels",10) == 0) {
            sim->max_levels = atoi(token[1]);
         } else if (strncmp(token[0],"max_ppc",7) == 0) {
            sim->max_parts_in_cell = atoi(token[1]);
         } else if (strncmp(token[0],"G",1) == 0) {
            sim->G = atof(token[1]);
         } else if (strncmp(token[0],"delta",5) == 0) {
            sim->delta = atof(token[1]);
         } else if (strncmp(token[0],"theta",5) == 0) {
            sim->theta = atof(token[1]);
         } else if (strncmp(token[0],"direct",6) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->compare_to_direct = TRUE;
            else
               sim->compare_to_direct = FALSE;
         } else if (strncmp(token[0],"comp_domain",11) == 0) {
            /* set initial size of first cell */
            for (d=0;d<DIM;d++) {
               top->min[d] = atof(token[2*d+1]);
               top->max[d] = atof(token[2*d+2]);
            }
         } else if (strncmp(token[0],"boundary",8) == 0) {
            if (strncmp(token[1],"wall",4) == 0) {
               for (d=0;d<DIM;d++) {
                  sim->bdry[d][0] = WALL;
                  sim->bdry[d][1] = WALL;
               }
            } else if (strncmp(token[1],"open",4) == 0) {
               for (d=0;d<DIM;d++) {
                  sim->bdry[d][0] = OPEN;
                  sim->bdry[d][1] = OPEN;
               }
            }
         } else if (strncmp(token[0],"image_size",10) == 0) {
            file->out_img_size = atoi(token[1]);
            for (d=0;d<DIM;d++) file->out->n[d] = atoi(token[1]);
            // fprintf(stdout,"img size is %d\n",file->out_img_size);
         } else if (strncmp(token[0],"image_depth",11) == 0) {
            file->image_depth = atoi(token[1]);
         } else if (strncmp(token[0],"use_self_grav",13) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->use_self_grav = TRUE;
            else
               sim->use_self_grav = FALSE;
         } else if (strncmp(token[0],"use_uniform_grav",16) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->use_uniform_grav = TRUE;
            else
               sim->use_uniform_grav = FALSE;
         } else if (strncmp(token[0],"use_contact",11) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->use_contact = TRUE;
            else
               sim->use_contact = FALSE;
         } else if (strncmp(token[0],"galaxy",6) == 0) {
            for (d=0;d<2*DIM+1;d++)
               sim->galaxy[sim->num_galaxy][d] = atof(token[d+1]);
            sim->num_galaxy++;
         } else if (strncmp(token[0],"g",1) == 0) {
            for (d=0;d<DIM;d++) sim->g[d] = atof(token[d+1]);
         } else if (strncmp(token[0],"rk",2) == 0) {
            sim->rk = atof(token[1]);
         } else if (strncmp(token[0],"rc",2) == 0) {
            sim->rc = atof(token[1]);
         } else if (strncmp(token[0],"rad",3) == 0) {
            sim->new_part_rad = atof(token[1]);
         } else if (strncmp(token[0],"mass",4) == 0) {
            sim->new_part_mass = atof(token[1]);
         } else if (strncmp(token[0],"contact_per_grav",4) == 0) {
            sim->contact_per_grav = atoi(token[1]);
         } else if (strncmp(token[0],"block",5) == 0) {
            for (d=0;d<2*DIM+1;d++)
               sim->block[sim->num_blocks][d] = atof(token[d+1]);
            sim->num_blocks++;
         } else if (strncmp(token[0],"niblock",7) == 0) {
            for (d=0;d<2*DIM+1;d++)
               sim->niblock[sim->num_niblocks][d] = atof(token[d+1]);
            sim->num_niblocks++;
         } else if (strncmp(token[0],"sphere",6) == 0) {
            for (d=0;d<DIM+2;d++)
               sim->sphere[sim->num_spheres][d] = atof(token[d+1]);
            sim->num_spheres++;
         } else if (strncmp(token[0],"nisphere",8) == 0) {
            for (d=0;d<DIM+2;d++)
               sim->nisphere[sim->num_nispheres][d] = atof(token[d+1]);
            sim->num_nispheres++;
         } else if (strncmp(token[0],"add_part",8) == 0) {
            for (d=1;d<2*DIM+1;d++)
               sim->indiv_part[sim->num_indiv_parts][d] = atof(token[d+1]);
            sim->num_indiv_parts++;
         } else if (strncmp(token[0],"read_part",9) == 0) {
            strcpy(sim->read_part[sim->num_read_part],token[1]);
            for (d=0;d<2*DIM;d++)
               sim->part_xform[sim->num_read_part][d] = atof(token[d+2]);
            sim->num_read_part++;
         } else if (strncmp(token[0],"read_stat",9) == 0) {
            strcpy(sim->read_stat[sim->num_read_stat],token[1]);
            for (d=0;d<DIM;d++)
               sim->stat_xform[sim->num_read_stat][d] = atof(token[d+2]);
            sim->num_read_stat++;
         } else if (strncmp(token[0],"density",7) == 0) {
            // sim->ff2->n[0] = atoi(token[1]);
            // sim->ff2->n[1] = atoi(token[1]);
            // sim->ff3->n[0] = atoi(token[1]);
            // sim->ff3->n[1] = atoi(token[1]);
            // sim->ff3->n[2] = atoi(token[1]);
            sim->use_density_field = TRUE;
#ifndef GRAV_ONLY
         } else if (strncmp(token[0],"strand",6) == 0) {
            sim->block[sim->num_strands][0] = atof(token[1]);
            sim->block[sim->num_strands][1] = atof(token[2]);
            sim->block[sim->num_strands][2] = atof(token[3]);
            sim->block[sim->num_strands][3] = atof(token[4]);
            sim->block[sim->num_strands][4] = atof(token[5]);
            sim->block[sim->num_strands][5] = atof(token[6]);
            sim->block[sim->num_strands][6] = atof(token[7]);
            sim->block[sim->num_strands][7] = atof(token[8]);
            sim->block[sim->num_strands][8] = atof(token[9]);
            // and more for growing strands!
            fprintf(stderr,"strlen(%s) is %d\n",token[10],strlen(token[10]));
            if (strlen(token[10]) > 1) {
               sim->block[sim->num_strands][9] = atof(token[10]);
               sim->block[sim->num_strands][10] = atof(token[11]);
               sim->block[sim->num_strands][11] = atof(token[12]);
            } else {
               sim->block[sim->num_strands][9] = 1./0.;
            }
            sim->num_strands++;
#endif
         }

         fscanf(infile,"%[\n]",twochar);    /* read newline */
      }
   }

   fclose(infile);

   /* ret_val is 0 is file is OK and read, 1 if not */
   return(0);
}


/*
 * Parse the command-line arguments
 */
int parse_args(int argc,char **argv,fileprop_ptr file) {

   int i,some_logical;
   FLOAT some_float;

   (void) strcpy(file->exectuable_fn,argv[0]);
   if (argc < 2) (void) Usage(file->exectuable_fn,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(file->exectuable_fn,0);
   if (strncmp(argv[1], "help", 4) == 0)
      (void) Usage(file->exectuable_fn,0);
   (void) strcpy(file->input_fn,argv[1]);
   for (i=2; i<argc; i++) {
      (void) Usage(file->exectuable_fn,0);
   }

   return(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for part-nd */
   static char **cpp, *help_message[] =
   {
       "  The input file should be a raw ASCII file containing 2 columns, the first",
       "  with a parameter name (dt,particles,G,write_rad) and the second with either",
       "  an integer or floating-point value or with a \'yes\' or \'no\'.",
       " ",
       "  There is running output to stdout, and PGM and Radiance data can be written",
       "  to disk.",
       " ",
       NULL
   };

   fprintf(stderr,"\n  Usage:  %s infile\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}


/*
 * allocate memory for a three-dimensional array of FLOATs
 */
FLOAT*** allocate_3d_array_F(int nx, int ny, int nz) {

   int i,j;
   FLOAT ***array = (FLOAT ***)malloc(nx * sizeof(FLOAT **));

   array[0] = (FLOAT **)malloc(nx * ny * sizeof(FLOAT *));
   array[0][0] = (FLOAT *)malloc(nx * ny * nz * sizeof(FLOAT));

   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   for (i=0; i<nx; i++) {
      if (i!=0)
         array[i][0] = array[0][0] + i * ny * nz;
      for (j=1; j<ny; j++)
         array[i][j] = array[i][0] + j * nz;
   }

   return(array);
}


/*
 * allocate memory for a two-dimensional array of FLOATs
 */
FLOAT** allocate_2d_array_F(int nx, int ny) {

   int i,j;
   FLOAT **array = (FLOAT **)malloc(nx * sizeof(FLOAT *));

   array[0] = (FLOAT *)malloc(nx * ny * sizeof(FLOAT));

   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}


/*
 * allocate memory for a two-dimensional array of png_byte
 */
png_byte** allocate_2d_array_pb(int nx, int ny, int depth) {

   int i,j,bytesperpixel;
   png_byte **array;

   if (depth <= 8) bytesperpixel = 1;
   else bytesperpixel = 2;
   array = (png_byte **)malloc(ny * sizeof(png_byte *));
   array[0] = (png_byte *)malloc(bytesperpixel * nx * ny * sizeof(png_byte));

   for (i=1; i<ny; i++)
      array[i] = array[0] + i * bytesperpixel * nx;

   return(array);
}


