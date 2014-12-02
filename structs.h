/***********************************************************
 *
 *  structs.h - data structures supporting part-nd
 *
 *  Copyright (C) 2000-2005  Mark J. Stock
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
 ***********************************************************/


/* 
 * Here are the includes
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <malloc.h>
#include "png.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Here are the defines:
 */

// probably keep this how it is
#ifndef FLOAT
//#define FLOAT double
#define FLOAT float
#endif

// use float or double-precision arithmetic
// make a gravitation-optimized code
//#define GRAV_ONLY

// are we using mass as the last dimension in the tree?
//#define MASS_IS_DIMENSION

//
// IMPORTANT: everything here and below shouldn't be changed
//

#if DIM==1
// number of spatial dimensions --- must set NCHILD manually (sorry!)
// NCHILD must be 2^(DIM+1) to accomodate mass as a tree variable or 2^DIM if not
#define NCHILD 2
#ifdef MASS_IS_DIMENSION
#define NCHILD 4
#else
#define NCHILD 2
#endif
#endif

#if DIM==2
// number of spatial dimensions --- must set NCHILD manually (sorry!)
// NCHILD must be 2^(DIM+1) to accomodate mass as a tree variable or 2^DIM if not
#ifdef MASS_IS_DIMENSION
#define NCHILD 8
#else
#define NCHILD 4
#endif
#endif

#if DIM==3
// number of spatial dimensions --- must set NCHILD manually (sorry!)
// NCHILD must be 2^(DIM+1) to accomodate mass as a tree variable or 2^DIM if not
#ifdef MASS_IS_DIMENSION
#define NCHILD 16
#else
#define NCHILD 8
#endif
#endif

#if DIM==4
// number of spatial dimensions --- must set NCHILD manually (sorry!)
// NCHILD must be 2^(DIM+1) to accomodate mass as a tree variable or 2^DIM if not
#ifdef MASS_IS_DIMENSION
#define NCHILD 32
#else
#define NCHILD 16
#endif
#endif

// all of these should not need changing, either
#define TRUE 1
#define FALSE 0
#define PI 3.1415926535
#define TWOPI 6.28318531
#define PIOT 1.57079633
#define OPEN 0
#define WALL 1
#define PERIODIC 2
#define LARGE 9.9e+99
#define mod(a,b) ((a)-(b)*floor((a)/(b)))


/*
 * Pointers to the two types of data structures below
 */
typedef struct particle_record *particle_ptr;
typedef struct cell_record *cell_ptr;
typedef struct strand_record *strand_ptr;


/*
 * A data structure definition of a participating particle.
 * The 1st dimension in x and u can be extended to 2 if a 
 * predictor-corrector method is to be used in the advection routine.
 */
typedef struct particle_record {
   unsigned int index;		// index of the node
   FLOAT mass;			// particle's mass
   FLOAT logmass;		// natural logarithm of particle's mass
   FLOAT rad;			// particle's radius
#ifndef GRAV_ONLY
   FLOAT s_rad;			// particle's effective radius for strand-neighbor spring force
#endif
   FLOAT x[1][DIM];		// node's location in space
   FLOAT u[1][DIM];		// node's velocity in space
#ifndef GRAV_ONLY
   FLOAT a[DIM];		// node's acceleration in space
#endif
   FLOAT ga[DIM];		// node's acceleration due to gravitational
 				// interactions alone
   char stationary;		// is the particle static (immobile)
   char flag;			// generic flag
   particle_ptr next;		// pointer to the next node in the list
#ifndef GRAV_ONLY
   particle_ptr s_last;		// the previous particle in this strand
   particle_ptr s_next;		// the next particle in this strand
   strand_ptr parent;		// pointer to the parent strand
#endif
} PARTICLE;

/*
 * A data structure definition of a participating strand
 */
typedef struct strand_record {
   unsigned int index;          // index of the strand
   FLOAT dens;			// density of the strand
   FLOAT rad;			// radius of the strand
   FLOAT E;			// Young's modulus of the strand
   FLOAT k;			// stiffness coefficient of the segments
   FLOAT c;			// damping coefficient of the segs
   FLOAT b;			// bending coefficient of the segs
   int nseg;			// number of segments (NOT nodes)
   particle_ptr start;		// pointer to the first node
   FLOAT dl;			// length per segment
   char grows;			// does the strand grow from one end?
   FLOAT growdir[DIM];		// here's the direction it grows!
} STRAND;

/*
 * A data structure definition of a cell in the nd-tree
 */
typedef struct cell_record {
   char level;			// recursion level
   char has_subcells;		// generic flag

   // for now, include both sets of pointers
   particle_ptr first;		// pointer to the three nodes, CCW
   // or
   cell_ptr s[NCHILD];		// pointers to the eight subcells

   FLOAT min[DIM+1];		// minimum coordinate of cube
   FLOAT mid[DIM+1];		// midpoint of cube
   FLOAT max[DIM+1];		// minimum coordinate of cube

   unsigned int num;		// number of particles in this cell
   FLOAT mass;			// total mass of particles in this cell
   FLOAT cm[DIM];		// physical center of mass
} CELL;


/*
 * Field structure and pointer
 */
typedef struct field3_record *field3_ptr;
typedef struct field3_record {
   int n[3];                    // integer dimensions of the field values
   FLOAT d[3];                  // size of each cell
   FLOAT o[3];                  // origin
   FLOAT ***rho;                // density array
} FIELD3;

typedef struct field2_record *field2_ptr;
typedef struct field2_record {
   int n[2];                    // integer dimensions of the field values
   FLOAT d[2];                  // size of each cell
   FLOAT o[2];                  // origin
   FLOAT **rho;                 // density array
} FIELD2;


typedef struct file_properties *fileprop_ptr;
typedef struct file_properties {
   int use_input_file;		// TRUE if input file is given on command line
   char input_fn[80];		// name of the input file, if not stdin
   char exectuable_fn[80];	// name of the executable
   char out_fn_root[80];	// root name of the output files
   //int write_gif;		// write a 2D GIF file of the scene
   int write_pgm;		// write a 2D ASCII PGM file of the scene
   int write_png;		// write PNG files
   int write_dot;		// write a simple 1-pixel-per-particle image
   int write_rad;		// write a 3D Radiance description of the scene
   int write_obj;		// write a scene description in obj format
   int write_part;		// write a part-nd description of the scene
   int write_seg;		// write a segmented file of the motion of the particles
   int out_img_size;		// the height and width of the raster output
   FIELD2 out_field;		// initialize struct for output array, re-usable, now
   field2_ptr out;
   int image_depth;		// bits per pixel: 8 or 16
   png_byte **image;		// memory for the png image
} FILEP;


typedef struct simulation_properties *sim_ptr;
typedef struct simulation_properties {
   int step;			// index of current integration step
   int maxstep;			// maximum number of steps, useless, really
   double time;			// current simulation time
   double start_time;		// time at the start of the simulation
   double end_time;		// time at the end of the simulation
   double dt;			// the dt for the current time step
   double max_dt;		// maximum dt to allow
   double output_dt;		// dt between output steps
   double next_output_time;	// time at which the next set of output file will be written
   int next_output_index;	// integer index for next set of output files
   int max_levels;		// maximum number of levels in nd-tree
   int max_parts_in_cell;	// maximum number of particles in a cell
   FLOAT vmax;			// current maximum speed of fastest particle
   int cell_count[20];		// count of number of cells at various levels
   int particle_cnt;		// initial particle count
   FLOAT new_part_rad;		// default particle radius
   FLOAT new_part_mass;		// default particle mass

   char use_self_grav;
   FLOAT G;			// gravitation constant
   FLOAT delta;			// added to distance checks before computing force
   FLOAT theta;			// accuracy parameter, Barnes&Hut criteria
   int compare_to_direct;	// run a direct simulation and find error
   FLOAT sserr,ssval;		// sum of square errors, sum square values
   FLOAT smerr,smval;		// sum of square errors, sum square values, by mass
   int scnt;			// number of entries in error analysis
   FLOAT rmserror;		// RMS acceleration error over all particles
   FLOAT rmsmerror;		// RMS acceleration error over all particles, by mass
   int contact_per_grav;	// number of contact steps per gravitational step

   char use_uniform_grav;
   FLOAT g[DIM];		// additional uniform gravitation vector

   char use_contact;
   FLOAT rk;			// repulsive force constant
   int rn;			// repulsive force exponent
   FLOAT rc;			// repulsive force damping constant
   FLOAT ak;			// attractive force constant
   int an;			// attractive force exponent

   int num_strands;		// number of strands
   // particle_ptr strand_start[1000];	// starting particles for the strands
   STRAND strand[1000];		// entries for the strands

   char bdry[DIM][2];		// boundary types, 0=OPEN, 1=WALL, 2=PERIODIC

   char use_density_field;      // flags use of density field calculations
   char construct_dens_surface; // flags construction of an isosurface of density
   FIELD2 density_field2;	// initialize struct for density field
   field2_ptr ff2;
   FIELD3 density_field3;	// initialize struct for density field
   field3_ptr ff3;

   // construction variables
   FLOAT block[100][2*DIM+1];	// store data for placing blocks of particles
   int num_blocks;
   FLOAT niblock[100][2*DIM+1];	// same, but for non-intersecting blocks
   int num_niblocks;

   FLOAT sphere[100][DIM+2];	// store data for placing blocks of particles
   int num_spheres;
   FLOAT nisphere[100][DIM+2];	// same, but for non-intersecting spheres
   int num_nispheres;

   FLOAT galaxy[100][2*DIM+1];	// store data for placing galaxies
   int num_galaxy;

   FLOAT indiv_part[100][2*DIM+1];// store data for one individual particle
   int num_indiv_parts;

   char read_part[10][80];	// part files to be read in and included
   FLOAT part_xform[10][2*DIM];	// translation and velocity of particles to be
 				// read in from said file
   int num_read_part;

   char read_stat[10][80];	// part files to be read in and included as stationary particles
   FLOAT stat_xform[10][DIM];	// translation of particles to be read in from said file
   int num_read_stat;

} SIMP;

/*
 * A timer, flags when stuff should be done
 */
typedef struct event_timer {
   double start_time;		// time at the start of the simulation
   double dt;			// time step for what the process is timing
   double next_time;		// next flagged time
   int next_index;		// integer index for next flagged time
} TIMER;


// from findvel.c
extern int find_new_vels2(sim_ptr,cell_ptr,int);
extern int find_new_vels(sim_ptr,cell_ptr,cell_ptr,int);
extern int find_acc_on_this_cells_parts(sim_ptr,cell_ptr,cell_ptr,int);
extern int find_grav_acc_on_this_part(sim_ptr,cell_ptr,particle_ptr);
extern int find_cont_acc_on_this_part(sim_ptr,cell_ptr,particle_ptr);
extern int find_acc_on_this_cell(sim_ptr,cell_ptr,cell_ptr,FLOAT*);
extern FLOAT find_vmax(cell_ptr,cell_ptr,FLOAT);

// from ndtree.c
extern cell_ptr build_new_tree(fileprop_ptr,sim_ptr,cell_ptr);
extern int add_particle_to_cell(sim_ptr,particle_ptr,cell_ptr);
extern int split_cell(sim_ptr,cell_ptr);
extern int count_cells(sim_ptr,cell_ptr);
extern int replace_all_particles(sim_ptr,cell_ptr,cell_ptr);
extern int clean_up_all_cells(sim_ptr,cell_ptr);
extern int clean_up_cells(sim_ptr,cell_ptr,int);
extern int merge_cell(cell_ptr);
extern int find_all_cells_cm(sim_ptr,cell_ptr);
extern int find_cell_cm(cell_ptr,int);

// from setup.c
extern int initialize_system(fileprop_ptr,sim_ptr,cell_ptr);
extern int add_particles_from_file(sim_ptr,cell_ptr,int,char);
extern int add_box_of_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT*);
extern int add_box_of_ni_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT*);
extern int add_sphere_of_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT);
extern int add_sphere_of_ni_particles(sim_ptr,cell_ptr,int,FLOAT*,FLOAT);
extern int is_this_space_open(cell_ptr,FLOAT*,FLOAT);
int add_galaxy (sim_ptr,cell_ptr,int,FLOAT*,FLOAT*,FLOAT,FLOAT*);
void pick_point_on_sphere (FLOAT*,FLOAT);
extern particle_ptr new_particle(int,FLOAT,FLOAT,FLOAT*,FLOAT*);
extern particle_ptr new_stationary_particle(int,FLOAT,FLOAT*);
extern int read_input_file(fileprop_ptr,sim_ptr,cell_ptr);
extern int parse_args(int,char**,fileprop_ptr);
extern int Usage(char[80],int);
extern FLOAT*** allocate_3d_array_F(int,int,int);
extern FLOAT** allocate_2d_array_F(int,int);

// from writeout.c
extern int write_output(fileprop_ptr,sim_ptr,cell_ptr);
extern int write_particle_count(cell_ptr);
extern int write_pgm_plane(fileprop_ptr,cell_ptr,char*);
extern int write_pgm_density_3d(fileprop_ptr,cell_ptr,field3_ptr,char*);
extern int write_pgm_density_2d(fileprop_ptr,cell_ptr,field2_ptr,char*);
extern int put_part_into_array(cell_ptr,cell_ptr,FLOAT**,int,int);
extern int write_rad_strand(fileprop_ptr,sim_ptr,cell_ptr);
extern int write_rad(fileprop_ptr,sim_ptr,cell_ptr);
extern int write_rad_cell(cell_ptr,FILE*,int);
extern int write_part(fileprop_ptr,sim_ptr,cell_ptr);
extern int write_part_cell(cell_ptr,FILE*);

// from density.c
extern int create_density_field_3d(cell_ptr,cell_ptr,field3_ptr);
extern int add_particle_to_field_3d(cell_ptr,particle_ptr,field3_ptr);
extern int create_density_field_2d(cell_ptr,cell_ptr,field2_ptr);
extern int add_particle_to_field_2d(cell_ptr,particle_ptr,field2_ptr);

// end

