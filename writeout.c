/*************************************************************
 *
 *  writeout.c - output subroutines for part-nd
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

int write_output(fileprop_ptr,sim_ptr,cell_ptr);
int write_particle_count(cell_ptr);
int write_2d_dots(fileprop_ptr,cell_ptr,int);
int write_pgm_density_3d(fileprop_ptr,cell_ptr,field3_ptr,char*);
int write_2d_density(fileprop_ptr,cell_ptr,field2_ptr,int);
int put_part_into_array(cell_ptr,cell_ptr,FLOAT**,int,int);
int write_rad_strand(fileprop_ptr,sim_ptr,cell_ptr);
int write_rad(fileprop_ptr,sim_ptr,cell_ptr);
int write_rad_cell(cell_ptr,FILE*,int);
int write_part(fileprop_ptr,sim_ptr,cell_ptr);
int write_part_cell(cell_ptr,FILE*);
int write_reg_part(fileprop_ptr,sim_ptr,cell_ptr);
int register_this_cells_parts(cell_ptr,particle_ptr*,int);
int write_png(char*,png_byte**,int,int,int);


/*
 * Controls writing of output, once decision to write has been made
 */
int write_output(fileprop_ptr file,sim_ptr sim,cell_ptr top) {

   /* actually write the output here */
   if (file->write_dot) {
      fprintf(stdout,"  writing dots, step %d\n",sim->next_output_index);
      write_2d_dots(file,top,sim->next_output_index);
   }

   if (file->write_rad) {
      fprintf(stdout,"  writing rad, step %d\n",sim->next_output_index);
#ifndef GRAV_ONLY
      if (sim->num_strands > 0) write_rad_strand(file,sim,top);
#endif
      write_rad(file,sim,top);
   }
   if (file->write_part) {
      fprintf(stdout,"  writing part, step %d\n",sim->next_output_index);
      // write_part(file,sim,top);
      write_reg_part(file,sim,top);
   }


   // create the density field, hopefully more later
   if (sim->use_density_field) {
      fprintf(stdout,"  writing dens, step %d\n",sim->next_output_index);
      create_density_field_2d(top,top,sim->ff2);
      write_2d_density(file,top,sim->ff2,sim->next_output_index);
   }

   return(0);
}


/*
 * write a heirarchical description of the particles in cells
 */
int write_particle_count(cell_ptr cell){

   int i,j,k;

   if (cell->level == 0) fprintf(stdout,"\n");

   if (cell->num > 0) {
      for (i=-1;i<cell->level;i++) fprintf(stdout,"  ");
      fprintf(stdout,"cell at level %d has %d particles\n",cell->level,cell->num);
   }

   if (cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         write_particle_count(cell->s[i]);

   return(0);
}


/*
 * Write a PGM image of the points projected onto the xy-plane
 */
int write_2d_dots(fileprop_ptr file,cell_ptr cell,int index){

   int do_fireworks = FALSE;
   int i,j;
   int nx = file->out_img_size;
   int ny = file->out_img_size;
   char filename[80];
   FLOAT scale = 1.0;
   FIELD2 my_array;		// initialize struct for density field
   field2_ptr array = &my_array;
   int printval;
   FILE *outfile;

   // fprintf(stderr,"I'm here %d\n",file->out_img_size); fflush(stderr);

   // make the full filename, default to png
   if (file->write_pgm) {
      sprintf(filename,"%sdots_%04d.pgm",file->out_fn_root,index);
   } else {
      sprintf(filename,"%sdots_%04d.png",file->out_fn_root,index);
   }

   /* zero the array */
   /* Not zeroing this array must use the pervious numbers, and this
      thing just adds the new locations to the old, looks like the
      burst of a firework! */
   if (do_fireworks) {
      /* do not clear memory */
   } else {
      for (j=ny-1; j>=0; j--) {
         for (i=0; i<nx; i++) {
            file->out->rho[i][j] = 0.0;
         }
      }
   }

   /* fill in the array */
   put_part_into_array(cell,cell,file->out->rho,nx,ny);

   /* if it's a PGM, write it here, else, write a PNG */
   if (file->write_pgm) {

      /* open file for writing */
      outfile = fopen(filename,"w");
      if (outfile==NULL) {
         fprintf(stderr,"Could not open output file %s\n",filename);
         fflush(stderr);
         exit(0);
      }

      /* plot a z-plane */
      if (file->image_depth == 8) {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*file->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               fprintf(outfile,"%d\n",printval);
            }
         }
      } else {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,65535);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*file->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               fprintf(outfile,"%d\n",printval);
            }
         }
      }

      /* close file */
      fclose(outfile);

   } else {

      // convert the floating pt "rho" into "png_byte"
      if (file->image_depth == 8) {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*file->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               file->image[ny-1-j][i] = (png_byte)printval;
            }
         }
         write_png(filename,file->image,nx,ny,8);
      } else {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*file->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               file->image[ny-1-j][2*i] = (png_byte)(printval/256);
               file->image[ny-1-j][2*i+1] = (png_byte)(printval%256);
            }
         }
         write_png(filename,file->image,nx,ny,16);
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a PGM image of the density field
 */
int write_pgm_density_3d(fileprop_ptr file,cell_ptr cell,field3_ptr ff,char *filename){

   int do_middle_slice = FALSE;
   int i,j,k;
   int kindex = ff->n[1]/2;
   int nx = ff->n[0];
   int ny = ff->n[2];
   FLOAT scale = 2.0;
   FLOAT **array = (FLOAT **)malloc(nx * sizeof(FLOAT *));
   int printval;
   FILE *outfile;

   /* allocate space for 2D array */
   array[0] = (FLOAT *)malloc(nx * ny * sizeof(FLOAT));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   /* zero the array */
   for (j=ny-1; j>=0; j--)
      for (i=0; i<nx; i++)
         array[i][j] = 0.0;

   /* fill the array */
   for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) {
         if (do_middle_slice) array[i][j] = ff->rho[i][kindex][j];
         else {
            for (k=0; k<ff->n[1]; k++)
               array[i][j] += ff->rho[i][k][j];
            array[i][j] *= 0.0025;
         }
      }

   // scale the array?


   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   /* plot a y-plane */
   fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
   for (j=ny-1; j>=0; j--) {
      for (i=0; i<nx; i++) {
         printval = (int)(256.0*array[i][j]*scale);
         if (printval<0) printval = 0;
         if (printval>255) printval = 255;
         fprintf(outfile,"%d\n",printval);
      }
   }

   /* close file */
   fclose(outfile);

   /* free memory from 2D array! */
   free(array[0]);
   free(array);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a PGM image of the density field
 */
//int write_2d_density(fileprop_ptr file,cell_ptr cell,field2_ptr ff,char *filename){
int write_2d_density(fileprop_ptr file,cell_ptr cell,field2_ptr ff,int index){

   int contrast_enhance = TRUE;
   int i,j;
   int nx = ff->n[0];
   int ny = ff->n[1];
   static FLOAT maxval = -1.0;
   FLOAT scale = 1.0;
   int printval;
   char filename[80];
   FILE *outfile;

   // make the full filename, default to png
   if (file->write_pgm) {
      sprintf(filename,"%sdens_%04d.pgm",file->out_fn_root,index);
   } else {
      sprintf(filename,"%sdens_%04d.png",file->out_fn_root,index);
   }

   // scale ff->rho to unit-maximum
   if (contrast_enhance || maxval < 0.0) {
      for (j=ny-1; j>=0; j--) {
         for (i=0; i<nx; i++) {
             if (ff->rho[i][j] > maxval) maxval = ff->rho[i][j];
         }
      }
      maxval *= 0.8;	// allow larger areas to be white
   }
   // fprintf(stdout,"maxval is %g\n",maxval);


   /* if it's a PGM, write it here, else, write a PNG */
   if (file->write_pgm) {

      /* open file for writing */
      outfile = fopen(filename,"w");
      if (outfile==NULL) {
         fprintf(stderr,"Could not open output file %s\n",filename);
         fflush(stderr);
         exit(0);
      }

      /* plot a y-plane */
      if (file->image_depth == 8) {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256.0*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               fprintf(outfile,"%d\n",printval);
            }
         }
      } else {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,65535);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536.0*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               fprintf(outfile,"%d\n",printval);
            }
         }
      }

      /* close file */
      fclose(outfile);

   } else {

      // convert the floating pt "rho" into "png_byte"
      if (file->image_depth == 8) {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               file->image[ny-1-j][i] = (png_byte)printval;
            }
         }
         write_png(filename,file->image,nx,ny,8);
      } else {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               file->image[ny-1-j][2*i] = (png_byte)(printval/256);
               file->image[ny-1-j][2*i+1] = (png_byte)(printval%256);
            }
         }
         write_png(filename,file->image,nx,ny,16);
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Put all of the particles in this cell onto the array
 */
int put_part_into_array(cell_ptr top,cell_ptr cell,FLOAT** array,int nx,int ny){

   int i,j,k;
   int xdim,ydim;
   FLOAT xstart,xsize,ystart,ysize,mult;
   particle_ptr curr;

   xdim = 0;	/* use x-values for x-dim of image */
   // ydim = 2;	/* use z-values for y-dim of image */
   ydim = 1;	/* use z-values for y-dim of image */

   if (cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         /* put_part_into_array(top,cell,&array[0][0],nx,ny); */
         put_part_into_array(top,cell->s[i],array,nx,ny);
   } else {
      /* run through the list of elems and place each on the grid */
      // mult = nx*ny*0.05;
      mult = 255.0*sqrt(nx*ny);
      xstart = top->min[xdim];
      xsize = top->max[xdim]-xstart;
      ystart = top->min[ydim];
      ysize = top->max[ydim]-ystart;
      curr = cell->first;
      while (curr) {
         i = (int)(nx*(curr->x[0][xdim]-xstart)/xsize);
         j = (int)(ny*(curr->x[0][ydim]-ystart)/ysize);
         if (i>-1 && j>-1 && i<nx && j<ny)
            array[i][j] += mult * curr->rad;
            // array[i][j] += mult * curr->mass;
         curr = curr->next;
      }
   }
   return(0);
}


/*
 * Write a Radiance description of the particles and cell edges
 */
int write_rad(fileprop_ptr file,sim_ptr sim,cell_ptr top){

   int type = 3;
   char filename[80];
   FILE *outfile;

   // make the filename
   sprintf(filename,"%s%04d.rad",file->out_fn_root,sim->next_output_index);

   /* open file for writing */
   outfile = fopen(filename,"w");
   // outfile = fopen(filename,"a");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# Radiance scene description of part-nd run\n");
   if (type == 0) {
      /* write all cell edges, and particles as lights */
      fprintf(outfile,"void plastic edgec 0 0 5 0.2 0.2 0.2 0.0 0.0\n");
      fprintf(outfile,"void light partc 0 0 3 1000 1000 1000\n");
   } else if (type == 1) {
      /* write stars as plastic objects, and a light source overhead */
      /* fprintf(outfile,"void plastic wallc 0 0 5  0.2 0.2 0.2 0.0 0.0\n"); */
      /* fprintf(outfile,"wallc polygon floor 0 0 12  0 0 0  1 0 0  1 1 0  0 1 0\n"); */
      fprintf(outfile,"void light lightc 0 0 3  10 10 10\n");
      fprintf(outfile,"lightc polygon l1 0 0 12  0 0 2  0 1 2  1 1 2  1 0 2\n");
      fprintf(outfile,"void plastic partc 0 0 5  0.5 0.5 0.5 0.0 0.0\n");
      if (sim->bdry[2][0] == WALL) fprintf(outfile,"partc polygon floor 0 0 12 %g %g %g %g %g %g %g %g %g %g %g %g\n", top->min[0],top->min[1],top->min[2],top->max[0],top->min[1],top->min[2], top->max[0],top->max[1],top->min[2],top->min[0],top->max[1],top->min[2]);
   } else if (type == 2) {
      /* write stars as glow objects */
      fprintf(outfile,"void glow partc 0 0 4  1.0 1.0 1.0 0.0\n");
   }

   write_rad_cell(top,outfile,type);

   /* close file */
   fclose(outfile);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a Radiance description of just this cell
 */
int write_rad_cell(cell_ptr c,FILE* outfile,int type){

   int i,j,k,cnt;
   FLOAT zp1,zp2,rad;
   particle_ptr curr;

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         write_rad_cell(c->s[i],outfile,type);
   } else {
      if (type == 0) {
         /* write out the edges */
         rad = 0.001;
         fprintf(outfile,"edgec cylinder e1 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->max[0],c->min[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e2 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->max[2],c->max[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e3 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->min[2],c->max[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e4 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->max[2],c->max[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e5 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->min[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e6 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->max[2],c->min[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e7 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->min[2],c->max[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e8 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->max[2],c->max[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e9 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->min[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e10 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->min[2],c->min[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e11 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->min[2],c->max[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e12 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->max[1],c->min[2],c->max[0],c->max[1],c->max[2],rad);
         fflush(outfile);
      } else if (type >= 1) {
         /* don't write out the edges */
      }

      /* write out the particles */
      curr = c->first;
      cnt = 0;
      while (curr) {
         /* if (type == 0) rad = 0.001*sqrt(curr->mass); */
         /* else if (type == 1) rad = 0.5*sqrt(curr->mass); */
         /* else if (type >= 1) rad = pow(curr->mass,1./3.)/50.; */
         /* rad = pow(curr->mass,1./3.)/50.; */
         // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->rad);

         // allow printing radiance files for 2-dimensional objects
         if (DIM < 3) {
            zp1 = 0.0;
         } else {
            zp1 = curr->x[0][2];
         }

         // if we're using cylinders instead of spheres:
         // for 50,000 particles+, use 0.05; for less, try 0.1
         rad = 0.009;
         if (type==3) {
            // this is leftover from dla-nd!
            // rad = 0.05*curr->rad*pow(curr->mass,1./3.);
            // write the *actual* radius, fool.
            rad = curr->rad;
            // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0][0],curr->x[0][1],curr->x[0][2],rad);
            fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0][0],curr->x[0][1],zp1,rad);
         } else {
            // just write the spheres
            rad = curr->rad;
            // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->rad);
            fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0][0],curr->x[0][1],zp1,curr->rad);
         }

         // fflush(outfile);
         cnt++;
         // moved this outside of the loop bc of changes made 2002-09-08
         curr = curr->next;
      }
   }

   /* if all went well, return zero */
   return(0);
}


#ifndef GRAV_ONLY
/*
 * Write a Radiance description of the strands
 */
int write_rad_strand(fileprop_ptr file,sim_ptr sim,cell_ptr top){

   int write_genworm = TRUE;
   int write_genworm_endcaps = TRUE;
   int i,j,cnt,numseg;
   FLOAT dir[2][3],ddss[2][3],len,radcurv[2];
   char filename[80];
   FILE *outfile;
   particle_ptr this,next,nextnext,last;

   // make the filename
   sprintf(filename,"%sstrand_%04d.rad",file->out_fn_root,sim->next_output_index);

   /* open file for writing */
   outfile = fopen(filename,"w");
   // outfile = fopen(filename,"a");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# Radiance scene description of part3d run\n");
   // fprintf(outfile,"void light lightc 0 0 3  10 10 10\n");
   // fprintf(outfile,"lightc polygon l1 0 0 12  0 0 2  0 1 2  1 1 2  1 0 2\n");
   // fprintf(outfile,"void plastic partc 0 0 5  0.5 0.5 0.5 0.0 0.0\n");
   // if (sim->bdry[2][0] == WALL) fprintf(outfile,"partc polygon floor 0 0 12 %g %g %g %g %g %g %g %g %g %g %g %g\n", top->min[0],top->min[1],top->min[2],top->max[0],top->min[1],top->min[2], top->max[0],top->max[1],top->min[2],top->min[0],top->max[1],top->min[2]);

   for (i=0; i<sim->num_strands; i++) {
      fprintf(outfile,"\n# strand %d\n",i);
      // this = sim->strand_start[i];
      this = sim->strand[i].start;
      cnt = 0;
      while (this) {
         next = this->s_next;
         if (next) nextnext = next->s_next;
         last = this->s_last;

         if (write_genworm && next) {

            // the directions are merely the lines joining i-1 and i+1
            if (last) {
               for (j=0;j<3;j++) dir[0][j] = next->x[0][j]-last->x[0][j];
               for (j=0;j<3;j++) ddss[0][j] = next->x[0][j]-2*this->x[0][j]+last->x[0][j];
            }
            if (nextnext) {
               for (j=0;j<3;j++) dir[1][j] = nextnext->x[0][j]-this->x[0][j];
               for (j=0;j<3;j++) ddss[1][j] = nextnext->x[0][j]-2*next->x[0][j]+this->x[0][j];
            }

            // but, if this is the first particle/segment, compute start direction
            if (!last) {
               for (j=0;j<3;j++) dir[0][j] = -3*this->x[0][j] + 4*next->x[0][j] - nextnext->x[0][j];
               for (j=0;j<3;j++) ddss[0][j] = ddss[1][j];
            }

            // or, if this is the last particle/segment, compute end direction
            if (!nextnext) {
               for (j=0;j<3;j++) dir[1][j] = 3*next->x[0][j] - 4*this->x[0][j] + last->x[0][j];
               for (j=0;j<3;j++) ddss[1][j] = ddss[0][j];
            }

            // normalize both
            len = sqrt(dir[0][0]*dir[0][0]+dir[0][1]*dir[0][1]+dir[0][2]*dir[0][2]);
            for (j=0;j<3;j++) dir[0][j] /= len;
            len = sqrt(dir[1][0]*dir[1][0]+dir[1][1]*dir[1][1]+dir[1][2]*dir[1][2]);
            for (j=0;j<3;j++) dir[1][j] /= len;
            // find the segment length, and scale each direction by it
            // (using the halflength still shows the segment bounds)
            len = sqrt(pow(next->x[0][0]-this->x[0][0],2)+pow(next->x[0][1]-this->x[0][1],2)+pow(next->x[0][2]-this->x[0][2],2));
            for (j=0;j<3;j++) dir[0][j] *= len;
            for (j=0;j<3;j++) dir[1][j] *= len;

            fprintf(outfile,"!genworm partc w_%d_%d ",i,cnt);
            fprintf(outfile,"'hermite (%g,%g,%g,%g,t)' ",this->x[0][0],next->x[0][0],dir[0][0],dir[1][0]);
            fprintf(outfile,"'hermite (%g,%g,%g,%g,t)' ",this->x[0][1],next->x[0][1],dir[0][1],dir[1][1]);
            fprintf(outfile,"'hermite (%g,%g,%g,%g,t)' ",this->x[0][2],next->x[0][2],dir[0][2],dir[1][2]);

            // find the radius of curvature (at both ends)

            for (j=0;j<3;j++) ddss[0][j] /= len*len;
            for (j=0;j<3;j++) ddss[1][j] /= len*len;
            // fprintf(stdout,"\n%g %g %g\n",ddss[0][0],ddss[0][1],ddss[0][2]);
            // fprintf(stdout,"%g %g %g\n",ddss[1][0],ddss[1][1],ddss[1][2]);

            radcurv[0] = 0.0;
            for (j=0;j<3;j++) radcurv[0] += ddss[0][j]*ddss[0][j];
            radcurv[0] = 1.0/sqrt(radcurv[0]);

            radcurv[1] = 0.0;
            for (j=0;j<3;j++) radcurv[1] += ddss[1][j]*ddss[1][j];
            radcurv[1] = 1.0/sqrt(radcurv[1]);

            len = sqrt(radcurv[0]*radcurv[1]);

            // write the worm radius and number of segments
            if (sim->num_strands > 100) numseg = (int)(1.5 + 10.0*atan(this->s_rad/len));
            else if (sim->num_strands > 10) numseg = (int)(1.5 + 20.0*atan(this->s_rad/len));
            else numseg = (int)(1.5 + 30.0*atan(this->s_rad/len));
            // fprintf(stdout,"%g %d\n",len,numseg);
            if (numseg > 15) numseg = 15;
            fprintf(outfile,"%g %d\n",this->rad,numseg);

            // genworm ends in a sphere, we want a cylinder
            if (write_genworm_endcaps) {

               // if this is the first particle
               if (!last) {

                  // compute the direction vector
                  len = sqrt(dir[0][0]*dir[0][0]+dir[0][1]*dir[0][1]+dir[0][2]*dir[0][2]);
                  // normalize it, and reverse the direction
                  for (j=0;j<3;j++) dir[0][j] /= -1.0*len;

                  // use the other dir, dir[1][0:2] as the end point
                  for (j=0;j<3;j++) dir[1][j] = this->x[0][j]+this->rad*dir[0][j];

                  // write the cylinder
                  fprintf(outfile,"partc cylinder c_%d_%d 0 0 7 ",i,cnt);
                  fprintf(outfile,"%g %g %g  %g %g %g  %g\n",this->x[0][0],this->x[0][1],this->x[0][2],dir[1][0],dir[1][1],dir[1][2],this->rad);

                  // write the ring
                  fprintf(outfile,"partc ring r_%d_%d 0 0 8 ",i,cnt);
                  fprintf(outfile,"%g %g %g  %g %g %g  0.0 %g\n",dir[1][0],dir[1][1],dir[1][2],dir[0][0],dir[0][1],dir[0][2],this->rad);

               // or, if this is the last particle
               } else if (!nextnext) {

                  // compute the direction vector at "next"
                  len = sqrt(dir[1][0]*dir[1][0]+dir[1][1]*dir[1][1]+dir[1][2]*dir[1][2]);
                  // normalize it, do not reverse the direction
                  for (j=0;j<3;j++) dir[1][j] /= len;

                  // use the other dir, dir[0][0:2] as the end point
                  for (j=0;j<3;j++) dir[0][j] = next->x[0][j]+next->rad*dir[1][j];

                  // write the cylinder
                  fprintf(outfile,"partc cylinder c_%d_%d 0 0 7 ",i,cnt);
                  fprintf(outfile,"%g %g %g  %g %g %g  %g\n",next->x[0][0],next->x[0][1],next->x[0][2],dir[0][0],dir[0][1],dir[0][2],next->rad);

                  // write the ring
                  fprintf(outfile,"partc ring r_%d_%d 0 0 8 ",i,cnt);
                  fprintf(outfile,"%g %g %g  %g %g %g  0.0 %g\n",dir[0][0],dir[0][1],dir[0][2],dir[1][0],dir[1][1],dir[1][2],next->rad);

               }
            }

         } else if (!write_genworm) {

            // write using raw segments, capped with rings

            // if this is the last particle, cap the strand
            if (!next) fprintf(outfile,"partc ring re_%d\n0 0 8 %g %g %g %g %g %g %g %g\n",i,this->x[0][0],this->x[0][1],this->x[0][2],this->x[0][0]-this->s_last->x[0][0],this->x[0][1]-this->s_last->x[0][1],this->x[0][2]-this->s_last->x[0][2],0.0,this->rad);

            // if this is the first particle, cap the strand
            if (!this->s_last) fprintf(outfile,"partc ring rs_%d\n0 0 8 %g %g %g %g %g %g %g %g\n",i,this->x[0][0],this->x[0][1],this->x[0][2],this->x[0][0]-next->x[0][0],this->x[0][1]-next->x[0][1],this->x[0][2]-next->x[0][2],0.0,this->rad);

            // write the sphere unless it's an end
            if (next && this->s_last) fprintf(outfile,"partc sphere s_%d_%d\n0 0 4 %g %g %g %g\n",i,cnt,this->x[0][0],this->x[0][1],this->x[0][2],this->rad);

            // write the segment only if this isn't the last one
            if (next) fprintf(outfile,"partc cylinder c_%d_%d\n0 0 7 %g %g %g %g %g %g %g\n",i,cnt,this->x[0][0],this->x[0][1],this->x[0][2],next->x[0][0],next->x[0][1],next->x[0][2],this->rad);

         }

         cnt++;
         this = this->s_next;
         fflush(outfile);
      }
   }

   /* close file */
   fclose(outfile);

   // fprintf(stderr,"Done writing strands.\n");

   // wait until a keystroke
   // getchar();

   /* if all went well, return zero */
   return(0);
}
#endif


/*
 * Write a re-readable native particle3d description of the particles
 */
int write_part(fileprop_ptr file,sim_ptr sim,cell_ptr top){

   char filename[80];
   FILE *outfile;
   particle_ptr curr;

   // make the filename
   sprintf(filename,"%s%04d.part",file->out_fn_root,sim->next_output_index);

   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# scene description of part-nd run\n");

   write_part_cell(top,outfile);

   /* close file */
   fclose(outfile);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a native particle3d description of just this cell
 */
int write_part_cell(cell_ptr c,FILE* outfile){

   int i,cnt;
   FLOAT rad;
   particle_ptr curr;

   /* fprintf(outfile,"\n# cell at level %d\n",c->level); */

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         write_part_cell(c->s[i],outfile);
   } else {
      /* write out the particles */
      curr = c->first;
      cnt = 0;
      while (curr) {
         fprintf(outfile,"%g %g %g %g %g %g %g\n",curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->rad,curr->u[0][0],curr->u[0][1],curr->u[0][2]);
         // fprintf(outfile,"%g %g %g %g\n",curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->rad);
         /* fflush(outfile); */
         cnt++;
         curr = curr->next;
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a re-readable native particle3d description of the particles, but
 * keep the particles "registered" (row N is always the same particle)
 */
int write_reg_part(fileprop_ptr file,sim_ptr sim,cell_ptr top){

   int i,d;
   static int first_time = TRUE;
   static particle_ptr *reglist;
   particle_ptr curr;
   char filename[80];
   FILE *outfile;

   // make the filename
   sprintf(filename,"%s%04d.part",file->out_fn_root,sim->next_output_index);

   // if this is the first time through, make a static array of particle
   // pointers
   if (first_time) {
      reglist = (particle_ptr *)malloc(top->num * sizeof(particle_ptr));
      first_time = FALSE;
      // fill the list with particle pointers
      register_this_cells_parts(top,reglist,0);
   }

   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# scene description of part-nd run\n");

   for (i=0;i<top->num;i++) {
      curr = reglist[i];
      // location
      for (d=0;d<DIM;d++)
         fprintf(outfile,"%g ",curr->x[0][d]);
      // radius
      fprintf(outfile,"%g ",curr->rad);
      // all but one velocity
      for (d=0;d<DIM-1;d++)
         fprintf(outfile,"%g ",curr->u[0][d]);
      // the last velocity
      fprintf(outfile,"%g\n",curr->u[0][DIM-1]);

      //fprintf(outfile,"%g %g %g %g %g %g %g\n",curr->x[0][0],curr->x[0][1],curr->x[0][2],curr->rad,curr->u[0][0],curr->u[0][1],curr->u[0][2]);
   }

   /* close file */
   fclose(outfile);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a native particle3d description of just this cell
 */
int register_this_cells_parts(cell_ptr c,particle_ptr *list,int total){

   int i,cnt;
   FLOAT rad;
   particle_ptr curr;

   /* fprintf(outfile,"\n# cell at level %d\n",c->level); */

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         total = register_this_cells_parts(c->s[i],list,total);
   } else {
      /* write out the particles */
      curr = c->first;
      while (curr) {
         list[total] = curr;
         total++;
         curr = curr->next;
      }
   }

   /* if all went well, return zero */
   return(total);
}


/*
 * write a png file
 */
int write_png(char *file_name,png_byte** image,int xres,int yres,int depth) {

   png_uint_32 k,height,width;
   // png_byte image[height][width*bytes_per_pixel];
   // png_bytep row_pointers[yres];
   int bytes_per_pixel=1;
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
   png_colorp palette;
   png_voidp user_error_ptr;
   // char *file_name = "out.png";

   height=yres;
   width=xres;

   /* open the file */
   fp = fopen(file_name, "wb");
   // fp = stdout;
   if (fp == NULL)
      return (-1);

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);
      // user_error_ptr, user_error_fn, user_warning_fn);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return (-1);
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      return (-1);
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return (-1);
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, height, width, depth, PNG_COLOR_TYPE_GRAY,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   /* Optional gamma chunk is strongly suggested if you have any guess
    * as to the correct gamma of the image.
    */
   png_set_gAMA(png_ptr, info_ptr, 2.2);

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* The easiest way to write the image (you may have a different memory
    * layout, however, so choose what fits your needs best).  You need to
    * use the first method if you aren't handling interlacing yourself.
    */
   // for (k = 0; k < height; k++)
     // row_pointers[k] = image + k*width*bytes_per_pixel;

   /* One of the following output methods is REQUIRED */
   // png_write_image(png_ptr, row_pointers);
   png_write_image(png_ptr, image);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   /* close the file */
   fclose(fp);

   /* that's it */
   return (0);
}


