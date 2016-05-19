#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
//Constants
#define NA 125
#define KB 1.38064852E-23
#define NUM_BINS 100
#define NUM_EQUIL 50
#define NUM_STATS 100
#define SAMPLE_FREQ 1
#define SIGMA 3.4E-10
#define MC_MOVE 1E-10
#define VOL_MOVE 0.1
#define EPSILON 120.0*KB
#define TEMPERATURE 90.0
#define PRESSURE 101325.0
#define DENSITY 1.3954 //g/cm^3 at bp
#define ABS2ANG 1.0e10
//Global Variables
typedef struct atoms{
    float x;
    float y;
    float z;
    uint8_t atomic_number;
} atom;
atom old_argons[NA];
atom new_argons[NA];
float L = 0;
float first_L = 0;
float energies[NA * NA];
float avg_E = 0.0;
float avg_L = 0.0;
float accepted_trans_moves = 0.0;
float accepted_vol_moves = 0.0;
float total_trans_moves = 0.0;
float total_vol_moves = 0.0;
float bins[NUM_BINS];
float currdensity = 0.0;
//CODE
float random() {
    return ((float)rand() / (float)RAND_MAX);
}
void setup_system_new(atom a[]){
    float width = ceil(cbrtf((float) NA));
    first_L = 6*SIGMA;
    L = first_L;
    float spacing =  SIGMA;
    float offset = spacing/2.0;
    uint8_t count = 0;
    for(float xpos = offset; xpos <= spacing*width; xpos +=spacing){
        for(float ypos = offset; ypos <= spacing*width; ypos +=spacing){
            for(float zpos = offset; zpos <= spacing*width; zpos +=spacing){
                a[count].atomic_number = 18;
                a[count].x = xpos;
                a[count].y = ypos;
                a[count].z = zpos;
                count++;
            }
        }
    }
}
void setup_system(atom a[]){
  // initialize energy array
  //for () {
  //}
    float width = ceil(cbrtf((float) NA));
    first_L = cbrtf((float) NA/(DENSITY*100*100*100*6.022E23/39.948));
    L = first_L;
    float spacing =  L/width;
    float offset = spacing/2.0;
    uint8_t count = 0;
    for(float xpos = offset; xpos <= spacing*width; xpos +=spacing){
        for(float ypos = offset; ypos <= spacing*width; ypos +=spacing){
            for(float zpos = offset; zpos <= spacing*width; zpos +=spacing){
                a[count].atomic_number = 18;
                a[count].x = xpos;
                a[count].y = ypos;
                a[count].z = zpos;
                count++;
            }
        }
    }
}
float minimum_image_argon(float x1,float y1,float z1,float x2,float y2,float z2){
    float dx = (x2-x1) - L*round((x2-x1)/(L*0.5));
    float dy = (y2-y1) - L*round((y2-y1)/(L*0.5));
    float dz = (z2-z1) - L*round((z2-z1)/(L*0.5));
    float r = sqrtf(powf(dx,2)+powf(dy,2)+powf(dz,2));
    return r;
}
float total_energy(atom a[]){
  // only call me once!
    float total = 0.0;
    float r;
    for(int8_t i = NA; i>0; i--){
        for(int8_t j = i-1;j>0;j--){
            r= minimum_image_argon(a[i-1].x,a[i-1].y,a[i-1].z,a[j-1].x,a[j-1].y,a[j-1].z);
            // energies[i] = 4*EPSILON*(powf(SIGMA/r, 12)-powf(SIGMA/r,6));
            total += 4*EPSILON*(powf(SIGMA/r, 12)-powf(SIGMA/r,6));
        }
    }
    return total;
}
void print_locations(atom a[], atom b[]){
    for(int8_t i = NA; i>0; i--){
        printf("(x,y,z) : (%g,%g,%g) \t (%g,%g,%g) \n", a[i-1].x,a[i-1].y,a[i-1].z,b[i-1].x,b[i-1].y,b[i-1].z);
    }
}
void print_XYZ(atom a[], FILE* xyzfile){
    if (xyzfile == NULL) {
      printf("%d\nARGON BOX\n", NA);
      for(int8_t i = NA; i>0; i--){
        printf("Ar\t%g\t%g\t%g\n", a[i-1].x,a[i-1].y,a[i-1].z);
      }
    } else {
      fprintf(xyzfile, "%d\nARGON BOX\n", NA);
      for(int8_t i = NA; i>0; i--){
        fprintf(xyzfile, "Ar\t%g\t%g\t%g\n", a[i-1].x * ABS2ANG,a[i-1].y  * ABS2ANG, a[i-1].z * ABS2ANG);
      }
    }
}
void translational_move_all(){
    total_trans_moves+=(float)NA;
    float old_e = total_energy(old_argons);
    for(int8_t i = NA; i>0; i--){
      float random_x = (random()-0.5)*MC_MOVE;
      float random_y= (random()-0.5)*MC_MOVE;
      float random_z = (random()-0.5)*MC_MOVE;
      new_argons[i].x += random_x;
      new_argons[i].y += random_y;
      new_argons[i].z += random_z;
      if(new_argons[i].x<0){
          new_argons[i].x+=L;
      }
      if(new_argons[i].y<0){
          new_argons[i].y+=L;
      }
      if(new_argons[i].z<0){
          new_argons[i].z+=L;
      }
      if(new_argons[i].x>L){
          new_argons[i].x-=L;
      }
      if(new_argons[i].y>L){
          new_argons[i].y-=L;
      }
      if(new_argons[i].z>L){
          new_argons[i].z-=L;
      }
      float new_e = total_energy(new_argons);
      if (random()<expf(-(1.0/(KB *TEMPERATURE))*(new_e - old_e))){
          accepted_trans_moves +=1.0;
          old_e = new_e;
          memcpy(old_argons, new_argons, sizeof(new_argons));
      } else{
          memcpy(new_argons,old_argons,sizeof(new_argons));
      }
    }
}
void translational_move(){
    total_trans_moves+=1.0;
    float old_e = total_energy(old_argons);
    uint8_t random_particle = (int) round(random()*NA);
    float random_x = (random()-0.5)*MC_MOVE;
    float random_y= (random()-0.5)*MC_MOVE;
    float random_z = (random()-0.5)*MC_MOVE;
    new_argons[random_particle].x += random_x;
    new_argons[random_particle].y += random_y;
    new_argons[random_particle].z += random_z;
    if(new_argons[random_particle].x<0){
        new_argons[random_particle].x+=L;
    }
    if(new_argons[random_particle].y<0){
        new_argons[random_particle].y+=L;
    }
    if(new_argons[random_particle].z<0){
        new_argons[random_particle].z+=L;
    }
    if(new_argons[random_particle].x>L){
        new_argons[random_particle].x-=L;
    }
    if(new_argons[random_particle].y>L){
        new_argons[random_particle].y-=L;
    }
    if(new_argons[random_particle].z>L){
        new_argons[random_particle].z-=L;
    }
    float new_e = total_energy(new_argons);
    if (random()<expf(-(1.0/(KB *TEMPERATURE))*(new_e - old_e))){
        accepted_trans_moves +=1.0;
        memcpy(old_argons, new_argons, sizeof(new_argons));
    } else{
        memcpy(new_argons,old_argons,sizeof(new_argons));
    }
}
void volume_move(){
    total_vol_moves+=1.0;
    float old_e = total_energy(old_argons);
    float old_vol = powf(L,3);
    float new_vol = expf(logf(old_vol) + (random() - 0.5)*VOL_MOVE);
    float new_L = cbrtf(new_vol);
    for(int8_t i = NA; i>0;i--){
        new_argons[i].x*=(new_L/L);
        new_argons[i].y*=(new_L/L);
        new_argons[i].z*=(new_L/L);
    }
    float new_e = total_energy(new_argons);
    if (random()<expf(-(1.0/(KB *TEMPERATURE))*((new_e - old_e) + (PRESSURE*(new_vol- old_vol)) - ((NA+1)*(KB*TEMPERATURE)*logf(new_vol/old_vol))))){
        accepted_vol_moves +=1.0;
        memcpy(old_argons, new_argons, sizeof(new_argons));
        L=new_L;
    } else{
        memcpy(new_argons,old_argons,sizeof(new_argons));
    }
}
void setup_histogram(){
  for(int16_t i = 0; i <NUM_BINS; i++){
      bins[i] = 0.0;
  }
}
void make_histogram(){
  float r = 0.0;
  int index = 0;
  for(int8_t i = NA; i>0; i--){
      for(int8_t j = i-1;j>0;j--){
          r=minimum_image_argon(old_argons[i-1].x,old_argons[i-1].y,old_argons[i-1].z,old_argons[j-1].x,old_argons[j-1].y,old_argons[j-1].z);
          index = (int)(floor(r*NUM_BINS/(first_L*0.5)));
          if(index<NUM_BINS){
              bins[index]+=2.0;
          }
      }
  }
}
void save_histogram(){
    FILE *hist;
    hist = fopen("hist.data","w");
    float increment  = first_L/NUM_BINS;
    float rho = NA/powf(avg_L ,3);
    float curr_r = 0.0;
    float next_r = 0.0;
    for(int16_t i = 0; i < NUM_BINS; i++){
        curr_r = i*increment;
        next_r = (i+1)*increment;
        //*increment
        bins[i] = bins[i]/((NUM_STATS / SAMPLE_FREQ)*(4.0/3.0)*rho*M_PI*(powf(next_r,3)-powf(curr_r,3)));
        fprintf(hist,"%g\n",bins[i]);
    }
    fclose(hist);
}
void print_startup_info(){
  printf("*****************\n");
  printf("Argon Monte Carlo\nBy: Shiv Upadhyay\n");
  printf("*****************\n\n");
  printf("*****************\n");
  printf("Simulation Details\n");
  printf("*****************\n\n");
  printf("Number of argons:\t%d\n",NA);
  printf("Temperature:\t%g\n",TEMPERATURE);
  printf("Equilibration steps:\t%d\n",NUM_EQUIL);
  printf("Stats steps:\t%d\n",NUM_STATS);
  printf("Sampling frequency:\t%d\n",SAMPLE_FREQ);
  printf("Sigma:\t%g\n",SIGMA);
  printf("Epsilon:\t%g\n",EPSILON);
  printf("Max translational move:\t%g\n",MC_MOVE);
  printf("Max Vol move:\t%g\n\n",VOL_MOVE);
  printf("First L:\t%g\n",first_L);
}

int main(){
    FILE *output, *xyzfile;
    output = fopen("output.data","w");
    xyzfile = fopen("trajectory.xyz", "w");
    setup_histogram();
    setup_system(new_argons);
    // Comment out one of the two lines, they should be equivalent to
    // each other.
    setup_system(old_argons);
    // memcpy(&old_argons, &new_argons, sizeof(new_argons));
    print_startup_info();
    /* Don't attempt a vol move every n steps
        frenkel smit pg 119 */
    srand(time(NULL));
    print_XYZ(old_argons, NULL);
    printf("\n*****************\n");
    printf("Equilibration run\n");
    printf("*****************\n\n");
    fprintf(output,"Energy \t Vol \t L\n");
    for(int32_t n_equil = NUM_EQUIL; n_equil>=0; n_equil--){
        if(random()*(NA+1)+1<= NA){
            translational_move();
        }
        else{
            volume_move();
        }
    }
    printf("\nStats run\n");
    printf("*****************\n\n");
    float temp_E = 0.0;
    float count = 0.0;
    float percent1 = 0.0;
    float percent2 = 0.0;
    accepted_vol_moves = 0.0;
    accepted_trans_moves = 0.0;
    for(int32_t n_stats = NUM_STATS; n_stats>=0; n_stats--){
        temp_E = total_energy(old_argons);
        if (temp_E >1E6){
          print_XYZ(old_argons, NULL);
        }
        count +=1.0;
        percent1 = accepted_vol_moves/total_vol_moves;
        percent2 = accepted_trans_moves/total_trans_moves;
        printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n",L,total_energy(old_argons)*6.022E23/125,currdensity,total_vol_moves,total_trans_moves,percent1,percent2);
        if(random()*(NA+1)+1<= NA){
            translational_move();
        }
        else{
            volume_move();
        }
        if(n_stats % SAMPLE_FREQ == 0){
            make_histogram();
            temp_E = total_energy(old_argons);
            currdensity = (NA*39.948/6.022E23)/powf(L*100.0,3);
            fprintf(output,"%g\t%g\t%g\n", temp_E, powf(L,3),L);
            avg_E += temp_E/((float)NUM_STATS / SAMPLE_FREQ);
            avg_L += L/((float)NUM_STATS/ SAMPLE_FREQ);
        }
        print_XYZ(old_argons, xyzfile);
    }
    save_histogram();
    printf("Avg E %g\n", avg_E);
    printf("Avg L %g\n", avg_L);
    printf("\nXYZ OF ARGS\n");
    printf("*****************\n\n");
    print_XYZ(old_argons, NULL);
    fclose(output);
    fclose(xyzfile);
    printf("\n\n\nSuccessful Termination\n" );
    return 0;
}
