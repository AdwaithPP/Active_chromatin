 
/*
===========================================================
Brownian Dynamics Simulation of Active Chromatin Polymer
===========================================================

Description:
- Simulates a coarse-grained chromatin polymer (N beads)
- Includes:
    * Harmonic spring forces (bonded)
    * Soft repulsive interactions (non-bonded)
    * Active forces 


Output:
- XYZ trajectory file

===========================================================
*/

// ==============================
// SIMULATION PARAMETERS
// ==============================
// Nbeads  : number of beads in polymer
// N_STEPS : total simulation steps
// DT      : time step
// K       : spring constant
// L       : equilibrium bond length
// Teff    : effective temperature (noise strength)
//a_ext    : Extensile active force parameter
//b_ang    : Angular active force parameter

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define sqr(xx)     ((xx) * (xx))
#define Cube(xx)    ((xx) * (xx) * (xx))
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))


// Random number generator (Numerical Recipes' ran2)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum) {
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

// Gaussian noise 
double gaussian(long *idum) {
    double u1 = ran2(idum);
    double u2 = ran2(idum);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main() {
    int Nbeads = 250;
    int N_STEPS = 2000000;
    double DT = 0.001;
     double K = 5.9;
    double L = 1.0;
    double x[Nbeads], y[Nbeads], z[Nbeads];
    FILE *op;
    op = fopen("Output.xyz", "w");
    int i, j, bead;
   
double *FactiveR = malloc(Nbeads * sizeof(double));
double *FactiveL = malloc(Nbeads * sizeof(double));
double *FangR    = malloc(Nbeads * sizeof(double));
double *FangL    = malloc(Nbeads * sizeof(double));
if (!FactiveR || !FactiveL || !FangR || !FangL) {
    fprintf(stderr, "Allocation failed\n"); exit(1);
}

/* Active forces initialize to zero */
for (i = 0; i < Nbeads; ++i) {
    FactiveR[i] = FactiveL[i] = FangR[i] = FangL[i] = 0.0;
}

double sigma = L;   


//soft-potential parameters

double inv_r2,inv_r6, f_scalar ;

double epsilon = 1.0;
double sigma6 = pow(sigma, 6);
double sigma12 = sigma6 * sigma6;


//Adding Soft-potential Parameters

double V0 = 2.75;       
double rm = 1.35;    
double eta1 = 2.25;       
double eta2 = 3.0;       
double rm_eta1 = pow(rm, eta1);   


//activity parameter

double a_ext=28;
double b_ang=9.5;
double Teff=15.5;




double rcut = pow(2.0, 1.0/6.0) * sigma;
 double rcut2 = rcut * rcut; 




long idum = -(time(NULL)+3133);
  
   
    double dx,dy,dz,rr;
    // Initialize SAW
x[0] = y[0] = z[0] = 0.0;


double dx9,dy9,dz9,dr9,fx,fy,fz; 
    
double r0 = sigma + 0.2;
const double EPS_R = 1e-8;
double *xnew = malloc(Nbeads * sizeof(double));
double *ynew = malloc(Nbeads * sizeof(double));
double *znew = malloc(Nbeads * sizeof(double));



for(i = 1; i < Nbeads; i++) {
    int valid = 0;
    while (valid<1) {
        // generate random direction
        dx = 2.0 * ran2(&idum) - 1.0;
        dy = 2.0 * ran2(&idum) - 1.0;
        dz = 2.0 * ran2(&idum) - 1.0;
        rr = sqrt(sqr(dx) + sqr(dy) + sqr(dz));

        dx = r0 * dx / rr;
        dy = r0 * dy / rr;
        dz = r0 * dz / rr;

        // new tentative position
        x[i] = x[i - 1] + dx;
        y[i] = y[i - 1] + dy;
        z[i] = z[i - 1] + dz;

        // check against all previous beads
        valid = 1; 
        for(j = 0; j < i; j++) {
            double xr = x[i] - x[j];
            double yr = y[i] - y[j];
            double zr = z[i] - z[j];
            double dist2 = sqr(xr) + sqr(yr) + sqr(zr);
            if (dist2 < sigma * sigma) {
                valid = 0;
                break;
            }
        }
    }
}


    printf("Initial self-avoiding configuration generated with %d beads.\n", Nbeads);






int current_step = 0;

int start = (int)(ran2(&idum) * 3);
int num_groups = 0;
int groups[100][3]; 
double Factive_group[100] = {0}; // Store Factive for each group
double Fang_group[100] = {0};    // Store Fang for each group




int time_periods[100];  // Store time periods for each group
int current_event[100]; // Track the active event for each group

// Define groups and assign initial forces
while (start + 2 < Nbeads) {
    groups[num_groups][0] = start;
    groups[num_groups][1] = start + 1;
    groups[num_groups][2] = start + 2;

    // Select a random event and n
double n = ran2(&idum) * 4.5+0.5;
    double rand_val = ran2(&idum);
               


   if (rand_val < 0.5) {  // 50% chance
    Fang_group[num_groups] = b_ang * n;
    Factive_group[num_groups] = 0.0;
} else {  // 50% chance
    Factive_group[num_groups] = a_ext * n;
    Fang_group[num_groups] = 0.0;
}
    time_periods[num_groups] = (int)((2.0 / n) * 200);
    current_event[num_groups] = 0; 

    num_groups++;
    start += 3  ;
}


int last_update[100] = {0};  // track last update time for each group


for (current_step = 0; current_step < N_STEPS; current_step++) {

  
    // Check if we need to update forces for each group
    for (i = 0; i < num_groups; i++) {
        if (current_step - last_update[i] >= time_periods[i]) {
    // update this group
    last_update[i] = current_step;

            // Select a new event for this group
            // Select a random event and n
   double n = ran2(&idum) * 4.5+0.5;
           double rand_val = ran2(&idum);
           


if (rand_val < 0.5) {  // 50% chance
    Fang_group[i] = b_ang * n;
    Factive_group[i] = 0.0;
} else {  // 50% chance
    Factive_group[i] = a_ext * n;
    Fang_group[i] = 0.0;
}

            time_periods[i] = (int)((2.0 / n) * 200);
            current_event[i] = current_step;  // Track last update
        }



        // Assign forces based on current step
        for (j = 0; j < 3; j++) {
            int bead = groups[i][j];
            if (j == 0) {
                FactiveR[bead] = Factive_group[i];
                FangR[bead] = Fang_group[i];
            }
            if (j == 1) {
                FactiveR[bead] = Factive_group[i];
                FactiveL[bead] = Factive_group[i];
            }
            if (j == 2) {
                FactiveL[bead] = Factive_group[i];
                FangL[bead] = Fang_group[i];
            }
        }
    }



   
        for (bead = 0; bead < Nbeads; bead++) {
        
            double dxl = 0, dyl = 0, dzl = 0, drl = 1;
            double dxr = 0, dyr = 0, dzr = 0, drr = 1;
            double dx13l = 0, dy13l = 0, dz13l = 0, dr13l = 1;
            double dx13r = 0, dy13r = 0, dz13r = 0, dr13r = 1;

            if (bead > 0) {
                dxl = x[bead] - x[bead - 1];
                dyl = y[bead] - y[bead - 1];
                dzl = z[bead] - z[bead - 1];
                drl = sqrt(dxl * dxl + dyl * dyl + dzl * dzl);
                if (drl<EPS_R)
                {dxl=dyl=dzl=0;
                drl=1;
                }
            }
            if (bead < Nbeads - 1) {
                dxr = x[bead + 1] - x[bead];
                dyr = y[bead + 1] - y[bead];
                dzr = z[bead + 1] - z[bead];
                drr = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
                if (drr<EPS_R)
                {dxr=dyr=dzr=0;
                drr=1;
                }
            }
            if (bead > 1) {
                dx13l = x[bead] - x[bead - 2];
                dy13l = y[bead] - y[bead - 2];
                dz13l = z[bead] - z[bead - 2];
                dr13l = sqrt(dx13l * dx13l + dy13l * dy13l + dz13l * dz13l);
                if (dr13l<EPS_R)
                {dx13l=dy13l=dz13l=0;
                dr13l=1;
                }
            }
            if (bead < Nbeads - 2) {
                dx13r = x[bead + 2] - x[bead];
                dy13r = y[bead + 2] - y[bead];
                dz13r = z[bead + 2] - z[bead];
                dr13r = sqrt(dx13r * dx13r + dy13r * dy13r + dz13r * dz13r);
                if (dr13r<EPS_R)
                {dx13r=dy13r=dz13r=0;
                dr13r=1;
                }
            }
            
           fx=fy=fz=0;


           for (j = 0; j < Nbeads; j++) {
    if (abs(bead - j) > 1) {  
        dx9 = x[bead] - x[j];
        dy9 = y[bead] - y[j];
        dz9 = z[bead] - z[j];
        dr9 = dx9*dx9 + dy9*dy9 + dz9*dz9;  /* r^2 */

        if (dr9 > 1e-12) {   /* protect against r ~ 0 */
            double r = sqrt(dr9);

            if (r < rm) {   /* APPLY soft force only for r < r_m */
                   if (r < 1e-8) r = 1e-8;     /* floor r for stability */
               
                double r_eta1 = pow(r, eta1);    /* r^{n1} */
                double u = 1.0 - r_eta1 / rm_eta1;   /* u > 0 for r < rm */

               
                if (u < 1e-16) u = 1e-16;

                /* u^{n2-1} */
                double u_pow = pow(u, eta2 - 1.0);

                /* radial magnitude F(r) = V0*n1*n2 * u^{n2-1} * r^{n1-1} / rm^{n1} */
                double Fsoft = V0 * eta1 * eta2 * u_pow * pow(r, eta1-2) / rm_eta1;
                
                fx += Fsoft * dx9;
                fy += Fsoft * dy9;
                fz += Fsoft * dz9;
            }
        }
    } 
} 






         
          xnew[bead] = x[bead] - DT * K * (drl - L) * (dxl / drl) + DT * K * (drr - L) * (dxr / drr) + sqrt(2 * DT*Teff) * gaussian(&idum) + DT * FactiveL[bead] * (dxl / drl) - DT * FactiveR[bead] * (dxr / drr) - DT * FangL[bead] * (dx13l / dr13l) + DT * FangR[bead] * (dx13r / dr13r)+fx*DT;
                ynew[bead] = y[bead] - DT * K * (drl - L) * (dyl / drl) + DT * K * (drr - L) * (dyr / drr) + sqrt(2 * DT*Teff) * gaussian(&idum) + DT * FactiveL[bead]  * (dyl / drl) - DT * FactiveR[bead] * (dyr / drr) - DT * FangL[bead] * (dy13l / dr13l) + DT * FangR[bead] * (dy13r / dr13r)+fy*DT;
                znew[bead] = z[bead] - DT * K * (drl - L) * (dzl / drl) + DT * K * (drr - L) * (dzr / drr) + sqrt(2 * DT*Teff) * gaussian(&idum) + DT * FactiveL[bead] * (dzl / drl) - DT * FactiveR[bead] * (dzr / drr) - DT * FangL[bead] * (dz13l / dr13l) + DT * FangR[bead] * (dz13r / dr13r)+fz*DT;
         

                    
            }
        
      //Updating positions
        
        for (int bead = 0; bead < Nbeads; bead++) {
    x[bead] = xnew[bead];
    y[bead] = ynew[bead];
    z[bead] = znew[bead];
}   
        
    if (current_step % 100 == 0 && current_step > 600000) {
            fprintf(op, "%d\n", Nbeads);
            fprintf(op, "%d\n", current_step);
            for (j = 0; j < Nbeads; j++) {
                fprintf(op, "%d %f %f %f\n", 1, x[j], y[j], z[j]);

            }
        }  







}



        
   
    


    printf("Simulation completed.\n");
    
    /* cleanup */
    fclose(op);
   

    free(xnew); free(ynew); free(znew);
    free(FactiveR); free(FactiveL); free(FangR); free(FangL);


    return 0;
}

