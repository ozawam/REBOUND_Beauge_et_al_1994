#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* r);

double tmax = 3.0e4 * 2.0 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();

//Setup constants.
    r->dt                   = 0.005*2.0*M_PI;                // initial timestep

//Setup boundary.
//   r->boundary    = REB_BOUNDARY_OPEN;

// Setup integrator.
//   r->integrator           = REB_INTEGRATOR_IAS15;
//   r->integrator    = REB_INTEGRATOR_MERCURIUS;
//   r->integrator           = REB_INTEGRATOR_LEAPFROG;
   r->integrator    = REB_INTEGRATOR_WHFAST;
   r->ri_whfast.corrector = 11;    // Turn on symplectic correctors (11th order).
//   r->ri_whfast.safe_mode = 0;     // Turn of safe mode. Need to call integrator_synchronize() before outputs.

//Setup gravity solvers.
//    r->gravity    = REB_GRAVITY_TREE;
//    r->opening_angle2    = 0.1;        // This constant determines the accuracy of the tree code gravity estimate.
//    const double boxsize = 4.0;
//    reb_configure_box(r,boxsize,1,1,1);

// Setup collision.
 //   r->collision            = REB_COLLISION_NONE;
    r->collision            = REB_COLLISION_DIRECT;
 //   r->collision    = REB_COLLISION_TREE;
    r->collision_resolve    = reb_collision_resolve_hardsphere;        // Choose merger collision routine.
 //  r->softening     = 1.0e-6;        // Gravitational softening length

// Setup callback function for velocity dependent forces.
   r->additional_forces     = additional_forces; 
   r->force_is_velocity_dependent = 1;

// Setup callback function for outputs.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = 0;
//    r->usleep        = 10000;        // Slow down integration (for visualization only)

 //   double boxsize = 3.0;
 //   reb_configure_box(r,boxsize,1,1,1);


   // Star
    struct reb_particle star = {0};
    star.m = 1.0;
    star.r = 0.1;
    reb_add(r, star);

    //Jupiter
    struct reb_particle Jupiter = {0};
        double aJ    = 5.20260 ;
        double eJ    = 0.04851;// 
        double incJ  = 0.0;//two-dimention
        double OmegaJ = 0.0;//reb_random_uniform(0,2.*M_PI);
        double apsisJ = reb_random_uniform(0,2.*M_PI);
        double phiJ     = reb_random_uniform(0,2.*M_PI);     

	double M_sun = 1.989*pow(10.0,33.0);//g
    Jupiter.m = 1.898*pow(10.0,30.0)/M_sun;
    Jupiter = reb_tools_orbit_to_particle(r->G, star, Jupiter.m, aJ, eJ,incJ, OmegaJ, apsisJ, phiJ);
         reb_add(r, Jupiter);


   // Planetesimal disk parameters
    double M_earth = 5.9724*pow(10.0,27.0);//g
    double one_au = 1.49597870*pow(10.0,13.0);//cm
    double total_disk_mass = 15.0*M_earth/M_sun;
    int N_planetesimals = 10;
    double _rho_planet = 3.0;//g/cm^3
    double rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
    double planetesimal_mass = 0.0; //4.0/3.0 * pow(1000.0/one_au,3.0) *M_PI * rho_planet;//r=10m 2019/05/20 ;//total_disk_mass/N_planetesimals;
//    double delta_a = 0.02;
//    double central_a = 1.0;
//    double amin = central_a - (delta_a/2.0), amax = central_a + (delta_a/2.0);   //planet at inner edge of disk
    double amin = 6.5, amax = 15.0;
//    double powerlaw = 1.0;
    
    r->N_active =  2;

    // Generate Planetesimal Disk
     for (int i=1;i<= N_planetesimals;i++){
        struct reb_particle pt = {0};
        double a    = 9.5;//reb_random_uniform(amin,amax);
       
//        double e    = reb_random_rayleigh(0.0043);
//        double inc  = reb_random_rayleigh(0.00215);
        double e    = 0.01;// reb_random_rayleigh(0.003972877/sqrt(2.0));//circular orbit
        double inc  = 0.0;//reb_random_rayleigh(0.003972877/2.0/sqrt(2.0));//two-dimention

        double Omega = 0.0;//reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi     = reb_random_uniform(0,2.*M_PI);

        double rho = 2.0;
        double one_au = 1.49597870*pow(10.0,13.0);//cm
  
        double radius_f = 1.0;//fold enlargement(radius enhancement)

        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis, phi); //change 3d to 2d in 2019/05/20
        pt.m = planetesimal_mass;
        pt.r = pow((3.0*planetesimal_mass*M_sun)/(4.0*M_PI*rho) ,1.0/3.0)/ one_au * radius_f; //unit is AU.
        pt.lastcollision = 0;
        reb_add(r, pt);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

  //  reb_integrate(r, INFINITY);
    reb_integrate(r, tmax);

}


void additional_forces(struct reb_simulation* const r){
    // velocity dependent drag force.
	struct reb_particle* const particles = r->particles;
	//Setup constants. Units AU,Msun,T G=1
	double one_au = 1.49597870*pow(10.0,13.0);//cm
        double M_sun = 1.989*pow(10.0,33.0);//g
	double _rho_planet = 2.0;//g/cm^3
	double rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
	double mass_sun = 1.0;// particles[0].m;
	double lum_sun = 1.0;
	//posとvelはstar. を代用すること
	double PI = 4.0*atan(1.0);
	double alpha = 11.0 / 4.0; // exponent of rho
        double beta = 0.5; // exponent of temperature
	double Cd = 73.0; // coeficient aero drag
	double f_gas = 1.0/sqrt(100.0/280.0) ; // scaling factor of gas density //2019/05/21 considering the  effect of assumption of T =100
	double tau_gas = 1e6 * 2.0 * PI; // 1e6[yr] * 2PI[T]/[yr]

	//Calculation of Aero Drag
	//coef
//	double coef_rho_gas = exp(-1.0* r->t /tau_gas) * 2400.0 * f_gas * 1.49597871e13 * 1.49597871e13 / (0.047*2 * 1.9884e33) * pow(lum_sun, -0.125 )* pow(mass_sun, 0.5);
	double coef_rho_gas =  2400.0 * f_gas * 1.49597871e13 * 1.49597871e13 / (0.047*2 * 1.9884e33) * pow(lum_sun, -0.125 )    * pow(mass_sun, 0.5); //no gas dissipation
	double coef_cs_vk = 1.0/29.78 * pow(lum_sun, 0.125) * pow(mass_sun, -0.5); // about 0.0018 for solar luminosity and solar mass
	// 1.0/29.78 is cs/vk at 1AU (cs=1km/sec, vkep=29.78km/sec)
	double coef_acc_gd = 0.5*Cd*PI;

	//double dragcoefficient = 1e-10;

	const int N = r->N; //Sun + planetesimals
	for (int i=2;i<N;i++){
        //variable
	double r_sq = particles[i].x*particles[i].x + particles[i].y*particles[i].y;
	double inv_r = 1.0 / sqrt(r_sq);
	double r = r_sq * inv_r;
	double r_sqrt = sqrt(r);
	double r_4rt = sqrt(r_sqrt);
	double ev[3] = {-particles[i].y*inv_r, particles[i].x*inv_r, 0.0}; // unit vector of kepler velocity
	double vkep[3] = { sqrt(mass_sun * inv_r) * ev[0],sqrt(mass_sun * inv_r) * ev[1],sqrt(mass_sun * inv_r) * ev[2],}; // kepler velocity //0 x, 1 y, 2 z
	double cs_vk = coef_cs_vk * r_4rt;
	double eta = 0.5*(alpha+beta)*cs_vk*cs_vk;
	double vgas[3] = {0.995*vkep[0],0.995*vkep[1],0.995*vkep[2]};//{(1.0 - eta)*vkep[0],(1.0 - eta)*vkep[1],(1.0 - eta)*vkep[2]};
	double u[3] = {particles[i].vx - vgas[0],particles[i].vy - vgas[1],particles[i].vz - vgas[2]};
	double m_50 = 4.0/3.0 * pow(1000.0/one_au,3.0) *M_PI * rho_planet;//r=10m 2019/05/20
	double rate_50 = 8.958600e+25/M_sun/m_50;
	double rplanet = 1000/one_au;//r=10m 2019/05/20  //cbrt(3.0*particles[i].m/rate_50/(4.0*PI*rho_planet));//cbrt(3.0*particles[i].m/(4.0*PI*rho_planet));
	double rho_gas = coef_rho_gas * inv_r * inv_r * inv_r * r_4rt;

        //double sys_acc_gd[3] = { (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[0]*u[0]) * u[0]) / particles[i].m, (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[1]*u[1]) * u[1]) / particles[i].m, (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[2]*u[2]) * u[2]) / particles[i].m } ;  
	double sys_acc_gd[3] = {(-11.1*Cd/(rplanet*one_au * rplanet*one_au) *u[0]),(-11.1*Cd/(rplanet*one_au * rplanet*one_au) *u[1]) ,(-11.1*Cd/(rplanet*one_au * rplanet*one_au) *u[2]) };
	particles[i].ax += sys_acc_gd[0];
        particles[i].ay += sys_acc_gd[1];
        particles[i].az += sys_acc_gd[2];
    }
}


void heartbeat(struct reb_simulation* r){
     //   int snap_n = (int)(r->t)/(100.0*2.0*M_PI);
        int snap_n = (r->hash_ctr);
        char SNAP[30];
    if (reb_output_check(r, 1.0*2.0*M_PI)){  
        reb_output_timing(r, tmax);
           }

    if (reb_output_check(r, 100.0*2.0*M_PI)){
   //Synchronize particles manually at end of timestep
   //    reb_integrator_synchronize(r);
       sprintf(SNAP,"output/snap01/snap%05d.dat",snap_n);
       reb_output_orbits(r,SNAP);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }
}

