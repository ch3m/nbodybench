#include <mkl_vsl.h>
#include <omp.h>
#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include <iostream>

using namespace boost::python;
using namespace std;

typedef struct ParticleSystemType ParticleSystemType ;

// Class contains methods to interface with Python
class nbody {
  int Nparticles;
  int nSteps;
  float dt;
  int nProme;
  struct ParticleSystemType { float *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz; };
  float box_sizeX, box_sizeY, box_sizeZ;
  float m_hiX,m_hiY,m_hiZ;
  float m_loX,m_loY,m_loZ;

	// Demonstrate extracting a derived object type that is built-in to Boost
	public: 
    void setSimulation(int N, int nSteps, float dt, int nProme) {
      this->Nparticles = N;
      this->nSteps = nSteps;
      this->dt = dt;
      this->nProme = nProme;
    }
  
    void setBox(float Lx, float Ly, float Lz) {
      this->box_sizeX = Lx;
      this->box_sizeY = Ly; 
      this->box_sizeZ = Lz;
      this->m_hiX = Lx/2.0;
      this->m_loX = -(this->m_hiX);
      this->m_hiY = Ly/2.0;
      this->m_loY = -(this->m_hiY);
      this->m_hiZ = Lz/2.0;
      this->m_loZ = -(this->m_hiZ);
    }
    
    //realiza benchamark de unistride
	  int unitstride(){
      int NP = this->Nparticles;
      float boxLx = this->box_sizeX;
      float boxLy = this->box_sizeY;
      float boxLz = this->box_sizeZ;

      // Particle data stored as an array of structures
      ParticleSystemType p; // Particle system stored as a structure of arrays
      float *buf = (float*) malloc(9*NP*sizeof(float)); // Malloc all data
      p.x = buf+0*NP; p.y = buf+1*NP; p.z = buf+2*NP;
      p.vx = buf+3*NP; p.vy = buf+4*NP; p.vz = buf+5*NP;
      p.fx = buf+6*NP; p.fy = buf+7*NP; p.fz = buf+8*NP;

      // Initialize particles
      VSLStreamStatePtr rnStream;
      vslNewStream( &rnStream, VSL_BRNG_MT19937, 1 );
      vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, 9*NP, buf, -this->box_sizeX, this->box_sizeX);

      // Propagate particles
      printf("Propagating particles using %d threads...\n", omp_get_max_threads());
      cout << "Steps | Total Energy | Kinetic Energy | Potential Energy | InstantTemperature" << endl;

      for (int step = 1; step <= nSteps; step++) {
        float Pot_energy = 0.0f;

        #pragma omp parallel for schedule(static) private (Pot_energy)
        for (int i = 0; i < NP; i++) { 
          // float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f; // Components of force on particle i
          float fxi=0.0f; 
          float fyi=0.0f; 
          float fzi=0.0f;
          #pragma vector always
          for (int j = 0; j < NP; j++) { // Vectorized loop over particles that exert force law of universal gravity calculation.
          // for (int j = i+1; j < NP-1; j++){
            if (j != i) {

              // Newton’slaw of universal gravity calculation.
              float dx = p.x[j] - p.x[i];
              float dy = p.y[j] - p.y[i];
              float dz = p.z[j] - p.z[i];

              //Baundary conditions
              if (dx >= this->m_hiX) dx -= boxLx;
              else if (dx < this->m_loX) dx += boxLx;
              if (dy >= this->m_hiY) dy -= boxLy;
              else if (dy < this->m_loY) dy += boxLy;
              if (dz >= this->m_hiZ) dz -= boxLz;
              else if (dz < this->m_loZ) dz += boxLz;

              const float rij2 = dx*dx + dy*dy + dz*dz;
              const float invrij2 = 1.0f/rij2;
              const float invrij6 = invrij2*invrij2*invrij2;   
              const float pot = 4.0*invrij6*(invrij6 - 1.0);
              Pot_energy =  pot + Pot_energy;
              const float force = 24.0*invrij6*( 2.0*invrij6 - 1.0 )*invrij2;
             // const float drSquared = dx*dx + dy*dy + dz*dz;
             // const float drPowerN32 = 1.0f/(drSquared*sqrtf(drSquared));

              fxi += force*dx;
              fyi += force*dy;
              fzi += force*dz;
              // Reduction to calculate the net force
              //Fx += dx * drPowerN32; Fy += dy * drPowerN32; Fz += dz * drPowerN32;
            }             
          }
          if(i==1000) cout << "pot: " << Pot_energy << endl;
          p.vx[i] += fxi*0.5*this->dt;
          p.vy[i] += fyi*0.5*this->dt;
          p.vz[i] += fzi*0.5*this->dt;
         // pBox.VelX[i] = pBox.VelX[i] + pBox.ForceX[i]*0.5*DeltaT;
           
          // Move particles in response to the gravitational force
//          p.vx[i] += this->dt*Fx; p.vy[i] += this->dt*Fy; p.vz[i] += this->dt*Fz;
        }
        float kinetic_energy = 0.0f;
        #pragma unroll
        for (int i = 0 ; i < NP; i++) { // Not much work, serial loop
          p.x[i] += p.vx[i]*this->dt; 
          p.y[i] += p.vy[i]*this->dt; 
          p.z[i] += p.vz[i]*this->dt;

          if(p.x[i] > this->box_sizeX){
            p.x[i]=p.x[i] - this->box_sizeX;
          }else if(p.x[i] < 0.0){
            p.x[i]=p.x[i] + this->box_sizeX;
          }
          if(p.y[i] > this->box_sizeY){
            p.y[i]=p.y[i] - this->box_sizeY;
          }else if(p.y[i] < 0.0){
            p.y[i]=p.y[i] + this->box_sizeY;
          }
          if(p.z[i] > this->box_sizeZ){
            p.z[i]=p.z[i] - this->box_sizeZ;
          }else if(p.z[i] < 0.0){
            p.z[i]=p.z[i] + this->box_sizeZ;
          }

          kinetic_energy += p.vx[i]*p.vx[i]+p.y[i]*p.vy[i]+p.vz[i]*p.vz[i];
        }
        if(step%this->nProme == 0){          
          cout.precision(6);
          cout << step << "\t"
            << (kinetic_energy+Pot_energy)/NP << "\t"
            << kinetic_energy/NP << "\t"
            << Pot_energy/NP << "\n";
          // printf("%d  %.4f\n",step,kinetic_energy);
          // fflush(stdout);
        }
      }
      free(buf);
      return 0;
    }
};

// Expose classes and methods to Python
BOOST_PYTHON_MODULE(NbodyBench) {
	class_<nbody> ("interface")
		.def("setSimulation", &nbody::setSimulation)
		.def("setBox", &nbody::setBox)
		.def("unitstride", &nbody::unitstride)
	;
	
}
