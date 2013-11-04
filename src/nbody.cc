#include <mkl_vsl.h>
#include <omp.h>
#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include <iostream>
#include <fstream>

using namespace boost::python;
using namespace std;

typedef struct ParticleSystemType ParticleSystemType ;

// Class contains methods to interface with Python
class nbody {
  int Nparticles;
  int nSteps;
  float dt;
  float Temperature;
  int nProme;
  //struct ParticleSystemType { float *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz; };
  struct ParticleSystemType { float *x, *y, *z, *vx, *vy, *vz; };
  float box_sizeX, box_sizeY, box_sizeZ;
  float m_hiX,m_hiY,m_hiZ;
  float m_loX,m_loY,m_loZ;
  string boxfilename;
  float *buf;


	// Demonstrate extracting a derived object type that is built-in to Boost
	public: 
    void setSimulation(int N, int nSteps, float dt, int nProme, float T) {
      this->Nparticles = N;
      this->nSteps = nSteps;
      this->dt = dt;
      this->nProme = nProme;
      this->Temperature = T;
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

    void printParams() {
      cout << "+++++++++++++++++++++++++++++++++++++++++++++"<< endl;
      cout << "Numero de particulas: " << this->Nparticles << endl;
      cout << "Numero de pasos: " << this->nSteps << endl;
      cout << "Delta T: " << this->dt << endl;
      cout << "Numero de promedios: " << this->nProme << endl;
      cout << "Caja (Lx,Ly,Lz): " << this->box_sizeX << "," 
        << this->box_sizeY << "," << this->box_sizeZ << endl;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }

    void setBoxFilename(string s) {
      boxfilename = s;
    }

    void read_box() {
      cout << "Reading box from file: " << boxfilename << endl;
      ifstream BoxConfigFile;
      BoxConfigFile.open(boxfilename.c_str(), ifstream::in); // input
      if(!BoxConfigFile) {
        printf("Cannot open Box Config file: %s\n",boxfilename.c_str());
        exit(15);
      }
      BoxConfigFile >> this->Nparticles;
      BoxConfigFile >> this->box_sizeX;
      BoxConfigFile >> this->box_sizeY;
      BoxConfigFile >> this->box_sizeZ;
      int NP = this->Nparticles;
      long int size_mem = 9*(NP)*sizeof(float);
      ParticleSystemType p;
      buf = (float*) malloc(size_mem);
      if(buf==NULL){
        cout << "Error en el alojamiento" << endl;
        exit(15);
      } 
      p.x = buf + 0*NP; p.y = buf + 1*NP; p.z = buf + 2*NP;
      p.vx = buf + 3*NP; p.vy = buf + 4*NP; p.vz = buf + 5*NP;
      //p.fx = buf + 6*NP; p.vy = buf + 7*NP; p.vz = buf + 8*NP;
      cout << "Alojando: " <<  size_mem/1024 << " kb" << " para " << NP << " particulas" << endl;
      for( int i = 0; i < NP; i++ ){
        BoxConfigFile >> p.x[i];
        BoxConfigFile >> p.y[i];
        BoxConfigFile >> p.z[i];
        BoxConfigFile >> p.vx[i];
        BoxConfigFile >> p.vy[i];
        BoxConfigFile >> p.vz[i];
      }
      BoxConfigFile.close();
    }

    
    
    //realiza benchamark de unistride
	  int unitstride(){
      int NP = this->Nparticles;
      float boxLx = this->box_sizeX;
      float boxLy = this->box_sizeY;
      float boxLz = this->box_sizeZ; 
      float scale, TemperatureInstant;     

      // Particle data stored as an array of structures
      ParticleSystemType p;
      p.x = buf+0*NP; p.y = buf+1*NP; p.z = buf+2*NP;
      p.vx = buf+3*NP; p.vy = buf+4*NP; p.vz = buf+5*NP;
      //p.fx = buf+6*NP; p.fy = buf+7*NP; p.fz = buf+8*NP;
      cout << "ultima particula: " << p.vx[NP -1] << "," << p.vy[NP -1] << "," << p.vz[NP -1] << "\n" ;


      // Propagate particles
      printf("Propagating particles using %d threads...\n", omp_get_max_threads());
      cout << "Steps | Total Energy | Kinetic Energy | Potential Energy | forceX | ForceY | ForceZ" << endl;

      for (int step = 1; step <= nSteps; step++) {
        float Pot_energy = 0.0f;

        // #pragma omp parallel for schedule(static) private (Pot_energy)
        for (int i = 0; i < NP; i++) { 
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

              fxi += force*dx;
              fyi += force*dy;
              fzi += force*dz;
            }             
          }
          // cout << "fxi: " << fxi<< "\n";
          p.vx[i] = fxi*0.5*this->dt;
          p.vy[i] = fyi*0.5*this->dt;
          p.vz[i] = fzi*0.5*this->dt;
           
        }

        float kinetic_energy = 0.0f;
        float sumaX = 0.0;
        float sumaY = 0.0;
        float sumaZ = 0.0;
        #pragma unroll
        for (int i = 0 ; i < NP; i++) { // Not much work, serial loop
          p.x[i] += p.vx[i]*this->dt; 
          p.y[i] += p.vy[i]*this->dt; 
          p.z[i] += p.vz[i]*this->dt;

          kinetic_energy += p.vx[i]*p.vx[i]+p.y[i]*p.vy[i]+p.vz[i]*p.vz[i];
        }

        kinetic_energy = 0.5*kinetic_energy;
        TemperatureInstant = (2.0*kinetic_energy) / (3.0*(NP-3));
        scale = sqrt((double)(Temperature/TemperatureInstant));
        
        #pragma unroll  
        for(int i=0; i < NP; i++ ){
          p.vx[i] = scale*p.vx[i];
          p.vy[i] = scale*p.vy[i];
          p.vz[i] = scale*p.vz[i];
        }

        if(step%this->nProme == 0){          
          cout.precision(6);
          cout << step << "\t"
            // << (kinetic_energy+Pot_energy)/NP << "\t"
            // << kinetic_energy/NP << "\t"
            // << Pot_energy/NP << "\t"
            << (kinetic_energy+Pot_energy) << "\t"
            << kinetic_energy << "\t"
            << Pot_energy << "\t"
            << sumaX << "\t"
            << sumaY << "\t"
            << sumaZ << "\n";
          // printf("%d  %.4f\n",step,kinetic_energy);
          // fflush(stdout);
        }
      }
      //free(buf);
      return 0;
    }
};

// Expose classes and methods to Python
BOOST_PYTHON_MODULE(NbodyBench) {
	class_<nbody> ("interface")
		.def("setSimulation", &nbody::setSimulation)
		.def("setBox", &nbody::setBox)
		.def("unitstride", &nbody::unitstride)
    .def("printParams", &nbody::printParams)
    .def("setBoxFilename", &nbody::setBoxFilename)
    .def("read_box", &nbody::read_box)
	;
	
}
