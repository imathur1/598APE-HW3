#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gperftools/profiler.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include <vector>

#ifndef PROFILE
   #define PROFILE 0
#endif

#ifndef PRINT_STEPS
   #define PRINT_STEPS 0
#endif

float tdiff(struct timeval *start, struct timeval *end) {
  return (end->tv_sec-start->tv_sec) + 1e-6*(end->tv_usec-start->tv_usec);
}

struct Planet {
   double mass;
   double x;
   double y;
   double vx;
   double vy;
};

unsigned long long seed = 100;

const int N = 200;
const int PARALLEL_THRESHOLD = 50;

unsigned long long randomU64() {
  seed ^= (seed << 21);
  seed ^= (seed >> 35);
  seed ^= (seed << 4);
  return seed;
}

double randomDouble()
{
   unsigned long long next = randomU64();
   next >>= (64 - 26);
   unsigned long long next2 = randomU64();
   next2 >>= (64 - 26);
   return ((next << 27) + next2) / (double)(1LL << 53);
}

void simulation_serial(Planet* planets[2], int nplanets, int timesteps, double dt) {
   // Precompute the mass products since mass does not change over time
   double* massProducts = (double*)malloc(sizeof(double) * nplanets * nplanets);

   for (int i = 0; i < nplanets; i++) {
      for (int j = i; j < nplanets; j++) {
         double massProduct = planets[0][i].mass * planets[0][j].mass;
         massProducts[i * nplanets + j] = massProduct * massProduct * massProduct;
      }
   }

   for (int t = 0; t < timesteps; t++) {
      int current = t % 2;
      int next = (t + 1) % 2;
      memcpy(planets[next], planets[current], nplanets * sizeof(Planet));
      for (int i = 0; i < nplanets; i++) {
         for (int j = i; j < nplanets; j++) {
               double dx = planets[current][j].x - planets[current][i].x;
               double dy = planets[current][j].y - planets[current][i].y;
               double distSqr = dx*dx + dy*dy + 0.0001;
               double invDist3 = massProducts[i * nplanets + j] / (distSqr * sqrt(distSqr));
               planets[next][i].vx += dt * dx * invDist3;
               planets[next][i].vy += dt * dy * invDist3;
               if (j != i) {
                  planets[next][j].vx += dt * -dx * invDist3;
                  planets[next][j].vy += dt * -dy * invDist3;
               }
         }
         planets[next][i].x += dt * planets[next][i].vx;
         planets[next][i].y += dt * planets[next][i].vy;
      }

      #if PRINT_STEPS
          if (t < N) {
             printf("x: %0.32g y: %0.32g vx: %0.32g vy: %0.32g\n", 
                    planets[next][nplanets-1].x, 
                    planets[next][nplanets-1].y, 
                    planets[next][nplanets-1].vx, 
                    planets[next][nplanets-1].vy);
          }
      #endif
   }

   free(massProducts);
}

void simulation_parallel(Planet* planets[2], int nplanets, int timesteps, double dt) {
   // Precompute the mass products since mass does not change over time
   double* massProducts = (double*)malloc(sizeof(double) * nplanets * nplanets);
   
   #pragma omp parallel for schedule(dynamic)
   for (int i = 0; i < nplanets; i++) {
      for (int j = i; j < nplanets; j++) {
         double massProduct = planets[0][i].mass * planets[0][j].mass;
         massProducts[i * nplanets + j] = massProduct * massProduct * massProduct;
      }
   }

   int nthreads = omp_get_max_threads();
   std::vector<double> X(nplanets, 0.0);
   std::vector<double> Y(nplanets, 0.0);
   for (int t = 0; t < timesteps; ++t) {
      int current = t % 2;
      int next = (t + 1) % 2;
      memcpy(planets[next], planets[current], nplanets * sizeof(Planet));

      // Allocate thread-local accumulators
      std::vector< std::vector<double> > VX(nthreads, std::vector<double>(nplanets, 0.0));
      std::vector< std::vector<double> > VY(nthreads, std::vector<double>(nplanets, 0.0));

      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         #pragma omp for schedule(dynamic)
         for (int i = 0; i < nplanets; ++i) {
            for (int j = i; j < nplanets; ++j) {
               double dx = planets[current][j].x - planets[current][i].x;
               double dy = planets[current][j].y - planets[current][i].y;
               double distSqr = dx * dx + dy * dy + 0.0001;
               double invDist3 = massProducts[i * nplanets + j] / (distSqr * sqrt(distSqr));
               
               VX[tid][i] += dt * dx * invDist3;
               VY[tid][i] += dt * dy * invDist3;
               if (j != i) {
                  VX[tid][j] += dt * -dx * invDist3;
                  VY[tid][j] += dt * -dy * invDist3;
               }
            }
         }
      }

      // Combine thread-local results
      for (int i = 0; i < nplanets; ++i) {
         X[i] = 0;
         Y[i] = 0;
         for (int tid = 0; tid < nthreads; ++tid) {
            X[i] += VX[tid][i];
            Y[i] += VY[tid][i];
         }
      }

      // Update positions and velocities
      for (int i = 0; i < nplanets; ++i) {
         planets[next][i].vx += X[i];
         planets[next][i].vy += Y[i];
         planets[next][i].x += dt * planets[next][i].vx;
         planets[next][i].y += dt * planets[next][i].vy;
      }

      #if PRINT_STEPS
          if (t < N) {
             printf("x: %0.32g y: %0.32g vx: %0.32g vy: %0.32g\n", 
                    planets[next][nplanets-1].x, 
                    planets[next][nplanets-1].y, 
                    planets[next][nplanets-1].vx, 
                    planets[next][nplanets-1].vy);
          }
      #endif
   }

   free(massProducts);
}


int nplanets;
int timesteps;
double dt;
double G;

int main(int argc, const char** argv){
   if (argc < 2) {
      printf("Usage: %s <nplanets> <timesteps>\n", argv[0]);
      return 1;
   }
   nplanets = atoi(argv[1]);
   timesteps = atoi(argv[2]);
   dt = 0.001;
   G = 6.6743;

   // Allocate two blocks of memory for double buffering
   Planet* planets[2];
   planets[0] = (Planet*)malloc(sizeof(Planet) * nplanets);
   planets[1] = (Planet*)malloc(sizeof(Planet) * nplanets);
   
   // Initialize first buffer
   double scale = pow(1 + nplanets, 0.4);
   for (int i=0; i<nplanets; i++) {
      planets[0][i].mass = randomDouble() * 10 + 0.2;
      planets[0][i].x = ( randomDouble() - 0.5 ) * 100 * scale;
      planets[0][i].y = ( randomDouble() - 0.5 ) * 100 * scale;
      planets[0][i].vx = randomDouble() * 5 - 2.5;
      planets[0][i].vy = randomDouble() * 5 - 2.5;
   }

   #if PROFILE
      ProfilerStart("my_profile.prof");
   #endif

   struct timeval start, end;
   gettimeofday(&start, NULL);

   if (nplanets >= PARALLEL_THRESHOLD) {
      simulation_parallel(planets, nplanets, timesteps, dt);
   } else {
      simulation_serial(planets, nplanets, timesteps, dt);
   }

   gettimeofday(&end, NULL);

   #if PROFILE
      ProfilerStop();
   #endif

   int final = timesteps % 2;
   printf("Total time to run simulation %0.6f seconds, final location %f %f\n", 
          tdiff(&start, &end), 
          planets[final][nplanets-1].x, 
          planets[final][nplanets-1].y);

   free(planets[0]);
   free(planets[1]);

   return 0;   
}