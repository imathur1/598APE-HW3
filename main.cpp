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

void simulation_serial(Planet* planets, int nplanets, int timesteps, double dt) {
   // Precompute the mass products since mass does not change over time
   double* massProducts = (double*)malloc(sizeof(double) * nplanets * nplanets);
   for (int i = 0; i < nplanets; ++i) {
      for (int j = i + 1; j < nplanets; ++j) {
         double massProduct = planets[i].mass * planets[j].mass;
         massProducts[i * nplanets + j] = massProduct * massProduct * massProduct;
      }
   }

   for (int t = 0; t < timesteps; ++t) {
      for (int i = 0; i < nplanets; ++i) {
         for (int j = i + 1; j < nplanets; ++j) {
            double dx = planets[j].x - planets[i].x;
            double dy = planets[j].y - planets[i].y;
            double distSqr = dx*dx + dy*dy + 0.0001;
            double invDist3 = dt * massProducts[i * nplanets + j] / (distSqr * sqrt(distSqr));
            double vxDiff = dx * invDist3;
            double vyDiff = dy * invDist3;
            planets[i].vx += vxDiff;
            planets[i].vy += vyDiff;
            planets[j].vx -= vxDiff;
            planets[j].vy -= vyDiff;
         }
         planets[i].x += dt * planets[i].vx;
         planets[i].y += dt * planets[i].vy;
      }

      #if PRINT_STEPS
          if (t < N) {
             printf("x: %0.32g y: %0.32g vx: %0.32g vy: %0.32g\n", 
                    planets[nplanets-1].x, 
                    planets[nplanets-1].y, 
                    planets[nplanets-1].vx, 
                    planets[nplanets-1].vy);
          }
      #endif
   }

   free(massProducts);
}

void simulation_parallel(Planet* planets, int nplanets, int timesteps, double dt) {
   int chunk_size = sqrt(nplanets) + 1; 
   int nthreads = omp_get_max_threads();

   // Precompute the mass products since mass does not change over time
   double* massProducts = (double*)malloc(sizeof(double) * nplanets * nplanets);
   #pragma omp parallel for schedule(dynamic, chunk_size)
   for (int i = 0; i < nplanets; ++i) {
      for (int j = i + 1; j < nplanets; ++j) {
         double massProduct = planets[i].mass * planets[j].mass;
         massProducts[i * nplanets + j] = massProduct * massProduct * massProduct;
      }
   }

   double** thread_VX = new double*[nthreads];
   double** thread_VY = new double*[nthreads];
   for (int tid = 0; tid < nthreads; tid++) {
      thread_VX[tid] = (double*)malloc(nplanets * sizeof(double));
      thread_VY[tid] = (double*)malloc(nplanets * sizeof(double));
   }

   for (int t = 0; t < timesteps; ++t) {
      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         memset(thread_VX[tid], 0, nplanets * sizeof(double));
         memset(thread_VY[tid], 0, nplanets * sizeof(double));
      }

      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         #pragma omp for schedule(dynamic, chunk_size)
         for (int i = 0; i < nplanets; ++i) {
            for (int j = i + 1; j < nplanets; ++j) {
               double dx = planets[j].x - planets[i].x;
               double dy = planets[j].y - planets[i].y;
               double distSqr = dx * dx + dy * dy + 0.0001;
               double invDist3 = massProducts[i * nplanets + j] / (distSqr * sqrt(distSqr));
               double vxDiff = dx * invDist3;
               double vyDiff = dy * invDist3;
               thread_VX[tid][i] += vxDiff;
               thread_VY[tid][i] += vyDiff;
               thread_VX[tid][j] -= vxDiff;
               thread_VY[tid][j] -= vyDiff;
            }
         }
      }

      // Combine thread-local results
      #pragma omp parallel for schedule(static)
      for (int i = 0; i < nplanets; ++i) {
         double vxTotal = 0.0, vyTotal = 0.0;
         for (int tid = 0; tid < nthreads; ++tid) {
            vxTotal += thread_VX[tid][i];
            vyTotal += thread_VY[tid][i];
         }
         planets[i].vx += dt * vxTotal;
         planets[i].vy += dt * vyTotal;
         planets[i].x += dt * planets[i].vx;
         planets[i].y += dt * planets[i].vy;
      }

      #if PRINT_STEPS
          if (t < N) {
             printf("x: %0.32g y: %0.32g vx: %0.32g vy: %0.32g\n", 
                    planets[nplanets-1].x, 
                    planets[nplanets-1].y, 
                    planets[nplanets-1].vx, 
                    planets[nplanets-1].vy);
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

   Planet* planets = (Planet*)malloc(sizeof(Planet) * nplanets);   
   double scale = pow(1 + nplanets, 0.4);
   for (int i=0; i<nplanets; ++i) {
      planets[i].mass = randomDouble() * 10 + 0.2;
      planets[i].x = ( randomDouble() - 0.5 ) * 100 * scale;
      planets[i].y = ( randomDouble() - 0.5 ) * 100 * scale;
      planets[i].vx = randomDouble() * 5 - 2.5;
      planets[i].vy = randomDouble() * 5 - 2.5;
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

   printf("Total time to run simulation %0.6f seconds, final location %f %f\n", 
         tdiff(&start, &end), 
         planets[nplanets-1].x, 
         planets[nplanets-1].y);
   free(planets);

   return 0;   
}