#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "timer.h"

#define pi 3.14159265359
#define THREADS 16

int main(void) {
	int   MaxCount = 1000;
	int   xMin = -2;
	int   xMax = 2;
	int   yMin = -2;
	int   yMax = 2;
	int   steps = 5000; //NOTE: Runs steps^2 steps, represents steps done in real and imaginary axes

	float Tol = .0001;
	double complex  r1 = 1 + 0*I;
	double complex  r2 = -0.5 + sin(2 * pi / 3)*I;
	double complex  r3 = -0.5 - sin(2 * pi / 3)*I;
	int points[4] = { 0 };
	StartTimer();
	
	omp_set_num_threads(THREADS);
	//#pragma omp parallel for default(shared) num_threads(THREADS) reduction(+:points) schedule(static,steps/THREADS)
	#pragma omp parallel for default(shared) num_threads(THREADS) reduction(+:points) schedule(dynamic)
	for (int y = 0; y < steps; y++) {
		for (int x = 0; x < steps; x++) {
			double complex z = (xMin + (xMax - xMin) * 1.0 * x / (steps - 1)) + (yMin + (yMax - yMin) * 1.0 * y / (steps - 1)) * I;
			int  count = 0;
			while ((count < MaxCount) && cabs(z - r1) >= Tol && cabs(z - r2) >= Tol && cabs(z - r3) >= Tol) {
				if (cabs(z) > 0) {
					z = z - (z*z*z - 1.0) / (z*z*3.0); //change fraction to desired function divided by its derivative to change fractal function
				}
				count++;
			}

			if (cabs(z - r1) < Tol && abs(cimag(z)) < Tol) {
				points[1]++;
			} else if (cabs(z - r2) <= Tol && cimag(z) > -Tol) {
				points[2]++;
			} else if (cabs(z - r3) <= Tol && cimag(z) < Tol) {
				points[3]++;
			} else {
				points[0]++;	
			}
		}
	}
	

	// Code below is a version of the above loop without the nested for loops. Ultimately ran worse for me on bridges2. 
	/*
	#pragma omp parallel for default(shared) num_threads(THREADS) reduction(+:points) schedule(dynamic)
	for (int k = 0; k < (steps * steps); k++) {
		int y = k / steps;
		int x = k % steps;
		double complex z = (xMin + (xMax - xMin) * 1.0 * x / (steps - 1)) + (yMin + (yMax - yMin) * 1.0 * y / (steps - 1)) * I;
		int count = 0;
		while ((count < MaxCount) && cabs(z - r1) >= Tol && cabs(z - r2) >= Tol && cabs(z - r3) >= Tol) {
			if (cabs(z) > 0) {
				z = z - (z*z*z - 1.0) / (z*z*3.0);
			}
			count++;
		}

		if (cabs(z - r1) < Tol && abs(cimag(z)) < Tol) {
			points[1]++;
		} else if (cabs(z - r2) <= Tol && cimag(z) > -Tol) {
			points[2]++;
		} else if (cabs(z - r3) <= Tol && cimag(z) < Tol) {
			points[3]++;
		} else {
			points[0]++;
		}
	}
	*/

	double runtime = GetTimer();
	printf("Newton Fractal for %d points:\n", steps * steps);
	printf("Points that converged to no roots : %d (%.2f%%)\n", points[0], 100.0 * points[0] / (steps * steps));
	printf("Points that converged to root %.2f + %.2fi: %d (%.2f%%)\n", creal(r1), cimag(r1), points[1], 100.0 * points[1] / (steps * steps));
	printf("Points that converged to root %.2f + %.2fi: %d (%.2f%%)\n", creal(r2), cimag(r2), points[2], 100.0 * points[2] / (steps * steps));
	printf("Points that converged to root %.2f + %.2fi: %d (%.2f%%)\n", creal(r3), cimag(r3), points[3], 100.0 * points[3] / (steps * steps));
	printf("Time taken: %f s\n", runtime / 1000);
}
