#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define timeTotal 100.0


// Run Velocity Verlet algorithm 
void velocityVerlet(double posStart, double velStart, double timeStep,
                    double time[], double pos[], double vel[],
                    double totalEnergy[], int size, bool useTime) {

  double currPos;
  double currVel;
  double currAcc;

  double newAcc;

  double currTotalEnergy;

  // Set initial conditions
  pos[0] = posStart;
  vel[0] = velStart;
  totalEnergy[0] = 0.5 * vel[0] * vel[0] + 0.5 * pos[0] * pos[0];
 
  currPos = pos[0];
  currVel = vel[0];
  if (useTime) {
    currAcc = -1.0 * (pos[0] * cos(time[0]) + vel[0] * sin(time[0]));
  }
  else {
    currAcc = -1.0 * currPos;
  }
  // Propagate Velocity Verlet over time
  for (int i = 1; i < size; i++) {
    // Obtain position x(t)
    currPos = currPos + currVel * timeStep + 0.5 * currAcc * timeStep * timeStep;

    // Obtain acceleration x''(t + dt)
    if (useTime) {
      newAcc = -1.0 * (pos[0] * cos(time[i - 1]) + vel[0] * sin(time[i - 1]));
    }
    else {
      newAcc = -1.0 * currPos;
    }
    

    // Obtain velocty x'(t + dt)
    if (useTime) {
      currVel = currVel + currAcc * timeStep;
    }
    else {
      currVel = currVel + 0.5 * (currAcc + newAcc) * timeStep;
    }

    // Update position, velocity arrays
    pos[i] = currPos;
    vel[i] = currVel;
    currAcc = newAcc;

    // Obtain total energy of system
    currTotalEnergy = 0.5 * currVel * currVel + 0.5 * currPos * currPos;
    totalEnergy[i] = currTotalEnergy;
  }
}

void gnuplotData(char * dataFilename, char * dataFilenameT, char * dataFilenameR,
                 double velStart, double meanE, double maxP, double maxE) {

  char plotName[500];
  char * extPos;

  strcpy(plotName, dataFilename);  
  extPos = strrchr(plotName, '.');
  strcpy(extPos, ".png");

  char ylimP[100];
  if (maxP > 10.0) {
    strcpy(ylimP, "unset yrange; ");
  }
  else {
    snprintf(ylimP, sizeof(ylimP), "set yrange[%g:%g]; ", -10.0, 10.0);
  }

  char ylimE[100];
  if (maxE > 32.5) {
    strcpy(ylimE, "unset yrange; ");
  }
  else {
    snprintf(ylimE, sizeof(ylimE), "set yrange[%g:%g]; ", meanE - 0.5, meanE + 0.5);
  }

  char plotCmd[1000] = "gnuplot -p -e \"set terminal png size 1500,1200;";
  char exactFunction[100];
  strcat(plotCmd, " set datafile separator '\\t';");
  strcat(plotCmd, " set output '");
  strcat(plotCmd, plotName);
  strcat(plotCmd, "'; set multiplot layout 2,1; ");
  strcat(plotCmd, "set samples 1000;");
  strcat(plotCmd, ylimP);
  strcat(plotCmd, "set xlabel 'Time'; set ylabel 'x(t)'; ");

  strcat(plotCmd, "plot '");
  strcat(plotCmd, dataFilename);
  strcat(plotCmd, "' using 1:2 title 'Position' with lines lc 'blue' lw 6, '");
  strcat(plotCmd, dataFilenameR);
  strcat(plotCmd, "' using 1:2 title 'Time-reversed' with lines lc 'goldenrod' lw 4, ");
  snprintf(exactFunction, sizeof(exactFunction), "%g * sin(x) ", velStart);
  strcat(plotCmd, exactFunction);
  strcat(plotCmd, "title 'Exact Solution' with lines lc 'black' lw 1, '");
  strcat(plotCmd, dataFilenameT);
  strcat(plotCmd, "' using 1:2 title 'Second-order Taylor' with lines lc 'dark-red' lw 1; ");
  strcat(plotCmd, ylimE);
  strcat(plotCmd, "set ylabel 'E'; ");
  strcat(plotCmd, "plot '");
  strcat(plotCmd, dataFilename);
  strcat(plotCmd, "' using 1:4 title 'Total Energy' with lines lc 'red' lw 6, '");
  strcat(plotCmd, dataFilenameT);
  strcat(plotCmd, "' using 1:4 title 'Second-order Taylor' with lines lc 'dark-red' lw 2;\"\n");

  printf("Plotting data ...\n");

  FILE * gnuplot = popen(plotCmd, "w");
  fflush(gnuplot);
  pclose(gnuplot);
}

void dumpData(double * time, double * pos, double * vel, double * totalEnergy, char * fname, int size) {
  FILE * dataFile = fopen(fname, "w");
  fprintf(dataFile, "# %s\t%s\t%s\t%s\n", "Time", "Position", "Velocity",
                                          "Total Energy");
  for (int i = 0; i < size; i++) {
    fprintf(dataFile, "%0.4f\t%0.4f\t%0.4f\t%0.4f\n", time[i], pos[i], vel[i],
                                          totalEnergy[i]);
  }
  fclose(dataFile);
}

void makeFilename(char * fname, bool reverse, bool useTime,
                  double timeStep, double velocity) {
  snprintf(fname, 34, "dt%g", timeStep);
  strcat(fname, "_");
  char tmp[33];
  snprintf(tmp, 33, "v%g", velocity);
  strcat(fname, tmp);
  if (reverse) {
    strcat(fname, "_rev");
  }
  if (useTime) {
    strcat(fname, "_T");
  }
  strcat(fname, ".dat");
} 

int main(int argc, char ** argv) {

  int numSteps;
  int numStepsR;
  double currStep;
  double maxP;
  double maxE;
  double meanE;

  double posStart = 0.0;
  double velStart[4] = {1, 2, 4, 8};
  const double timeStep[6] = {0.0002, 0.001, 0.01, 1, 2, 4};

  double posStartRev;
  double velStartRev;

  double * time;
  double * pos;
  double * vel;
  double * totalEnergy;

  double * timeR;
  double * posR;
  double * velR;
  double * totalEnergyR;

  double * timeT;
  double * posT;
  double * velT;
  double * totalEnergyT;

  // Loop over set of time step sizes to simulate
  for (int i = 0; i < 6; i++) {

    currStep = timeStep[i];
    numSteps = timeTotal / currStep;
    numStepsR = numSteps / 2;

    time = malloc(sizeof(double)*numSteps);
    pos  = malloc(sizeof(double)*numSteps);
    vel  = malloc(sizeof(double)*numSteps);
    totalEnergy = malloc(sizeof(double)*numSteps);

    posT = malloc(sizeof(double)*numSteps);
    velT = malloc(sizeof(double)*numSteps);
    totalEnergyT = malloc(sizeof(double)*numSteps);

    timeR = malloc(sizeof(double)*numStepsR);
    posR = malloc(sizeof(double)*numStepsR);
    velR = malloc(sizeof(double)*numStepsR);
    totalEnergyR = malloc(sizeof(double)*numStepsR);

    
    for (int k = 0; k < numSteps; k++) {
      time[k]  = k * currStep;
    }

    for (int k = 0; k < numSteps / 2; k++) {
      timeR[k] = (numStepsR - k) * currStep;
    }

    // Loop over set of initial velocities to simulate
    for (int j = 0; j < 4; j++) {
      posStart = 0.0;

      // Run forward Velocity Verlet simulation
      velocityVerlet(posStart, velStart[j], currStep, time, pos, vel, totalEnergy,
                     numSteps, false);

      // Run forward Velocity Verlet simulation based on Taylor series wrt time
      velocityVerlet(posStart, velStart[j], currStep, time, posT, velT, totalEnergyT,
                     numSteps, true);

      // Set initial conditions for reverse-time trajectory of Velocity Verlet simulation
      posStartRev = pos[numStepsR];
      velStartRev = -1.0 * vel[numStepsR];
      
      // Run reverse-time trajectory of Velocity Verlet simulation
      velocityVerlet(posStartRev, velStartRev, currStep, timeR, posR, velR, totalEnergyR,
                     numStepsR, false);

      maxP = 0.0;
      maxE = 0.0;
      meanE = 0.0;
      for (int k = 0; k < numSteps; k++) {
        if (totalEnergyT[k] > maxE) {
          maxE = totalEnergyT[k];
        }
        if (posT[k] > maxP) {
          maxP = posT[k];
        }
        meanE = meanE + totalEnergyT[k];
      }
      meanE = meanE / numSteps;

      
      char fname[500];
      makeFilename(fname, false, false, currStep, velStart[j]);
      dumpData(time, pos, vel, totalEnergy, fname, numSteps);

      char fnameT[500];
      makeFilename(fnameT, false, true, currStep, velStart[j]);
      dumpData(time, posT, velT, totalEnergyT, fnameT, numSteps);

      char fnameR[500];
      makeFilename(fnameR, true, false, currStep, velStart[j]);
      dumpData(timeR, posR, velR, totalEnergyR, fnameR, numStepsR);

      gnuplotData(fname, fnameT, fnameR, velStart[j], meanE, maxP, maxE);
    }
  }
}
