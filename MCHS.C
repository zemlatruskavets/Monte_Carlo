//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////                                                                  //////
//////                      HARD SPHERE MONTE CARLO                     //////
//////                          Alexander Tkalych                       //////
//////                             05 12 2013                           //////
//////                                                                  //////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//                Include headers and assign global variables               //

//////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679


//////////////////////////////////////////////////////////////////////////////

//                          Initialize all variables                        //

//////////////////////////////////////////////////////////////////////////////


// ARRAYS
double     *particle;
double     *coords;
double     *correlation;
double     *distance;


// CONSTANTS
int        iteration;
int        I;
double     N;
double     L;
double     D;
double     R;
double     n;
double     B;
double     b;
double     S;
double     success;
double     failure;


// FUNCTIONS
void     MonteCarlo(void);
void     WriteCoords(FILE *fp, double *coords);


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//                                                                          //
//                            START THE PROGRAM                             //
//                                                                          //


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc != 7)
    {
        printf("Usage: MC <geometry.xyz> <iterations> <cube length> <number of particles> <reduced density> <number of bins>\n");
        exit(1);
    }
    
    int i;
    
    // Time how long the whole process takes
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
    // Open the geometry file
    FILE *fp = fopen(argv[1], "r");
    
    // Read in the number of iterations
    sscanf(argv[2], "%i", &I);
    
    // Read in the cube length
    sscanf(argv[3], "%lf", &L);
    
    // Read in the number of particles
    sscanf(argv[4], "%lf", &N);
    
    // Read in the reduced density
    sscanf(argv[5], "%lf", &n);
    
    // Read in the number of bins
    sscanf(argv[6], "%lf", &B);
    
    // Calculate the diameter
    D = pow(((6 * n) / (N * pi)), (1 / 3.));
    
    // Calculate the bin size
    b = (3.5 / B);
    
    
    
    
    
    //////////////////////////////////////////////////////////////////////////////
    
    //      Create various arrays using heap memory in which to store data      //
    
    //////////////////////////////////////////////////////////////////////////////
    
    
    // Allocate particle, coordinate, encounter, and distance arrays
    particle    = (double *)  malloc(sizeof(double)     * N);
    coords      = (double *)  malloc(sizeof(double) * 3 * N);
    correlation = (double *)  malloc(sizeof(double)     * B);
    distance    = (double *)  malloc(sizeof(double)     * N);
    
    //////////////////////////////////////////////////////////////////////////////
    
    //          Read data into the various arrays from the geometry file        //
    
    //////////////////////////////////////////////////////////////////////////////
    
    
    // Open the input file and fill out coordinate array
    int EndOfFile;
    EndOfFile = 1;
    
    while (EndOfFile != -1)
    {
        EndOfFile = fscanf(fp, "%lf %lf %lf", &coords[3 * i], &coords[3 * i + 1], &coords[3 * i + 2]);
        i++;
    }
    
    fclose(fp);
    
    
    // Print information about the simulation
    printf("\n\nSIMULATION PARAMETERS\n"
           "--------------------------\n"
           "Number of iterations = %i\n"
           "Box length = %g\n"
           "Number of particles = %g\n"
           "Reduced density = %g\n"
           "Diameter = %lf\n"
           "Number of bins = %g\n"
           "Bin size = %g\n"
           "--------------------------\n\n", I, L, N, n, D, B, b);
    
    
    // Index every particle in the system
    for (i = 0; i < N; i++)
    {
        particle[i] = i;
    }
    
    
    // Print the initial coordinates
    printf("\n\nINITIAL COORDINATES\n"
           "-------------------------------\n");
    
    for (i = 0; i < N; i++)
    {
        printf("%g %g %g %g\n", particle[i] + 1, coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
    }
    
    printf("-------------------------------\n\n");
    
    
    
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
    //  So far, the program has determined the particle number, the diameter,   //
    //  the bin size, and the geometry.  The individual particle coordinates    //
    //  have been placed into arrays based on the particles' indicies and       //
    //  coordinates.                                                            //
    
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
    
    
    //////////////////////////////////////////////////////////////////////////////
    
    //           Now, a file is opened in which the data is output              //
    //           and the simulation begins                                      //
    
    //////////////////////////////////////////////////////////////////////////////
    
    
    
    printf("\n\n..............................................................................................................\n");
    printf("..............................................................................................................\n");
    printf("\n\n                                              S T A R T I N G                                 \n\n");
    printf("..............................................................................................................\n");
    printf("..............................................................................................................\n\n\n");
    
    
    
    // Define the file where the VMD trajectory is printed
    fp = fopen("trajectory.xyz", "w");
    FILE *fpoutput;
    
    
    
    // Move the particles and calculate the pair correlation function
    for (iteration = 0; iteration < I; iteration++)
    {
        // Move the particles randomly and evaluate overlap criteria
        MonteCarlo();
        
        // Record the particle positions at the end of every iteration
        WriteCoords(fp, coords);
    }
    
    
    // Display the number of successful and unsuccessful random moves, along with the acceptance rate
    double rate = (success / (success + failure));
    printf("SUCCESS RATE\n"
           "-----------------------\n"
           "Successes = %g\n"
           "Failures = %g\n"
           "Success rate = %g\n"
           "-----------------------\n\n", success, failure, rate);
    
    
    // Print the bin population in a separate file
    fpoutput = fopen("bins.dat", "w");
    
    for (i = 0; i < B; i++)
    {
        // Calculate the volume of the shell of thickness b
        double S;
        S = 4 * pi * (i * b) * (i * b) * b;
        
        // This is the normalized bin population
        fprintf(fpoutput, "%lf %lf\n", (i * b), (correlation[i] / (N * S * I)));
    }
    
    fclose(fpoutput);
    
    
    // Calculate the amount of time the program took
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("This run took %g minutes (i.e. %g seconds).\n\n", (cpu_time_used / 60.), cpu_time_used);
    
    
    
    printf("\n\n..............................................................................................................\n");
    printf("..............................................................................................................\n");
    printf("\n\n                                                F I N I S H E D                                 \n\n");
    printf("..............................................................................................................\n");
    printf("..............................................................................................................\n\n\n\a");
    
    
    fclose(fp);
    
    
    
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//      This sections contains the various functions called in              //
//                      the body of the program                             //

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////

//          This function selectively moves each particle randomly          //

//////////////////////////////////////////////////////////////////////////////


void MonteCarlo()
{
    int i;
    int j;
    
    // Move every particle in the system
    for (i = 0; i < N; i++)
    {
        int IJ_overlap = 0;
        double Dx, Dy, Dz;
        
        // Set the distance array elements to 0 at the beginning of each loop
        int y;
        for (y = 0; y < N; y++)
        {
            distance[y] = 0;
        }
        
        
        // Move one of the particles by a small, random amount
        Dx = coords[3 * i    ] + ((double)rand() / (double)RAND_MAX * (3) * D);
        Dy = coords[3 * i + 1] + ((double)rand() / (double)RAND_MAX * (3) * D);
        Dz = coords[3 * i + 2] + ((double)rand() / (double)RAND_MAX * (3) * D);
        
        // Enforce periodic boundary conditions
        if (Dx > L)
            Dx = L - Dx;
        
        if (Dy > L)
            Dy = L - Dy;
        
        if (Dz > L)
            Dz = L - Dz;
        
        // Calculate whether or not overlap occurs
        for (j = 0; j < N; j++)
        {
            double dx, dy, dz, r2;
            
            // Do not test whether the particle overlaps with itself
            if (i == j) continue;
            
            // Take account of the periodic boundary conditions when calculating the distances
            dx  = Dx - coords[3 * j    ];
            dy  = Dy - coords[3 * j + 1];
            dz  = Dz - coords[3 * j + 2];
            
            // Ensure the interparticle distance is the smallest possible
            dx = dx - floor(dx + (L / 2.));
            dy = dy - floor(dy + (L / 2.));
            dz = dz - floor(dz + (L / 2.));
            
            // Calculate and store the interparticle distances
            r2 = ((dx * dx) + (dy * dy) + (dz * dz));
            distance[j]  = (sqrt(r2) / D);
            
            // If two particles overlap, set IJ_overlap to a nonzero value
            if (distance[j] < 1)
                ++IJ_overlap;
        }
        
        
        // If the particle does not overlap with any of its neighbours, keep the new coordinates
        if (IJ_overlap == 0)
        {
            ++success;
            
            // Reassign coordinates to the successful random displacement
            coords[3 * i    ] = Dx;
            coords[3 * i + 1] = Dy;
            coords[3 * i + 2] = Dz;
            
            // Put each distance into its appropriate bin
            int t;
            for (t = 0; t < N; t++)
            {
                // Do not store zero distances or the distance between a particle and itself
                if ((i == t) || (distance[t] == 0)) continue;
                int k;
                for (k = 0; k < B; k++)
                {
                    // If the distance is greater than a given value and smaller than the next, it goes in this bin
                    if ((distance[t] > (k * b)) && (distance[t] < ((k + 1) * b)))
                        ++correlation[k];
                }
                
            }
        }
        
        // If IJ_overlap is nonzero, at least one pair of particles overlapped
        else
        {
            ++failure;
            
            for (j = 0; j < N; j++)
            {
                double dx, dy, dz, r2;
                
                // Do not test whether the particle overlaps with itself
                if (i == j) continue;
                
                // Take account of the periodic boundary conditions when calculating the distances
                dx  = coords[3 * i    ] - coords[3 * j    ];
                dy  = coords[3 * i + 1] - coords[3 * j + 1];
                dz  = coords[3 * i + 2] - coords[3 * j + 2];
                
                // Ensure the interparticle distance is the smallest possible
                dx = dx - floor(dx + (L / 2.));
                dy = dy - floor(dy + (L / 2.));
                dz = dz - floor(dz + (L / 2.));
                
                // Calculate and store the interparticle distances
                r2 = ((dx * dx) + (dy * dy) + (dz * dz));
                distance[j]  = (sqrt(r2) / D);
            }
            
            // Put each distance into its appropriate bin
            int t;
            for (t = 0; t < N; t++)
            {
                // Do not store zero distances or the distance between a particle and itself
                if ((i == t) || (distance[t] == 0)) continue;
                int k;
                for (k = 0; k < B; k++)
                {
                    // If the distance is greater than a given value and smaller than the next, it goes in this bin
                    if ((distance[t] > (k * b)) && (distance[t] < ((k + 1) * b)))
                        ++correlation[k];
                }
                
            }
            
        }
    }
}




//////////////////////////////////////////////////////////////////////////////

//       This function stores the particle coordinates after each step      //

//////////////////////////////////////////////////////////////////////////////


void WriteCoords(FILE *fp, double *coords)
{
    // VMD requires the number of particles in each frame be listed at the beginning
    fprintf(fp, "%g\n\n", N);
    int i;
    for (i = 0; i < N; i++)
    {
        // VMD requires an atom ID; in this case, each particle is designated as an oxygen
        fprintf(fp, "O %f %f %f\n", coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
    }
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////                                                                  //////
//////                               DONE                               //////
//////                                                                  //////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
