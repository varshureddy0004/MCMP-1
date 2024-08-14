#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

/*This code is for reading and writing to files for the 2023-24 COMP528 CA1*/

/*Use the functions in this file to read from the input file, and write to the output file*/

/*You should use this file when compiling your code*/

/*Declare these functions at the top of each 'main' file*/

/*If there are any issues with this code, please contact: h.j.forbes@liverpool.ac.uk*/

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
double euclideanDistance(double x1, double y1, double x2, double y2);
double **createDistanceMatrix(double **coords, int numOfCoords);
void printMatrix(double **matrix, int size);
int *insertAtLowestCost(const int *tour, int tourSize, int newPoint, double **distanceMatrix, int *newTourSize);
int findNextVertexMinimax(const int *tour, int tourSize, double **distanceMatrix, int numPoints);

int *insertAtLowestCost(const int *tour, int tourSize, int newPoint, double **distanceMatrix, int *newTourSize)
{
    *newTourSize = tourSize + 1; // New size of the tour
    int *newTour = malloc(*newTourSize * sizeof(int));
    if (!newTour)
        return NULL; // Memory allocation failed

    int bestPosition = 0;
    double bestCostIncrease = INFINITY;

    // Calculate the best position to insert the new point
    for (int i = 0; i < tourSize - 1; i++)
    {
        double costIfInserted = distanceMatrix[tour[i]][newPoint] + distanceMatrix[newPoint][tour[i + 1]] - distanceMatrix[tour[i]][tour[i + 1]];
        if (costIfInserted < bestCostIncrease)
        {
            bestCostIncrease = costIfInserted;
            bestPosition = i + 1;
        }
    }

    // Copy the tour to the new tour up to the best position
    for (int i = 0; i < bestPosition; i++)
    {
        newTour[i] = tour[i];
    }

    // Insert the new point
    newTour[bestPosition] = newPoint;

    // Copy the remaining part of the tour
    for (int i = bestPosition; i < tourSize; i++)
    {
        newTour[i + 1] = tour[i];
    }

    return newTour;
}

int findNextVertexMinimax(const int *tour, int tourSize, double **distanceMatrix, int numPoints)
{
    double minMaxDistance = INFINITY;
    int nextVertex = -1;

#pragma omp parallel
    {
        double localMinMaxDistance = INFINITY;
        int localNextVertex = -1;

#pragma omp for nowait
        for (int i = 0; i < numPoints; i++)
        {
            bool isInTour = false;
            for (int j = 0; j < tourSize; j++)
            {
                if (i == tour[j])
                {
                    isInTour = true;
                    break;
                }
            }

            if (!isInTour)
            {
                double maxDistance = -INFINITY;

                // Find the maximum distance from this vertex to any vertex in the tour
                for (int j = 0; j < tourSize; j++)
                {
                    double dist = distanceMatrix[i][tour[j]];
                    if (dist > maxDistance)
                    {
                        maxDistance = dist;
                    }
                }

                // Update the minimax candidate if this one is better
                if (maxDistance < localMinMaxDistance)
                {
                    localMinMaxDistance = maxDistance;
                    localNextVertex = i;
                }
            }
        }

// Critical section to update the global minimum of the maximized distances
#pragma omp critical
        {
            if (localMinMaxDistance < minMaxDistance)
            {
                minMaxDistance = localMinMaxDistance;
                nextVertex = localNextVertex;
            }
        }
    }

    return nextVertex;
}

double euclideanDistance(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double **createDistanceMatrix(double **coords, int numOfCoords)
{
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++)
    {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    for (int i = 0; i < numOfCoords; i++)
    {
        for (int j = 0; j < numOfCoords; j++)
        {
            matrix[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
        }
    }

    return matrix;
}

void printMatrix(double **matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfCoords(char *filename)
{
    FILE *file = fopen(filename, "r");
    int numOfCoords = 0;

    if (file == NULL)
    {
        return -1;
    }

    char line[100];

    while (fgets(line, sizeof(line), file) != NULL)
    {
        numOfCoords++;
    }

    return numOfCoords;
}

/*Gets the data from the file. Returns as an array of doubles, Ignores the first integer*/
double **readCoords(char *filename, int numOfCoords)
{
    FILE *file = fopen(filename, "r");
    int i;

    char line[100];

    if (file == NULL)
    {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

    double **coords = (double **)malloc(numOfCoords * sizeof(double *));

    for (i = 0; i < numOfCoords; i++)
    {
        coords[i] = (double *)malloc(2 * sizeof(int));
        if (coords[i] == NULL)
        {
            perror("Memory Allocation Failed");
        }
    }

    int lineNum = 0;
    while (fgets(line, sizeof(line), file) != NULL)
    {
        double x, y;
        if (sscanf(line, "%lf,%lf", &x, &y) == 2)
        {
            coords[lineNum][0] = x;
            coords[lineNum][1] = y;
            lineNum++;
        }
    }

    return coords;
}

void *writeTourToFile(int *tour, int tourLength, char *filename)
{

    FILE *file = fopen(filename, "w");
    int i;

    if (file == NULL)
    {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

    fprintf(file, "%d \n", tourLength);

    printf("Writing output data\n");
    for (i = 0; i < tourLength; i++)
    {
        fprintf(file, "%d ", tour[i]);
    }
}

/*int main(int argc, char *argv[])
{

    char *filename = argv[1]; // Get the filename from command line arguments
    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);
    int tour[] = {0, 1,2,3,4,5,6,7,8};
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int closestPoint = findClosestPointToTour(tour, tourSize, distanceMatrix, numOfCoords);
    printf("Closest point not in the tour is: %d\n", closestPoint);
    int newPoint = 4; // The new point to insert

    int newTourSize;
    int *newTour = insertAtLowestCost(tour, tourSize, newPoint, distanceMatrix, &newTourSize);

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++)
    {
        printf("%d ", newTour[i]);
    }
    printf("\n");

    free(newTour);
    // Free allocated memory

    for (int i = 0; i < numOfCoords; i++)
    {
        free(coords[i]);
        free(distanceMatrix[i]);
    }

    free(coords);
    free(distanceMatrix);

    return 0;
}*/

int main(int argc, char *argv[])
{

    char *filename = argv[1];       // Get the filename from command line arguments
    char *outputFilename = argv[2]; // Second argument: output filename

    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);

    // Initial tour starts and ends at point 0 (assumed to be the starting point)
    int tour[] = {0, 0};
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int *currentTour = malloc(tourSize * sizeof(int));
    memcpy(currentTour, tour, tourSize * sizeof(int));

    int closestPoint;
    int newTourSize = tourSize;
    int *newTour;

    double totalExecutionStart = omp_get_wtime();
    double totalFindVertexTime = 0.0, totalInsertCostTime = 0.0;

    while (1)
    {
        double funcStart = omp_get_wtime();
        closestPoint = findNextVertexMinimax(currentTour, newTourSize, distanceMatrix, numOfCoords);
        double funcEnd = omp_get_wtime();
        totalFindVertexTime += funcEnd - funcStart;

        if (closestPoint == -1)
        {
            break;
        }

        funcStart = omp_get_wtime();
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        funcEnd = omp_get_wtime();
        totalInsertCostTime += funcEnd - funcStart;

        free(currentTour);     // Free the old tour
        currentTour = newTour; // Update currentTour to point to the new tour
    }
    double totalExecutionEnd = omp_get_wtime();

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++)
    {
        printf("%d ", currentTour[i]);
    }
    printf("\n");
    writeTourToFile(currentTour, newTourSize, outputFilename);

    printf("Total time for finding next vertex: %f seconds\n", totalFindVertexTime);
    printf("Total time for inserting at lowest cost: %f seconds\n", totalInsertCostTime);
    printf("Total execution time: %f seconds\n", totalExecutionEnd - totalExecutionStart);

    free(currentTour); // Clean up dynamically allocated memory for the tour
    // Clean up distance matrix and coords
    for (int i = 0; i < numOfCoords; i++)
    {
        free(distanceMatrix[i]);
        free(coords[i]);
    }
    free(distanceMatrix);
    free(coords);

    return 0;
}
