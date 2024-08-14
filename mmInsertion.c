#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>


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
int findNextVertexMinimax(const int *tour, int tourSize, const double **distanceMatrix, int numPoints);


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

int findNextVertexMinimax(const int *tour, int tourSize, const double **distanceMatrix, int numPoints) {
    double minMaxDistance = INFINITY;
    int nextVertex = -1;

    // Iterate through all possible vertices to find the minimax vertex
    for (int i = 0; i < numPoints; i++) {
        bool isInTour = false;
        for (int j = 0; j < tourSize; j++) {
            if (i == tour[j]) {
                isInTour = true;
                break;
            }
        }

        if (!isInTour) {
            double maxDistance = -INFINITY;

            // Find the maximum distance from this vertex to any vertex in the tour
            for (int j = 0; j < tourSize; j++) {
                double dist = distanceMatrix[i][tour[j]];
                if (dist > maxDistance) {
                    maxDistance = dist;
                }
            }

            // Update the minimax candidate if this one is better
            if (maxDistance < minMaxDistance) {
                minMaxDistance = maxDistance;
                nextVertex = i;
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

/*int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    char *filename = argv[1]; // Get the filename from command line arguments
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

        while (1)
    {
        closestPoint = findNextVertexMinimax(currentTour, newTourSize, distanceMatrix, numOfCoords);
        if (closestPoint == -1)
        {
            break;
        }

        int *newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        free(currentTour);     // Free the old tour
        currentTour = newTour; // Update currentTour to point to the new tour
    }
    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++)
    {
        printf("%d ", currentTour[i]);
    }
    printf("\n");


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
*/
/*int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    clock_t startTotal, endTotal, startFunc, endFunc;
    double totalTime, findVertexTotalTime = 0.0, insertCostTotalTime = 0.0;

    startTotal = clock(); // Start the total execution timer

    char *filename = argv[1]; // Get the filename from command line arguments
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

    while (1)
    {
        startFunc = clock(); // Start timing for findNextVertexMinimax
        closestPoint = findNextVertexMinimax(currentTour, newTourSize, (const double **)distanceMatrix, numOfCoords);
        endFunc = clock(); // End timing for findNextVertexMinimax
        findVertexTotalTime += ((double) (endFunc - startFunc)) / CLOCKS_PER_SEC;

        if (closestPoint == -1)
        {
            break;
        }

        startFunc = clock(); // Start timing for insertAtLowestCost
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        endFunc = clock(); // End timing for insertAtLowestCost
        insertCostTotalTime += ((double) (endFunc - startFunc)) / CLOCKS_PER_SEC;

        free(currentTour);     // Free the old tour
        currentTour = newTour; // Update currentTour to point to the new tour
    }

    endTotal = clock(); // End the total execution timer
    totalTime = ((double) (endTotal - startTotal)) / CLOCKS_PER_SEC;

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++)
    {
        printf("%d ", currentTour[i]);
    }
    printf("\n");

    printf("Total time for finding next vertex: %f seconds\n", findVertexTotalTime);
    printf("Total time for inserting at lowest cost: %f seconds\n", insertCostTotalTime);
    printf("Total execution time: %f seconds\n", totalTime);

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
*/

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <input filename> <output filename>\n", argv[0]);
        return 1;
    }

    clock_t startTotal, endTotal, startFunc, endFunc;
    double totalTime, findVertexTotalTime = 0.0, insertCostTotalTime = 0.0;

    startTotal = clock(); // Start the total execution timer

    char *inputFilename = argv[1]; // Get the input filename from command line arguments
    char *outputFilename = argv[2]; // Get the output filename from command line arguments
    int numOfCoords = readNumOfCoords(inputFilename);
    double **coords = readCoords(inputFilename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);

    // Initial tour starts and ends at point 0 (assumed to be the starting point)
    int tour[] = {0, 0};
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int *currentTour = malloc(tourSize * sizeof(int));
    memcpy(currentTour, tour, tourSize * sizeof(int));

    int closestPoint;
    int newTourSize = tourSize;
    int *newTour;

    while (1)
    {
        startFunc = clock(); // Start timing for findNextVertexMinimax
        closestPoint = findNextVertexMinimax(currentTour, newTourSize, (const double **)distanceMatrix, numOfCoords);
        endFunc = clock(); // End timing for findNextVertexMinimax
        findVertexTotalTime += ((double) (endFunc - startFunc)) / CLOCKS_PER_SEC;

        if (closestPoint == -1)
        {
            break;
        }

        startFunc = clock(); // Start timing for insertAtLowestCost
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        endFunc = clock(); // End timing for insertAtLowestCost
        insertCostTotalTime += ((double) (endFunc - startFunc)) / CLOCKS_PER_SEC;

        free(currentTour);     // Free the old tour
        currentTour = newTour; // Update currentTour to point to the new tour
    }

    writeTourToFile(currentTour, newTourSize, outputFilename); // Write tour to the output file

    endTotal = clock(); // End the total execution timer
    totalTime = ((double) (endTotal - startTotal)) / CLOCKS_PER_SEC;

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
