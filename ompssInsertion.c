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
int findClosestPointToTour(int *tour, int tourSize, double **distanceMatrix, int numPoints);
int *insertAtLowestCost(const int *tour, int tourSize, int newPoint, double **distanceMatrix, int *newTourSize);
bool isPointInTour(int *tour, int tourSize, int point); 


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

int findClosestPointToTour(int *tour, int tourSize, double **distanceMatrix, int numPoints) {
    double minDistance = INFINITY; // Initialize as the maximum possible distance
    int closestPoint = -1;         // Initialize to -1 to indicate no valid point found

    // OpenMP parallel block
    #pragma omp parallel
    {
        double localMinDistance = INFINITY;
        int localClosestPoint = -1;

        // Parallel for loop with reduction on localMinDistance
        #pragma omp for nowait // No implicit barrier at the end of the for loop
        for (int i = 0; i < numPoints; i++) {
            if (!isPointInTour(tour, tourSize, i)) {  // Helper function to check if point is in the tour
                double totalDistance = 0;
                for (int j = 0; j < tourSize - 1; j++) {
                    totalDistance += distanceMatrix[tour[j]][i];
                }

                if (totalDistance < localMinDistance) {
                    localMinDistance = totalDistance;
                    localClosestPoint = i;
                }
            }
        }

        // Critical section to update the global minimum and closest point
        #pragma omp critical
        {
            if (localMinDistance < minDistance) {
                minDistance = localMinDistance;
                closestPoint = localClosestPoint;
            }
        }
    }

    return closestPoint;
}

bool isPointInTour(int *tour, int tourSize, int point) {
    for (int i = 0; i < tourSize; i++) {
        if (tour[i] == point) {
            return true;
        }
    }
    return false;
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

void freeResources(int **tour, double ***distanceMatrix, double ***coords, int numOfCoords) {
    free(*tour);  // Free the tour array

    // Free the distance matrix
    for (int i = 0; i < numOfCoords; i++) {
        free((*distanceMatrix)[i]);
    }
    free(*distanceMatrix);

    // Free the coordinates array
    for (int i = 0; i < numOfCoords; i++) {
        free((*coords)[i]);
    }
    free(*coords);
}



/* int main(int argc, char *argv[]) {

    char *filename = argv[1];
    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);

    int tour[] = {0, 0};
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int *currentTour = malloc(tourSize * sizeof(int));
    memcpy(currentTour, tour, tourSize * sizeof(int));

    int closestPoint;
    int newTourSize = tourSize;
    int *newTour;

    while (1) {

        closestPoint = findClosestPointToTour(currentTour, newTourSize, distanceMatrix, numOfCoords);

        if (closestPoint == -1) {
            break;
        }

        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        free(currentTour);
        currentTour = newTour;
    }

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++) {
        printf("%d ", currentTour[i]);
    }
    printf("\n");


    freeResources(&currentTour, &distanceMatrix, &coords, numOfCoords);

    return 0;
}
*/

/* int main(int argc, char *argv[]) {
    double start, end; // For total execution time
    double total_time_used;

    start = omp_get_wtime(); // Start timer for total execution

    char *filename = argv[1];
    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);

    int tour[] = {0, 0};
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int *currentTour = malloc(tourSize * sizeof(int));
    memcpy(currentTour, tour, tourSize * sizeof(int));

    int closestPoint;
    int newTourSize = tourSize;
    int *newTour;

    double findClosestPointToTourTotalTime = 0.0;
    double insertAtLowestCostTotalTime = 0.0;
    double funcTimeStart, funcTimeEnd; // For function execution times

    while (1) {
        funcTimeStart = omp_get_wtime(); // Start timing for finding closest point
        closestPoint = findClosestPointToTour(currentTour, newTourSize, distanceMatrix, numOfCoords);
        funcTimeEnd = omp_get_wtime(); // End timing
        findClosestPointToTourTotalTime += funcTimeEnd - funcTimeStart;

        if (closestPoint == -1) {
            break;
        }

        funcTimeStart = omp_get_wtime(); // Start timing for insertAtLowestCost
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        funcTimeEnd = omp_get_wtime(); // End timing
        insertAtLowestCostTotalTime += funcTimeEnd - funcTimeStart;

        free(currentTour);
        currentTour = newTour;
    }

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++) {
        printf("%d ", currentTour[i]);
    }
    printf("\n");

    freeResources(&currentTour, &distanceMatrix, &coords, numOfCoords);

    end = omp_get_wtime(); // End timer for total execution
    total_time_used = end - start;

    // Print the total time for each function after the while loop
    printf("Total time for findClosestPointToTour: %f seconds\n", findClosestPointToTourTotalTime);
    printf("Total time for insertAtLowestCost: %f seconds\n", insertAtLowestCostTotalTime);
    printf("Total execution time: %f seconds\n", total_time_used);

    return 0;
}
*/

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input filename> <output filename>\n", argv[0]);
        return 1;
    }

    char *inputFilename = argv[1];   // First argument: input filename
    char *outputFilename = argv[2];  // Second argument: output filename

    double start = omp_get_wtime();  // Start timer for total execution

    int numOfCoords = readNumOfCoords(inputFilename);
    double **coords = readCoords(inputFilename, numOfCoords);
    double **distanceMatrix = createDistanceMatrix(coords, numOfCoords);

    int tour[] = {0, 0};  // Initial tour
    int tourSize = sizeof(tour) / sizeof(tour[0]);

    int *currentTour = malloc(tourSize * sizeof(int));
    memcpy(currentTour, tour, tourSize);

    int closestPoint;
    int newTourSize = tourSize;
    int *newTour;

    double funcTimeStart, funcTimeEnd;
    double findClosestPointTotalTime = 0.0;
    double insertAtLowestCostTotalTime = 0.0;

    while (1) {
        funcTimeStart = omp_get_wtime();  // Start timing for finding closest point
        closestPoint = findClosestPointToTour(currentTour, newTourSize, distanceMatrix, numOfCoords);
        funcTimeEnd = omp_get_wtime();  // End timing
        findClosestPointTotalTime += funcTimeEnd - funcTimeStart;

        if (closestPoint == -1) {
            break;
        }

        funcTimeStart = omp_get_wtime();  // Start timing for insertAtLowestCost
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        funcTimeEnd = omp_get_wtime();  // End timing
        insertAtLowestCostTotalTime += funcTimeEnd - funcTimeStart;

        free(currentTour);
        currentTour = newTour;
    }

    // Writing the tour to the output file instead of printing it to the console
    writeTourToFile(currentTour, newTourSize, outputFilename);

    freeResources(&currentTour, &distanceMatrix, &coords, numOfCoords);

    double totalExecutionTime = omp_get_wtime() - start;  // Calculate total execution time

    printf("Total time for finding closest point: %f seconds\n", findClosestPointTotalTime);
    printf("Total time for inserting at lowest cost: %f seconds\n", insertAtLowestCostTotalTime);
    printf("Total execution time: %f seconds\n", totalExecutionTime);

    return 0;
}