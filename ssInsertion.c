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
int findClosestPointToTour(int *tour, int tourSize, double **distanceMatrix, int numPoints);
int *insertAtLowestCost(const int *tour, int tourSize, int newPoint, double **distanceMatrix, int *newTourSize);

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

int findClosestPointToTour(int *tour, int tourSize, double **distanceMatrix, int numPoints)
{
    double minDistance = INFINITY; // Use INFINITY defined in math.h to initialize as the maximum possible
    int closestPoint = -1;         // Initialize to -1 to indicate no point found if that remains the case

    for (int i = 0; i < numPoints; i++)
    {
        // Check if this point is already in the tour
        bool inTour = false;
        for (int j = 0; j < tourSize; j++)
        {
            if (tour[j] == i)
            {
                inTour = true;
                break;
            }
        }

        if (!inTour)
        {
            double totalDistance = 0;
            for (int j = 0; j < (tourSize - 1); j++)
            {
                totalDistance += distanceMatrix[tour[j]][i];
            }

            if (totalDistance < minDistance)
            {
                minDistance = totalDistance;
                closestPoint = i;
            }
        }
    }

    return closestPoint;
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
/*int main(int argc, char *argv[]) {
    clock_t start, end; // For total execution time
    double cpu_time_used;
    
    start = clock(); // Start timer for total execution

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
    clock_t funcTimeStart, funcTimeEnd; // For function execution times

    while (1) {
        funcTimeStart = clock(); // Start timing for finding closest point
        closestPoint = findClosestPointToTour(currentTour, newTourSize, distanceMatrix, numOfCoords);
        funcTimeEnd = clock(); // End timing
        findClosestPointToTourTotalTime += (double) (funcTimeEnd - funcTimeStart) / CLOCKS_PER_SEC;

        if (closestPoint == -1) {
            break;
        }

        funcTimeStart = clock(); // Start timing for insertAtLowestCost
        newTour = insertAtLowestCost(currentTour, newTourSize, closestPoint, distanceMatrix, &newTourSize);
        funcTimeEnd = clock(); // End timing
        insertAtLowestCostTotalTime += (double) (funcTimeEnd - funcTimeStart) / CLOCKS_PER_SEC;

        free(currentTour);
        currentTour = newTour;
    }

    printf("New Tour: ");
    for (int i = 0; i < newTourSize; i++) {
        printf("%d ", currentTour[i]);
    }
    printf("\n");

    freeResources(&currentTour, &distanceMatrix, &coords, numOfCoords);

    end = clock(); // End timer for total execution
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    // Print the total time for each function after the while loop
    printf("Total time for findClosestPointToTour: %f seconds\n", findClosestPointToTourTotalTime);
    printf("Total time for insertAtLowestCost: %f seconds\n", insertAtLowestCostTotalTime);
    printf("Total execution time: %f seconds\n", cpu_time_used);

    return 0;
}
*/

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input filename> <output filename>\n", argv[0]);
        return 1;
    }

    char *inputFilename = argv[1];
    char *outputFilename = argv[2];

    int numOfCoords = readNumOfCoords(inputFilename);
    double **coords = readCoords(inputFilename, numOfCoords);
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

    // Write the final tour to the output file instead of printing it
    writeTourToFile(currentTour, newTourSize, outputFilename);

    // Free all resources
    freeResources(&currentTour, &distanceMatrix, &coords, numOfCoords);

    return 0;
}