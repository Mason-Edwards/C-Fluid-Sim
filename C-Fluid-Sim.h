#pragma once

// Grid macros
#define SIZE 100
#define PIXELSIZE 10

// Macro for going from an x,y index to a 1d array index.
// x and y are in brackets so that when a calulation gets passed in it doesnt mess with the order of operations.
// E.g INDEX(1-1, 1) turning into (1-1*size+1) which turns into 1-(1*size)+1 which can cause out of bounds.
#define INDEX(x, y) ((x)+(SIZE * y))

// Sim macros
#define DELTATIME 1.0f
#define DIFFUSION 0.00000000001f
#define VISCOSITY 0.0000001f
#define ITTERATION 16

typedef enum ArrayType {
    ArrayType_OTHER = 0,
    ArrayType_Vx = 1,
    ArrayType_Vy = 2
} ArrayType;

typedef struct Grid {
    int size;
    int area;
    float deltaTime;
    float diffusion;
    float viscosity;

    float* scratchSpace;
    float* density;

    float* Vx;
    float* Vy;

    float* prev_Vx;
    float* prev_Vy;
} Grid;

Grid* createGrid(int size, float diffusion, float viscosity, float deltaTime);
void printGrid(Grid* grid);
void freeGrid(Grid* grid);
void printDensities(Grid* grid);
void addDensity(Grid* grid, int x, int y, float amount);
void addVelocity(Grid* grid, int x, int y, float amountX, float amountY);
int showGUI(Grid* grid);
void drawDensity(Grid* grid, SDL_Renderer* renderer, int i, int j, int cellSize);
void drawVelocityX(Grid* grid, SDL_Renderer* renderer, int i, int j, int cellSize);

static void setBoundary(ArrayType arrayType, float* array, int size);
static void solve(ArrayType arrayType, float* array, float* prevArray, float flow, float deltaTime, int iteration, int size);
void diffuse(ArrayType arrayType, float* array, float* prevArray, float flow, float deltaTime, int iteration, int size);
static void project(float* velocX, float* velocY, float* p, float* div, int iter, int N);
static void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt, int N);
void simStep(Grid* grid);