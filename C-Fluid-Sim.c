#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "C-Fluid-Sim.h"

// TODO Fix warnings

int main(void)
{
    Grid* grid = createGrid(SIZE, DIFFUSION, VISCOSITY, DELTATIME);



    //runSim(grid);

    freeGrid(grid);
    return 0;
}

Grid* createGrid(int size, float diffusion, float viscosity, float deltaTime)
{
    Grid* grid = malloc(sizeof(*grid));
    int area = size * size;

    grid->size = size;
    grid->area = area;
    grid->diffusion = diffusion;
    grid->viscosity = viscosity;
    grid->deltaTime = deltaTime;

    grid->scratchSpace = calloc(area, sizeof(*grid->scratchSpace));
    grid->density = calloc(area, sizeof(*grid->density));

    grid->Vx = calloc(area, sizeof(*grid->Vx));
    grid->Vy = calloc(area, sizeof(*grid->Vy));

    grid->prev_Vx = calloc(area, sizeof(*grid->prev_Vx));
    grid->prev_Vy = calloc(area, sizeof(*grid->prev_Vy));

    return grid;
}

void freeGrid(Grid* grid)
{
    free(grid->scratchSpace);
    free(grid->density);
    free(grid->Vx);
    free(grid->Vy);
    free(grid->prev_Vx);
    free(grid->prev_Vy);

    free(grid);
}

void addDensity(Grid* grid, int x, int y, float amount)
{
    grid->density[INDEX(x, y)] += amount;
}

void addVelocity(Grid* grid, int x, int y, float amountX, float amountY)
{
    grid->Vx[INDEX(x, y)] += amountX;
    grid->Vy[INDEX(x, y)] += amountY;
}

void printGrid(Grid* grid)
{
    for (int i = 0; i < grid->area; i++)
    {
        if (i % grid->size == 0)
        {
            printf("\n");
        }
        else if (grid->density[i] > 0.0f)
        {
            printf("#");
        }
        else
        {
            printf(".");
        }
    }
}

void printDensities(Grid* grid)
{
    for (int i = 1; i < grid->size; i++)
    {
        for (int j = 1; j < grid->size; j++)
        {
            printf("[%i, %i] = %f\n", i, j, grid->density[INDEX(j, i)]);            
        }
    }
}

/*
* Defuse a force though the grid.
* @param arrayType An enum representing which array has been passed in (E.g. Vx array, Vy array, etc.).
* @param array Array representing the current array.
* @param prevArray Array representing the previous values of the array passed in.
* @param flow Grid diffusion or viscosity value.
* @param deltaTime Grid deltaTime value.
* @param iteration Grids iteration value.
* @param size Grids size value.
*/
void diffuse(ArrayType arrayType, float* array, float* prevArray, float flow, float deltaTime, int iteration, int size)
{
    // TODO Figure out this?
    float f = deltaTime * flow * (size - 2) * (size - 2);

    // TODO 7?
    solve(arrayType, array, prevArray, f, 1 + 6 * f, iteration, size);
}

static void solve(ArrayType arrayType, float* array, float* prevArray, float flow, float deltaTime, int iteration, int size)
{
    // TODO Figure out this?
    float cRecip = 1.0f / deltaTime;

    for (int j = 1; j < size - 1; j++) {
        for (int i = 1; i < size - 1; i++) {
            array[INDEX(i,j)] =
                (prevArray[INDEX(i, j)]
                    + flow * (array[INDEX(i + 1, j)]
                        + array[INDEX(i - 1, j)]
                        + array[INDEX(i, j + 1)]
                        + array[INDEX(i, j - 1)]
                        + array[INDEX(i, j)]
                        + array[INDEX(i, j)]
                        )) * cRecip;
        }
    }
    setBoundary(arrayType, array, size);
}

/*
* Sets the boundaries of the simulation. This is to ensure that fluid doesn't "leak" out of the grid.
* @param arrayType An enum representing which array has been passed in (E.g. Vx array, Vy array, etc.).
* @param grid Array representing the grid for a certain force (E.g. Vx array, Vy array, etc.).
* @param size The size of the array on a given axis (E.g. If the grid is 20x20, then the size is 20).
*/
static void setBoundary(ArrayType arrayType, float* array, int size)
{
    // Loop though the top and bottom row cells of the grid, excluding the corners.
    for(int i = 1; i < size - 1; i++)
    {
        array[INDEX(i, 0)] = arrayType == ArrayType_Vy ? -array[INDEX(i, 1)] : array[INDEX(i, 1)];
        array[INDEX(i, size - 1)] = arrayType == ArrayType_Vy ? -array[INDEX(i, size - 2)] : array[INDEX(i, size - 2)];
    }

    // Loop though the left and right columns cells of the grid, excluding the corners.
    for(int i = 1; i < size - 1; i++)
    {
        array[INDEX(0, i)] = arrayType == ArrayType_Vx ? -array[INDEX(1, i)] : array[INDEX(1, i)];
        array[INDEX(size - 1, i)] = arrayType == ArrayType_Vx ? -array[INDEX(size - 2, i)] : array[INDEX(size - 2, i)];
    }

    // Set the corners manually as they behave differently.
    // Top left
    array[INDEX(0, 0)] = 0.33f * (array[INDEX(1, 0)]
        + array[INDEX(0, 1)]);
    // Bottom left
    array[INDEX(0, size - 1)] = 0.33f * (array[INDEX(1, size - 1)]
        + array[INDEX(0, size - 2)]);
    
    // Top right
    array[INDEX(size - 1, 0)] = 0.33f * (array[INDEX(size - 2, 0)]
        + array[INDEX(size - 1, 1)]);

    // Bottom right
    array[INDEX(size - 1, size - 1)] = 0.33f * (array[INDEX(size - 2, size - 1)]
        + array[INDEX(size - 1, size - 2)]);
    
}

static void project(float* velocX, float* velocY, float* p, float* div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[INDEX(i, j)] = -0.5f * (
                velocX[INDEX(i + 1, j)]
                - velocX[INDEX(i - 1, j)]
                + velocY[INDEX(i, j + 1)]
                - velocY[INDEX(i, j - 1)]
                ) / N;
            p[INDEX(i, j)] = 0;
        }
    }
    setBoundary(0, div, N);
    setBoundary(0, p, N);
    solve(0, p, div, 1, 6, iter, N);


    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[INDEX(i, j)] -= 0.5f * (p[INDEX(i + 1, j)]
                - p[INDEX(i - 1, j)]) * N;
            velocY[INDEX(i, j)] -= 0.5f * (p[INDEX(i, j + 1)]
                - p[INDEX(i, j - 1)]) * N;
        }
    }
    
    setBoundary(1, velocX, N);
    setBoundary(2, velocY, N);
}

static void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt, int N)
{
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = (float)N;
    float ifloat, jfloat;
    int i, j;


    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            // TODO Why is this returning shit
            tmp1 = dtx * velocX[INDEX(i, j)];
            tmp2 = dty * velocY[INDEX(i, j)];

            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = floorf(x);
            i1 = i0 + 1.0f;

            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;

            t1 = y - j0;
            t0 = 1.0f - t1;


            int i0i = (int)i0;
            int i1i = (int)i1;

            int j0i = (int)j0;
            int j1i = (int)j1;

            /*d[INDEX(i, j)] =
                s0 * (t0 * (d0[INDEX(i0i, j0i)]
                    + d0[INDEX(i0i, j0i)])
                    + (t1 * (d0[INDEX(i0i, j1i)]
                        + d0[INDEX(i0i, j1i)])))
                + s1 * (t0 * (d0[INDEX(i1i, j0i)]
                    + d0[INDEX(i1i, j0i)])
                    + (t1 * (d0[INDEX(i1i, j1i)]
                        + d0[INDEX(i1i, j1i)])));*/

            d[INDEX(i, j)] =
                s0 * (t0 * d0[INDEX(i0i, j0i)] + t1 * d0[INDEX(i0i, j1i)]) +
                s1 * (t0 * d0[INDEX(i1i, j0i)] + t1 * d0[INDEX(i1i, j1i)]);
        }
    }
    
    setBoundary(b, d, N);
}

// Sim steps:
//Diffuse all three velocity components.
//Fix up velocities so they keep things incompressible
//Move the velocities around according to the velocities of the fluid(confused yet ? ).
//Fix up the velocities again
//Diffuse the dye.
//Move the dye around according to the velocities.
// 
void runSim(Grid* grid)
{
    int n = grid->size;
    float visc = grid->viscosity;
    float diff = grid->diffusion;
    float dt = grid->deltaTime;
    float* vx = grid->Vx;
    float* vy = grid->Vy;
    float* vx0 = grid->prev_Vx;
    float* vy0 = grid->prev_Vy;
    float* s = grid->scratchSpace;
    float* density = grid->density;

    diffuse(1, vx, vx0, visc, dt, 4, n);
    diffuse(2, vy, vy0, visc, dt, 4, n);

    project(vx, vy, vx0, vy0, 4, n);

    advect(1, vx, vx0, vx0, vy0, dt, n);
    advect(2, vy, vy0, vx0, vy0, dt, n);

    project(vx, vy, vx0, vy0, 4, n);

    diffuse(0, s, density, diff, dt, 4, n);
    advect(0, density, s, vx, vy, dt, n);
}