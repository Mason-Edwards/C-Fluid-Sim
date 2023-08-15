#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <time.h>
// Try CMAKE next time around.
#include <SDL.h>

#include "C-Fluid-Sim.h"

// TODO Issues
// 1.Velocity should be in the direction that the mouse is being moved.
// 2.When velocity of squares has dissipated they are no longer reached by density, causing you to need to
// move the mouse to an area to put some velocity there.
//
// As a fix, pressing w or s will add velocity to the whole grid meaning the density acts somewhat normal.

int main(int argc, char *argv[])
{
    Grid* grid = createGrid(SIZE, DIFFUSION, VISCOSITY, DELTATIME);

    for (int i = 0; i < grid->area; i++)
    {
        grid->Vx[i] += 0.5f;
        grid->Vy[i] += 0.5f;
    }

    showGUI(grid);

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
    float density = grid->density[INDEX(x, y)];
    if (density + amount > 50.f)
    {
        grid->density[INDEX(x, y)] = 50.f;
    }
    else
    {
        grid->density[INDEX(x, y)] += amount;
    }
}

void addVelocity(Grid* grid, int x, int y, float amountX, float amountY)
{
    float vx = grid->Vx[INDEX(x, y)];
    float vy = grid->Vx[INDEX(x, y)];

    if (vx + amountX > 0.7f)
    {
        grid->Vx[INDEX(x, y)] = 0.5f;
    }
    else {
        grid->Vx[INDEX(x, y)] += amountX;
    }

    if (vy + amountY > 0.7f)
    {
        grid->Vy[INDEX(x, y)] = 0.7f;
    }
    else {
        grid->Vy[INDEX(x, y)] += amountY;
    }
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
            printf("[%i, %i] = %f\n", i, j, grid->density[INDEX(i, j)]);            
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
    for (int i = 1; i < size - 1; i++) {
        array[INDEX(i, 0)] = arrayType == 2 ? -array[INDEX(i, 1)] : array[INDEX(i, 1)];
        array[INDEX(i, size - 1)] = arrayType == 2 ? -array[INDEX(i, size - 2)] : array[INDEX(i, size - 2)];
    }
    for (int j = 1; j < size - 1; j++) {
        array[INDEX(0, j)] = arrayType == 1 ? -array[INDEX(1, j)] : array[INDEX(1, j)];
        array[INDEX(size - 1, j)] = arrayType == 1 ? -array[INDEX(size - 2, j)] : array[INDEX(size - 2, j)];
    }

    array[INDEX(0, 0)] = 0.5f * (array[INDEX(1, 0)] + array[INDEX(0, 1)]);
    array[INDEX(0, size - 1)] = 0.5f * (array[INDEX(1, size - 1)] + array[INDEX(0, size - 2)]);
    array[INDEX(size - 1, 0)] = 0.5f * (array[INDEX(size - 2, 0)] + array[INDEX(size - 1, 1)]);
    array[INDEX(size - 1, size - 1)] = 0.5f * (array[INDEX(size - 2, size - 1)] + array[INDEX(size - 1, size - 2)]);
    
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
void simStep(Grid* grid)
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

    diffuse(1, vx0, vx, visc, dt, 4, n);
    diffuse(2, vy0, vy, visc, dt, 4, n);

    project(vx0, vy0, vx, vy, 4, n);

    advect(1, vx, vx0, vx0, vy0, dt, n);
    advect(2, vy, vy0, vx0, vy0, dt, n);

    project(vx, vy, vx0, vy0, 4, n);
    diffuse(0, s, density, diff, dt, 4, n);
    advect(0, density, s, vx, vy, dt, n);
}

int showGUI(Grid *grid)
{

    static SDL_MouseMotionEvent prevMouseEvent;
    prevMouseEvent.x = 0;
    prevMouseEvent.y = 0;

    int cellSize = PIXELSIZE;
    int gridWidth = grid->size;
    int gridHeight = grid->size;

    // + 1 so that the last grid lines fit in the screen.
    int windowWidth = (gridWidth * cellSize) + 1;
    int windowHeight = (gridHeight * cellSize) + 1;

    SDL_Rect square = {
        .x = (gridWidth - 1) / 2 * cellSize,
        .y = (gridHeight - 1) / 2 * cellSize,
        .w = cellSize,
        .h = cellSize,
    };

    SDL_Color background = { 22, 22, 22, 255 }; // Barely Black
    SDL_Color lineColour = { 44, 44, 44, 255 }; // Dark grey

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Initialize SDL: %s",
            SDL_GetError());
        return EXIT_FAILURE;
    }

    SDL_Window* window;
    SDL_Renderer* renderer;
    if (SDL_CreateWindowAndRenderer(windowWidth, windowHeight, 0, &window,
        &renderer) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION,
            "Create window and renderer: %s", SDL_GetError());
        return EXIT_FAILURE;
    }

    SDL_SetWindowTitle(window, "Fluid Simulation");
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    SDL_bool quit = SDL_FALSE;
    SDL_bool mouseActive = SDL_FALSE;
    SDL_bool mouseHover = SDL_FALSE;

    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_w:
                    for (int i = 0; i < grid->area; i++)
                    {
                        grid->Vx[i] = 0.0005f;
                        grid->Vy[i] = 0.0005f;
                    }
                    break;
                case SDLK_UP:
                case SDLK_s:
                    for (int i = 0; i < grid->area; i++)
                    {
                        grid->Vx[i] = -0.001f;
                        grid->Vy[i] = -0.000001f;
                    }
                case SDLK_DOWN:
                case SDLK_a:
                case SDLK_LEFT:
                case SDLK_d:
                case SDLK_RIGHT:
                    break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                square.x = (event.motion.x / cellSize);
                square.y = (event.motion.y / cellSize);
                printf("x: %i | y: %i | INDEX: %i : density %f | Vx %f | Vy %f\n",
                    square.x, square.y, INDEX(square.x, square.y),
                    grid->density[INDEX(square.x, square.y)], grid->Vx[INDEX(square.x, square.y)],
                    grid->Vy[INDEX(square.x, square.y)]
                );
                break;
            case SDL_MOUSEMOTION:
                if (!mouseActive)
                    mouseActive = SDL_TRUE;
                    int x = (event.motion.x / cellSize);
                    int y = (event.motion.y / cellSize);
                    square.x = x;
                    square.y = y;

                    /*Sint32 diffX = event.motion.x - prevMouseEvent.x;
                    Sint32 diffY = event.motion.y - prevMouseEvent.y;

                    prevMouseEvent.x = event.motion.x;
                    prevMouseEvent.y = event.motion.y;*/

                    // Randomise velocity?
                    //srand((unsigned)time(NULL));
                    //float vx = -1 + 2 * ((float)rand()) / RAND_MAX;
                    //float vy = -1 + 2 * ((float)rand()) / RAND_MAX;

                    addDensity(grid, square.x, square.y, 0.7f);
                    addVelocity(grid, square.x, square.y, 0.5f, 0.5f);
                break;
            case SDL_WINDOWEVENT:
                if (event.window.event == SDL_WINDOWEVENT_ENTER && !mouseHover)
                    mouseHover = SDL_TRUE;
                else if (event.window.event == SDL_WINDOWEVENT_LEAVE && mouseHover)
                    mouseHover = SDL_FALSE;
                break;
            case SDL_QUIT:
                quit = SDL_TRUE;
                break;
            }
        }

        // Draw grid background.
        SDL_SetRenderDrawColor(renderer, background.r, background.g,
            background.b, background.a);
        SDL_RenderClear(renderer);

        // Draw grid lines.
        SDL_SetRenderDrawColor(renderer, lineColour.r, lineColour.g,
            lineColour.b, lineColour.a);

        for (int x = 0; x < 1 + gridWidth * cellSize;
            x += cellSize) {
            SDL_RenderDrawLine(renderer, x, 0, x, windowWidth);
        }

        for (int y = 0; y < 1 + gridHeight * cellSize;
            y += cellSize) {
            SDL_RenderDrawLine(renderer, 0, y, windowHeight, y);
        }

        // Draw densities
        for (int i = 0; i < grid->size; i++)
        {
            for (int j = 0; j < grid->size; j++)
            {
                //// Fade density
                //grid->density[INDEX(i, j)] -= 0.0001f;
                //if (grid->density[INDEX(i, j)] < 0.0f) 
                //    grid->density[INDEX(i, j)] = 0.0001f;

                drawDensity(grid, renderer, i, j, cellSize);
                //drawVelocityX(grid, renderer, i, j, cellSize);
            }
        }

        SDL_RenderPresent(renderer);
        simStep(grid);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return EXIT_SUCCESS;
}

void drawDensity(Grid *grid, SDL_Renderer *renderer, int i, int j, int cellSize)
{
    SDL_Rect cell = {
        .x = i * cellSize,
        .y = j * cellSize,
        .w = cellSize,
        .h = cellSize,
    };

    Uint8 density = 0;
    float gridDensity = grid->density[INDEX(i, j)];
    if (gridDensity > 0.0f)
    {
        //printf("[%i, %i] = %.9g\n", j, i, gridDensity);
        density = (int)(gridDensity * 255.0f);
        //density = 255;
    }

    SDL_Color colour = { 255, 255, 255, density }; // White
    SDL_SetRenderDrawColor(renderer, colour.r, colour.g, colour.b, colour.a);
    SDL_RenderFillRect(renderer, &cell);
}

// TODO Test velocities
void drawVelocityX(Grid *grid, SDL_Renderer* renderer, int i, int j, int cellSize)
{
    SDL_Rect cell = {
        .x = i * cellSize,
        .y = j * cellSize,
        .w = cellSize,
        .h = cellSize,
    };

    Uint8 velocity = 0;
    float gridVelocity = grid->Vx[INDEX(i, j)];
    if (gridVelocity > 0.0f)
    {
        //printf("[%i, %i] = %.9g\n", j, i, gridVelocity);
        //velocity = (Uint8)(gridVelocity * 1000);
        velocity = 255;
    }

    SDL_Color colour = { 255, 255, 255, velocity }; // White
    SDL_SetRenderDrawColor(renderer, colour.r, colour.g, colour.b, colour.a);
    SDL_RenderFillRect(renderer, &cell);
}