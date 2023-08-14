#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Try CMAKE next time around.
#include <SDL.h>

#include "C-Fluid-Sim.h"

int main(int argc, char *argv[])
{
    Grid* grid = createGrid(SIZE, DIFFUSION, VISCOSITY, DELTATIME);

    showGUI();

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

int showGUI()
{
    int grid_cell_size = 36;
    int grid_width = 29;
    int grid_height = 23;

    // + 1 so that the last grid lines fit in the screen.
    int window_width = (grid_width * grid_cell_size) + 1;
    int window_height = (grid_height * grid_cell_size) + 1;

    // Place the grid cursor in the middle of the screen.
    SDL_Rect grid_cursor = {
        .x = (grid_width - 1) / 2 * grid_cell_size,
        .y = (grid_height - 1) / 2 * grid_cell_size,
        .w = grid_cell_size,
        .h = grid_cell_size,
    };

    // The cursor ghost is a cursor that always shows in the cell below the
    // mouse cursor.
    SDL_Rect grid_cursor_ghost = { grid_cursor.x, grid_cursor.y, grid_cell_size,
                                  grid_cell_size };

    // Dark theme.
    SDL_Color grid_background = { 22, 22, 22, 255 }; // Barely Black
    SDL_Color grid_line_color = { 44, 44, 44, 255 }; // Dark grey
    SDL_Color grid_cursor_ghost_color = { 44, 44, 44, 255 };
    SDL_Color grid_cursor_color = { 255, 255, 255, 255 }; // White

    // Light Theme.
    // SDL_Color grid_background = {233, 233, 233, 255}; // Barely white
    // SDL_Color grid_line_color = {200, 200, 200, 255}; // Very light grey
    // SDL_Color grid_cursor_ghost_color = {200, 200, 200, 255};
    // SDL_Color grid_cursor_color = {160, 160, 160, 255}; // Grey

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Initialize SDL: %s",
            SDL_GetError());
        return EXIT_FAILURE;
    }

    SDL_Window* window;
    SDL_Renderer* renderer;
    if (SDL_CreateWindowAndRenderer(window_width, window_height, 0, &window,
        &renderer) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION,
            "Create window and renderer: %s", SDL_GetError());
        return EXIT_FAILURE;
    }

    SDL_SetWindowTitle(window, "SDL Grid");

    SDL_bool quit = SDL_FALSE;
    SDL_bool mouse_active = SDL_FALSE;
    SDL_bool mouse_hover = SDL_FALSE;

    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                case SDLK_w:
                case SDLK_UP:
                    grid_cursor.y -= grid_cell_size;
                    break;
                case SDLK_s:
                case SDLK_DOWN:
                    grid_cursor.y += grid_cell_size;
                    break;
                case SDLK_a:
                case SDLK_LEFT:
                    grid_cursor.x -= grid_cell_size;
                    break;
                case SDLK_d:
                case SDLK_RIGHT:
                    grid_cursor.x += grid_cell_size;
                    break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                grid_cursor.x = (event.motion.x / grid_cell_size) * grid_cell_size;
                grid_cursor.y = (event.motion.y / grid_cell_size) * grid_cell_size;
                break;
            case SDL_MOUSEMOTION:
                grid_cursor_ghost.x = (event.motion.x / grid_cell_size) * grid_cell_size;
                grid_cursor_ghost.y = (event.motion.y / grid_cell_size) * grid_cell_size;

                if (!mouse_active)
                    mouse_active = SDL_TRUE;
                break;
            case SDL_WINDOWEVENT:
                if (event.window.event == SDL_WINDOWEVENT_ENTER && !mouse_hover)
                    mouse_hover = SDL_TRUE;
                else if (event.window.event == SDL_WINDOWEVENT_LEAVE && mouse_hover)
                    mouse_hover = SDL_FALSE;
                break;
            case SDL_QUIT:
                quit = SDL_TRUE;
                break;
            }
        }

        // Draw grid background.
        SDL_SetRenderDrawColor(renderer, grid_background.r, grid_background.g,
            grid_background.b, grid_background.a);
        SDL_RenderClear(renderer);

        // Draw grid lines.
        SDL_SetRenderDrawColor(renderer, grid_line_color.r, grid_line_color.g,
            grid_line_color.b, grid_line_color.a);

        for (int x = 0; x < 1 + grid_width * grid_cell_size;
            x += grid_cell_size) {
            SDL_RenderDrawLine(renderer, x, 0, x, window_height);
        }

        for (int y = 0; y < 1 + grid_height * grid_cell_size;
            y += grid_cell_size) {
            SDL_RenderDrawLine(renderer, 0, y, window_width, y);
        }

        // Draw grid ghost cursor.
        if (mouse_active && mouse_hover) {
            SDL_SetRenderDrawColor(renderer, grid_cursor_ghost_color.r,
                grid_cursor_ghost_color.g,
                grid_cursor_ghost_color.b,
                grid_cursor_ghost_color.a);
            SDL_RenderFillRect(renderer, &grid_cursor_ghost);
        }

        // Draw grid cursor.
        SDL_SetRenderDrawColor(renderer, grid_cursor_color.r,
            grid_cursor_color.g, grid_cursor_color.b,
            grid_cursor_color.a);
        SDL_RenderFillRect(renderer, &grid_cursor);

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return EXIT_SUCCESS;
}