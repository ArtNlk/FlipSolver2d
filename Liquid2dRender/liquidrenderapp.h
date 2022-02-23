#ifndef LIQUIDRENDERAPP_H
#define LIQUIDRENDERAPP_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "flipsolver2d.h"

class LiquidRenderApp
{
public:
    LiquidRenderApp();

    void init();
    void run();

protected:
    static void resizeCallback(GLFWwindow* window, int width, int height);

    GLFWwindow* window;
    FlipSolver solver;
};

#endif // LIQUIDRENDERAPP_H
