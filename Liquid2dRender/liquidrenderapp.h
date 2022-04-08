#ifndef LIQUIDRENDERAPP_H
#define LIQUIDRENDERAPP_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <memory>

#include "flipsolver2d.h"
#include "fluidrenderer.h"

class LiquidRenderApp
{
public:
    LiquidRenderApp();

    void init();
    void run();

protected:
    static void resizeCallback(GLFWwindow* window, int width, int height);

    GLFWwindow* m_window;
    std::shared_ptr<FlipSolver> m_solver;
    FluidRenderer m_fluidRenderer;
};

#endif // LIQUIDRENDERAPP_H
