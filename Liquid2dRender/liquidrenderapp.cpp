#include "liquidrenderapp.h"
#include <stdexcept>

#include "linearindexable2d.h"

LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(new FlipSolver(120,250,1,0.01,2)),
    m_fluidRenderer(m_solver)
{

}

void LiquidRenderApp::init()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    m_window = glfwCreateWindow(800, 600, "Flip fluid 2d", NULL, NULL);
    if (m_window == NULL)
    {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(m_window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize glad");
    }

    glViewport(0, 0, 800, 600);

    glfwSetFramebufferSizeCallback(m_window,resizeCallback);

    m_fluidRenderer.init();
}

void LiquidRenderApp::run()
{
    while(!glfwWindowShouldClose(m_window))
    {
        m_fluidRenderer.render();
        glfwSwapBuffers(m_window);
        glfwPollEvents();
    }

    glfwTerminate();
}

void LiquidRenderApp::resizeCallback(GLFWwindow *window, int width, int height)
{
    glViewport(0, 0, width, height);
}
