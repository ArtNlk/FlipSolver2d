#include "liquidrenderapp.h"
#include <stdexcept>
#include <cstdlib>

#include "linearindexable2d.h"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;

LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(new FlipSolver(50,75,1,0.01,2,1,false)),
    m_fluidRenderer(m_solver)
{
    m_solver->grid().setMaterial(10,10,FluidCellMaterial::SOLID);
    m_solver->grid().setMaterial(10,11,FluidCellMaterial::FLUID);
    for(int i = 0; i < 1000; i++)
    {
        m_solver->grid().setU(rand() % (m_solver->grid().sizeI() / 2), rand() % m_solver->grid().sizeJ(), rand() % 20 - 10,true);
        m_solver->grid().setV(rand() % (m_solver->grid().sizeI() / 2), rand() % m_solver->grid().sizeJ(), rand() % 20 - 10,true);
    }
    LiquidRenderApp::GLFWCallbackWrapper::SetApplication(this);
    m_fluidRenderer.setRenderMode(FluidRenderMode::RENDER_V);
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

    glfwSetFramebufferSizeCallback(m_window,LiquidRenderApp::GLFWCallbackWrapper::ResizeCallback);
    glfwSetKeyCallback(m_window, LiquidRenderApp::GLFWCallbackWrapper::KeyboardCallback);

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

void LiquidRenderApp::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    switch(action)
    {
        case GLFW_PRESS:
            switch(key)
            {
                case GLFW_KEY_P:
                    m_solver->extrapolateVelocityField();
                    m_fluidRenderer.updateGrid();
                break;

                case GLFW_KEY_T:
                    m_fluidRenderer.renderMode()++;
                    m_fluidRenderer.updateGrid();
                break;
            }
        break;
    }
}
