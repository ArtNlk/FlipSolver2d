#include "liquidrenderapp.h"
#include <stdexcept>
#include <cstdlib>

#include "linearindexable2d.h"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;

LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(new FlipSolver(m_gridSizeI,m_gridSizeJ,1,1,1,1,false)),
    m_fluidRenderer(m_solver)
{
    LiquidRenderApp::GLFWCallbackWrapper::SetApplication(this);
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
                case GLFW_KEY_M:
                    m_fluidRenderer.renderMode()++;
                    m_fluidRenderer.updateGrid();
                break;

                case GLFW_KEY_U:
                    m_fluidRenderer.updateGrid();
                break;

                case GLFW_KEY_E:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        initGridForExtrapolation();
                    }
                    else
                    {
                        m_solver->extrapolateVelocityField();
                        m_fluidRenderer.updateGrid();
                    }
                break;

                case GLFW_KEY_P:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        initGridForProjection();
                    }
                    else
                    {
                        m_solver->project();
                        m_fluidRenderer.updateGrid();
                    }
                break;
            }
        break;
    }
}

void LiquidRenderApp::resetGrid()
{
    m_solver->grid().fill();
}

void LiquidRenderApp::initGridForExtrapolation()
{
    resetGrid();
    for(int i = 0; i < 1000; i++)
    {
        m_solver->grid().setU(rand() % (m_solver->grid().sizeI() / 2), rand() % m_solver->grid().sizeJ(), static_cast<float>(rand() % 20 - 10)/10,true);
        m_solver->grid().setV(rand() % (m_solver->grid().sizeI() / 2), rand() % m_solver->grid().sizeJ(), static_cast<float>(rand() % 20 - 10)/10,true);
    }
    m_fluidRenderer.updateGrid();
}

void LiquidRenderApp::initGridForProjection()
{
    resetGrid();
    srand(0);
    Index2d fluidTopLeft(1,1);
    Index2d fluidBottomRight(9,9);
    m_solver->grid().fillMaterialRect(FluidCellMaterial::SOLID,0,0,9,9);
    m_solver->grid().fillMaterialRect(FluidCellMaterial::FLUID,fluidTopLeft,fluidBottomRight);
    for (int i = 0; i < m_solver->grid().sizeI(); i++)
    {
        for (int j = 0; j < m_solver->grid().sizeJ(); j++)
        {
            m_solver->grid().setU(i,j,static_cast<float>(rand() % 20 - 10) / 10);
            m_solver->grid().setV(i,j,static_cast<float>(rand() % 20 - 10) / 10);
            //m_solver->grid().setU(i,j,1.0);
            //m_solver->grid().setV(i,j,0.0);
        }
    }
    m_solver->grid().fillKnownFlagsU(true);
    m_solver->grid().fillKnownFlagsV(true);
    //m_solver->grid().velocityGridU().fillRect(10.f,fluidTopLeft,fluidBottomRight);
    //m_solver->grid().knownFlagsGridU().fillRect(true,fluidTopLeft,fluidBottomRight);
    //m_solver->grid().velocityGridV().fillRect(23.f,fluidTopLeft,fluidBottomRight);
    //m_solver->grid().knownFlagsGridV().fillRect(true,fluidTopLeft,fluidBottomRight);
    m_fluidRenderer.updateGrid();
}
