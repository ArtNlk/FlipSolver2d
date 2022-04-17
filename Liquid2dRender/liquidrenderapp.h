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
    void resizeCallback(GLFWwindow* window, int width, int height);
    void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

    GLFWwindow* m_window;
    std::shared_ptr<FlipSolver> m_solver;
    FluidRenderer m_fluidRenderer;

    //Taken from https://stackoverflow.com/questions/7676971/pointing-to-a-function-that-is-a-class-member-glfw-setkeycallback
    class GLFWCallbackWrapper
    {
    public:
        GLFWCallbackWrapper() = delete;
        GLFWCallbackWrapper(const GLFWCallbackWrapper&) = delete;
        GLFWCallbackWrapper(GLFWCallbackWrapper&&) = delete;
        ~GLFWCallbackWrapper() = delete;

        static void ResizeCallback(GLFWwindow *window, int width, int height)
        {
            s_application->resizeCallback(window,width,height);
        }
        static void KeyboardCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
        {
            s_application->keyCallback(window,key,scancode,action,mods);
        }
        static void SetApplication(LiquidRenderApp *application)
        {
            GLFWCallbackWrapper::s_application = application;
        }
    private:
        static LiquidRenderApp* s_application;
    };
};

#endif // LIQUIDRENDERAPP_H
