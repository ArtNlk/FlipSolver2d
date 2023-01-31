#ifndef LIQUIDRENDERAPP_H
#define LIQUIDRENDERAPP_H

#include "flipsolver2d.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <ft2build.h>
#include FT_FREETYPE_H
#include <memory>

#include "nlohmann/json.hpp"

#include "flipsmokesolver.h"
#include "fluidrenderer.h"
#include "geometry2d.h"
#include "textmenurenderer.h"
#include "emitter.h"

using json = nlohmann::json;

class LiquidRenderApp
{
public:
    LiquidRenderApp();

    void init();
    void run();
    void requestRender();

protected:
    void resizeCallback(GLFWwindow* window, int width, int height);
    void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

    void loadJson(std::string fileName);

    template<class T>
    T tryGetValue(json input, std::string key, T defaultValue);

    void settingsFromJson(json settingsJson);
    void solverFromJson(json solverJson);
    Emitter emitterFromJson(json emitterJson);
    Obstacle obstacleFromJson(json obstacleJson);
    Sink sinkFromJson(json sinkJson);
    void addObjectFromJson(json objectJson);
    void setupFluidrender();
    void setupFluidrenderQuad();
    void addVert(std::vector<float> &vertexVector, float x, float y, float u, float v);
    void formQuad(std::vector<unsigned int> &indexVector, std::vector<float> &vertexVector);
    void updateFluidrenderQuad();
    void updateFluidrenderBuffers();
    void updateFluidrenderQuadVertex(Vertex v, int vertexIndex);
    void render();
    void resizeFluidrenderQuad();

    GLFWwindow* m_window;
    std::shared_ptr<FlipSolver> m_solver;
    FluidRenderer m_fluidRenderer;
    TextMenuRenderer m_textMenuRenderer;

    unsigned int m_fluidgrid_vbo;
    unsigned int m_fluidgrid_vao;
    unsigned int m_fluidgrid_ebo;
    unsigned int m_textureQuadShaderProgram;

    std::vector<float> m_fluidgridQuadVerts;
    std::vector<unsigned int> m_fluidgridQuadIndices;

    int m_windowWidth;
    int m_windowHeight;

    bool m_renderRequested;

    static const int m_startWindowWidth = 900;
    static const int m_startWindowHeight = 900;

    static constexpr float m_gridDrawFraction = 0.85;

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
