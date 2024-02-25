#ifndef LIQUIDRENDERAPP_H
#define LIQUIDRENDERAPP_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <memory>
#include <chrono>

#include "nlohmann/json.hpp"

#include "imgui.h"

#include "flipsolver2d.h"
#include "nbflipsolver.h"
#include "flipsmokesolver.h"
#include "flipfiresolver.h"

#include "fluidrenderer.h"
#include "geometry2d.h"
#include "emitter.h"

using json = nlohmann::json;

class LiquidRenderApp
{
public:
    LiquidRenderApp();
    ~LiquidRenderApp();

    void init();
    void run();
    void requestRender();

protected:
    void resizeCallback(GLFWwindow* window, int width, int height);

    void loadJson(std::string fileName);

    void populateFlipSolverParamsFromJson(FlipSolverParameters* p, json settingsJson);
    void populateNBFlipSolverParamsFromJson(NBFlipParameters* p, json settingsJson);
    void populateSmokeSolverParamsFromJson(SmokeSolverParameters* p, json settingsJson);
    void populateFireSolverParamsFromJson(FireSolverParameters* p, json settingsJson);

    SimulationMethod simMethodFromName(const std::string& name);


    template<class T>
    T tryGetValue(json input, std::string key, T defaultValue);

    void objectsFromJson(json solverJson);
    Emitter emitterFromJson(json emitterJson);
    Obstacle obstacleFromJson(json obstacleJson);
    Sink sinkFromJson(json sinkJson);
    void addObjectFromJson(json objectJson);

    void render();
    void renderSceneViewPanel(bool update);
    bool renderControlsPanel();
    void renderStatsPanel();

    bool gridRenderCombo();
    bool vectorRenderCombo();
    bool particleRenderCombo();

    const std::string stepStageToString(SolverStage stage) const;

    GLFWwindow* m_window;
    std::shared_ptr<FlipSolver> m_solver;
    FluidRenderer m_fluidRenderer;
    static const char* m_configFilePath;

    unsigned int m_fluidgrid_vbo;
    unsigned int m_fluidgrid_vao;
    unsigned int m_fluidgrid_ebo;
    unsigned int m_textureQuadShaderProgram;

    std::vector<float> m_fluidgridQuadVerts;
    std::vector<unsigned int> m_fluidgridQuadIndices;

    int m_windowWidth;
    int m_windowHeight;

    bool m_renderRequested;

    int m_simStepsLeft;
    bool m_recording;

    SolverTimeStats m_lastFrameStats;

    static const int m_startWindowWidth = 900;
    static const int m_startWindowHeight = 900;

    static constexpr float m_gridDrawFraction = 0.85;

    std::chrono::time_point<std::chrono::high_resolution_clock> m_lastFrameTime;

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
        static void SetApplication(LiquidRenderApp *application)
        {
            GLFWCallbackWrapper::s_application = application;
        }
    private:
        static LiquidRenderApp* s_application;
    };

    ImGuiIO* m_io;
};

#endif // LIQUIDRENDERAPP_H
