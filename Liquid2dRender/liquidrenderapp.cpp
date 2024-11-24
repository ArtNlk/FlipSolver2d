#include "liquidrenderapp.h"
#include <chrono>
#include <memory>
#include <ratio>
#include <stdexcept>
#include <filesystem>

#include "GLFW/glfw3.h"
#include "fluidrenderer.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "flipfiresolver.h"
#include "flipsmokesolver.h"
#include "flipsolver2d.h"
#include "nbflipsolver.h"
#include "logger.h"
#include "jsonscenereader.h"

#include "nlohmann/json.hpp"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;
const char* LiquidRenderApp::m_configFilePath = "./config.json";


LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(nullptr),
    m_fluidRenderer(m_startWindowWidth,m_startWindowHeight),
    m_renderRequested(false),
    m_simStepsLeft(0),
    m_recording(false)
{
    m_windowWidth = m_startWindowWidth;
    m_windowHeight = m_startWindowHeight;
}

LiquidRenderApp::~LiquidRenderApp()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(m_window);
    glfwTerminate();
}

void LiquidRenderApp::init()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    {
        std::ifstream configFile(m_configFilePath);
        if(!configFile.is_open())
        {
            std::cout << "error opening config file " << m_configFilePath;
        }
        nlohmann::json configJson;
        configFile >> configJson;
        m_solver = JsonSceneReader::loadJson(configJson["scene"]);
    }

    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());

    m_windowWidth = mode->width * 0.8f;
    m_windowHeight = mode->height * 0.8f;

    m_window = glfwCreateWindow(m_windowWidth, m_windowHeight, "Flip fluid 2d", nullptr, nullptr);
    if (m_window == nullptr)
    {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(m_window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize glad");
    }

    std::cout << "Render: " << glGetString(GL_RENDERER);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    m_io = &ImGui::GetIO(); (void)m_io;
    m_io->ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    m_io->ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGuiIO& io = ImGui::GetIO();
    // ImFont* fontDefault = io.Fonts->AddFontDefault();
    // fontDefault->Scale = 14.f/13.f;
    //ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(m_window, true);
    ImGui_ImplOpenGL3_Init();

    //glViewport(0, 0, m_windowWidth, m_windowHeight);

    m_fluidRenderer.init(m_solver);
    m_fluidRenderer.particlesEnabled() = true;

    //setupFluidrender();
    //resizeFluidrenderQuad();
    m_solver->updateSinks();
    m_solver->updateSources();
    m_solver->updateSolids();
    m_fluidRenderer.update();
    m_renderRequested = true;

    m_lastFrameTime = std::chrono::high_resolution_clock::now();
}

void LiquidRenderApp::run()
{
    while(!glfwWindowShouldClose(m_window))
    {
        glfwPollEvents();
        render();
    }
}

void LiquidRenderApp::resizeCallback(GLFWwindow *window, int width, int height)
{
    m_windowWidth = width;
    m_windowHeight = height;
    //resizeFluidrenderQuad();
    //m_textMenuRenderer.resize(width,height);
    m_renderRequested = true;
}

void LiquidRenderApp::requestRender()
{
    m_renderRequested = true;
}

void LiquidRenderApp::render()
{
    std::chrono::duration<float, std::milli> dur =
        std::chrono::high_resolution_clock::now() - m_lastFrameTime;

    if(dur.count() < 33.3f)
    {
        return;
    }
    m_lastFrameTime = std::chrono::high_resolution_clock::now();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    //ImGui::ShowDemoWindow();
    static bool firstUpdate = true;
    if(firstUpdate)
    {
        m_fluidRenderer.update();
        m_fluidRenderer.render();
        firstUpdate = false;
    }

    const bool update = renderControlsPanel();
    renderSceneViewPanel(update);
    renderStatsPanel();

    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(m_window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);

    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(m_window);
}

void LiquidRenderApp::renderSceneViewPanel(bool update)
{
    ImGui::Begin("Scene view");
    ImGui::SetWindowSize(ImVec2(0,0));

    if(update)
    {
        std::cout << "update true" << std::endl;
        m_fluidRenderer.update();
        m_fluidRenderer.render();
    }

    ImVec2 availSize = ImGui::GetContentRegionAvail();
    ImVec2 fluidSize = ImVec2(m_solver->gridSizeJ(),m_solver->gridSizeI());
    float scaleX = availSize.x/fluidSize.x;
    float scaleY = availSize.y/fluidSize.y;
    float scale = std::min(scaleX, scaleY);
    ImVec2 renderSize(fluidSize.x*scale,fluidSize.y*scale);

    m_fluidRenderer.resizeTexture(renderSize.x,renderSize.y);

    ImGui::Image((void*)(intptr_t)m_fluidRenderer.renderTexture(),
                 renderSize,
                 ImVec2(0,1),ImVec2(1,0));

    ImGui::End();
}

bool LiquidRenderApp::renderControlsPanel()
{
    ImGui::Begin("Controls");
    ImGui::SeparatorText("Run");

    bool update = false;
    if(ImGui::Button("Step one frame"))
    {
        m_simStepsLeft = 1;
    }

    static int steps = 30;
    ImGui::InputInt("Steps to run",&steps,1,10);
    steps = std::clamp(steps,0,1000000000);

    ImGui::BeginDisabled(m_simStepsLeft > 0);
    if(ImGui::Button("Step multiple frames"))
    {
        m_simStepsLeft = steps;
    }

    if(ImGui::Button("Record multiple frames"))
    {
        m_simStepsLeft = steps;
        m_recording = true;
        static bool isFirst = true;
        if(isFirst)
        {
            std::filesystem::create_directory("./output");
            for (const auto& entry : std::filesystem::directory_iterator("./output"))
                std::filesystem::remove_all(entry.path());
            isFirst = false;
        }
    }
    ImGui::EndDisabled();

    ImGui::BeginDisabled(m_simStepsLeft == 0);
    if(ImGui::Button("Stop"))
    {
        m_simStepsLeft = 0;
    }
    ImGui::EndDisabled();

    if(m_simStepsLeft > 0)
    {
        m_solver->stepFrame();
        m_lastFrameStats = m_solver->timeStats();
        std::cout << "Stat: " << m_lastFrameStats.timings().at(SolverStage::PRESSURE) / m_lastFrameStats.substepCount() << std::endl;
        update = true;
        m_simStepsLeft--;
        if(m_recording)
        {
            m_fluidRenderer.dumpToTga(std::to_string(m_solver->fps()) + "_fps_" + std::to_string(m_solver->frameNumber()) + ".tga");
        }
    }
    else if(m_recording)
    {
        m_recording = false;
    }
    ImGui::Value("Steps left",m_simStepsLeft);

    ImGui::SeparatorText("Render control");

    update |= gridRenderCombo();
    update |= vectorRenderCombo();
    update |= particleRenderCombo();

    update |= ImGui::Checkbox("Vectors enabled",&m_fluidRenderer.vectorRenderEnabled());
    update |= ImGui::Checkbox("Particles enabled",&m_fluidRenderer.particlesEnabled());
    update |= ImGui::Checkbox("Geometry outline enabled",&m_fluidRenderer.geometryEnabled());
    update |= ImGui::Checkbox("Sinks and sources enabled",&m_fluidRenderer.extrasEnabled());

    update |= ImGui::InputInt("Particle size",&m_fluidRenderer.particleSize(),1,10);
    m_fluidRenderer.particleSize() = std::clamp(m_fluidRenderer.particleSize(),1,1000000000);

    ImGui::End();
    return update;
}

void LiquidRenderApp::renderStatsPanel()
{
    ImGui::Begin("Stats");
    ImGui::SeparatorText("Application stats");
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / m_io->Framerate, m_io->Framerate);
    ImGui::SeparatorText("Solver stats");
    ImGui::Text("Total frame time: %.3f ms", m_lastFrameStats.frameTime());
    ImGui::Text("Substeps taken: %d", m_lastFrameStats.substepCount());

    SolverStage stage = static_cast<SolverStage>(0);
    do
    {
        ImGui::Text(stepStageToString(stage).c_str(), m_lastFrameStats.timings().at(stage),
                    (m_lastFrameStats.timings().at(stage)/m_lastFrameStats.frameTime()) * 100.f);
    }while(nextEnum<SolverStage,SOLVER_STAGE_COUNT>(stage));

    ImGui::SeparatorText("Max PCG iterations over substeps:");
    ImGui::Text("Pressure: %d/%d", m_lastFrameStats.pressureIterations(), m_solver->pcgIterationLimit());
    ImGui::Text("Density: %d/%d", m_lastFrameStats.densityIterations(), m_solver->pcgIterationLimit());
    ImGui::Text("Viscosity: %d/%d", m_lastFrameStats.viscosityIterations(), m_solver->pcgIterationLimit());

    ImGui::End();
}

bool LiquidRenderApp::gridRenderCombo()
{
    bool updated = false;
    if (ImGui::BeginCombo("Grid render mode", m_fluidRenderer.currentFluidRenderModeName().c_str()))
    {
        FluidRenderMode fluidMode = static_cast<FluidRenderMode>(0);

        do
        {
            if (ImGui::Selectable(m_fluidRenderer.fluidModeToName(fluidMode).c_str(),
                                  m_fluidRenderer.gridRenderMode() == fluidMode))
            {
                m_fluidRenderer.setRenderMode(fluidMode);
                updated = true;
            }
        }while(nextEnum<FluidRenderMode,GRID_RENDER_ITER_END>(fluidMode));

        ImGui::EndCombo();
    }
    return updated;
}

bool LiquidRenderApp::vectorRenderCombo()
{
    bool updated = false;
    if (ImGui::BeginCombo("Vector render mode", m_fluidRenderer.currentVectorRenderModeName().c_str()))
    {
        VectorRenderMode vectorMode = static_cast<VectorRenderMode>(0);

        do
        {
            if (ImGui::Selectable(m_fluidRenderer.vectorModeToName(vectorMode).c_str(),
                                  m_fluidRenderer.vectorRenderMode() == vectorMode))
            {
                m_fluidRenderer.vectorRenderMode() = vectorMode;
                updated = true;
            }
        }while(nextEnum<VectorRenderMode,VECTOR_RENDER_ITER_END>(vectorMode));

        ImGui::EndCombo();
    }
    return updated;
}

bool LiquidRenderApp::particleRenderCombo()
{
    bool updated = false;
    if (ImGui::BeginCombo("Particle render mode", m_fluidRenderer.currentParticleRenderModeName().c_str()))
    {
        ParticleRenderMode particleMode = static_cast<ParticleRenderMode>(0);

        do
        {
            if (ImGui::Selectable(m_fluidRenderer.particleModeToName(particleMode).c_str(),
                                  m_fluidRenderer.particleRenderMode() == particleMode))
            {
                m_fluidRenderer.particleRenderMode() = particleMode;
                updated = true;
            }
        }while(nextEnum<ParticleRenderMode,PARTICLE_RENDER_ITER_END>(particleMode));

        ImGui::EndCombo();
    }
    return updated;
}

const std::string LiquidRenderApp::stepStageToString(SolverStage stage) const
{
    switch(stage)
    {
    case ADVECTION:
        return "Advection: %.3f ms %.1f%%";
        break;
    case DECOMPOSITION:
        return "Decomposition: %.3f ms %.1f%%";
        break;
    case DENSITY:
        return "Density correction: %.3f ms %.1f%%";
        break;
    case PARTICLE_REBIN:
        return "Particle rebinning: %.3f ms %.1f%%";
        break;
    case PARTICLE_TO_GRID:
        return "Particle to grid transfer: %.3f ms %.1f%%";
        break;
    case GRID_UPDATE:
        return "Grid update: %.3f ms %.1f%%";
        break;
    case AFTER_TRANSFER:
        return "After particle transfer: %.3f ms %.1f%%";
        break;
    case PRESSURE:
        return "Pressure update: %.3f ms %.1f%%";
        break;
    case VISCOSITY:
        return "Viscosity update: %.3f ms %.1f%%";
        break;
    case REPRESSURE:
        return "Secondary pressure update: %.3f ms %.1f%%";
        break;
    case PARTICLE_UPDATE:
        return "Particle values update: %.3f ms %.1f%%";
        break;
    case PARTICLE_RESEED:
        return "Particle reseeding: %.3f ms %.1f%%";
        break;
    case SOLVER_STAGE_COUNT:
    default:
        return "Invalid timing stage value!";
        break;
    }
}
