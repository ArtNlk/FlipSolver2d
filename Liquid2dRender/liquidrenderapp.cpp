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
        json configJson;
        configFile >> configJson;
        loadJson(configJson["scene"]);
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

void LiquidRenderApp::loadJson(std::string fileName)
{
    try
    {
        std::ifstream sceneFile(fileName);
        if(!sceneFile.is_open())
        {
            std::cout << "errorOpening scene file " << fileName;
        }
        json sceneJson;
        sceneFile >> sceneJson;
        json settingsJson = sceneJson["settings"];

        std::string simTypeName = settingsJson["simType"].get<std::string>();
        SimulationMethod simMethod = simMethodFromName(simTypeName);


        switch(simMethod)
        {
            case SIMULATION_LIQUID:
            {
                FlipSolverParameters p;
                populateFlipSolverParamsFromJson(&p, settingsJson);
                m_solver.reset(new FlipSolver(&p));
            }
            break;

            case SIMULATION_SMOKE:
            {
                SmokeSolverParameters p;
                populateFlipSolverParamsFromJson(&p, settingsJson);
                populateSmokeSolverParamsFromJson(&p, settingsJson);
                m_solver.reset(new FlipSmokeSolver(&p));
            }
            break;

            case SIMULATION_FIRE:
            {
                FireSolverParameters p;
                populateFlipSolverParamsFromJson(&p, settingsJson);
                populateFireSolverParamsFromJson(&p, settingsJson);
                m_solver.reset(new FlipFireSolver(&p));
            }
            break;

            case SIMULATION_NBFLIP:
            {
                NBFlipParameters p;
                populateFlipSolverParamsFromJson(&p, settingsJson);
                populateNBFlipSolverParamsFromJson(&p, settingsJson);
                m_solver.reset(new NBFlipSolver(&p));
            }
            break;
        }

        m_solver->initAdditionalParameters();
        objectsFromJson(sceneJson["solver"]);
    }
    catch (std::exception &e)
    {
        debug() << e.what();
        std::cout << e.what();
    }
}

void LiquidRenderApp::populateFlipSolverParamsFromJson(FlipSolverParameters *p, json settingsJson)
{
    float s = tryGetValue(settingsJson,"scale",1.f);
    p->fluidDensity = tryGetValue(settingsJson,"density",1.f);
    p->seed = tryGetValue(settingsJson,"seed",0);
    p->dx =
    p->particlesPerCell = settingsJson["particlesPerCell"].get<int>();
    std::pair<float,float> v = tryGetValue(settingsJson,"globalAcceleration",std::pair(9.8,0.f));
    // v.first *= s;
    // v.second *= s;
    p->globalAcceleration = v;
    p->resolution = settingsJson["resolution"].get<int>();
    p->fps = settingsJson["fps"].get<int>();
    p->maxSubsteps = tryGetValue(settingsJson,"maxSubsteps",30);
    p->picRatio = tryGetValue(settingsJson,"picRatio",0.03);
    p->cflNumber = tryGetValue(settingsJson,"cflNumber",10.f);
    p->particleScale = tryGetValue(settingsJson,"particleScale",0.8);
    p->pcgIterLimit = tryGetValue(settingsJson,"pcgIterLimit",200);
    p->viscosityEnabled = tryGetValue(settingsJson,"viscosityEnabled",false);
    p->domainSizeI = settingsJson["domainSizeI"].get<float>() * s;
    p->domainSizeJ = settingsJson["domainSizeJ"].get<float>() * s;
    p->useHeavyViscosity = tryGetValue(settingsJson, "heavyViscosity", false);
    p->sceneScale = s;

    if(p->domainSizeI > p->domainSizeJ)
    {
        p->dx = static_cast<float>(p->domainSizeI) /
                                                p->resolution;
        p->gridSizeI = p->resolution;
        p->gridSizeJ = (static_cast<float>(p->domainSizeJ) /
                    static_cast<float>(p->domainSizeI)) * p->resolution;
    }
    else
    {
        p->dx = static_cast<float>(p->domainSizeJ) /
                                                p->resolution;
        p->gridSizeJ = p->resolution;
        p->gridSizeI = (static_cast<float>(p->domainSizeI) /
                    static_cast<float>(p->domainSizeJ)) * p->resolution;
    }

    std::string parameterHandling = tryGetValue(settingsJson,"parameterHandlingMethod",std::string("particle"));
    if(parameterHandling == "particle")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::PARTICLE;
    }
    if(parameterHandling == "hybrid")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::HYBRID;
    }
    if(parameterHandling == "grid")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::GRID;
    }

    std::string simTypeName = settingsJson["simType"].get<std::string>();
    p->simulationMethod = simMethodFromName(simTypeName);
}

void LiquidRenderApp::populateNBFlipSolverParamsFromJson(NBFlipParameters *p, json settingsJson)
{
    return;
}

void LiquidRenderApp::populateSmokeSolverParamsFromJson(SmokeSolverParameters *p, json settingsJson)
{
    p->ambientTemperature = tryGetValue(settingsJson,"ambientTemperature",273.0f);
    p->temperatureDecayRate = tryGetValue(settingsJson,"temperatureDecayRate",0.0);
    p->concentrationDecayRate = tryGetValue(settingsJson,"concentrationDecayRate",0.0);
    p->buoyancyFactor = tryGetValue(settingsJson,"buoyancyFactor",1.0);
    p->sootFactor = tryGetValue(settingsJson,"sootFactor",1.0);
}

void LiquidRenderApp::populateFireSolverParamsFromJson(FireSolverParameters *p, json settingsJson)
{
    populateSmokeSolverParamsFromJson(p,settingsJson);
    p->ignitionTemperature = tryGetValue(settingsJson,"ignitionTemp",250.f);
    p->burnRate = tryGetValue(settingsJson,"burnRate",0.05f);
    p->smokeProportion = tryGetValue(settingsJson,"smokeEmission",1.f);
    p->heatProportion = tryGetValue(settingsJson,"heatEmission",1.f);
    p->divergenceProportion = tryGetValue(settingsJson,"billowing",0.1f);
}

SimulationMethod LiquidRenderApp::simMethodFromName(const std::string &name)
{
    if(name == "fluid" || name == "flip")
    {
        return SimulationMethod::SIMULATION_LIQUID;
    }
    if(name == "smoke")
    {
        return SimulationMethod::SIMULATION_SMOKE;
    }
    if(name == "fire")
    {
        return SimulationMethod::SIMULATION_FIRE;
    }
    if(name == "nbflip")
    {
        return SimulationMethod::SIMULATION_NBFLIP;
    }

    return SimulationMethod::SIMULATION_LIQUID;
}

//    SimSettings::airDensity() = tryGetValue(settingsJson,"airDensity",SimSettings::fluidDensity() * 0.001f);
//    SimSettings::surfaceTensionFactor() = tryGetValue(settingsJson,"surfaceTensionFactor",0.0);

void LiquidRenderApp::objectsFromJson(json solverJson)
{
    std::vector<json> objects = solverJson["objects"]
                                        .get<std::vector<json>>();
    for(json &geo : objects)
    {
        addObjectFromJson(geo);
    }
}

Emitter LiquidRenderApp::emitterFromJson(json emitterJson)
{
    std::vector<std::pair<float,float>> verts = emitterJson["verts"]
                                                .get<std::vector<std::pair<float,float>>>();
    float temp = tryGetValue(emitterJson,"temperature",273.f);
    float conc = tryGetValue(emitterJson,"concentration",1.f);
    float viscosity = tryGetValue(emitterJson,"viscosity",0.f);
    float fuel = tryGetValue(emitterJson,"fuel",1.f);
    float div = tryGetValue(emitterJson,"divergence",0.f);
    bool transfersVelocity = tryGetValue(emitterJson,"transferVelocity",false);
    std::pair<float,float> velocity = tryGetValue(emitterJson,"velocity",std::pair<float,float>(0.f, 0.f));
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(m_solver->sceneScale() * Vec3(v.first,v.second));
    }

    Emitter output(geo);

    output.setTemperature(temp);
    output.setConcentrartion(conc);
    output.setViscosity(viscosity);
    output.setFuel(fuel);
    output.setDivergence(div);
    output.setVelocity(velocity);
    output.setVelocityTransfer(transfersVelocity);

    return output;
}

Obstacle LiquidRenderApp::obstacleFromJson(json obstacleJson)
{
    float friction = tryGetValue(obstacleJson,"friction",0);
    std::vector<std::pair<float,float>> verts = obstacleJson["verts"]
                                                .get<std::vector<std::pair<float,float>>>();
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(m_solver->sceneScale() * Vec3(v.first,v.second));
    }

    return Obstacle(friction,geo);
}

Sink LiquidRenderApp::sinkFromJson(json sinkJson)
{
    std::vector<std::pair<float,float>> verts = sinkJson["verts"]
                                                .get<std::vector<std::pair<float,float>>>();
    float div = tryGetValue(sinkJson,"divergence",0.f);
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(m_solver->sceneScale() * Vec3(v.first,v.second));
    }

    return Sink(div,geo);
}

void LiquidRenderApp::addObjectFromJson(json objectJson)
{
    std::string geoType = objectJson["type"].get<std::string>();

    bool enabled = tryGetValue(objectJson,"enabled",true);
    if(!enabled)
    {
        return;
    }

    if(geoType == "solid")
    {
        Obstacle o = obstacleFromJson(objectJson);
        m_solver->addGeometry(o);
    }
    else if(geoType == "source")
    {
        Emitter e = emitterFromJson(objectJson);
        m_solver->addSource(e);
    }
    else if(geoType == "sink")
    {
        Sink s = sinkFromJson(objectJson);
        m_solver->addSink(s);
    }
    else if(geoType == "fluid")
    {
        Emitter e = emitterFromJson(objectJson);
        m_solver->addInitialFluid(e);
    }
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

template<class T>
T LiquidRenderApp::tryGetValue(json input, std::string key, T defaultValue)
{
    return input.contains(key) ? input[key].get<T>() : defaultValue;
}
