#include "liquidrenderapp.h"
#include <memory>
#include <stdexcept>
#include <cstdlib>
#include <filesystem>

#include "flipfiresolver.h"
#include "flipsmokesolver.h"
#include "flipsolver2d.h"
#include "globalcallbackhandler.h"
#include "logger.h"
#include "simsettings.h"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;
const char* LiquidRenderApp::m_configFilePath = "./config.json";


LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(nullptr),
    m_fluidRenderer(m_startWindowWidth,m_startWindowHeight),
    m_textMenuRenderer(0,0,m_startWindowWidth,m_startWindowHeight,m_fluidRenderer),
    m_renderRequested(false)
{
    m_windowWidth = m_startWindowWidth;
    m_windowHeight = m_startWindowHeight;
    LiquidRenderApp::GLFWCallbackWrapper::SetApplication(this);
}

void LiquidRenderApp::init()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GlobalCallbackHandler::instance().init(this,
                                           &m_fluidRenderer,
                                           &m_textMenuRenderer);

    {
        std::ifstream configFile(m_configFilePath);
        if(!configFile.is_open())
        {
            std::cout << "errorOpening config file " << m_configFilePath;
        }
        json configJson;
        configFile >> configJson;
        loadJson(configJson["scene"]);
    }

    m_window = glfwCreateWindow(m_windowWidth, m_windowHeight, "Flip fluid 2d", NULL, NULL);
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

    glViewport(0, 0, m_windowWidth, m_windowHeight);

    glfwSetFramebufferSizeCallback(m_window,LiquidRenderApp::GLFWCallbackWrapper::ResizeCallback);
    glfwSetKeyCallback(m_window, LiquidRenderApp::GLFWCallbackWrapper::KeyboardCallback);

    m_fluidRenderer.init(m_solver);
    m_textMenuRenderer.init();

    setupFluidrender();
    resizeFluidrenderQuad();
    m_solver->updateSinks();
    m_solver->updateSources();
    m_solver->updateSolids();
    m_fluidRenderer.updateGrid();
    m_renderRequested = true;
}

void LiquidRenderApp::run()
{
    while(!glfwWindowShouldClose(m_window))
    {
        render();
        glfwWaitEvents();
    }

    glfwTerminate();
}

void LiquidRenderApp::resizeCallback(GLFWwindow *window, int width, int height)
{
    m_windowWidth = width;
    m_windowHeight = height;
    resizeFluidrenderQuad();
    m_textMenuRenderer.resize(width,height);
    m_renderRequested = true;
}

void LiquidRenderApp::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    switch(action)
    {
        case GLFW_PRESS:
            switch(key)
            {
                case GLFW_KEY_M:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        m_fluidRenderer.gridRenderMode()--;
                        m_fluidRenderer.updateGrid();
                    }
                    else
                    {
                        m_fluidRenderer.gridRenderMode()++;
                        m_fluidRenderer.updateGrid();
                    }
                break;

                case GLFW_KEY_V:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        m_fluidRenderer.toggleVectors();
                    }
                    else
                    {
                        m_fluidRenderer.vectorRenderMode()++;
                    }
                    m_fluidRenderer.update();
                break;

                case GLFW_KEY_G:
                    m_fluidRenderer.toggleGeometry();
                    m_fluidRenderer.update();
                break;

                case GLFW_KEY_U:
                    m_fluidRenderer.update();
                break;

                case GLFW_KEY_P:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        m_fluidRenderer.toggleParticles();
                    }
                    else
                    {
                        m_fluidRenderer.particleRenderMode()++;
                    }
                    m_fluidRenderer.update();
                break;

                case GLFW_KEY_R:
                {
                    int frame = 0;
                    static bool isFirst = true;
                    if(isFirst)
                    {
                        std::filesystem::create_directory("./output");
                        for (const auto& entry : std::filesystem::directory_iterator("./output"))
                                std::filesystem::remove_all(entry.path());
                        isFirst = false;
                    }
                    bool run = true;
                    for(int second = 0; second < 30; second++)
                    {
                        for(int i = 0; i < SimSettings::fps(); i++)
                        {
                            m_solver->stepFrame();
                            m_fluidRenderer.update();
                            render();
                            m_fluidRenderer.dumpToPng(std::to_string(SimSettings::fps()) + "fps_" + std::to_string(m_solver->frameNumber()) + ".png");
                            frame++;
                            glfwPollEvents();
                            if(glfwWindowShouldClose(m_window))
                            {
                                glfwTerminate();
                                run = false;
                            }
                            if(!run) break;
                        }
                        if(!run) break;
                    }
                }
                break;

                case GLFW_KEY_F:
                    m_fluidRenderer.dumpToPng("frame.png");
                break;

                case GLFW_KEY_SPACE:
                {
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        m_solver->stepFrame();
                        m_fluidRenderer.update();
                        render();
                    }
                    else
                    {
                        bool run = true;
                        for(int i = 0; i < SimSettings::fps(); i++)
                        {
                            m_solver->stepFrame();
                            m_fluidRenderer.update();
                            render();
                            glfwPollEvents();
                            if(glfwWindowShouldClose(m_window))
                            {
                                glfwTerminate();
                                run = false;
                            }
                            if(!run) break;
                        }
                    }

                }
                break;
            }
        break;
    }
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
        settingsFromJson(sceneJson["settings"]);

        switch(SimSettings::simMethod())
        {
            case SIMULATION_LIQUID:
                m_solver.reset(new FlipSolver(SimSettings::gridSizeI(),SimSettings::gridSizeJ(),1,true));
            break;

            case SIMULATION_SMOKE:
                m_solver.reset(new FlipSmokeSolver(SimSettings::gridSizeI(),SimSettings::gridSizeJ(),1,true));
            break;

            case SIMULATION_FIRE:
                m_solver.reset(new FlipFireSolver(SimSettings::gridSizeI(),SimSettings::gridSizeJ(),1, true));
            break;
        }

        solverFromJson(sceneJson["solver"]);
    }
    catch (std::exception &e)
    {
        debug() << e.what();
        std::cout << e.what();
    }
}

void LiquidRenderApp::settingsFromJson(json settingsJson)
{
    SimSettings::sceneScale() = tryGetValue(settingsJson,"scale",1.f);
    float s = SimSettings::sceneScale();
    SimSettings::domainSizeI() = std::ceil(s*settingsJson["domainSizeI"].get<int>());
    SimSettings::domainSizeJ() = std::ceil(s*settingsJson["domainSizeJ"].get<int>());
    SimSettings::resolution() = settingsJson["resolution"].get<int>();
    SimSettings::fps() = settingsJson["fps"].get<int>();
    SimSettings::frameDt() = 1.f / SimSettings::fps();
    SimSettings::maxSubsteps() = tryGetValue(settingsJson,"maxSubsteps",30);
    SimSettings::stepDt() = SimSettings::frameDt() / SimSettings::maxSubsteps();
    SimSettings::fluidDensity() = tryGetValue(settingsJson,"density",1.f);
    SimSettings::burnRate() = tryGetValue(settingsJson,"burnRate",0.05f);
    SimSettings::ignitionTemp() = tryGetValue(settingsJson,"ignitionTemp",250.f);
    SimSettings::smokeProportion() = tryGetValue(settingsJson,"smokeEmission",1.f);
    SimSettings::heatProportion() = tryGetValue(settingsJson,"heatEmission",1.f);
    SimSettings::divergenceProportion() = tryGetValue(settingsJson,"billowing",0.1f);
    SimSettings::airDensity() = tryGetValue(settingsJson,"airDensity",SimSettings::fluidDensity() * 0.001f);
    SimSettings::randomSeed() = tryGetValue(settingsJson,"seed",0);
    SimSettings::particlesPerCell() = settingsJson["particlesPerCell"].get<int>();
    SimSettings::cflNumber() = tryGetValue(settingsJson,"cflNumber",10.f);
    SimSettings::picRatio() = tryGetValue(settingsJson,"picRatio",0.03);
    SimSettings::ambientTemp() = tryGetValue(settingsJson,"ambientTemperature",273.0f);
    std::pair<float,float> v = tryGetValue(settingsJson,"globalAcceleration",std::pair(9.8,0.f));
    SimSettings::globalAcceleration() = Vertex(v.first,v.second);
    SimSettings::tempDecayRate() = tryGetValue(settingsJson,"temperatureDecayRate",0.0);
    SimSettings::concentrartionDecayRate() = tryGetValue(settingsJson,"concentrationDecayRate",0.0);
    SimSettings::particleScale() = tryGetValue(settingsJson,"particleScale",0.8);
    SimSettings::surfaceTensionFactor() = tryGetValue(settingsJson,"surfaceTensionFactor",0.0);
    SimSettings::pcgIterLimit() = tryGetValue(settingsJson,"pcgIterLimit",200);

    std::string simTypeName = settingsJson["simType"].get<std::string>();
    if(simTypeName == "fluid")
    {
        SimSettings::simMethod() = SimulationMethod::SIMULATION_LIQUID;
    }
    else if(simTypeName == "smoke")
    {
        SimSettings::simMethod() = SimulationMethod::SIMULATION_SMOKE;
    }
    else if(simTypeName == "fire")
    {
        SimSettings::simMethod() = SimulationMethod::SIMULATION_FIRE;
    }

    if(SimSettings::domainSizeI() > SimSettings::domainSizeJ())
    {
        SimSettings::dx() = static_cast<float>(SimSettings::domainSizeI()) /
                                                SimSettings::resolution();
        SimSettings::gridSizeI() = SimSettings::resolution();
        SimSettings::gridSizeJ() = (static_cast<float>(SimSettings::domainSizeJ()) /
                    static_cast<float>(SimSettings::domainSizeI())) * SimSettings::resolution();
    }
    else
    {
        SimSettings::dx() = static_cast<float>(SimSettings::domainSizeJ()) /
                                                SimSettings::resolution();
        SimSettings::gridSizeJ() = SimSettings::resolution();
        SimSettings::gridSizeI() = (static_cast<float>(SimSettings::domainSizeI()) /
                    static_cast<float>(SimSettings::domainSizeJ())) * SimSettings::resolution();
    }
}

void LiquidRenderApp::solverFromJson(json solverJson)
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
        geo.addVertex(SimSettings::sceneScale() * Vertex(v.first,v.second));
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
        geo.addVertex(SimSettings::sceneScale() * Vertex(v.first,v.second));
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
        geo.addVertex(SimSettings::sceneScale() * Vertex(v.first,v.second));
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

void LiquidRenderApp::setupFluidrender()
{
    static const char* vertexShaderSource = \
            "#version 330 core\n\
            layout (location = 0) in vec2 aPos;\n\
            layout (location = 1) in vec2 aTexCoords;\n\
            \n\
            out vec2 TexCoords;\n\
            \n\
            void main()\n\
            {\n\
                gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);\n\
                TexCoords = aTexCoords;\n\
            }\n";

    static const char* fragmentShaderSource = \
            "#version 330 core\n\
            out vec4 FragColor;\n\
            \n\
            in vec2 TexCoords;\n\
            \n\
            uniform sampler2D screenTexture;\n\
            \n\
            void main()\n\
            { \n\
                FragColor = texture(screenTexture, TexCoords);\n\
                //FragColor = vec4(TexCoords.x,TexCoords.y,0,1);\n\
            }";
    int success;
    unsigned int vShader;
    unsigned int fShader;
    vShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vShader,1,&vertexShaderSource,NULL);
    glCompileShader(vShader);
    glGetShaderiv(vShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Render app: Fluid grid quad vertex shader compilation error";
    }

    fShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fShader);
    glGetShaderiv(fShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Render app: Fluid grid quad fragment shader compilation error";
    }

    m_textureQuadShaderProgram = glCreateProgram();
    glAttachShader(m_textureQuadShaderProgram,vShader);
    glAttachShader(m_textureQuadShaderProgram,fShader);
    glLinkProgram(m_textureQuadShaderProgram);
    glGetProgramiv(m_textureQuadShaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        std::cout << "Render app: Fluid grid quad shader program linking error!";
    }

    glDeleteShader(vShader);
    glDeleteShader(fShader);

    glGenVertexArrays(1, &m_fluidgrid_vao);
    glGenBuffers(1,&m_fluidgrid_vbo);
    glGenBuffers(1,&m_fluidgrid_ebo);
    setupFluidrenderQuad();
}

void LiquidRenderApp::setupFluidrenderQuad()
{
    m_fluidgridQuadVerts.reserve(4*5);//4 verts, 5 floats each (xyz,uv)
    m_fluidgridQuadIndices.reserve(6);//two tris for one quad

    addVert(m_fluidgridQuadVerts,-1,1,0,1);
    addVert(m_fluidgridQuadVerts,0,1,1,1);
    addVert(m_fluidgridQuadVerts,0,0,1,0);
    addVert(m_fluidgridQuadVerts,-1,0,0,0);

    formQuad(m_fluidgridQuadIndices,m_fluidgridQuadVerts);
    updateFluidrenderBuffers();
}

void LiquidRenderApp::addVert(std::vector<float> &vertexVector, float x, float y, float u, float v)
{
    vertexVector.push_back(x);
    vertexVector.push_back(y);
    vertexVector.push_back(0);
    vertexVector.push_back(u);
    vertexVector.push_back(v);
}

void LiquidRenderApp::formQuad(std::vector<unsigned int> &indexVector, std::vector<float> &vertexVector)
{
    if(vertexVector.size() < 4*5)
    {
        std::cout << "Unable to form quad, vertex vector has < 4 verts!";
        debug() << "Unable to form quad, vertex vector has < 4 verts!";
        return;
    }

    if(vertexVector.size() % 5 != 0)
    {
        std::cout << "Unable to form quad, vertex vector length % 5 != 0!";
        debug() << "Unable to form quad, vertex vector length % 5 != 0!";
        return;
    }

    int lastVertexIdx = (vertexVector.size() / 5) - 1;
    int topLeftIdx = lastVertexIdx - 3;
    int topRightIdx = lastVertexIdx - 2;
    int bottomRightIdx = lastVertexIdx - 1;
    int bottomLeftIdx = lastVertexIdx - 0;

    indexVector.push_back(topLeftIdx);
    indexVector.push_back(bottomRightIdx);
    indexVector.push_back(topRightIdx);
    indexVector.push_back(topLeftIdx);
    indexVector.push_back(bottomLeftIdx);
    indexVector.push_back(bottomRightIdx);
}

void LiquidRenderApp::updateFluidrenderQuad()
{
    glBindVertexArray(m_fluidgrid_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_fluidgrid_vbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_fluidgridQuadVerts.size(),m_fluidgridQuadVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void LiquidRenderApp::updateFluidrenderBuffers()
{
    updateFluidrenderQuad();
    glBindVertexArray(m_fluidgrid_vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_fluidgrid_ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int) * m_fluidgridQuadIndices.size(), m_fluidgridQuadIndices.data(),GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 5, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(float) * 5, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

void LiquidRenderApp::updateFluidrenderQuadVertex(Vertex v, int vertexIndex)
{
    int vertexParamStartIndex = vertexIndex*5;
    m_fluidgridQuadVerts[vertexParamStartIndex+0] = v.x();
    m_fluidgridQuadVerts[vertexParamStartIndex+1] = v.y();
}

void LiquidRenderApp::requestRender()
{
    m_renderRequested = true;
}

void LiquidRenderApp::render()
{
    if(!m_renderRequested) return;
    m_fluidRenderer.render();

    glViewport(0,0,m_windowWidth,m_windowHeight);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(.0f, .0f, .0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUseProgram(m_textureQuadShaderProgram);
    glBindVertexArray(m_fluidgrid_vao);
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, m_fluidRenderer.renderTexture());
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, m_fluidgridQuadIndices.size(),GL_UNSIGNED_INT,0);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBindVertexArray(0);

    m_textMenuRenderer.render();

    glfwSwapBuffers(m_window);
    m_renderRequested = false;
}

void LiquidRenderApp::resizeFluidrenderQuad()
{
    int workingWidth = m_windowWidth - m_textMenuRenderer.menuWidth();
    int maxGridDrawSize = std::min(workingWidth,m_windowHeight)*m_gridDrawFraction;
    float fgAsp = m_fluidRenderer.fluidGridAspect();
    int fluidDrawWidth = fgAsp > 1 ? maxGridDrawSize / fgAsp : maxGridDrawSize;
    int fluidDrawHeight = fgAsp > 1 ? maxGridDrawSize : maxGridDrawSize * fgAsp;
    updateFluidrenderQuadVertex(Vertex(static_cast<float>(fluidDrawWidth)*2 / m_windowWidth - 1, 1),1);
    updateFluidrenderQuadVertex(Vertex(static_cast<float>(fluidDrawWidth)*2 / m_windowWidth - 1, 1 - (static_cast<float>(fluidDrawHeight)*2 / m_windowHeight)),2);
    updateFluidrenderQuadVertex(Vertex(-1, 1 - (static_cast<float>(fluidDrawHeight)*2 / m_windowHeight)),3);
    updateFluidrenderQuad();
    m_fluidRenderer.resizeTexture(fluidDrawWidth,fluidDrawHeight);
}

template<class T>
T LiquidRenderApp::tryGetValue(json input, std::string key, T defaultValue)
{
    return input.contains(key) ? input[key].get<T>() : defaultValue;
}
