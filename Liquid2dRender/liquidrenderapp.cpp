#include "liquidrenderapp.h"
#include <stdexcept>
#include <cstdlib>
#include <filesystem>

#include "linearindexable2d.h"
#include "customassert.h"
#include "globalcallbackhandler.h"
#include "simsettings.h"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;


LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(new FlipSolver(1,true)),
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
    //loadJson("./scenes/waterfall.json");
    //loadJson("./scenes/dam_break.json");
    //loadJson("./scenes/test_scene.json");
    loadJson("./scenes/viscosity_test.json");
    int foo = 11;
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

    m_solver->init();
    m_fluidRenderer.init(m_solver);
    m_textMenuRenderer.init();

    setupFluidrender();
    resizeFluidrenderQuad();
    m_solver->updateSdf();
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
                    m_fluidRenderer.gridRenderMode()++;
                    m_fluidRenderer.updateGrid();
                break;

                case GLFW_KEY_A:
                    m_solver->advect();
                    m_fluidRenderer.update();
                break;

                case GLFW_KEY_T:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        initGridForProjection();
                    }
                    else
                    {
                        m_solver->project();
                        m_fluidRenderer.updateGrid();
                    }
                    m_fluidRenderer.update();
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

                case GLFW_KEY_E:
                    if(mods & GLFW_MOD_SHIFT)
                    {
                        initGridForExtrapolation();
                    }
                    else
                    {
                        m_solver->extrapolateVelocityField();
                    }
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

                case GLFW_KEY_S:
                if(mods & GLFW_MOD_SHIFT)
                {
                    for(int i = 0; i < SimulationStepStage::STAGE_ITER_END; i++)
                    {
                        m_solver->stagedStep();
                    }
                }
                else
                {
                    m_solver->stagedStep();
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
                        m_solver->step();
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
    std::ifstream sceneFile(fileName);
    if(!sceneFile.is_open())
    {
        std::cout << "errorOpening scene file " << fileName;
    }
    json sceneJson;
    sceneFile >> sceneJson;
    settingsFromJson(sceneJson["settings"]);
    solverFromJson(sceneJson["solver"]);
}

void LiquidRenderApp::settingsFromJson(json settingsJson)
{
    SimSettings::domainSizeI() = settingsJson["domainSizeI"].get<int>();
    SimSettings::domainSizeJ() = settingsJson["domainSizeJ"].get<int>();
    SimSettings::resolution() = settingsJson["resolution"].get<int>();
    SimSettings::fps() = settingsJson["fps"].get<int>();
    SimSettings::frameDt() = 1.f / SimSettings::fps();
    SimSettings::maxSubsteps() = settingsJson["maxSubsteps"].get<int>();
    SimSettings::stepDt() = SimSettings::frameDt() / SimSettings::maxSubsteps();
    SimSettings::density() = 1;
    SimSettings::randomSeed() = settingsJson["seed"].get<int>();
    SimSettings::particlesPerCell() = settingsJson["particlesPerCell"].get<int>();
    SimSettings::cflNumber() = settingsJson["cflNumber"].get<float>();
    SimSettings::picRatio() = settingsJson["picRatio"].get<float>();
    std::pair<float,float> v = settingsJson["globalAcceleration"]
                                    .get<std::pair<float,float>>();
    SimSettings::globalAcceleration() = Vertex(v.first,v.second);
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
        addGeometryFromJson(geo);
    }
}

void LiquidRenderApp::addGeometryFromJson(json geometryJson)
{
    std::vector<std::pair<float,float>> verts = geometryJson["verts"]
                                                .get<std::vector<std::pair<float,float>>>();
    std::string geoType = geometryJson["type"].get<std::string>();
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(Vertex(v.first,v.second));
    }

    if(geoType == "solid")
    {
        m_solver->addGeometry(geo);
    }
    else if(geoType == "source")
    {
        m_solver->addSource(geo);
    }
    else if(geoType == "sink")
    {
        m_solver->addSink(geo);
    }
    else if(geoType == "fluid")
    {
        m_solver->addInitialFluid(geo);
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
    Index2d fluidTopLeft(15,0);
    Index2d fluidBottomRight(99,99);
    static int seed = 0;
    srand(seed);
    rand();
    seed++;
    //m_solver->grid().fillMaterialRect(FluidCellMaterial::SOLID,0,0,99,99);
    //m_solver->grid().fillMaterialRect(FluidCellMaterial::FLUID,fluidTopLeft,fluidBottomRight);
    float u = (static_cast<float>(rand()) / RAND_MAX) - 0.5f;
    float v = (static_cast<float>(rand()) / RAND_MAX) - 0.5f;
    for (int i = 0; i < m_solver->grid().sizeI(); i++)
    {
        for (int j = 0; j < m_solver->grid().sizeJ(); j++)
        {
            //m_solver->grid().setU(i,j,static_cast<float>(rand()) / RAND_MAX);
            //m_solver->grid().setV(i,j,static_cast<float>(rand()) / RAND_MAX);
            m_solver->updateSolids();
            m_solver->grid().setU(i,j,0);
            m_solver->grid().setV(i,j,10);
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
