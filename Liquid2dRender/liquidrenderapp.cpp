#include "liquidrenderapp.h"
#include <stdexcept>
#include <cstdlib>

#include "linearindexable2d.h"
#include "customassert.h"

LiquidRenderApp* LiquidRenderApp::GLFWCallbackWrapper::s_application = nullptr;


LiquidRenderApp::LiquidRenderApp() :
    m_window(nullptr),
    m_solver(new FlipSolver(m_gridSizeI,m_gridSizeJ,1,1,1,1,false)),
    m_fluidRenderer(m_solver,800,600),
    m_textMenuRenderer(0,0,800,600,m_fluidRenderer)
{
    m_windowWidth = 800;
    m_windowHeight = 600;
    LiquidRenderApp::GLFWCallbackWrapper::SetApplication(this);
}

void LiquidRenderApp::init()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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

    m_fluidRenderer.init();
    m_textMenuRenderer.init();

    setupFluidrender();
    resizeFluidrenderQuad();
}

void LiquidRenderApp::run()
{
    while(!glfwWindowShouldClose(m_window))
    {
        render();
        glfwPollEvents();
    }

    glfwTerminate();
}

void LiquidRenderApp::resizeCallback(GLFWwindow *window, int width, int height)
{
    m_windowWidth = width;
    m_windowHeight = height;
    resizeFluidrenderQuad();
    m_textMenuRenderer.resize(width,height);
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
                        initGridForProjection();
                    }
                    else
                    {
                        m_solver->project();
                        m_fluidRenderer.updateGrid();
                    }
                    m_fluidRenderer.update();
                break;
            }
        break;
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

void LiquidRenderApp::render()
{
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
    Index2d fluidTopLeft(1,1);
    Index2d fluidBottomRight(99,99);
    m_solver->grid().fillMaterialRect(FluidCellMaterial::SOLID,0,0,99,99);
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
