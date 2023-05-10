#include "fluidrenderer.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <utility>

#include "flipsmokesolver.h"
#include "logger.h"
#include "materialgrid.h"
#include "glm/gtc/type_ptr.hpp"

#include "globalcallbackhandler.h"

#include "mathfuncs.h"

const Color FluidRenderer::m_emptyColor = Color(255,255,255);
const Color FluidRenderer::m_fluidLiquidColor = Color(44, 95, 150);
const Color FluidRenderer::m_solidColor = Color(94,94,94);
const Color FluidRenderer::m_sourceColor = Color(22, 196, 193);
const Color FluidRenderer::m_sinkColor = Color(11, 6, 64);
const Color FluidRenderer::m_velocityVectorColor = Color(255,255,255);
const Color FluidRenderer::m_geometryColor = Color(0,255,0);
const Color FluidRenderer::m_markerParticleFluidColor = Color(37,0,82);
const Color FluidRenderer::m_markerParticleAirColor = Color(77,190,255);

const char *FluidRenderer::m_vertexShaderSource =
        "#version 330 core\n"
        "layout (location = 0) in vec3 aPos;\n"
        "layout (location = 1) in vec3 vColor;\n"
        "uniform mat4 projection;"
        "out vec3 vertexColor;\n"
        "void main()\n"
        "{\n"
        "   vertexColor = vColor;\n"
        "   //vertexColor = aPos.xyz;\n"
        "   gl_Position = projection * vec4(aPos.xyz, 1.0);\n"
        "}\0";

const char *FluidRenderer::m_fragShaderSource =
        "#version 330 core\n"
        "out vec4 FragColor;\n"
        "in vec3 vertexColor;\n"
        "void main()\n"
        "{\n"
        "    FragColor = vec4(vertexColor.x, vertexColor.y, vertexColor.z, 1.0f);\n"
        "    //FragColor = vec4(vertexColor.x, vertexColor.y, vertexColor.z, 1.0f);\n"
        "}\0";

FluidRenderer::FluidRenderer(int textureWidth, int textureHeight) :
    m_gridVertexCount(0),
    m_textureWidth(textureWidth),
    m_textureHeight(textureHeight),
    m_gridRenderMode(FluidRenderMode::RENDER_MATERIAL),
    m_vectorRenderMode(VectorRenderMode::VECTOR_RENDER_CENTER),
    m_particleRenderMode(ParticleRenderMode::PARTICLE_RENDER_VELOCITY),
    m_vectorsEnabled(true),
    m_geometryEnabled(true),
    m_particlesEnabled(true),
    m_solver(nullptr)
{

}

void FluidRenderer::init(std::shared_ptr<FlipSolver> solver)
{
    m_solver = std::move(solver);
    float gridHeightF = m_solver->sizeI() * m_solver->dx();
    float gridWidthF = m_solver->sizeJ() * m_solver->dx();
    m_projection = glm::ortho(0.f,static_cast<float>(gridWidthF),
                                static_cast<float>(gridHeightF),0.f);
    m_gridVerts.reserve(m_solver.get()->cellCount()*m_vertexPerCell*m_vertexSize);
    m_gridIndices.reserve(m_solver.get()->cellCount()*m_vertexPerCell);

    m_vectorVerts.reserve(m_solver.get()->cellCount()*m_vectorVertsPerCell);

    int gridWidth = m_solver->gridSizeJ();
    int gridHeight = m_solver->gridSizeI();

    float dx = m_solver->dx();

    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float it = i * dx;
            float jt = j * dx;
            addGridQuad(Vertex(jt,it),Vertex(jt+dx,it+dx),Color(255,0,0));
            Vertex vStart = Vertex(jt + 0.5*dx,it + 0.5*dx);
            Vertex vEnd = Vertex(vStart.x() + 0.5*dx, vStart.y());
            addVector(vStart,vEnd,m_velocityVectorColor);
            setVectorColor(i,j,Color(255,0,0));
        }
    }
    loadGeometry();
    updateParticles();
    setupGl();
    updateGrid();
    updateVectors();
    setupGridVerts();
}

void FluidRenderer::render()
{
    glViewport(0,0,m_textureWidth,m_textureHeight);
    glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);
    glUseProgram(m_shaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "projection"),
                       1,
                       GL_FALSE,
                       glm::value_ptr(m_projection));
    //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glBindVertexArray(m_vao_grid);
    glDrawElements(GL_TRIANGLES,m_gridIndices.size(),GL_UNSIGNED_INT,0);
    if(m_vectorsEnabled)
    {
        glBindVertexArray(m_vao_vectors);
        glDrawArrays(GL_LINES,0,m_vectorVerts.size() / m_vertexSize);
    }
    if(m_geometryEnabled)
    {
        glBindVertexArray(m_vao_geometry);
        glDrawElements(GL_LINES,m_geometryIndices.size(),GL_UNSIGNED_INT,0);
    }
    if(m_particlesEnabled)
    {
        glPointSize(1);
        glBindVertexArray(m_vao_particles);
        glDrawArrays(GL_POINTS,0,m_particleVerts.size() / m_vertexSize);
    }
    glBindVertexArray(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FluidRenderer::update()
{
    updateGrid();
    updateVectors();
    updateParticles();
}

void FluidRenderer::addGridVertex(Vertex v, Color c)
{
    m_gridVerts.push_back(v.x());
    m_gridVerts.push_back(v.y());
    m_gridVerts.push_back(0.0f);
    m_gridVerts.push_back(c.rf());
    m_gridVerts.push_back(c.gf());
    m_gridVerts.push_back(c.bf());
    m_gridVertexCount++;
}

void FluidRenderer::addGridQuad(Vertex topLeft, Vertex bottomRight, Color c)
{
    int prevVertexTopIndex = m_gridVertexCount - 1;
    addGridVertex(topLeft,c);
    addGridVertex(Vertex(bottomRight.x(),topLeft.y()),c);
    addGridVertex(bottomRight,c);
    addGridVertex(Vertex(topLeft.x(),bottomRight.y()),c);
    m_gridIndices.push_back(prevVertexTopIndex+1);
    m_gridIndices.push_back(prevVertexTopIndex+3);
    m_gridIndices.push_back(prevVertexTopIndex+2);
    m_gridIndices.push_back(prevVertexTopIndex+1);
    m_gridIndices.push_back(prevVertexTopIndex+4);
    m_gridIndices.push_back(prevVertexTopIndex+3);
}

void FluidRenderer::addVectorVertex(Vertex v, Color c)
{
    m_vectorVerts.push_back(v.x());
    m_vectorVerts.push_back(v.y());
    m_vectorVerts.push_back(0.8f);//offset to draw vectors over the grid
    m_vectorVerts.push_back(c.rf());
    m_vectorVerts.push_back(c.gf());
    m_vectorVerts.push_back(c.bf());
}

void FluidRenderer::addVector(Vertex start, Vertex end, Color c)
{
    addVectorVertex(start,c);
    addVectorVertex(end,c);
}

void FluidRenderer::addParticle(Vertex particle, Color c)
{
    m_particleVerts.push_back(particle.y()*m_solver->dx());
    m_particleVerts.push_back(particle.x()*m_solver->dx());
    m_particleVerts.push_back(0.9f);
    m_particleVerts.push_back(c.rf());
    m_particleVerts.push_back(c.gf());
    m_particleVerts.push_back(c.bf());
}

void FluidRenderer::setupGl()
{
    int success;
    m_vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(m_vertexShader,1,&m_vertexShaderSource,NULL);
    glCompileShader(m_vertexShader);
    glGetShaderiv(m_vertexShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Vertex shader compilation error";
    }

    m_fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(m_fragShader, 1, &m_fragShaderSource, NULL);
    glCompileShader(m_fragShader);
    glGetShaderiv(m_vertexShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Fragment shader compilation error";
    }

    m_shaderProgram = glCreateProgram();
    glAttachShader(m_shaderProgram,m_vertexShader);
    glAttachShader(m_shaderProgram,m_fragShader);
    glLinkProgram(m_shaderProgram);
    glGetProgramiv(m_shaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        std::cout << "Shader program linking error!";
    }

    glDeleteShader(m_vertexShader);
    glDeleteShader(m_fragShader);

    //Grid buffers
    glGenVertexArrays(1, &m_vao_grid);
    glGenBuffers(1,&m_vbo_grid);
    glGenBuffers(1,&m_ebo_grid);

    //Vector buffers
    glGenVertexArrays(1, &m_vao_vectors);
    glGenBuffers(1,&m_vbo_vectors);

    //Geometry buffers
    glGenVertexArrays(1, &m_vao_geometry);
    glGenBuffers(1,&m_vbo_geometry);
    glGenBuffers(1,&m_ebo_geometry);

    //Particle buffers
    glGenVertexArrays(1, &m_vao_particles);
    glGenBuffers(1,&m_vbo_particles);

    setupBuffers();
    //Frame buffers
    glGenFramebuffers(1,&m_framebuffer);

    //Textures
    glGenTextures(1, &m_renderTexture);

    //Render buffers
    glGenRenderbuffers(1, &m_renderbuffer);
    setupOffscreenBuffer();
}

void FluidRenderer::setupOffscreenBuffer()
{
    glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer);

    glBindTexture(GL_TEXTURE_2D, m_renderTexture);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,m_textureWidth,m_textureHeight,0,GL_RGB,GL_UNSIGNED_BYTE,NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_renderTexture, 0);

    glBindRenderbuffer(GL_RENDERBUFFER, m_renderbuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, m_textureWidth, m_textureHeight);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, m_renderbuffer);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete! Status = " << glCheckFramebufferStatus(GL_FRAMEBUFFER);
        debug() << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!\n" << glCheckFramebufferStatus(GL_FRAMEBUFFER);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FluidRenderer::setupBuffers()
{
    setupGridVerts();
    glBindVertexArray(m_vao_grid);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_grid);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int) * m_gridIndices.size(), m_gridIndices.data(),GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    setupVectorVerts();
    glBindVertexArray(m_vao_vectors);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    setupGeometryVerts();
    glBindVertexArray(m_vao_geometry);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_geometry);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int) * m_geometryIndices.size(), m_geometryIndices.data(),GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    setupParticleVerts();
    glBindVertexArray(m_vao_particles);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

void FluidRenderer::setupGridVerts()
{
    glBindVertexArray(m_vao_grid);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_grid);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_gridVerts.size(),m_gridVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void FluidRenderer::setupVectorVerts()
{
    glBindVertexArray(m_vao_vectors);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_vectors);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_vectorVerts.size(),m_vectorVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void FluidRenderer::setupGeometryVerts()
{
    glBindVertexArray(m_vao_geometry);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_geometry);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_geometryVerts.size(),m_geometryVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void FluidRenderer::setupParticleVerts()
{
    glBindVertexArray(m_vao_particles);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_particles);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_particleVerts.size(),m_particleVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void FluidRenderer::loadGeometry()
{
    for(Obstacle &g : m_solver->geometryObjects())
    {
        addGeometry(g.geometry());
    }
}

void FluidRenderer::addGeometry(Geometry2d &geometry)
{
    m_geometryVerts.reserve(m_geometryVerts.size() + geometry.vertextCount() * m_vertexSize);
    m_geometryIndices.reserve(m_geometryVerts.size() + geometry.vertextCount()*2);

    int startVertexIndex = m_geometryVerts.size() / m_vertexSize;
    int endVertexIndex = startVertexIndex;
    for(Vertex& v : geometry.verts())
    {
        //X and Y are swapped, because i is x and vertical, and j is y and horizontal, but in opengl x is horizontal and y is vertical
        m_geometryVerts.push_back(v.y());
        m_geometryVerts.push_back(v.x());
        m_geometryVerts.push_back(0.9f);
        m_geometryVerts.push_back(m_geometryColor.rf());
        m_geometryVerts.push_back(m_geometryColor.gf());
        m_geometryVerts.push_back(m_geometryColor.bf());
        endVertexIndex++;
    }

    for(int i = startVertexIndex; i < endVertexIndex; i++)
    {
        if(i == endVertexIndex - 1)
        {
            m_geometryIndices.push_back(i);
            m_geometryIndices.push_back(startVertexIndex);
            break;
        }
        m_geometryIndices.push_back(i);
        m_geometryIndices.push_back(i+1);
    }
}

void FluidRenderer::updateGridVerts()
{
    glBindVertexArray(m_vao_grid);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_grid);
    glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float) * m_gridVerts.size(),m_gridVerts.data());
    glBindVertexArray(0);
}

void FluidRenderer::updateVectorVerts()
{
    glBindVertexArray(m_vao_vectors);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_vectors);
    glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float) * m_vectorVerts.size(),m_vectorVerts.data());
    glBindVertexArray(0);
}

void FluidRenderer::updateGeometryVerts()
{
    glBindVertexArray(m_vao_geometry);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_geometry);
    glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float) * m_geometryVerts.size(),m_geometryVerts.data());
    glBindVertexArray(0);
}

void FluidRenderer::updateParticleVerts()
{
    glBindVertexArray(m_vao_particles);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_particles);
    glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(float) * m_particleVerts.size(),m_particleVerts.data());
    glBindVertexArray(0);
}

void FluidRenderer::updateGrid()
{
    m_solver->updateSinks();
    m_solver->updateSources();
    GlobalCallbackHandler::instance().registerRenderUpdate();
    switch(m_gridRenderMode)
    {
    case FluidRenderMode::RENDER_MATERIAL:
        updateGridFromMaterial();
        break;
    case FluidRenderMode::RENDER_VELOCITY:
        updateGridFromVelocity();
        break;
    case FluidRenderMode::RENDER_TEST:
        updateGridFromTestGrid();
        break;
    case FluidRenderMode::RENDER_U:
        updateGridFromUComponent();
        break;
    case FluidRenderMode::RENDER_V:
        updateGridFromVComponent();
        break;
    case FluidRenderMode::RENDER_OBSTACLE_SDF:
        updateGridFromObstacleSdf();
        break;
    case FluidRenderMode::RENDER_FLUID_SDF:
        updateGridFromFluidSdf();
        break;

    default:
        std::cout << "Bad render mode";
        break;
    }
    updateGridVerts();
}

unsigned int FluidRenderer::renderTexture()
{
    return m_renderTexture;
}

void FluidRenderer::resizeTexture(int width, int height)
{
    m_textureHeight = height;
    m_textureWidth = width;
    setupOffscreenBuffer();
}

float FluidRenderer::fluidGridAspect()
{
    return static_cast<float>(m_solver->gridSizeI()) / m_solver->gridSizeJ();
}

std::shared_ptr<FlipSolver> FluidRenderer::solver()
{
    return m_solver;
}

void FluidRenderer::updateGridFromMaterial()
{
    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            switch(m_solver->materialGrid().at(i,j))
            {
                case FluidMaterial::EMPTY:
                    setCellColor(i,j,m_emptyColor);
                break;

                case FluidMaterial::FLUID:
                    if(m_solver->simulationMethod() == SimulationMethod::SIMULATION_LIQUID
                    || m_solver->simulationMethod() == SimulationMethod::SIMULATION_NBFLIP)
                    {
                        setCellColor(i,j,m_fluidLiquidColor);
                    }
                    else
                    {
                        FlipSmokeSolver* solver = dynamic_cast<FlipSmokeSolver*>(m_solver.get());
                        float smokeConcentration = solver->smokeConcentration().at(i,j);
                        float opacity = simmath::lerp(0.01f,1.f,smokeConcentration);
                        float temp = solver->temperature().at(i,j);
                        float intensity = std::clamp((std::clamp(temp, 773.f,temp) - 773.f) / 277.f, 0.f, 1.f);
                        Color c = getBlackbodyColor(temp);
                        c = c*intensity;
                        c = Color::lerp(m_emptyColor,c,opacity);
                        setCellColor(i,j,c);
                    }
                break;

                case FluidMaterial::SOLID:
                    setCellColor(i,j,m_solidColor);
                break;

                case FluidMaterial::SOURCE:
                    setCellColor(i,j,m_sourceColor);
                break;

                case FluidMaterial::SINK:
                    setCellColor(i,j,m_sinkColor);
                break;

                default:
                    std::cout << "Unknown material " << m_solver->materialGrid().at(i,j) << "at i,j " << i << "," << j;
                break;
            }
        }
    }
}

void FluidRenderer::updateGridFromVelocity()
{
    //float min = -10;
    //float diff = max - min;

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            //float r = (m_solver->fluidVelocityGrid().getU(i,j) - min) / diff;
            //float g = (m_solver->fluidVelocityGrid().getV(i,j) - min) / diff;
            //float r = (std::abs(m_solver->fluidVelocityGrid().getU(i,j))) / m_velocityComponentRangeMax;
            //float g = (std::abs(m_solver->fluidVelocityGrid().getV(i,j))) / m_velocityComponentRangeMax;
            Vertex velocity = m_solver->fluidVelocityGrid().velocityAt(static_cast<float>(i)+0.5f,
                                                                       static_cast<float>(j)+0.5f);
            Color c = hueColorRamp(velocity.distFromZero() / m_velocityRangeMax * m_solver->dx());
            setCellColor(i,j,c);
        }
    }
}

void FluidRenderer::updateGridFromTestGrid()
{
    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float brightness = m_solver->testGrid().at(i,j);
            //setCellColor(i,j,Color(brightness,brightness,brightness));
            if(brightness > 0)
            {
                setCellColor(i,j,Color(255,255,255) * brightness);
            }
            else
            {
                setCellColor(i,j,Color(0,0,255) * -brightness);
            }
        }
    }
}

void FluidRenderer::updateGridFromUComponent()
{
//    auto minMaxValues = std::minmax_element(m_solver->grid().data().begin(),
//                                            m_solver->grid().data().end(),
//                                            [] (const FluidCell& a, const FluidCell& b)
//                                            {
//                                                return a.getFluidU() < b.getFluidU();
//                                            });
//    float min = minMaxValues.first->getFluidU();
//    float max = minMaxValues.second->getFluidU();

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float brightness = std::clamp(std::abs(m_solver->fluidVelocityGrid().getU(i,j))
                    / m_velocityComponentRangeMax,0.f,1.f);
            setCellColor(i,j,Color(brightness,brightness,brightness));
        }
    }
}

void FluidRenderer::updateGridFromVComponent()
{

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float brightness = std::clamp(std::abs(m_solver->fluidVelocityGrid().getV(i,j))
                    / m_velocityComponentRangeMax,0.f,1.f);
            setCellColor(i,j,Color(brightness,brightness,brightness));
        }
    }
}

void FluidRenderer::updateGridFromObstacleSdf()
{
    float max = std::max(m_solver->gridSizeI(),m_solver->gridSizeJ());

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float dist = m_solver->solidSdf().at(i,j);
            float brightness = std::pow(std::abs(dist) / max,0.4);
//            if (std::abs(dist - 10.f) < 1e-0f)
//            {
//                setCellColor(i,j,Color(255,255,255));
//                continue;
//            }
            if(dist <= 0.f)
            {
                setCellColor(i,j,Color(static_cast<int>(87 * brightness), static_cast<int>(202 * brightness), static_cast<int>(255 * brightness)));
            }
            else
            {
                setCellColor(i,j,Color(static_cast<int>(255 * brightness), static_cast<int>(168 * brightness), static_cast<int>(87 * brightness)));
            }
        }
    }
}

void FluidRenderer::updateGridFromFluidSdf()
{
    float max = std::max(m_solver->gridSizeI(),m_solver->gridSizeJ());

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float dist = m_solver->fluidSdf().at(i,j);
            float brightness = std::pow(std::abs(dist) / max,0.4);
            if(brightness > 1.f)
            {
                std::cout << "bad fluid sdf " << dist << " at " << i << ' ' << j << '\n';
            }
//            if(dist < -1.f)
//            {
//                setCellColor(i,j,Color(255,255,255));
//                continue;
//            }
            if (std::abs(dist) < 1.f)
            {
                if(dist <= 0.f)
                {
                    setCellColor(i,j,Color(static_cast<int>(20 * brightness), static_cast<int>(80 * brightness), static_cast<int>(255 * brightness)));
                }
                else
                {
                    setCellColor(i,j,Color(static_cast<int>(255 * brightness), static_cast<int>(80 * brightness), static_cast<int>(20 * brightness)));
                }
            }
            else
            {
                if(dist <= 0.f)
                {
                    setCellColor(i,j,Color(static_cast<int>(87 * brightness), static_cast<int>(202 * brightness), static_cast<int>(255 * brightness)));
                }
                else
                {
                    setCellColor(i,j,Color(static_cast<int>(255 * brightness), static_cast<int>(168 * brightness), static_cast<int>(87 * brightness)));
                }
            }
        }
    }
}



void FluidRenderer::updateVectors()
{
    GlobalCallbackHandler::instance().registerRenderUpdate();
    switch(m_vectorRenderMode)
    {
    case VECTOR_RENDER_CENTER:
        updateVectorsCentered();
        break;
    case VECTOR_RENDER_STAGGERED:
        updateVectorsStaggered();
        break;
    case VECTOR_RENDER_SOLID_SDF_GRADIENT:
        updateVectorsSolidSdfGrad();
        break;
    case VECTOR_RENDER_FLUID_SDF_GRADIENT:
        updateVectorsFluidSdfGrad();
        break;
    case VECTOR_RENDER_ITER_END:
        std::cout << "Invalid vector render mode!";
        break;
    }
    updateVectorVerts();
}

void FluidRenderer::updateParticles()
{
    switch (m_particleRenderMode) {
    case PARTICLE_RENDER_VELOCITY:
        reloadParticlesVelocity();
        break;
    case PARTICLE_RENDER_SOLID:
        reloadParticlesSolid();
        break;
    case PARTICLE_RENDER_TEST_VALUE:
        reloadParticlesFromTestValue();
        break;
    case PARTICLE_RENDER_ITER_END:
        std::cout << "Invalid particle render mode!";
        break;
    }
}

void FluidRenderer::updateVectorsStaggered()
{
    //float unitLengthVelocity = 0.75; //Velocity len(v) at which draw vector will be exactly length 1 grid cell side

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    Vertex vEnd;
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            Vertex gridspaceVelocity(m_solver->fluidVelocityGrid().getU(i,j),m_solver->fluidVelocityGrid().getV(i,j));
            //Vertex gridspaceVelocity(1,0);
            float scaleFactor = gridspaceVelocity.distFromZero() / m_solver->dx();
            //float scaleFactor = 1;
            Vertex newVector = Vertex(gridspaceVelocity.y() / scaleFactor,
                                       gridspaceVelocity.x() / scaleFactor);
            updateVector(i,j,newVector);
            setVectorColor(i,j,hueColorRamp(gridspaceVelocity.distFromZero() / m_velocityRangeMax));
        }
    }
}

void FluidRenderer::updateVectorsCentered()
{
    //float unitLengthVelocity = 0.75;

    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    Vertex vEnd;
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            Vertex gridspaceVelocity = m_solver->fluidVelocityGrid().velocityAt(static_cast<float>(i) + 0.5,
                                                        static_cast<float>(j) + 0.5);
            //Vertex gridspaceVelocity(1,0);
            float scaleFactor = gridspaceVelocity.distFromZero() / m_solver->dx();
            //float scaleFactor = 1;
            Vertex newVector = Vertex((gridspaceVelocity.y() / scaleFactor),
                                       (gridspaceVelocity.x() / scaleFactor));
            updateVector(i,j,newVector);
            setVectorColor(i,j,hueColorRamp(gridspaceVelocity.distFromZero() / m_velocityRangeMax));
        }
    }
}

void FluidRenderer::updateVectorsSolidSdfGrad()
{
    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    Vertex vEnd;
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            Vertex gridspaceSdf = simmath::gradCenteredGrid(i,j, m_solver->solidSdf());
            //Vertex gridspaceVelocity(1,0);
            float scaleFactor = gridspaceSdf.distFromZero() / m_solver->dx();
            //float scaleFactor = 1;
            Vertex newVector = Vertex((gridspaceSdf.y() / scaleFactor),
                                       (gridspaceSdf.x() / scaleFactor));
            updateVector(i,j,newVector);
        }
    }
}

void FluidRenderer::updateVectorsFluidSdfGrad()
{
    int gridHeight = m_solver->sizeI();
    int gridWidth = m_solver->sizeJ();
    Vertex vEnd;
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            Vertex gridspaceSdfGrad = simmath::gradCenteredGrid(i,j, m_solver->solidSdf());
            //Vertex gridspaceVelocity(1,0);
            float scaleFactor = gridspaceSdfGrad.distFromZero() / m_solver->dx();
            //float scaleFactor = 1;
            Vertex newVector = Vertex((gridspaceSdfGrad.y() / scaleFactor),
                                       (gridspaceSdfGrad.x() / scaleFactor));
            updateVector(i,j,newVector);
            setVectorColor(i,j,hueColorRamp(gridspaceSdfGrad.distFromZero() / 1.f));
        }
    }
}

void FluidRenderer::reloadParticlesSolid()
{
    int oldParticlesSize = m_particleVerts.size();
    m_particleVerts.clear();
    for(MarkerParticle &particle : m_solver->markerParticles())
    {
        Color particleColor = particle.material == FluidMaterial::FLUID? m_markerParticleFluidColor :
                                                                         m_markerParticleAirColor;
        addParticle(particle.position,particleColor);
    }
    if(m_particleVerts.size() > oldParticlesSize)
    {
        setupParticleVerts();
    }
    else
    {
        updateParticleVerts();
    }
}

void FluidRenderer::reloadParticlesVelocity()
{
    int oldParticlesSize = m_particleVerts.size();
    m_particleVerts.clear();
    for(MarkerParticle &particle : m_solver->markerParticles())
    {
        addParticle(particle.position,hueColorRamp(particle.velocity.distFromZero() / m_velocityRangeMax * m_solver->dx()));
    }
    if(m_particleVerts.size() > oldParticlesSize)
    {
        setupParticleVerts();
    }
    else
    {
        updateParticleVerts();
    }
}

void FluidRenderer::reloadParticlesFromTestValue()
{
    int oldParticlesSize = m_particleVerts.size();
    m_particleVerts.clear();
    for(MarkerParticle &particle : m_solver->markerParticles())
    {
        Color pColor = Color(particle.testValue,particle.testValue,particle.testValue);
        addParticle(particle.position,pColor);
    }
    if(m_particleVerts.size() > oldParticlesSize)
    {
        setupParticleVerts();
    }
    else
    {
        updateParticleVerts();
    }
}

void FluidRenderer::updateVector(int x, int y, Vertex newVector)
{
    int linearVectorStartVertexIndex = m_solver->linearIndex(x,y) * m_vectorVertsPerCell;
    int linearVectorEndVertexIndex = linearVectorStartVertexIndex + 1;
    int linearStartVertexCoordIndex = linearVectorStartVertexIndex * m_vertexSize;
    int linearEndVertexCoordIndex = linearVectorEndVertexIndex * m_vertexSize;

    m_vectorVerts[linearEndVertexCoordIndex] = m_vectorVerts[linearStartVertexCoordIndex] + newVector.x();
    m_vectorVerts[linearEndVertexCoordIndex + 1] = m_vectorVerts[linearStartVertexCoordIndex + 1] + newVector.y();
    m_vectorVerts[linearEndVertexCoordIndex + 2] = m_vectorVerts[linearStartVertexCoordIndex + 2] + newVector.z();
}

void FluidRenderer::setCellVertexColor(int vIndex, Color c)
{
    int linearIndex = vIndex * m_vertexSize + 3;
    m_gridVerts[linearIndex] = c.rf();
    m_gridVerts[linearIndex + 1] = c.gf();
    m_gridVerts[linearIndex + 2] = c.bf();
}

void FluidRenderer::setCellColor(int x, int y, Color c)
{
    int linearIndex = m_solver->linearIndex(x,y) * m_vertexPerCell;
    setCellVertexColor(linearIndex,c);
    setCellVertexColor(linearIndex+1,c);
    setCellVertexColor(linearIndex+2,c);
    setCellVertexColor(linearIndex+3,c);
}

void FluidRenderer::setVectorVertexColor(int vIndex, Color c)
{
    int linearIndex = vIndex * m_vertexSize + 3;
    m_vectorVerts[linearIndex + 0] = c.rf();
    m_vectorVerts[linearIndex + 1] = c.gf();
    m_vectorVerts[linearIndex + 2] = c.bf();
}

void FluidRenderer::setVectorColor(int x, int y, Color c)
{
    int linearIndex = m_solver->linearIndex(x,y) * m_vectorVertsPerCell;
    setVectorVertexColor(linearIndex, c);
    setVectorVertexColor(linearIndex + 1, c);
}

void FluidRenderer::setParticleColor(int idx, Color c)
{
    int linearIndex = idx*6 + 3;
    m_particleVerts[linearIndex + 0] = c.rf();
    m_particleVerts[linearIndex + 1] = c.gf();
    m_particleVerts[linearIndex + 2] = c.bf();
}

Color FluidRenderer::hueColorRamp(float val)
{
    float hue = std::clamp(val,0.f,1.f) * (m_hueMaxVelocity - m_hueMinVelocity) + m_hueMinVelocity;
    return Color::fromHSVA(hue,1.f,1.f);
}

Color FluidRenderer::velocityComponentColorRamp(float val)
{
    float hue = std::clamp(val,0.f,1.f) * (m_hueMaxVelocityComponent - m_hueMinVelocityComponent) + m_hueMinVelocityComponent;
    return Color::fromHSVA(hue,1.f,1.f);
}

//Taken from https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
Color FluidRenderer::getBlackbodyColor(float temp)
{
    int r = 0;
    int g = 0;
    int b = 0;
    temp /= 100;
    if(temp <= 66)
    {
        r = 255;
    }
    else
    {
        r = temp - 60;
        r = 329.698727446 * (std::pow(r,-0.1332047592));
        r = std::clamp(r,0,255);
    }

    if(temp <= 66)
    {
        g = temp;
        g = 99.4708025861 * std::log(g) - 161.1195681661;
        g = std::clamp(g,0,255);
    }
    else
    {
        g = temp - 60;
        g = 288.1221695283 * std::pow(g,-0.0755148492);
        g = std::clamp(g,0,255);
    }

    if(temp >= 66 )
    {
        b = 255;
    }
    else
    {

        if(temp <= 19)
        {
            b = 0;
        }
        else
        {
            b = temp - 10;
            b = 138.5177312231 * std::log(b) - 305.0447927307;
            b = std::clamp(b,0,255);
        }
    }

    return Color(r,g,b);
}

void FluidRenderer::dumpToTga(const std::string& fileName)
{
    FILE *fp = fopen(("./output/" + fileName).c_str(), "wb");
    if (!fp) {
        std::cout << "cannot open for tga frame dump: " << fileName << '\n';
        return;
    }

    std::vector<GLubyte> data(m_textureWidth*m_textureHeight*4,0);

    render();
    glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer);
    glReadPixels(0, 0, m_textureWidth, m_textureHeight, GL_BGRA, GL_UNSIGNED_BYTE, data.data());
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    std::array<GLubyte,18> header = {00,00,02, 00,00,00, 00,00,00, 00,00,00,
                        static_cast<GLubyte>(0xff & m_textureWidth),
                        static_cast<GLubyte>(0xff & m_textureWidth >> 8),
                        static_cast<GLubyte>(0xff & m_textureHeight),
                        static_cast<GLubyte>(0xff & m_textureHeight >> 8),
                        32, 0x00 }; // next-to-last = bit depth
    fwrite( header.data(), header.size(), 1, fp);
    fwrite(data.data(), data.size(), 1, fp);
    fclose(fp);
}
