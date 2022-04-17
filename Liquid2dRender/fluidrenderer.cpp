#include "fluidrenderer.h"

#include <iostream>
#include <cmath>
#include <algorithm>

const Color FluidRenderer::m_emptyColor = Color(255,255,255);
const Color FluidRenderer::m_fluidColor = Color(44, 95, 150);
const Color FluidRenderer::m_solidColor = Color(94,94,94);

const char *FluidRenderer::m_vertexShaderSource =
        "#version 330 core\n"
        "layout (location = 0) in vec3 aPos;\n"
        "layout (location = 1) in vec3 vColor;\n"
        "out vec3 vertexColor;\n"
        "void main()\n"
        "{\n"
        "   vertexColor = vColor;\n"
        "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
        "}\0";

const char *FluidRenderer::m_fragShaderSource =
        "#version 330 core\n"
        "out vec4 FragColor;\n"
        "in vec3 vertexColor;\n"
        "void main()\n"
        "{\n"
        "    FragColor = vec4(vertexColor.x, vertexColor.y, vertexColor.z, 1.0f);\n"
        "}\0";

FluidRenderer::FluidRenderer(std::shared_ptr<FlipSolver> solver) :
    m_solver(solver),
    m_vertexCount(0)
{
    m_verts.reserve(solver.get()->grid().cellCount()*m_vertexPerCell*m_vertexSize);

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    float cellSize = 2.f / std::max(gridWidth,gridHeight);

    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            addQuad(Vertex(cellSize*j - 1,1 - cellSize*i),Vertex(cellSize*(j+1) - 1,1 - cellSize*(i+1)),Color(255,0,0));
        }
    }
}

void FluidRenderer::init()
{
    setupGl();
    updateGrid();
    updateVerts();
}

void FluidRenderer::render()
{
    glUseProgram(m_shaderProgram);
    //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glBindVertexArray(m_vao);
    glDrawElements(GL_TRIANGLES,m_indices.size(),GL_UNSIGNED_INT,0);
}

void FluidRenderer::addVertex(Vertex v, Color c)
{
    m_verts.push_back(v.x());
    m_verts.push_back(v.y());
    m_verts.push_back(v.z());
    m_verts.push_back(c.rf());
    m_verts.push_back(c.gf());
    m_verts.push_back(c.bf());
    m_vertexCount++;
}

void FluidRenderer::addQuad(Vertex topLeft, Vertex bottomRight, Color c)
{
    int prevVertexTopIndex = m_vertexCount - 1;
    addVertex(topLeft,c);
    addVertex(Vertex(bottomRight.x(),topLeft.y()),c);
    addVertex(bottomRight,c);
    addVertex(Vertex(topLeft.x(),bottomRight.y()),c);
    m_indices.push_back(prevVertexTopIndex+1);
    m_indices.push_back(prevVertexTopIndex+3);
    m_indices.push_back(prevVertexTopIndex+2);
    m_indices.push_back(prevVertexTopIndex+1);
    m_indices.push_back(prevVertexTopIndex+4);
    m_indices.push_back(prevVertexTopIndex+3);
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

    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1,&m_vbo);
    glGenBuffers(1,&m_ebo);
    updateBuffers();
}

void FluidRenderer::updateBuffers()
{
    glBindVertexArray(m_vao);
    updateVerts();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int) * m_indices.size(), m_indices.data(),GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
}

void FluidRenderer::updateVerts()
{
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_verts.size(),m_verts.data(),GL_DYNAMIC_DRAW);
}

void FluidRenderer::updateGrid()
{
    switch(m_renderMode)
    {
    case FluidRenderMode::RENDER_MATERIAL:
        updateGridFromMaterial();
        break;
    case FluidRenderMode::RENDER_VELOCITY:
        updateGridFromVelocity();
        break;
    case FluidRenderMode::RENDER_U:
        updateGridFromUComponent();
        break;
    case FluidRenderMode::RENDER_V:
        updateGridFromVComponent();
        break;

    default:
        std::cout << "Bad render mode";
        break;
    }
    updateVerts();
}

void FluidRenderer::updateGridFromMaterial()
{
    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            switch(m_solver->grid().getMaterial(i,j))
            {
                case FluidCellMaterial::EMPTY:
                    setColor(i,j,m_emptyColor);
                break;

                case FluidCellMaterial::FLUID:
                    setColor(i,j,m_fluidColor);
                break;

                case FluidCellMaterial::SOLID:
                    setColor(i,j,m_solidColor);
                break;

                default:
                    std::cout << "Unknown material " << m_solver->grid().getMaterial(i,j) << "at i,j " << i << "," << j;
                break;
            }
        }
    }
}

void FluidRenderer::updateGridFromVelocity()
{
    float min = -10;
    float max = 10;
    float diff = max - min;

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float r = (m_solver->grid().getU(i,j) - min) / diff;
            float g = (m_solver->grid().getV(i,j) - min) / diff;
            setColor(i,j,Color(r,g,0.0));
        }
    }
}

void FluidRenderer::updateGridFromUComponent()
{
//    auto minMaxValues = std::minmax_element(m_solver->grid().data().begin(),
//                                            m_solver->grid().data().end(),
//                                            [] (const FluidCell& a, const FluidCell& b)
//                                            {
//                                                return a.getU() < b.getU();
//                                            });
//    float min = minMaxValues.first->getU();
//    float max = minMaxValues.second->getU();
    float min = -10;
    float max = 10;
    float diff = max - min;

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float brightness = (m_solver->grid().getU(i,j) - min) / diff;
            setColor(i,j,Color(brightness,brightness,brightness));
        }
    }
}

void FluidRenderer::updateGridFromVComponent()
{
    float min = -10;
    float max = 10;
    float diff = max - min;

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            float brightness = (m_solver->grid().getV(i,j) - min) / diff;
            setColor(i,j,Color(brightness,brightness,brightness));
        }
    }
}

void FluidRenderer::setVertexColor(int vIndex, Color c)
{
    int linearIndex = vIndex * 6 + 3;
    m_verts[linearIndex] = c.rf();
    m_verts[linearIndex + 1] = c.gf();
    m_verts[linearIndex + 2] = c.bf();
}

void FluidRenderer::setColor(int x, int y, Color c)
{
    int linearIndex = m_solver->grid().linearIndex(x,y) * 4;
    setVertexColor(linearIndex,c);
    setVertexColor(linearIndex+1,c);
    setVertexColor(linearIndex+2,c);
    setVertexColor(linearIndex+3,c);
}
