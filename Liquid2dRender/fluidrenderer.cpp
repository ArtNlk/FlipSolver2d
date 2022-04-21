#include "fluidrenderer.h"

#include <iostream>
#include <cmath>
#include <algorithm>

const Color FluidRenderer::m_emptyColor = Color(255,255,255);
const Color FluidRenderer::m_fluidColor = Color(44, 95, 150);
const Color FluidRenderer::m_solidColor = Color(94,94,94);
const Color FluidRenderer::m_velocityVectorColor = Color(255,255,255);

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
    m_gridVertexCount(0)
{
    m_gridVerts.reserve(solver.get()->grid().cellCount()*m_vertexPerCell*m_vertexSize);
    m_gridIndices.reserve(solver.get()->grid().cellCount()*m_vertexPerCell);

    m_vectorVerts.reserve(solver.get()->grid().cellCount()*m_vectorVertsPerCell);

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    m_cellSize = 2.f / std::max(gridWidth,gridHeight);

    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            addGridQuad(Vertex(m_cellSize*j - 1,1 - m_cellSize*i),Vertex(m_cellSize*(j+1) - 1,1 - m_cellSize*(i+1)),Color(255,0,0));
            Vertex vStart = Vertex(m_cellSize*j - 1 + m_cellSize / 2,1 - m_cellSize*i - m_cellSize/2);
            Vertex vEnd = Vertex(vStart.x() + m_cellSize / 2, vStart.y());
            addVector(vStart,vEnd,m_velocityVectorColor);
            setVectorColor(i,j,Color(255,0,0));
        }
    }
}

void FluidRenderer::init()
{
    setupGl();
    updateGrid();
    updateGridVerts();
}

void FluidRenderer::render()
{
    glUseProgram(m_shaderProgram);
    //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glBindVertexArray(m_vao_grid);
    glDrawElements(GL_TRIANGLES,m_gridIndices.size(),GL_UNSIGNED_INT,0);
    glBindVertexArray(m_vao_vectors);
    glDrawArrays(GL_LINES,0,m_vectorVerts.size());
}

void FluidRenderer::addGridVertex(Vertex v, Color c)
{
    m_gridVerts.push_back(v.x());
    m_gridVerts.push_back(v.y());
    m_gridVerts.push_back(v.z());
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
    m_vectorVerts.push_back(v.z());
    m_vectorVerts.push_back(c.rf());
    m_vectorVerts.push_back(c.gf());
    m_vectorVerts.push_back(c.bf());
}

void FluidRenderer::addVector(Vertex start, Vertex end, Color c)
{
    addVectorVertex(start,c);
    addVectorVertex(end,c);
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
    updateBuffers();
}

void FluidRenderer::updateBuffers()
{
    updateGridVerts();
    glBindVertexArray(m_vao_grid);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_grid);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(int) * m_gridIndices.size(), m_gridIndices.data(),GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    updateVectorVerts();
    glBindVertexArray(m_vao_vectors);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(float) * 6, reinterpret_cast<void*>(sizeof(float)*3));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

void FluidRenderer::updateGridVerts()
{
    glBindVertexArray(m_vao_grid);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_grid);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_gridVerts.size(),m_gridVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void FluidRenderer::updateVectorVerts()
{
    glBindVertexArray(m_vao_vectors);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo_vectors);
    glBufferData(GL_ARRAY_BUFFER,sizeof(float) * m_vectorVerts.size(),m_vectorVerts.data(),GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
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
    updateGridVerts();
    updateVectors();
    updateVectorVerts();
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
                    setCellColor(i,j,m_emptyColor);
                break;

                case FluidCellMaterial::FLUID:
                    setCellColor(i,j,m_fluidColor);
                break;

                case FluidCellMaterial::SOLID:
                    setCellColor(i,j,m_solidColor);
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
            setCellColor(i,j,Color(r,g,0.0));
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
            setCellColor(i,j,Color(brightness,brightness,brightness));
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
            setCellColor(i,j,Color(brightness,brightness,brightness));
        }
    }
}

void FluidRenderer::updateVectors()
{
    float maxLengthGridSizes = sqrt(2.f) / 2;
    float maxLength = 50;

    int gridHeight = m_solver->grid().sizeI();
    int gridWidth = m_solver->grid().sizeJ();
    Vertex vEnd;
    for (int i = 0; i < gridHeight; i++)
    {
        for (int j = 0; j < gridWidth; j++)
        {
            Vertex gridspaceVelocity(m_solver->grid().getU(i,j),m_solver->grid().getV(i,j));
            //Vertex gridspaceVelocity(1,0);
            float scaleFactor = maxLength / gridspaceVelocity.distFromZero();
            Vertex newVector = Vertex((gridspaceVelocity.x()/scaleFactor) * maxLengthGridSizes * m_cellSize,
                                       (gridspaceVelocity.y()/scaleFactor) * maxLengthGridSizes * m_cellSize);
            updateVector(i,j,newVector);
        }
    }
}

void FluidRenderer::updateVector(int x, int y, Vertex newVector)
{
    int linearVectorStartVertexIndex = m_solver->grid().linearIndex(x,y) * m_vectorVertsPerCell;
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
    int linearIndex = m_solver->grid().linearIndex(x,y) * m_vertexPerCell;
    setCellVertexColor(linearIndex,c);
    setCellVertexColor(linearIndex+1,c);
    setCellVertexColor(linearIndex+2,c);
    setCellVertexColor(linearIndex+3,c);
}

void FluidRenderer::setVectorVertexColor(int vIndex, Color c)
{
    int linearIndex = vIndex * m_vertexSize + 3;
    m_vectorVerts[linearIndex] = c.rf();
    m_vectorVerts[linearIndex + 1] = c.gf();
    m_vectorVerts[linearIndex + 2] = c.bf();
}

void FluidRenderer::setVectorColor(int x, int y, Color c)
{
    int linearIndex = m_solver->grid().linearIndex(x,y) * m_vectorVertsPerCell;
    setVectorVertexColor(linearIndex, c);
    //setVectorVertexColor(linearIndex + 1, c);
}
