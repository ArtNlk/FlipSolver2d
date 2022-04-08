#ifndef FLUIDRENDERER_H
#define FLUIDRENDERER_H
#include "glad/glad.h"

#include <memory>

#include "flipsolver2d.h"
#include "color.h"
#include "vertex.h"

class FluidRenderer
{
public:
    FluidRenderer(std::shared_ptr<FlipSolver> solver);
    void init();
    void render();

protected:
    void addVertex(Vertex v, Color c = Color());
    void addQuad(Vertex topLeft, Vertex bottomRight, Color C);
    void setupGl();
    void updateBuffers();
    void setVertexColor(int vIndex, Color c);
    void setColor(int x, int y, Color c);

    unsigned int m_vbo;
    unsigned int m_vao;
    unsigned int m_ebo;
    unsigned int m_vertexShader;
    unsigned int m_fragShader;
    unsigned int m_shaderProgram;

    int m_vertexCount;

    static const int m_vertexPerCell = 4;
    static const int m_vertexSize = 6;

    static const char *m_vertexShaderSource;
    static const char *m_fragShaderSource;

    std::vector<float> m_verts;
    std::vector<unsigned int> m_indices;

    std::shared_ptr<FlipSolver> m_solver;
};

#endif // FLUIDRENDERER_H
