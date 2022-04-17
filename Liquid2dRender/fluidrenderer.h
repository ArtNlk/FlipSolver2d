#ifndef FLUIDRENDERER_H
#define FLUIDRENDERER_H
#include "glad/glad.h"

#include <memory>

#include "flipsolver2d.h"
#include "color.h"
#include "vertex.h"

enum FluidRenderMode : int {RENDER_MATERIAL,RENDER_VELOCITY,RENDER_U,RENDER_V,RENDER_ITER_END};
inline FluidRenderMode& operator++(FluidRenderMode& state, int) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<FluidRenderMode>((i) % RENDER_ITER_END);
    return state;
}

class FluidRenderer
{
public:
    FluidRenderer(std::shared_ptr<FlipSolver> solver);
    void init();
    void render();
    void updateGrid();

    inline void setRenderMode(FluidRenderMode mode)
    {
        m_renderMode = mode;
    }

    inline FluidRenderMode getRenderMode()
    {
        return m_renderMode;
    }

    inline FluidRenderMode& renderMode()
    {
        return m_renderMode;
    }

protected:
    void addVertex(Vertex v, Color c = Color());
    void addQuad(Vertex topLeft, Vertex bottomRight, Color C);
    void setupGl();
    void updateBuffers();
    void updateVerts();
    void updateGridFromMaterial();
    void updateGridFromVelocity();
    void updateGridFromUComponent();
    void updateGridFromVComponent();
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

    static const Color m_emptyColor;
    static const Color m_fluidColor;
    static const Color m_solidColor;

    std::vector<float> m_verts;
    std::vector<unsigned int> m_indices;
    FluidRenderMode m_renderMode;

    std::shared_ptr<FlipSolver> m_solver;
};

#endif // FLUIDRENDERER_H
