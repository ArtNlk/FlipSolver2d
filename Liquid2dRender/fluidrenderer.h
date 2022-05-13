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
    FluidRenderer(std::shared_ptr<FlipSolver> solver, int textureWidth, int textureHeight);
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

    static const Color &velocityVectorColor();

    unsigned int renderTexture();

    void resizeTexture(int width, int height);

    float fluidGridAspect();

protected:
    void addGridVertex(Vertex v, Color c = Color());
    void addGridQuad(Vertex topLeft, Vertex bottomRight, Color c);
    void addVectorVertex(Vertex v, Color c = Color());
    void addVector(Vertex start, Vertex end, Color c);
    void setupGl();
    void setupOffscreenBuffer();
    void updateBuffers();
    void updateGridVerts();
    void updateVectorVerts();
    void updateGridFromMaterial();
    void updateGridFromVelocity();
    void updateGridFromUComponent();
    void updateGridFromVComponent();
    void updateVectors();
    void updateVector(int x, int y, Vertex newVector);
    void setCellVertexColor(int vIndex, Color c);
    void setCellColor(int x, int y, Color c);
    void setVectorVertexColor(int vIndex, Color c);
    void setVectorColor(int x, int y, Color c);

    unsigned int m_vbo_grid;
    unsigned int m_vao_grid;
    unsigned int m_ebo_grid;
    unsigned int m_vbo_vectors;
    unsigned int m_vao_vectors;
    unsigned int m_ebo_vectors;
    unsigned int m_vertexShader;
    unsigned int m_fragShader;
    unsigned int m_shaderProgram;
    unsigned int m_framebuffer;
    unsigned int m_renderbuffer;
    unsigned int m_renderTexture;

    int m_gridVertexCount;
    int m_textureWidth;
    int m_textureHeight;

    static const int m_vertexPerCell = 4;
    static const int m_vertexSize = 6;

    static const int m_vectorVertsPerCell = 2;

    static const char *m_vertexShaderSource;
    static const char *m_fragShaderSource;

    static const Color m_emptyColor;
    static const Color m_fluidColor;
    static const Color m_solidColor;
    static const Color m_velocityVectorColor;

    float m_cellWidth;
    float m_cellHeight;

    std::vector<float> m_gridVerts;
    std::vector<unsigned int> m_gridIndices;
    std::vector<float> m_vectorVerts;
    FluidRenderMode m_renderMode;

    std::shared_ptr<FlipSolver> m_solver;
};

#endif // FLUIDRENDERER_H
