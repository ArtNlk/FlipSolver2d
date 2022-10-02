#ifndef FLUIDRENDERER_H
#define FLUIDRENDERER_H
#include "flipsolverbase.h"
#include "glad/glad.h"
#include <glm/glm.hpp>

#include <memory>

#include "flipfluidsolver2d.h"
#include "color.h"
#include "geometry2d.h"

enum FluidRenderMode : int {RENDER_MATERIAL,
                            RENDER_VELOCITY,
                            RENDER_U,RENDER_V,
                            RENDER_SDF,
                            RENDER_KNOWN_FLAG_U,
                            RENDER_KNOWN_FLAG_V,
                            RENDER_KNOWN_FLAG_CENTERED,
                            RENDER_VISCOSITY,
                            RENDER_TEMPERATURE,
                            RENDER_SMOKE_CONCENTRATION,
                            RENDER_DIV_CONTROL,
                            GRID_RENDER_ITER_END};
inline FluidRenderMode& operator++(FluidRenderMode& state, int) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<FluidRenderMode>((i) % GRID_RENDER_ITER_END);
    return state;
}

enum VectorRenderMode : int {VECTOR_RENDER_CENTER,VECTOR_RENDER_STAGGERED,VECTOR_RENDER_SDF_GRADIENT,VECTOR_RENDER_ITER_END};
inline VectorRenderMode& operator++(VectorRenderMode& state, int) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<VectorRenderMode>((i) % VECTOR_RENDER_ITER_END);
    return state;
}

enum ParticleRenderMode : int {PARTICLE_RENDER_VELOCITY,PARTICLE_RENDER_SOLID,PARTICLE_RENDER_ITER_END};
inline ParticleRenderMode& operator++(ParticleRenderMode& state, int) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<ParticleRenderMode>((i) % PARTICLE_RENDER_ITER_END);
    return state;
}

class FluidRenderer
{
public:
    FluidRenderer(int textureWidth, int textureHeight);
    void init(std::shared_ptr<FlipSolverBase> solver);
    void render();
    void update();
    void updateGrid();
    void updateVectors();
    void updateParticles();

    inline void setRenderMode(FluidRenderMode mode)
    {
        m_gridRenderMode = mode;
    }

    inline FluidRenderMode& gridRenderMode()
    {
        return m_gridRenderMode;
    }

    inline VectorRenderMode& vectorRenderMode()
    {
        return m_vectorRenderMode;
    }

    inline ParticleRenderMode& particleRenderMode()
    {
        return m_particleRenderMode;
    }

    inline bool& vectorRenderEnabled()
    {
        return m_vectorsEnabled;
    }

    inline bool toggleVectors()
    {
        m_vectorsEnabled = !m_vectorsEnabled;
        return m_vectorsEnabled;
    }

    inline bool& geometryEnabled()
    {
        return m_geometryEnabled;
    }

    inline bool toggleGeometry()
    {
        m_geometryEnabled = !m_geometryEnabled;
        return m_geometryEnabled;
    }

    inline bool& particlesEnabled()
    {
        return m_particlesEnabled;
    }

    inline bool toggleParticles()
    {
        m_particlesEnabled = !m_particlesEnabled;
        return m_particlesEnabled;
    }

    static const Color &velocityVectorColor();

    unsigned int renderTexture();

    void dumpToPng(std::string fileName);

    void resizeTexture(int width, int height);

    float fluidGridAspect();

    std::shared_ptr<FlipSolverBase> solver();

protected:
    void addGridVertex(Vertex v, Color c = Color());
    void addGridQuad(Vertex topLeft, Vertex bottomRight, Color c);
    void addVectorVertex(Vertex v, Color c = Color());
    void addVector(Vertex start, Vertex end, Color c);
    void addParticle(Vertex particle, Color c);
    void setupGl();
    void setupOffscreenBuffer();
    void setupBuffers();
    void setupGridVerts();
    void setupVectorVerts();
    void setupGeometryVerts();
    void setupParticleVerts();
    void loadGeometry();
    void addGeometry(Geometry2d &geometry);
    void updateGridVerts();
    void updateVectorVerts();
    void updateGeometryVerts();
    void updateParticleVerts();
    void updateGridFromMaterial();
    void updateGridFromVelocity();
    void updateGridFromUComponent();
    void updateGridFromVComponent();
    void updateGridFromUKnownFlag();
    void updateGridFromVKnownFlag();
    void updateGridFromCenteredKnownFlag();
    void updateGridFromSdf();
    void updateGridFromViscosity();
    void updateGridFromTemperature();
    void updateGridFromConcentration();
    void updateGridFromDivergence();
    void updateVectorsStaggered();
    void updateVectorsCentered();
    void updateVectorsSdfGrad();
    void reloadParticlesSolid();
    void reloadParticlesVelocity();
    void updateVector(int x, int y, Vertex newVector);
    void setCellVertexColor(int vIndex, Color c);
    void setCellColor(int x, int y, Color c);
    void setVectorVertexColor(int vIndex, Color c);
    void setVectorColor(int x, int y, Color c);
    void setParticleColor(int idx, Color c);

    Color hueColorRamp(float val);
    Color velocityComponentColorRamp(float val);

    unsigned int m_vbo_grid;
    unsigned int m_vao_grid;
    unsigned int m_ebo_grid;
    unsigned int m_vbo_vectors;
    unsigned int m_vao_vectors;
    unsigned int m_ebo_vectors;
    unsigned int m_vbo_geometry;
    unsigned int m_vao_geometry;
    unsigned int m_ebo_geometry;
    unsigned int m_vbo_particles;
    unsigned int m_vao_particles;
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
    static const Color m_sourceColor;
    static const Color m_sinkColor;
    static const Color m_velocityVectorColor;
    static const Color m_markerParticleColor;
    static const Color m_geometryColor;

    static constexpr float m_velocityRangeMax = 50;
    static constexpr float m_velocityComponentRangeMax = 5;
    static constexpr float m_hueMinVelocity = 0.66f;
    static constexpr float m_hueMaxVelocity = 0.f;
    static constexpr float m_hueMinVelocityComponent = 0.66f;
    static constexpr float m_hueMaxVelocityComponent = 0.f;

    std::vector<float> m_gridVerts;
    std::vector<unsigned int> m_gridIndices;
    std::vector<float> m_vectorVerts;
    std::vector<float> m_geometryVerts;
    std::vector<float> m_particleVerts;
    std::vector<unsigned int> m_geometryIndices;
    glm::mat4 m_projection;
    FluidRenderMode m_gridRenderMode;
    VectorRenderMode m_vectorRenderMode;
    ParticleRenderMode m_particleRenderMode;
    bool m_vectorsEnabled;
    bool m_geometryEnabled;
    bool m_particlesEnabled;

    std::shared_ptr<FlipSolverBase> m_solver;
};

#endif // FLUIDRENDERER_H
