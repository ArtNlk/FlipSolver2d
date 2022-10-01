#include "textmenurenderer.h"
#include "simsettings.h"

TextMenuRenderer::TextMenuRenderer(int x, int y, int width, int height, FluidRenderer &fluidRenderer) :
    m_fluidRenderer(fluidRenderer),
    m_x(x),
    m_y(y),
    m_width(width),
    m_height(height),
    m_textRenderer(width,height)
{
    m_renderModeTexts.assign(FluidRenderMode::GRID_RENDER_ITER_END,"");
    m_renderModeTexts[FluidRenderMode::RENDER_MATERIAL] = "Render material";
    m_renderModeTexts[FluidRenderMode::RENDER_VELOCITY] = "Render velocity";
    m_renderModeTexts[FluidRenderMode::RENDER_U] = "Render U component";
    m_renderModeTexts[FluidRenderMode::RENDER_V] = "Render V component";
    m_renderModeTexts[FluidRenderMode::RENDER_SDF] = "Render SDF";
    m_renderModeTexts[FluidRenderMode::RENDER_KNOWN_FLAG_U] = "Render U flag";
    m_renderModeTexts[FluidRenderMode::RENDER_KNOWN_FLAG_V] = "Render V flag";
    m_renderModeTexts[FluidRenderMode::RENDER_KNOWN_FLAG_CENTERED] = "Render centered flag";
    m_renderModeTexts[FluidRenderMode::RENDER_VISCOSITY] = "Render viscosity";
    m_renderModeTexts[FluidRenderMode::RENDER_TEMPERATURE] = "Render temperature";
    m_renderModeTexts[FluidRenderMode::RENDER_SMOKE_CONCENTRATION] = "Render smoke concentrartion";

    m_vectorRenderModeTexts.assign(VectorRenderMode::VECTOR_RENDER_ITER_END + 1,"");
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_CENTER] = "Vector center avg velocity";
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_STAGGERED] = "Vector staggered velocity";
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_SDF_GRADIENT] = "Vector SDF gradient";
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_ITER_END] = "Vectors off";

    m_particleRenderModeTexts.assign(ParticleRenderMode::PARTICLE_RENDER_ITER_END + 1,"");
    m_particleRenderModeTexts[ParticleRenderMode::PARTICLE_RENDER_VELOCITY] = "Particle velocity";
    m_particleRenderModeTexts[ParticleRenderMode::PARTICLE_RENDER_SOLID] = "Particle render solid";
    m_particleRenderModeTexts[ParticleRenderMode::PARTICLE_RENDER_ITER_END] = "Particles off";
}

void TextMenuRenderer::init()
{
#if defined(_WIN32)
    const char* font = "/Consola.ttf";
#elif defined(__linux__)
    const char* font = "/open-sans/OpenSans-Regular.ttf";
#endif
    m_textRenderer.load(font,16);
    m_textRenderer.init();
}

void TextMenuRenderer::render()
{
    float widthBaseline = m_width - m_menuWidth;
    glm::vec2 currentTextPos = m_nextLineOffset;
    std::string temp;

    temp = "Frame: " + std::to_string(m_fluidRenderer.solver()->frameNumber());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Domain size (IxJ): " + std::to_string(SimSettings::domainSizeI()) + " x " + std::to_string(SimSettings::domainSizeJ());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Grid size (IxJ): " + std::to_string(SimSettings::gridSizeI()) + " x " + std::to_string(SimSettings::gridSizeJ());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Particles per cell: " + std::to_string(SimSettings::particlesPerCell());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "PIC ratio: " + std::to_string(SimSettings::picRatio());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "CFL: " + std::to_string(SimSettings::cflNumber());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Max substeps: " + std::to_string(SimSettings::maxSubsteps());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "dX: " + std::to_string(SimSettings::dx());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Step dT: " + std::to_string(SimSettings::stepDt());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Substeps: " + std::to_string(SimSettings::maxSubsteps());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Target FPS: " + std::to_string(SimSettings::fps());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "Global accel: ";
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "          I: " + std::to_string(SimSettings::globalAcceleration().x());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    temp = "          J: " + std::to_string(SimSettings::globalAcceleration().y());
    m_textRenderer.renderText(temp,
                              currentTextPos.x + widthBaseline,
                              m_height - currentTextPos.y,
                              1.0f,
                              Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    m_textRenderer.renderText(m_renderModeTexts[m_fluidRenderer.gridRenderMode()],
                                currentTextPos.x + widthBaseline,
                                m_height - currentTextPos.y,
                                1.0f,
                                Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    m_textRenderer.renderText(m_vectorRenderModeTexts[m_fluidRenderer.vectorRenderEnabled() ? m_fluidRenderer.vectorRenderMode() : VectorRenderMode::VECTOR_RENDER_ITER_END],
                                currentTextPos.x + widthBaseline,
                                m_height - currentTextPos.y,
                                1.0f,
                                Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    m_textRenderer.renderText(m_fluidRenderer.geometryEnabled() ? "Geometry on" : "Geometry off",
                                currentTextPos.x + widthBaseline,
                                m_height - currentTextPos.y,
                                1.0f,
                                Color(255,255,255));

    currentTextPos += m_nextLineOffset;
    m_textRenderer.renderText(m_particleRenderModeTexts[m_fluidRenderer.particlesEnabled() ? m_fluidRenderer.particleRenderMode() : ParticleRenderMode::PARTICLE_RENDER_ITER_END],
                                currentTextPos.x + widthBaseline,
                                m_height - currentTextPos.y,
                                1.0f,
                                Color(255,255,255));
}

void TextMenuRenderer::move(int x, int y)
{
    m_x = x;
    m_y = y;
}

void TextMenuRenderer::resize(int width, int height)
{
    m_width= width;
    m_height = height;
    m_textRenderer.resize(width,height);
}

float TextMenuRenderer::menuWidth()
{
    return m_menuWidth;
}
