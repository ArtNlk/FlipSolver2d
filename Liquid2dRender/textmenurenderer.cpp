#include "textmenurenderer.h"

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

    m_vectorRenderModeTexts.assign(VectorRenderMode::VECTOR_RENDER_ITER_END + 1,"");
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_CENTER] = "Vector center avg velocity";
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_STAGGERED] = "Vector staggered velocity";
    m_vectorRenderModeTexts[VectorRenderMode::VECTOR_RENDER_ITER_END] = "Vectors off";
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
    glm::vec2 currentTextPos = m_renderModeTextPosition;
    m_textRenderer.renderText(m_renderModeTexts[m_fluidRenderer.gridRenderMode()],
                                currentTextPos.x + widthBaseline,
                                m_height - currentTextPos.y,
                                1.0f,
                                Color(255,255,255));

    currentTextPos += m_vectorRenderModeOffset;
    m_textRenderer.renderText(m_vectorRenderModeTexts[m_fluidRenderer.vectorRenderEnabled() ? m_fluidRenderer.vectorRenderMode() : VectorRenderMode::VECTOR_RENDER_ITER_END],
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
