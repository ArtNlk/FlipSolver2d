#include "textmenurenderer.h"

TextMenuRenderer::TextMenuRenderer(int x, int y, int width, int height, FluidRenderer &fluidRenderer) :
    m_fluidRenderer(fluidRenderer),
    m_x(x),
    m_y(y),
    m_width(width),
    m_height(height),
    m_textRenderer(width,height)
{
    m_renderModeNames.assign(FluidRenderMode::RENDER_ITER_END,"");
    m_renderModeNames[FluidRenderMode::RENDER_MATERIAL] = "Render material";
    m_renderModeNames[FluidRenderMode::RENDER_VELOCITY] = "Render velocity";
    m_renderModeNames[FluidRenderMode::RENDER_U] = "Render U component";
    m_renderModeNames[FluidRenderMode::RENDER_V] = "Render V component";
}

void TextMenuRenderer::init()
{
    m_textRenderer.load("/open-sans/OpenSans-Regular.ttf",16);
    m_textRenderer.init();
}

void TextMenuRenderer::render()
{
    float widthBaseline = m_width - m_menuWidth;
    m_textRenderer.renderText(m_renderModeNames[m_fluidRenderer.renderMode()],
                                m_renderModeTextPosition.x + widthBaseline,
                                m_height - m_renderModeTextPosition.y,
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
