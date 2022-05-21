#ifndef TEXTMENURENDERER_H
#define TEXTMENURENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "textrenderer.h"
#include "fluidrenderer.h"

class TextMenuRenderer
{
public:
    TextMenuRenderer(int x, int y, int width, int height, FluidRenderer &fluidRenderer);

    void init();

    void render();

    void move(int x, int y);

    void resize(int width, int height);

    float menuWidth();

protected:
    FluidRenderer &m_fluidRenderer;
    TextRenderer m_textRenderer;
    std::vector<std::string> m_renderModeTexts;
    std::vector<std::string> m_vectorRenderModeTexts;

    int m_x;
    int m_y;
    int m_width;
    int m_height;

    glm::vec2 m_renderModeTextPosition = glm::vec2(0.f,20.f);
    glm::vec2 m_vectorRenderModeOffset = glm::vec2(0.f,20.f);
    glm::vec2 m_geometryRenderModeOffset = glm::vec2(0.f,20.f);

    float m_menuWidth = 300.f;
};

#endif // TEXTMENURENDERER_H
