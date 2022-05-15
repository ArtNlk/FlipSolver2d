#ifndef TEXTRENDERER_H
#define TEXTRENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <ft2build.h>
#include FT_FREETYPE_H
#include <string>
#include <map>

#include "color.h"

/// Holds all state information relevant to a character as loaded using FreeType
struct GraphicCharacter {
    unsigned int textureID; // ID handle of the glyph texture
    glm::ivec2   size;      // size of glyph
    glm::ivec2   bearing;   // offset from baseline to left/top of glyph
    unsigned int advance;   // horizontal offset to advance to next glyph
};

class TextRenderer
{
public:
    TextRenderer(int width, int height);

    void init();

    void load(std::string fontPath, int fontSize);

    void renderText(std::string text, float x, float y, float scale, Color color = Color(255,255,255));

    void resize(int width, int height);
protected:
    unsigned int m_vao;
    unsigned int m_vbo;
    unsigned int m_shaderProgram;
    int m_width;
    int m_height;

    std::map<char,GraphicCharacter> m_charmap;
    glm::mat4x4 m_projectionMatrix;
    static const std::string m_fontPath;
    static const char* m_vertexShaderSource;
    static const char* m_fragmentShaderSource;
};

#endif // TEXTRENDERER_H
