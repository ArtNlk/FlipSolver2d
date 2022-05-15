#include "textrenderer.h"

#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#if defined(_WIN32)
    const std::string TextRenderer::m_fontPath = "C:/Windows/Fonts";
#elif defined(__linux__)
    const std::string TextRenderer::m_fontPath = "/usr/share/fonts/truetype";
#endif

const char* TextRenderer::m_vertexShaderSource = \
        "#version 330 core\n\
        layout (location = 0) in vec4 vertex; // <vec2 pos, vec2 texCoords>\n\
        out vec2 TexCoords;\n\
        \n\
        uniform mat4 projection;\n\
        \n\
        void main()\n\
        {\n\
            gl_Position = projection * vec4(vertex.xy, 0.f, 1.f);\n\
            TexCoords = vertex.zw;\n\
        }  ";

const char* TextRenderer::m_fragmentShaderSource = \
        "#version 330 core\n\
        in vec2 TexCoords;\n\
        out vec4 color;\n\
        \n\
        uniform sampler2D text;\n\
        uniform vec3 textColor;\n\
        \n\
        void main()\n\
        {\n\
            vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);\n\
            color = vec4(textColor, 1.0) * sampled;\n\
        }  ";

TextRenderer::TextRenderer(int width, int height)
{
    m_width = width;
    m_height = height;
}

void TextRenderer::init()
{
    // load and configure shader

    int success;
    unsigned int vShader;
    unsigned int fShader;
    vShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vShader,1,&m_vertexShaderSource,NULL);
    glCompileShader(vShader);
    glGetShaderiv(vShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Text renderer: vertex shader compilation error";

        int maxLength = 0;
        glGetShaderiv(vShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::string errorLog(maxLength,'\0');
        glGetShaderInfoLog(vShader, maxLength, &maxLength, &errorLog[0]);

        std::cout << errorLog;
    }

    fShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fShader, 1, &m_fragmentShaderSource, NULL);
    glCompileShader(fShader);
    glGetShaderiv(fShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        std::cout << "Text renderer: fragment shader compilation error";
    }

    m_shaderProgram = glCreateProgram();
    glAttachShader(m_shaderProgram,vShader);
    glAttachShader(m_shaderProgram,fShader);
    glLinkProgram(m_shaderProgram);
    glGetProgramiv(m_shaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        std::cout << "Text renderer: shader program linking error!";
    }

    glDeleteShader(vShader);
    glDeleteShader(fShader);

    // configure VAO/VBO for texture quads
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    m_projectionMatrix = glm::ortho(0.0f, static_cast<float>(m_width),
                                    0.0f, static_cast<float>(m_height));
}

void TextRenderer::load(std::string fontPath, int fontSize)
{
    FT_Library ft;
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
        return;
    }

    FT_Face face;
    if (FT_New_Face(ft, (m_fontPath + fontPath).c_str(), 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
        return;
    }

    FT_Set_Pixel_Sizes(face, 0, fontSize);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // disable byte-alignment restriction

    for (unsigned char c = 0; c < 128; c++)
    {
        // load character glyph
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // generate texture
        unsigned int texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer
        );
        // set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // now store character for later use
        GraphicCharacter character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            static_cast<unsigned int>(face->glyph->advance.x)
        };
        m_charmap.insert(std::pair<char,GraphicCharacter>(static_cast<char>(c),character));
    }
    FT_Done_Face(face);
    FT_Done_FreeType(ft);
}

void TextRenderer::renderText(std::string text, float x, float y, float scale, Color color)
{
    // activate corresponding render state
    glUseProgram(m_shaderProgram);
    glUniform3f(glGetUniformLocation(m_shaderProgram, "textColor"), color.rf(),
                                                                    color.gf(),
                                                                    color.bf());
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "projection"),
                       1,
                       GL_FALSE,
                       glm::value_ptr(m_projectionMatrix));
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBindVertexArray(m_vao);

    // iterate through all characters
        std::string::const_iterator c;
        for (c = text.begin(); c != text.end(); c++)
        {
            GraphicCharacter ch = m_charmap[*c];

            float xpos = x + ch.bearing.x * scale;
            float ypos = y - (ch.size.y - ch.bearing.y) * scale;

            float w = ch.size.x * scale;
            float h = ch.size.y * scale;
            // update VBO for each character
            float vertices[6][4] = {
                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos,     ypos,       0.0f, 1.0f },
                { xpos + w, ypos,       1.0f, 1.0f },

                { xpos,     ypos + h,   0.0f, 0.0f },
                { xpos + w, ypos,       1.0f, 1.0f },
                { xpos + w, ypos + h,   1.0f, 0.0f }
            };
            // render glyph texture over quad
            glBindTexture(GL_TEXTURE_2D, ch.textureID);
            // update content of VBO memory
            glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            // render quad
            glDrawArrays(GL_TRIANGLES, 0, 6);
            // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
            x += (ch.advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64)
        }
        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D, 0);
}

void TextRenderer::resize(int width, int height)
{
    m_width = width;
    m_height = height;
    m_projectionMatrix = glm::ortho(0.0f, static_cast<float>(width),
                                    0.0f, static_cast<float>(height));
}
