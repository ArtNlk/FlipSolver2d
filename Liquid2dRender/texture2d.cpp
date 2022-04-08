#include "texture2d.h"

#include <iostream>

//Texture2d::Texture2d(int width, int height, int format) :
//    m_textureFormat(format),
//    m_width(width),
//    m_height(height),
//    m_pixelSize(0)
//{
//    switch(format)
//    {
//        case GL_RGBA:
//            m_pixelSize = 4;
//            break;

//        case GL_RGB:
//            m_pixelSize = 3;
//            break;

//        default:
//            std::cout << "Bad texture allocation format, defaulting to RGBA: " << format;
//            m_pixelSize = 4;
//            m_textureFormat = GL_RGBA;
//    }

//    m_data = new unsigned char[width*height];
//}
