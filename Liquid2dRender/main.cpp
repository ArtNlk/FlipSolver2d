#include <glad/glad.h>
#include <GLFW/glfw3.h>


#include <iostream>

#include "liquidrenderapp.h"

int main()
{
    LiquidRenderApp app;
    app.init();
    app.run();
    Logger::instance().close();
}
