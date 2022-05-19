#ifndef GLOBALCALLBACKHANDLER_H
#define GLOBALCALLBACKHANDLER_H

#include <memory>

class LiquidRenderApp;
class FluidRenderer;
class TextMenuRenderer;

class GlobalCallbackHandler
{
public:
    void init(LiquidRenderApp* app, FluidRenderer* liquidRenderer, TextMenuRenderer* textMenuRenderer);

    static GlobalCallbackHandler instance();

    void registerRenderUpdate();

protected:
    static GlobalCallbackHandler m_instance;
    LiquidRenderApp* m_renderapp;
    FluidRenderer* m_fluidRenderer;
    TextMenuRenderer* m_textMenuRenderer;
};

#endif // GLOBALCALLBACKHANDLER_H
