#include "globalcallbackhandler.h"

#include "liquidrenderapp.h"

GlobalCallbackHandler GlobalCallbackHandler::m_instance = GlobalCallbackHandler();

void GlobalCallbackHandler::init(LiquidRenderApp* app, FluidRenderer* liquidRenderer, TextMenuRenderer* textMenuRenderer)
{
    m_instance.m_renderapp = app;
    m_instance.m_fluidRenderer = liquidRenderer;
    m_instance.m_textMenuRenderer = textMenuRenderer;
}

GlobalCallbackHandler GlobalCallbackHandler::instance()
{
    return m_instance;
}

void GlobalCallbackHandler::registerRenderUpdate()
{
    m_renderapp->requestRender();
}
