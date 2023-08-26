#include <fenv.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#ifdef __linux__
#include <cfenv>
#endif

#include <iostream>

#include "liquidrenderapp.h"
#include "logger.h"

#ifdef __linux__
#include <execinfo.h>
#include <unistd.h>
void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}
#endif

int main()
{
#ifdef __linux__
    signal(SIGSEGV, handler);
    feenableexcept(FE_DIVBYZERO);
#endif
    LiquidRenderApp app;
    app.init();
    app.run();
    Logger::instance().close();
}
