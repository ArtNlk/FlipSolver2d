#ifndef JSONSCENEREADER_H
#define JSONSCENEREADER_H

#include <memory>

#include "nlohmann/json.hpp"

#include "solvers.h"

class JsonSceneReader
{
    using json = nlohmann::json;
public:
    JsonSceneReader() = default;

    static std::shared_ptr<FlipSolver> loadJson(std::string fileName);

protected:
    static void populateFlipSolverParamsFromJson(FlipSolverParameters* p, json settingsJson);
    static void populateNBFlipSolverParamsFromJson(NBFlipParameters* p, json settingsJson);
    static void populateSmokeSolverParamsFromJson(SmokeSolverParameters* p, json settingsJson);
    static void populateFireSolverParamsFromJson(FireSolverParameters* p, json settingsJson);

    static SimulationMethod simMethodFromName(const std::string& name);

    template<class T>
    static T tryGetValue(json input, std::string key, T defaultValue);

    static void objectsFromJson(json solverJson, std::shared_ptr<FlipSolver> solver);
    static Emitter emitterFromJson(json emitterJson, float sceneScale);
    static Obstacle obstacleFromJson(json obstacleJson, float sceneScale);
    static Sink sinkFromJson(json sinkJson, float sceneScale);
    static void addObjectFromJson(json objectJson, std::shared_ptr<FlipSolver> solver);
};

#endif // JSONSCENEREADER_H
