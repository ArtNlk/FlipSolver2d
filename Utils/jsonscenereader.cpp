#include "jsonscenereader.h"

#include <memory>
#include <fstream>

#include "logger.h"

std::shared_ptr<FlipSolver> JsonSceneReader::loadJson(std::string fileName)
{
    std::shared_ptr<FlipSolver> output;

    try
    {
        std::ifstream sceneFile(fileName);
        if(!sceneFile.is_open())
        {
            std::cout << "errorOpening scene file " << fileName;
        }
        json sceneJson;
        sceneFile >> sceneJson;
        json settingsJson = sceneJson["settings"];

        std::string simTypeName = settingsJson["simType"].get<std::string>();
        SimulationMethod simMethod = simMethodFromName(simTypeName);


        switch(simMethod)
        {
        case SIMULATION_LIQUID:
        {
            FlipSolverParameters p;
            populateFlipSolverParamsFromJson(&p, settingsJson);
            output.reset(new FlipSolver(&p));
        }
        break;

        case SIMULATION_SMOKE:
        {
            SmokeSolverParameters p;
            populateFlipSolverParamsFromJson(&p, settingsJson);
            populateSmokeSolverParamsFromJson(&p, settingsJson);
            output.reset(new FlipSmokeSolver(&p));
        }
        break;

        case SIMULATION_FIRE:
        {
            FireSolverParameters p;
            populateFlipSolverParamsFromJson(&p, settingsJson);
            populateFireSolverParamsFromJson(&p, settingsJson);
            output.reset(new FlipFireSolver(&p));
        }
        break;

        case SIMULATION_NBFLIP:
        {
            NBFlipParameters p;
            populateFlipSolverParamsFromJson(&p, settingsJson);
            populateNBFlipSolverParamsFromJson(&p, settingsJson);
            output.reset(new NBFlipSolver(&p));
        }
        break;
        }

        output->initAdditionalParameters();
        objectsFromJson(sceneJson["solver"], output);
    }
    catch (std::exception &e)
    {
        debug() << e.what();
        std::cout << e.what();

        return std::shared_ptr<FlipSolver>();
    }

    return output;
}

void JsonSceneReader::populateFlipSolverParamsFromJson(FlipSolverParameters *p, json settingsJson)
{
    float s = tryGetValue(settingsJson,"scale",1.f);
    p->fluidDensity = tryGetValue(settingsJson,"density",1.f);
    p->seed = tryGetValue(settingsJson,"seed",0);
    p->dx =
        p->particlesPerCell = settingsJson["particlesPerCell"].get<int>();
    std::pair<float,float> v = tryGetValue(settingsJson,"globalAcceleration",std::pair(9.8,0.f));
    // v.first *= s;
    // v.second *= s;
    p->globalAcceleration = v;
    p->resolution = settingsJson["resolution"].get<int>();
    p->fps = settingsJson["fps"].get<int>();
    p->maxSubsteps = tryGetValue(settingsJson,"maxSubsteps",30);
    p->picRatio = tryGetValue(settingsJson,"picRatio",0.03);
    p->cflNumber = tryGetValue(settingsJson,"cflNumber",10.f);
    p->particleScale = tryGetValue(settingsJson,"particleScale",0.8);
    p->pcgIterLimit = tryGetValue(settingsJson,"pcgIterLimit",200);
    p->viscosityEnabled = tryGetValue(settingsJson,"viscosityEnabled",false);
    p->domainSizeI = settingsJson["domainSizeI"].get<float>() * s;
    p->domainSizeJ = settingsJson["domainSizeJ"].get<float>() * s;
    p->useHeavyViscosity = tryGetValue(settingsJson, "heavyViscosity", false);
    p->sceneScale = s;

    if(p->domainSizeI > p->domainSizeJ)
    {
        p->dx = static_cast<float>(p->domainSizeI) /
                p->resolution;
        p->gridSizeI = p->resolution;
        p->gridSizeJ = (static_cast<float>(p->domainSizeJ) /
                        static_cast<float>(p->domainSizeI)) * p->resolution;
    }
    else
    {
        p->dx = static_cast<float>(p->domainSizeJ) /
                p->resolution;
        p->gridSizeJ = p->resolution;
        p->gridSizeI = (static_cast<float>(p->domainSizeI) /
                        static_cast<float>(p->domainSizeJ)) * p->resolution;
    }

    std::string parameterHandling = tryGetValue(settingsJson,"parameterHandlingMethod",std::string("particle"));
    if(parameterHandling == "particle")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::PARTICLE;
    }
    if(parameterHandling == "hybrid")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::HYBRID;
    }
    if(parameterHandling == "grid")
    {
        p->parameterHandlingMethod = ParameterHandlingMethod::GRID;
    }

    std::string simTypeName = settingsJson["simType"].get<std::string>();
    p->simulationMethod = simMethodFromName(simTypeName);
}

void JsonSceneReader::populateNBFlipSolverParamsFromJson(NBFlipParameters *p, json settingsJson)
{
    return;
}

void JsonSceneReader::populateSmokeSolverParamsFromJson(SmokeSolverParameters *p, json settingsJson)
{
    p->ambientTemperature = tryGetValue(settingsJson,"ambientTemperature",273.0f);
    p->temperatureDecayRate = tryGetValue(settingsJson,"temperatureDecayRate",0.0);
    p->concentrationDecayRate = tryGetValue(settingsJson,"concentrationDecayRate",0.0);
    p->buoyancyFactor = tryGetValue(settingsJson,"buoyancyFactor",1.0);
    p->sootFactor = tryGetValue(settingsJson,"sootFactor",1.0);
}

void JsonSceneReader::populateFireSolverParamsFromJson(FireSolverParameters *p, json settingsJson)
{
    populateSmokeSolverParamsFromJson(p,settingsJson);
    p->ignitionTemperature = tryGetValue(settingsJson,"ignitionTemp",250.f);
    p->burnRate = tryGetValue(settingsJson,"burnRate",0.05f);
    p->smokeProportion = tryGetValue(settingsJson,"smokeEmission",1.f);
    p->heatProportion = tryGetValue(settingsJson,"heatEmission",1.f);
    p->divergenceProportion = tryGetValue(settingsJson,"billowing",0.1f);
}

SimulationMethod JsonSceneReader::simMethodFromName(const std::string &name)
{
    if(name == "fluid" || name == "flip")
    {
        return SimulationMethod::SIMULATION_LIQUID;
    }
    if(name == "smoke")
    {
        return SimulationMethod::SIMULATION_SMOKE;
    }
    if(name == "fire")
    {
        return SimulationMethod::SIMULATION_FIRE;
    }
    if(name == "nbflip")
    {
        return SimulationMethod::SIMULATION_NBFLIP;
    }

    return SimulationMethod::SIMULATION_LIQUID;
}

//    SimSettings::airDensity() = tryGetValue(settingsJson,"airDensity",SimSettings::fluidDensity() * 0.001f);
//    SimSettings::surfaceTensionFactor() = tryGetValue(settingsJson,"surfaceTensionFactor",0.0);

void JsonSceneReader::objectsFromJson(json solverJson, std::shared_ptr<FlipSolver> solver)
{
    std::vector<json> objects = solverJson["objects"]
                                    .get<std::vector<json>>();
    for(json &geo : objects)
    {
        addObjectFromJson(geo, solver);
    }
}

Emitter JsonSceneReader::emitterFromJson(json emitterJson, float sceneScale)
{
    std::vector<std::pair<float,float>> verts = emitterJson["verts"]
                                                     .get<std::vector<std::pair<float,float>>>();
    float temp = tryGetValue(emitterJson,"temperature",273.f);
    float conc = tryGetValue(emitterJson,"concentration",1.f);
    float viscosity = tryGetValue(emitterJson,"viscosity",0.f);
    float fuel = tryGetValue(emitterJson,"fuel",1.f);
    float div = tryGetValue(emitterJson,"divergence",0.f);
    bool transfersVelocity = tryGetValue(emitterJson,"transferVelocity",false);
    std::pair<float,float> velocity = tryGetValue(emitterJson,"velocity",std::pair<float,float>(0.f, 0.f));
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(sceneScale * Vec3(v.first,v.second));
    }

    Emitter output(geo);

    output.setTemperature(temp);
    output.setConcentrartion(conc);
    output.setViscosity(viscosity);
    output.setFuel(fuel);
    output.setDivergence(div);
    output.setVelocity(velocity);
    output.setVelocityTransfer(transfersVelocity);

    return output;
}

Obstacle JsonSceneReader::obstacleFromJson(json obstacleJson, float sceneScale)
{
    float friction = tryGetValue(obstacleJson,"friction",0);
    std::vector<std::pair<float,float>> verts = obstacleJson["verts"]
                                                     .get<std::vector<std::pair<float,float>>>();
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(sceneScale * Vec3(v.first,v.second));
    }

    return Obstacle(friction,geo);
}

Sink JsonSceneReader::sinkFromJson(json sinkJson, float sceneScale)
{
    std::vector<std::pair<float,float>> verts = sinkJson["verts"]
                                                     .get<std::vector<std::pair<float,float>>>();
    float div = tryGetValue(sinkJson,"divergence",0.f);
    Geometry2d geo;
    for(auto v : verts)
    {
        geo.addVertex(sceneScale * Vec3(v.first,v.second));
    }

    return Sink(div,geo);
}

void JsonSceneReader::addObjectFromJson(json objectJson, std::shared_ptr<FlipSolver> solver)
{
    std::string geoType = objectJson["type"].get<std::string>();

    const float scale = solver->sceneScale();

    bool enabled = tryGetValue(objectJson,"enabled",true);
    if(!enabled)
    {
        return;
    }

    if(geoType == "solid")
    {
        Obstacle o = obstacleFromJson(objectJson,scale);
        solver->addGeometry(o);
    }
    else if(geoType == "source")
    {
        Emitter e = emitterFromJson(objectJson,scale);
        solver->addSource(e);
    }
    else if(geoType == "sink")
    {
        Sink s = sinkFromJson(objectJson,scale);
        solver->addSink(s);
    }
    else if(geoType == "fluid")
    {
        Emitter e = emitterFromJson(objectJson,scale);
        solver->addInitialFluid(e);
    }
}

template<class T>
T JsonSceneReader::tryGetValue(json input, std::string key, T defaultValue)
{
    return input.contains(key) ? input[key].get<T>() : defaultValue;
}
