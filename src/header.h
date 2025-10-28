#pragma once

#define kSirModelInitialTime	(0.0)
#define kSirModelFinalTime	(30.0)
#define kSirModelTimestepSize	(0.5)

typedef struct
{
	double	susceptibleToInfectedRate;
	double	infectedToRecoveredRate;
    double  birthRate;
    double  deathRate;
} vectorFieldParameters;

typedef struct
{
	double			susceptible;
	double			infected;
	double			recovered;
	vectorFieldParameters	interactionParameters;
} sirModelState;

typedef struct
{
	double	initTime;
	double	finalTime;
	double	integratorTimeStep;
    size_t  numberOfSteps;
} simulationParameters;

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
