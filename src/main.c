#include <stdio.h>
#include <uxhw.h>
#include "header.h"

/*
 *	This is a generic ode integrator function pointer, we can choose the specific integrator inside or just implement one
 */
sirModelState
odeIntegrationStep(sirModelState *  currentState, simulationParameters *  simParameters)
{
	sirModelState	outputState = {0};

	double	S0 = currentState->susceptible;
	double	I0 = currentState->infected;
	double	R0 = currentState->recovered;
	double	susceptibleToInfectedRate = currentState->interactionParameters.susceptibleToInfectedRate;
	double	infectedToRecoveredRate = currentState->interactionParameters.infectedToRecoveredRate;

	double	deltaT = simParameters->integratorTimeStep;

	outputState.susceptible = S0 + deltaT * (-susceptibleToInfectedRate * S0 * I0);
	outputState.infected = I0 + deltaT * (susceptibleToInfectedRate * S0 * I0 - infectedToRecoveredRate * I0);
	outputState.recovered = R0 + deltaT * (infectedToRecoveredRate * I0);

	return outputState;
};

int
main(int argc, char *  argv[])
{
	/*
	 *	Setup the interaction parameters
	 */
	vectorFieldParameters	sirModelInteractionParameters = {
		.susceptibleToInfectedRate = UxHwDoubleUniformDist(0.25, 0.35),
		.infectedToRecoveredRate = UxHwDoubleUniformDist(0.15, 0.25)
	};

	/*
	 *	Setup initial state and parameters
	 */
	sirModelState	sirState = {
		.susceptible = UxHwDoubleGaussDist(0.75, 0.01),
		.infected = UxHwDoubleGaussDist(0.15, 0.01),
		.recovered = UxHwDoubleGaussDist(0.1, 0.01),
		.interactionParameters = sirModelInteractionParameters
	};
	simulationParameters	simTimeParams = {
		.initTime = kSirModelInitialTime,
		.finalTime = kSirModelFinalTime,
		.integratorTimeStep = kSirModelTimestepSize
	};

	/*
	 *	Define an array of states to store all the information. Same for the timing info
	 */
	sirModelState	sirModelEvolution[(kSirModelFinalTime - kSirModelInitialTime) / kSirModelTimestepSize + 1];
	double		timeInstantArray[(kSirModelFinalTime - kSirModelInitialTime) / kSirModelTimestepSize + 1];

	sirModelEvolution[0] = sirState;
	timeInstantArray[0] = simulationParameters.initTime;
	for(size_t simIteration = 0; simIteration < (kSirModelFinalTime - kSirModelInitialTime) / kSirModelTimestepSize; simIteration++)
	{
		sirModelEvolution[simIteration + 1] = odeIntegrationStep(&sirModelEvolution[simIteration], &simTimeParams);
		timeInstantArray[simIteration + 1] = timeInstantArray[simIteration] + simTimeParams.integratorTimeStep;
	}

	sirModelState	finalState = sirModelEvolution[(kSirModelFinalTime - kSirModelInitialTime) / kSirModelTimestepSize + 1];
	fprintf(stdout, "Final susceptibles: %lf\nFinal infected: %lf\nFinal recovered: %lf\n", finalState.susceptible, finalState.infected, finalState.recovered);

	/*
	 *	Store the full simulation to a file in the cloud storage:
	 */
	FILE *	outputFile = fopen("sirModelOutput.txt", "w");
	if(!outputFile)
	{
		fprintf(stderr, "Error has occured. File cannot be opened.\n");
		return EXIT_FAILURE;
	}
	for(size_t simIteration = 0; simIteration < (kSirModelFinalTime - kSirModelInitialTime) / kSirModelTimestepSize + 1; simIteration++)
	{
		fprintf(outputFile, "%zu, %lf, %lf\n", simIteration, timeInstantArray[simIteration], sirModelEvolution[simIteration]);
	}
	fclose(outputFile);

	return EXIT_SUCCESS;
}
