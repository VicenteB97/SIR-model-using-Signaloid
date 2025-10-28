#include <stdio.h>
#include <stdlib.h>
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
	double	deathRate = currentState->interactionParameters.deathRate;
	double	birthRate = currentState->interactionParameters.birthRate;

	double	deltaT = simParameters->integratorTimeStep;

	outputState.susceptible = S0 + deltaT * (birthRate - deathRate*S0 - susceptibleToInfectedRate*S0*I0);
	outputState.infected = I0 + deltaT * (susceptibleToInfectedRate*S0*I0 - (infectedToRecoveredRate + deathRate)*I0);
	outputState.recovered = R0 + deltaT * (infectedToRecoveredRate*I0 - deathRate*R0);

	return outputState;
};

int
main(int argc, char *  argv[])
{
	/*
	 *	Setup the interaction parameters
	 */
	vectorFieldParameters	sirModelInteractionParameters = {
		.susceptibleToInfectedRate = UxHwDoubleGammaDist(900, 0.000333),
		.infectedToRecoveredRate = UxHwDoubleGammaDist(400, 0.0005),
		.birthRate = 0.025,
		.deathRate = 0.025
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
		.integratorTimeStep = kSirModelTimestepSize,
		.numberOfSteps = (kSirModelFinalTime - kSirModelInitialTime)/kSirModelTimestepSize + 1
	};

	/*
	 *	Define an array of states to store all the information. Same for the timing info
	 */
	sirModelState *	sirModelEvolution = (sirModelState *) malloc(sizeof(sirModelState) * simTimeParams.numberOfSteps);
	double *	timeInstantArray = (double *) malloc(sizeof(double) * simTimeParams.numberOfSteps);

	sirModelEvolution[0] = sirState;
	timeInstantArray[0] = simTimeParams.initTime;
	for(size_t simIteration = 0; simIteration < simTimeParams.numberOfSteps - 1; simIteration++)
	{
		sirModelEvolution[simIteration + 1] = odeIntegrationStep(&sirModelEvolution[simIteration], &simTimeParams);
		timeInstantArray[simIteration + 1] = timeInstantArray[simIteration] + simTimeParams.integratorTimeStep;
	}

	sirModelState	finalState = sirModelEvolution[simTimeParams.numberOfSteps - 1];
	fprintf(stdout, "Final susceptibles: %lf\nFinal infected: %lf\nFinal recovered: %lf\n", finalState.susceptible, finalState.infected, finalState.recovered);

	/*
	 *	Store the full simulation to a file in the cloud storage:
	 */
	FILE *	outputFile = fopen("sd0/Signaloid/sirExampleOutput/sirModelOutput.txt", "w");
	if(outputFile==NULL)
	{
		fprintf(stderr, "Error has occured. File cannot be opened.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t simIteration = 0; simIteration < simTimeParams.numberOfSteps; simIteration++)
	{
		fprintf(outputFile, "%lf, %lf, %lf, %lf\n", timeInstantArray[simIteration], sirModelEvolution[simIteration].susceptible, sirModelEvolution[simIteration].infected, sirModelEvolution[simIteration].recovered);
	}
	fclose(outputFile);

	free(sirModelEvolution);
	free(timeInstantArray);

	return EXIT_SUCCESS;
}
