package org.fog.scheduling;

import java.util.ArrayList;
import java.util.List;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.myLocalSearchAlgorithm.MyLocalSearchAlgorithm;
import org.fog.scheduling.myGAEntities.MyGeneticAlgorithm;
import org.fog.scheduling.myGAEntities.MyIndividual;
import org.fog.scheduling.myGAEntities.MyPopulation;
import org.fog.scheduling.myGAEntities.MyService;

public class MySchedulingAlgorithm {
	// Algorithm's Name
	public static final String GA = "Genetic Algorithm";
	public static final String HILL_CLIMBING = "Hill Climbing";
	public static final String TABU_SEARCH = "tabu search";
	public static final String SACRIFICED_HILL_CLIMBING = "Sacrified Hill Climbing";
	public static final String SIMULATED_ANNEALING = "Simulated Annealing";
	public static final String DEGRATED_CEILING = "Degrated Ceiling";
	public static final String COMBINATION_GA_LOCALSEARCH = "Combination of GA and Local Search";
	
	

	// Trade-off Between Time and Cost
	public static final double TIME_WEIGHT = 0.5;

	// Genetic Algorithm's Parameters
	public static final int POPULATION_SIZE = 3000; // Number of population's individuals
	public static final int OFFSPRING_SIZE = (int)(POPULATION_SIZE*0.9); // Number of offsprings
	public static final double MAX_TIME = 60.0; // Maximum executing time
	public static final int MAX_GENETIC_ITERATIONS = 400; // Maximum iterations
	public static final int MUTATION_SIZE = (int) 0.1 * OFFSPRING_SIZE;
	public static final double INITIAL_SELECTION_PRESSURE = 2.0;
//	public static final double ENDING_SELECTION_PRESSURE = 2.0;
	public static final double INITIAL_DIGITS_ONE_RATE = 0.9;
//	public static final double ENDING_DIGITS_ONE_RATE = 0.9;


	// Common parameters for Local Search
	public static final int LOCALSEARCH_ITERATIONS = 300; // 12000
	public static final int LOCALSEARCH_TIME = 120;
	
	// Parameters for Tabu Search
	public static final int TABU_LENGTH = 30;
	public static final int TABU_STABLE = 200;
	
	// Parameters for Sacrified Hill Climbing
//	public static final double INITIAL_SACRIFICE = 0.005;
	public static final double INITIAL_SACRIFICE = 0.005;
	public static final double DESCENDING_SPEED = 0.9;

	// Parameters for Simulated Annealing
	public static final double INITIAL_TEMPERATURE = 0.0005;
	public static final double ENDING_TEMPERATURE = 0.0005;

	
	// Combination of GA and Local Search
	public static void runGA_LocalSearch(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		// Create GA object
		MyGeneticAlgorithm myGA = new MyGeneticAlgorithm(POPULATION_SIZE, OFFSPRING_SIZE, MUTATION_SIZE);

		// Calculate the boundary of time and cost
		myGA.calcMinTimeCost(fogDevices, cloudletList);

		// Initialize the population
		MyPopulation population = myGA.initPopulation(cloudletList.size(), fogDevices.size() - 1);

		// Calculate fitness of each individuals and of population, 
		// sort the population descending after its fitness
		myGA.evalPopulation(population, fogDevices, cloudletList);

		population.printPopulation();

		// Save the offsprings of each generation
		MyPopulation offsprings;
		
		double startGA = System.currentTimeMillis();
		double timeGA;
		double fitnessGA;
		
		// Keep track of current generation
		int generationIndex = 0;
		while (generationIndex < MAX_GENETIC_ITERATIONS) {
			System.out.println("\n------------- Generation " + generationIndex + " --------------");

			// Select parents for offsprings
//			offsprings = myGA.selectBestOffsprings(population);
//			offsprings = myGA.selectOffspringsRandomlyIdentical(population);
//			offsprings = myGA.selectOffspringsRandomlyUniquely(population);
			offsprings = myGA.selectOffspringsPressure(population, INITIAL_SELECTION_PRESSURE);
			
			// Cross-over operation
			offsprings = myGA.crossoverOffspringsRandomTemplate(offsprings, INITIAL_DIGITS_ONE_RATE);
//			offsprings = myGA.crossoverOffsprings2Point(offsprings);
//			offsprings = myGA.crossoverOffsprings1Point(offsprings);
			
			// Mutation operation
			offsprings = myGA.mutateOffsprings(offsprings);
			
			// Replacement operation
			population = myGA.selectNextGeneration(population, offsprings, fogDevices, cloudletList);

			// Prints fittest individual from population
			System.out.println("\nBest solution of generation " + generationIndex + ": "
					+ population.getIndividual(0).getFitness());
			System.out.println("Makespan: (" + myGA.getMinTime() + ")--" + population.getIndividual(0).getTime());
			System.out.println("TotalCost: (" + myGA.getMinCost() + ")--" + population.getIndividual(0).getCost());
			// population.printPopulation();
			
			generationIndex++;
			
		}
		
		timeGA = (System.currentTimeMillis() - startGA)/1000;
		fitnessGA = population.getIndividual(0).getFitness();

		System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
		System.out.println("Found solution in " + generationIndex + " generations");
		population.getIndividual(0).printGene();
		System.out.println("\nBest solution: " + population.getIndividual(0).getFitness());
		

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();

		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);
		
		ArrayList<Double> timeRuns = new ArrayList<Double>();
		ArrayList<Double> localSearchFitness = new ArrayList<Double>();
		
		double maxFitness = Double.MIN_VALUE;
		for (int individualIndex = 0; individualIndex < 5; individualIndex++) {
			// Initiates an random individual
			MyIndividual individual = population.getIndividual(individualIndex);
			individual.printGene();
			
			double startTimeLocalSearch = System.currentTimeMillis();
			
//			individual = localSearch.hillClimbing(individual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
//			individual = localSearch.sacrificedHillClimbing(individual, INITIAL_SACRIFICE, DESCENDING_SPEED, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
//			individual = localSearch.simulatedAnnealing(individual, INITIAL_TEMPERATURE, ENDING_TEMPERATURE, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);		
			individual = localSearch.degradedCeiling(individual, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);		
			// Stable, iteration, time, tabuLength
//			individual = localSearch.tabuSearch(individual, fogDevices, cloudletList, 100, 10000, 20, 30);
			
			double timeRunLocalSearch = (System.currentTimeMillis() - startTimeLocalSearch)/1000;
			timeRuns.add(individualIndex, timeRunLocalSearch);
			localSearchFitness.add(individualIndex, individual.getFitness());
			
			if (individual.getFitness() > maxFitness)
				maxFitness = individual.getFitness();
		}
		
		// Calculate the average fitness
		double averageFitness = 0.0;
		for (int index = 0; index < localSearchFitness.size(); index++) {
			averageFitness += localSearchFitness.get(index);
		}
		averageFitness /= localSearchFitness.size();
		
		// Calculate the average time
		double averageTime = 0.0;
		for (int index = 0; index < timeRuns.size(); index++) {
			averageTime += timeRuns.get(index);
		}
		averageTime /= timeRuns.size();
		averageTime += timeGA;
		
		System.out.println("timeGA : " + timeGA);
		System.out.println("fitnessGA : " + fitnessGA);
		System.out.println("MaxIterationGA : " + MAX_GENETIC_ITERATIONS);
		System.out.println("timeRunsLocalSearch : " + timeRuns);
		System.out.println("localSearchFitness : " + localSearchFitness);
		System.out.println("LOCALSEARCH_ITERATIONS : " + LOCALSEARCH_ITERATIONS);
		System.out.println("averageFitness : " + averageFitness);
		System.out.println("BEST FITNESS : " + maxFitness);
		System.out.println("averageTime : " + averageTime);
		
		
	}


	// Genetic Algorithm
	public static void runGeneticAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		// Create GA object
				MyGeneticAlgorithm myGA = new MyGeneticAlgorithm(POPULATION_SIZE, OFFSPRING_SIZE, MUTATION_SIZE);

				// Calculate the boundary of time and cost
				myGA.calcMinTimeCost(fogDevices, cloudletList);

				// Initialize the population
				MyPopulation population = myGA.initPopulation(cloudletList.size(), fogDevices.size() - 1);

				// Calculate fitness of each individuals and of population, 
				// sort the population descending after its fitness
				myGA.evalPopulation(population, fogDevices, cloudletList);

				population.printPopulation();

				// Save the offsprings of each generation
				MyPopulation offsprings;
				
				double startGA = System.currentTimeMillis();
				double timeGA;
				double fitnessGA;
				
				// Keep track of current generation
				int generationIndex = 0;
				while (generationIndex < MAX_GENETIC_ITERATIONS) {
					System.out.println("\n------------- Generation " + generationIndex + " --------------");

					// Select parents for offsprings
//					offsprings = myGA.selectBestOffsprings(population);
//					offsprings = myGA.selectOffspringsRandomlyIdentical(population);
//					offsprings = myGA.selectOffspringsRandomlyUniquely(population);
					offsprings = myGA.selectOffspringsPressure(population, INITIAL_SELECTION_PRESSURE);
					
					// Cross-over operation
					offsprings = myGA.crossoverOffspringsRandomTemplate(offsprings, INITIAL_DIGITS_ONE_RATE);
//					offsprings = myGA.crossoverOffsprings2Point(offsprings);
//					offsprings = myGA.crossoverOffsprings1Point(offsprings);
					
					// Mutation operation
					offsprings = myGA.mutateOffsprings(offsprings);
					
					// Replacement operation
					population = myGA.selectNextGeneration(population, offsprings, fogDevices, cloudletList);

					// Prints fittest individual from population
					System.out.println("\nBest solution of generation " + generationIndex + ": "
							+ population.getIndividual(0).getFitness());
					System.out.println("Makespan: (" + myGA.getMinTime() + ")--" + population.getIndividual(0).getTime());
					System.out.println("TotalCost: (" + myGA.getMinCost() + ")--" + population.getIndividual(0).getCost());
					// population.printPopulation();
					
					generationIndex++;
					
				}
				
				timeGA = (System.currentTimeMillis() - startGA)/1000;
				fitnessGA = population.getIndividual(0).getFitness();

				System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
				System.out.println("Found solution in " + generationIndex + " generations");
				population.getIndividual(0).printGene();
				System.out.println("\nBest solution: " + population.getIndividual(0).getFitness());
				System.out.println("Time GA : " + timeGA);
				System.out.println("numIteration : " + MAX_GENETIC_ITERATIONS);
				
				
				
	}


	// Hill Climbing
	public static void runHillClimbing(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();

		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// Initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		
		double startTime = System.currentTimeMillis();
		
		individual = localSearch.hillClimbing(individual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
		
		double timeRun = (System.currentTimeMillis() - startTime)/1000;
		System.out.println("timeRun : " + timeRun);
		System.out.println("Best Individual : " + individual.getFitness());
		System.out.println("Iteration : " + LOCALSEARCH_ITERATIONS);
	}

	

	// Sacrificed Hill Climbing
	public static void runSacrificedHillClimbing(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		
		double startTime = System.currentTimeMillis();
		
		individual = localSearch.sacrificedHillClimbing(individual, INITIAL_SACRIFICE, DESCENDING_SPEED,
				LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
		
		double timeRun = (System.currentTimeMillis() - startTime)/1000;
		System.out.println("runSacrificedHillClimbing - timeRun : " + timeRun);
		System.out.println("Best Individual : " + individual.getFitness());
		System.out.println("Iteration : " + LOCALSEARCH_ITERATIONS);

	}

	
	

	
	// Simulated Annealing Using Many New Records
	public static void runSimulatedAnnealing(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		
		double startTime = System.currentTimeMillis();
		
		individual = localSearch.simulatedAnnealing(individual, INITIAL_TEMPERATURE, ENDING_TEMPERATURE,
				LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
		
		double timeRun = (System.currentTimeMillis() - startTime)/1000;
		System.out.println("runSimulatedAnnealing - timeRun : " + timeRun);
		System.out.println("Best Individual : " + individual.getFitness());
		System.out.println("Iteration : " + LOCALSEARCH_ITERATIONS);


	}
	

	

	// Degrated Ceiling
	public static void runDegratedCeiling(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		
		double startTime = System.currentTimeMillis();
		
		individual = localSearch.degradedCeiling(individual, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);

		
		double timeRun = (System.currentTimeMillis() - startTime)/1000;
		System.out.println("runDegratedCeiling - timeRun : " + timeRun);
		System.out.println("Best Individual : " + individual.getFitness());
		System.out.println("Iteration : " + LOCALSEARCH_ITERATIONS);
	}

	

	
	// Tabu Search algorithm
	public static void runTabuSearchAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		MyLocalSearchAlgorithm myLocalSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		myLocalSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiate an individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();

		// individual = myLocalSearch.tabuSearch(individual, fogDevices, cloudletList,
		// TABU_STABLE, LOCALSEARCH_ITERATIONS, LOCALSEARCH_TIME, TABU_LENGTH);
		// Stable, iteration, time, tabuLength
		individual = myLocalSearch.tabuSearch(individual, fogDevices, cloudletList, 100, LOCALSEARCH_ITERATIONS, 20, 30);
	}

}
