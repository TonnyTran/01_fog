package org.fog.scheduling;

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
	public static final String SIMULATED_ANNEALING_BEST_NEW_RECORD = "Simulated Annealing Best New Record";
	public static final String SIMULATED_ANNEALING_MANY_NEW_RECORDS = "Simulated Annealing Many New Records";
	public static final String DEGRATED_CEILING_BEST_NEW_RECORD = "Degrated Ceiling Best New Record";
	public static final String DEGRATED_CEILING_MANY_NEW_RECORDS = "Degrated Ceiling Many New Records";
	public static final String COMBINATION_GA_LOCALSEARCH = "Combination of GA and Local Search";
	
	

	// Trade-off Between Time and Cost
	public static final double TIME_WEIGHT = 0.5;

	// Genetic Algorithm's Parameters
	public static final int POPULATION_SIZE = 3000; // Number of population's individuals
	public static final int OFFSPRING_SIZE = (int)(POPULATION_SIZE*0.9); // Number of offsprings
	public static final double MAX_TIME = 60.0; // Maximum executing time
	public static final int MAX_GENETIC_ITERATIONS = 500; // Maximum iterations
	public static final int MUTATION_SIZE = (int) 0.1 * OFFSPRING_SIZE;
	public static final double INITIAL_SELECTION_PRESSURE = 2.0;
	public static final double ENDING_SELECTION_PRESSURE = 2.0;
	public static final double INITIAL_DIGITS_ONE_RATE = 0.9;
	public static final double ENDING_DIGITS_ONE_RATE = 0.9;


	// Common parameters for Local Search
	public static final int LOCALSEARCH_ITERATIONS = 12000;
	public static final int LOCALSEARCH_TIME = 120;
	
	// Parameters for Tabu Search
	public static final int TABU_LENGTH = 30;
	public static final int TABU_STABLE = 200;
	
	// Parameters for Sacrified Hill Climbing
//	public static final double INITIAL_SACRIFICE = 0.005;
	public static final double INITIAL_SACRIFICE = 0.008;
	public static final double DESCENDING_SPEED = 0.9;

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

		System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
		System.out.println("Found solution in " + generationIndex + " generations");
		population.getIndividual(0).printGene();
		System.out.println("\nBest solution: " + population.getIndividual(0).getFitness());
		

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();

		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);
		
		double maxFitness = Double.MIN_VALUE;
		for (int individualIndex = 0; individualIndex < 10; individualIndex++) {
			// Initiates an random individual
			MyIndividual individual = population.getIndividual(individualIndex);
			individual.printGene();
			individual = localSearch.hillClimbing(individual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
//			individual = localSearch.sacrificedHillClimbing(individual, INITIAL_SACRIFICE, DESCENDING_SPEED, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
					
//			individual = localSearch.simulatedAnnealing_ManyNewRecords(individual, INITIAL_TEMPERATURE, ENDING_TEMPERATURE, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);		
//			individual = localSearch.degradedCeiling_ManyNewRecords(individual, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);		
			
			// Stable, iteration, time, tabuLength
//			individual = localSearch.tabuSearch(individual, fogDevices, cloudletList, 100, 10000, 20, 30);

			
			if (individual.getFitness() > maxFitness)
				maxFitness = individual.getFitness();
		}
		
		System.out.println("BEST FITNESS : " + maxFitness);
	}


	// Genetic Algorithm
	public static void runGeneticAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		// Create GA object
		MyGeneticAlgorithm myGA = new MyGeneticAlgorithm(POPULATION_SIZE, OFFSPRING_SIZE, MUTATION_SIZE);

		// Calculate the boundary of time and cost
		myGA.calcMinTimeCost(fogDevices, cloudletList);

		// Initialize population
		MyPopulation population = myGA.initPopulation(cloudletList.size(), fogDevices.size() - 1);

		// Evaluate population
		myGA.evalPopulation(population, fogDevices, cloudletList);

		population.printPopulation();

		
		double currentDigitOneRate = INITIAL_DIGITS_ONE_RATE;
		double descendingSpeed = (INITIAL_DIGITS_ONE_RATE - ENDING_DIGITS_ONE_RATE)/MAX_GENETIC_ITERATIONS;			
		
		double currentPressure = INITIAL_SELECTION_PRESSURE;
		double pressureAscendingSpeed = (ENDING_SELECTION_PRESSURE - INITIAL_SELECTION_PRESSURE)/MAX_GENETIC_ITERATIONS;
		
		
		// Keep track of current generation
		int generationIndex = 0;

		MyPopulation offsprings;
		while (generationIndex < MAX_GENETIC_ITERATIONS) {
			System.out.println("\n------------- Generation " + generationIndex + " --------------");

//			offsprings = myGA.selectOffspringsRandomlyUniquely(population);
			offsprings = myGA.selectOffspringsPressure(population, currentPressure);
			offsprings = myGA.crossoverOffspringsRandomTemplate(offsprings, currentDigitOneRate);

//			 offsprings = myGA.crossoverOffsprings2Point(offsprings);
			 offsprings = myGA.mutateOffsprings(offsprings);
			population = myGA.selectNextGeneration(population, offsprings, fogDevices, cloudletList);

			// Prints fittest individual from population
			System.out.println("\nBest solution of generation " + generationIndex + ": "
					+ population.getIndividual(0).getFitness());
			System.out.println("Makespan: (" + myGA.getMinTime() + ")--" + population.getIndividual(0).getTime());
			System.out.println("TotalCost: (" + myGA.getMinCost() + ")--" + population.getIndividual(0).getCost());
			// Increment the current generation
			generationIndex++;
			// population.printPopulation();
			
			
//			currentDigitOneRate = currentDigitOneRate - descendingSpeed;
			
			currentPressure = currentPressure + pressureAscendingSpeed;
			
			System.out.println("CurrentDigitOneRate : " + currentDigitOneRate);
			System.out.println("currentPressure : " + currentPressure);
		}

		System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
		System.out.println("Found solution in " + generationIndex + " generations");
		population.getIndividual(0).printGene();
		System.out.println("\nBest solution: " + population.getIndividual(0).getFitness());
	}


	// Hill Climbing
	public static void runHillClimbing(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();

		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// Initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.hillClimbing(individual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
	}

	

	// Sacrificed Hill Climbing
	public static void runSacrificedHillClimbing(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.sacrificedHillClimbing(individual, INITIAL_SACRIFICE, DESCENDING_SPEED,
				LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);

	}

	
	
	
	
	
	// Simulated Annealing Using Best New Record
	public static void runSimulatedAnnealing_BestNewRecord(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.simulatedAnnealing_BestNewRecord(individual, INITIAL_TEMPERATURE, ENDING_TEMPERATURE,
				LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
		
		MyIndividual clonedIndividual;
		double bestFitness = individual.getFitness();
		
		for (int index = 0; index < 10; index++)
		{
			clonedIndividual = MyService.clonedIndividual(individual);
			clonedIndividual = localSearch.hillClimbing(clonedIndividual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
			if (bestFitness < clonedIndividual.getFitness())
				bestFitness = clonedIndividual.getFitness();
		}
			
		System.out.println("Fitness of Annealing : " + individual.getFitness());
		System.out.println("Best Fitness : " + bestFitness);

	}
	
	
	
	
	
	// Simulated Annealing Using Many New Records
	public static void runSimulatedAnnealing_ManyNewRecords(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.simulatedAnnealing_ManyNewRecords(individual, INITIAL_TEMPERATURE, ENDING_TEMPERATURE,
				LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
		
		MyIndividual clonedIndividual;
		double bestFitness = individual.getFitness();
		
		for (int index = 0; index < 10; index++)
		{
			clonedIndividual = MyService.clonedIndividual(individual);
			clonedIndividual = localSearch.hillClimbing(clonedIndividual,LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);
			if (bestFitness < clonedIndividual.getFitness())
				bestFitness = clonedIndividual.getFitness();
		}
			
		System.out.println("Fitness of Annealing : " + individual.getFitness());
		System.out.println("Best Fitness : " + bestFitness);

	}
	

	
	
	
	
	// Degrated Ceiling
	public static void runDegratedCeiling_BestNewRecord(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.degradedCeiling_BestNewRecord(individual, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);

	}
	
	
	
	
	
	// Degrated Ceiling
	public static void runDegratedCeiling_ManyNewRecords(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		MyLocalSearchAlgorithm localSearch = new MyLocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiates an random individual
		MyIndividual individual = new MyIndividual(cloudletList.size(), fogDevices.size() - 1, true);
		individual.printGene();
		individual = localSearch.degradedCeiling_ManyNewRecords(individual, LOCALSEARCH_ITERATIONS, fogDevices, cloudletList);

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
		individual = myLocalSearch.tabuSearch(individual, fogDevices, cloudletList, 100, 10000, 20, 30);
	}

}
