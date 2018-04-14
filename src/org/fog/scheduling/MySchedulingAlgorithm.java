package org.fog.scheduling;

import java.util.List;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.myGAEntities.MyGeneticAlgorithm;
import org.fog.scheduling.myGAEntities.MyPopulation;

public class MySchedulingAlgorithm {
	// Algorithm's Name
	public static final String GA = "Genetic Algorithm";
	public static final String LOCAL_SEARCH = "local search";
	public static final String TABU_SEARCH = "tabu search";

	// Trade-off Between Time and Cost
	public static final double TIME_WEIGHT = 0.5;

	// Genetic Algorithm's Parameters
	public static final int POPULATION_SIZE = 2000; // Number of population's individuals
	public static final int OFFSPRING_SIZE = 1800; // Number of offsprings
	public static final double MAX_TIME = 60.0; // Maximum executing time
	public static final int MAX_ITERATIONS = 1500; // Maximum iterations
	public static final int MUTATION_SIZE = (int) 0.1 * OFFSPRING_SIZE;
	public static final double SELECTION_PRESSURE = 2.0;
	public static final double DIGITS_ONE_RATE = 1.1;

	// Tabu Search's Parameters
	public static final int TABU_CONSTANT = 10;
	
	
	
	
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

		// Keep track of current generation
		int generationIndex = 0;

		MyPopulation offsprings;
		while (generationIndex < MAX_ITERATIONS) {
			System.out.println("\n------------- Generation " + generationIndex + " --------------");
			
//			offsprings = myGA.selectOffspringsRandomly2(population);
			offsprings = myGA.selectOffspringsPressure(population, SELECTION_PRESSURE);
			offsprings = myGA.crossoverOffspringsRandomTemplate(offsprings, DIGITS_ONE_RATE);
			
//			offsprings = myGA.crossoverOffsprings2Point(offsprings);
//			offsprings = myGA.mutateOffsprings(offsprings);
			population = myGA.selectNextGeneration(population, offsprings, fogDevices, cloudletList);
			
			
		

			// Print fittest individual from population
			System.out.println("\nBest solution of generation " + generationIndex + ": " + population.getIndividual(0).getFitness());
			System.out.println("Makespan: (" + myGA.getMinTime() + ")--" + population.getIndividual(0).getTime());
			System.out.println("TotalCost: (" + myGA.getMinCost() + ")--" + population.getIndividual(0).getCost());
			// Increment the current generation
			generationIndex++;
			// population.printPopulation();
		}

		/**
		 * We're out of the loop now, which means we have a perfect solution on our
		 * hands. Let's print it out to confirm that it is actually all ones, as
		 * promised.
		 */
		// population.printPopulation();

		// LocalSearchAlgorithm localSearch = new LocalSearchAlgorithm();
		// // Calculate the boundary of time and cost
		// localSearch.calcMinTimeCost(fogDevices, cloudletList);
		// localSearch.searchBestOne(population.getFittest(0), fogDevices,
		// cloudletList);

		System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
		System.out.println("Found solution in " + generationIndex + " generations");
		population.getIndividual(0).printGene();
		System.out.println("\nBest solution: " + population.getIndividual(0).getFitness());
	}
}
