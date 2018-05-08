package org.fog.scheduling;

import java.util.ArrayList;
import java.util.List;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.gaEntities.BeeAlgorithm;
import org.fog.scheduling.gaEntities.GeneticAlgorithm;
import org.fog.scheduling.gaEntities.Individual;
import org.fog.scheduling.gaEntities.Population;
import org.fog.scheduling.gaEntities.Service;
import org.fog.scheduling.localSearchAlgorithm.LocalSearchAlgorithm;
import org.fog.scheduling.localSearchAlgorithm.Pair;

import java.io.*;

public class SchedulingAlgorithm {
	public static final String output = "test.txt";
	// Algorithm name
	public static final String GA = "Genetic Algorithm";
	public static final String LOCAL_SEARCH = "local search";
	public static final String TABU_SEARCH = "tabu search";
	public static final String BEE = "Bee Algorithm";
	// the weight value defines the trade-off between time and cost
	public static final double TIME_WEIGHT = 0.5;

	// GA parameters
	public static final int NUMBER_INDIVIDUAL = 100;
	public static final int NUMBER_ITERATION = 1500;

	public static final double MUTATION_RATE = 0.1;
	public static final double CROSSOVER_RATE = 0.95;
	public static final int NUMBER_ELITISM_INDIVIDUAL = 1;

	// Bee parameter
	public static final int NUMBER_INDIVIDUAL2 = 100;
	public static final int NUMBER_ITERATION2 = 1500;

	public static final double MUTATION_RATE2 = 0.15;
	public static final double CROSSOVER_RATE2 = 1;
	public static final int NUMBER_ELITISM_INDIVIDUAL2 = 10;

	// Tabu Search parameters
	public static final int TABU_CONSTANT = 10;

	// GA run
	public static void runGeneticAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		// Create GA object
		GeneticAlgorithm ga = new GeneticAlgorithm(NUMBER_INDIVIDUAL, MUTATION_RATE, CROSSOVER_RATE,
				NUMBER_ELITISM_INDIVIDUAL);

		// Calculate the boundary of time and cost
		ga.calcMinTimeCost(fogDevices, cloudletList);

		// Initialize population
		Population population = ga.initPopulation(cloudletList.size(), fogDevices.size() - 1);

		// Evaluate population
		ga.evalPopulation(population, fogDevices, cloudletList);

		population.printPopulation();

		// Keep track of current generation
		int generation = 0;

		/**
		 * Start the evolution loop
		 * 
		 * Every genetic algorithm problem has different criteria for finishing. In this
		 * case, we know what a perfect solution looks like (we don't always!), so our
		 * isTerminationConditionMet method is very straightforward: if there's a member
		 * of the population whose chromosome is all ones, we're done!
		 */
		while (generation < NUMBER_ITERATION) {
			System.out.println("\n------------- Generation " + generation + " --------------");
			// population.printPopulation();
			// Apply crossover
			population = ga.crossoverPopulation(population, fogDevices, cloudletList);

			// Evaluate population
			population = ga.evalPopulation(population, fogDevices, cloudletList);

			// //select individuals to form the population in current generation.
			// ga.selectPopulation(population);
			//
			// Apply mutation
			population = ga.mutatePopulation(population, fogDevices, cloudletList);
			//
			// Evaluate population
			ga.evalPopulation(population, fogDevices, cloudletList);

			// select individuals to form the population in current generation.
			// ga.selectPopulation(population);

			population.getFittest(0).printGene();

			// Print fittest individual from population
			System.out.println(
					"\nBest solution of generation " + generation + ": " + population.getFittest(0).getFitness());
			System.out.println("Makespan: (" + ga.getMinTime() + ")--" + population.getFittest(0).getTime());
			System.out.println("TotalCost: (" + ga.getMinCost() + ")--" + population.getFittest(0).getCost());
			// Increment the current generation
			generation++;
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
		System.out.println("Found solution in " + generation + " generations");
		population.getFittest(0).printGene();
		System.out.println("\nBest solution: " + population.getFittest(0).getFitness());
		writeOutput(output, population.getFittest(0).getFitness() + "\t");
	}

	// bee
	public static void runBeeAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		// Create GA object
		double currentFitness = 0;
		int count = 0;
		int count2 = 0;
		int maxIterBeforeSearch = 600;
		boolean searchFood = false;
		BeeAlgorithm bee = new BeeAlgorithm(NUMBER_INDIVIDUAL2, MUTATION_RATE2, CROSSOVER_RATE2,
				NUMBER_ELITISM_INDIVIDUAL2);

		// Calculate the boundary of time and cost
		bee.calcMinTimeCost(fogDevices, cloudletList);

		// Initialize population
		Population population = bee.initPopulation(cloudletList.size(), fogDevices.size() - 1);

		// Evaluate population
		bee.evalPopulation(population, fogDevices, cloudletList);

		population.printPopulation();

		// Keep track of current generation
		int generation = 0;
		/**
		 * Start the evolution loop
		 * 
		 * Every genetic algorithm problem has different criteria for finishing. In this
		 * case, we know what a perfect solution looks like (we don't always!), so our
		 * isTerminationConditionMet method is very straightforward: if there's a member
		 * of the population whose chromosome is all ones, we're done!
		 */

		while (generation < NUMBER_ITERATION2) {
			System.out.println("\n------------- Generation " + generation + " --------------");

			bee.incNumIteration();

			//Remove same individual after certain interation with same max fitness
			if (currentFitness == population.getIndividual(0).getFitness()) {
				if (count == 10) {
					population = bee.changeClone(population, fogDevices, cloudletList);
					population = bee.evalPopulation(population, fogDevices, cloudletList);
					count = 0;
				} else
					count++;

			} else {
				currentFitness = population.getIndividual(0).getFitness();
				count = 0;
			}
			
			
			
			
			//Crossover population
			population = bee.crossoverPopulation(population, fogDevices, cloudletList);
			population = bee.evalPopulation(population, fogDevices, cloudletList);
			
			//Mutate population
			population = bee.mutatePopulation(population, fogDevices, cloudletList);
			bee.evalPopulation(population, fogDevices, cloudletList);
			
			//Search food after certain interation
			if (count2 == maxIterBeforeSearch) {
				searchFood = true;
				count2 = 0;
				maxIterBeforeSearch = (int) (maxIterBeforeSearch * 0.5);
				if (maxIterBeforeSearch == 0)
					maxIterBeforeSearch = 1;
			} else
				count2++;
			
			if (searchFood == true) {
				population = bee.searchFood(population, fogDevices, cloudletList);
				population = bee.evalPopulation(population, fogDevices, cloudletList);
				searchFood = false;
			}

			population.getFittest(0).printGene();

			// Print fittest individual from population
			System.out.println(
					"\nBest solution of generation " + generation + ": " + population.getFittest(0).getFitness());
			System.out.println("Makespan: (" + bee.getMinTime() + ")--" + population.getFittest(0).getTime());
			System.out.println("TotalCost: (" + bee.getMinCost() + ")--" + population.getFittest(0).getCost());
			// Increment the current generation
			generation++;
			// population.printPopulation();

		}

		
		System.out.println(">>>>>>>>>>>>>>>>>>>RESULTS<<<<<<<<<<<<<<<<<<<<<");
		System.out.println("Found solution in " + generation + " generations");
		population.getFittest(0).printGene();
		System.out.println("\nBest solution: " + population.getFittest(0).getFitness());
		writeOutput(output, population.getFittest(0).getFitness() + "\t");
	}

	
	// local search algorithm
	public static void runLocalSearchAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {

		LocalSearchAlgorithm localSearch = new LocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiate an individual
		Individual individual = new Individual(cloudletList.size(), fogDevices.size() - 1, 1);
		individual.printGene();
		individual = localSearch.hillCliming(individual, fogDevices, cloudletList);
	}

	// Tabu Search algorithm
	public static void runTabuSearchAlgorithm(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		LocalSearchAlgorithm localSearch = new LocalSearchAlgorithm();
		// Calculate the boundary of time and cost
		localSearch.calcMinTimeCost(fogDevices, cloudletList);

		// initiate an individual
		Individual individual = new Individual(cloudletList.size(), fogDevices.size() - 1, 1);
		individual.printGene();
		individual = localSearch.tabuSearch(individual, fogDevices, cloudletList, 100, 100000, 300, 30);
	}

	// write output to file
	public static void writeOutput(String output, String data) {

		try {

			File file = new File(output);
			FileWriter fw;
			fw = new FileWriter(file.getAbsoluteFile(), true);
			BufferedWriter bw = new BufferedWriter(fw);
			bw.append(data);
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
