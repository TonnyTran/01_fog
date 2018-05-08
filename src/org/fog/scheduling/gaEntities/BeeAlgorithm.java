package org.fog.scheduling.gaEntities;

import java.util.ArrayList;
import java.util.List;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.SchedulingAlgorithm;
import org.fog.scheduling.localSearchAlgorithm.Pair;

/**
 * The GeneticAlgorithm class is our main abstraction for managing the
 * operations of the genetic algorithm. This class is meant to be
 * problem-specific, meaning that (for instance) the "calcFitness" method may
 * need to change from problem to problem.
 * 
 * This class concerns itself mostly with population-level operations, but also
 * problem-specific operations such as calculating fitness, testing for
 * termination criteria, and managing mutation and crossover operations (which
 * generally need to be problem-specific as well).
 * 
 * Generally, GeneticAlgorithm might be better suited as an abstract class or an
 * interface, rather than a concrete class as below. A GeneticAlgorithm
 * interface would require implementation of methods such as
 * "isTerminationConditionMet", "calcFitness", "mutatePopulation", etc, and a
 * concrete class would be defined to solve a particular problem domain. For
 * instance, the concrete class "TravelingSalesmanGeneticAlgorithm" would
 * implement the "GeneticAlgorithm" interface. This is not the approach we've
 * chosen, however, so that we can keep each chapter's examples as simple and
 * concrete as possible.
 * 
 * @author bkanber
 *
 */
public class BeeAlgorithm {
	private int populationSize;

	/**
	 * Mutation rate is the fractional probability than an individual gene will
	 * mutate randomly in a given generation. The range is 0.0-1.0, but is generally
	 * small (on the order of 0.1 or less).
	 */
	private double mutationRate;

	/**
	 * Crossover rate is the fractional probability that two individuals will "mate"
	 * with each other, sharing genetic information, and creating offspring with
	 * traits of each of the parents. Like mutation rate the rance is 0.0-1.0 but
	 * small.
	 */
	private double crossoverRate;

	/**
	 * Elitism is the concept that the strongest members of the population should be
	 * preserved from generation to generation. If an individual is one of the
	 * elite, it will not be mutated or crossover.
	 */
	private int elitismCount;

	private int numIteration;

	private double minTime;
	private double minCost;

	public BeeAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount) {
		this.populationSize = populationSize;
		this.mutationRate = mutationRate;
		this.crossoverRate = crossoverRate;
		this.elitismCount = elitismCount;
		this.numIteration = 0;
	}

	/**
	 * calculate the lower boundary of time and cost
	 *
	 */
	public void calcMinTimeCost(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		this.minTime = calcMinTime(fogDevices, cloudletList);
		this.minCost = calcMinCost(fogDevices, cloudletList);
	}

	private double calcMinCost(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		double minCost = 0;
		for (Cloudlet cloudlet : cloudletList) {
			double minCloudletCost = Double.MAX_VALUE;
			for (FogDevice fogDevice : fogDevices) {
				double cost = calcCost(cloudlet, fogDevice);
				if (minCloudletCost > cost) {
					minCloudletCost = cost;
				}
			}
			// the minCost is defined as the sum of all minCloudletCost
			minCost += minCloudletCost;
		}
		return minCost;
	}

	// the method calculates the cost (G$) when a fogDevice executes a cloudlet
	private double calcCost(Cloudlet cloudlet, FogDevice fogDevice) {
		double cost = 0;
		// cost includes the processing cost
		cost += fogDevice.getCharacteristics().getCostPerSecond() * cloudlet.getCloudletLength()
				/ fogDevice.getHost().getTotalMips();
		// cost includes the memory cost
		cost += fogDevice.getCharacteristics().getCostPerMem() * cloudlet.getMemRequired();
		// cost includes the bandwidth cost
		cost += fogDevice.getCharacteristics().getCostPerBw()
				* (cloudlet.getCloudletFileSize() + cloudlet.getCloudletOutputSize());
		return cost;
	}

	// the function calculate the lower bound of the solution about time execution
	private double calcMinTime(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		double minTime = 0;
		double totalLength = 0;
		double totalMips = 0;
		for (Cloudlet cloudlet : cloudletList) {
			totalLength += cloudlet.getCloudletLength();
		}
		for (FogDevice fogDevice : fogDevices) {
			totalMips += fogDevice.getHost().getTotalMips();
		}
		minTime = totalLength / totalMips;
		return minTime;
	}

	/**
	 * Initialize population
	 * 
	 * @param chromosomeLength
	 *            The length of the individuals chromosome
	 * @return population The initial population generated
	 */
	public Population initPopulation(int chromosomeLength, int maxValue) {
		// Initialize population
		Population population = new Population(this.populationSize, chromosomeLength, maxValue);
		return population;
	}

	/**
	 * Calculate fitness for an individual.
	 * 
	 * In this case, the fitness score is very simple: it's the number of ones in
	 * the chromosome. Don't forget that this method, and this whole
	 * GeneticAlgorithm class, is meant to solve the problem in the "AllOnesGA"
	 * class and example. For different problems, you'll need to create a different
	 * version of this method to appropriately calculate the fitness of an
	 * individual.
	 * 
	 * @param individual
	 *            the individual to evaluate
	 * @return double The fitness value for individual
	 */
	public double calcFitness(Individual individual, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// clear the fogDevice - task list before calculate
		for (FogDevice fogDevice : fogDevices) {
			fogDevice.getCloudletListAssignment().clear();
		}

		// Loop over individual's genes to all the task assigned to the fogDevice
		for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
			// add current cloudlet to fog device respectively
			fogDevices.get(individual.getGene(geneIndex)).getCloudletListAssignment().add(cloudletList.get(geneIndex));
		}

		// Calculate makespan and cost
		double makespan = 0;
		double execTime = 0;
		double totalCost = 0;
		for (FogDevice fogDevice : fogDevices) {
			double totalLength = 0;
			for (Cloudlet cloudlet : fogDevice.getCloudletListAssignment()) {
				totalLength += cloudlet.getCloudletLength();
				// the total cost is sum of the cost execution of each cloudlet
				totalCost += calcCost(cloudlet, fogDevice);
			}
			// execTime is the time that fogDevice finishes its list cloudlet assignment
			execTime = totalLength / fogDevice.getHostList().get(0).getTotalMips();
			// makespan is defined as when the last cloudlet finished or when all fogDevices
			// finish its work.
			if (execTime > makespan) {
				makespan = execTime;
			}
		}

		// store makespan
		individual.setTime(makespan);
		// store cost
		individual.setCost(totalCost);

		// Calculate fitness
		double fitness = SchedulingAlgorithm.TIME_WEIGHT * minTime / makespan
				+ (1 - SchedulingAlgorithm.TIME_WEIGHT) * minCost / totalCost;

		// Store fitness
		individual.setFitness(fitness);
		return fitness;
	}

	/**
	 * Evaluate the whole population
	 * 
	 * Essentially, loop over the individuals in the population, calculate the
	 * fitness for each, and then calculate the entire population's fitness. The
	 * population's fitness may or may not be important, but what is important here
	 * is making sure that each individual gets evaluated.
	 * 
	 * @param population
	 *            the population to evaluate
	 */
	public Population evalPopulation(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// System.out.println("Min time evaluated: " + minTime);
		// System.out.println("Min cost evaluated: " + minCost);
		double populationFitness = 0;

		// Loop over population evaluating individuals and summing population fitness
		for (Individual individual : population.getPopulation()) {
			populationFitness += calcFitness(individual, fogDevices, cloudletList);
		}

		// sort population with increasing fitness value
		population.sortPopulation();

		population.setPopulationFitness(populationFitness);
		return population;
	}

	/**
	 * Check if population has met termination condition
	 * 
	 * For this simple problem, we know what a perfect solution looks like, so we
	 * can simply stop evolving once we've reached a fitness of one.
	 * 
	 * @param population
	 * @return boolean True if termination condition met, otherwise, false
	 */
	public boolean isTerminationConditionMet(Population population) {
		for (Individual individual : population.getPopulation()) {
			if (individual.getFitness() == 1) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Select parent for crossover
	 * 
	 * @param population
	 *            The population to select parent from
	 * @return The individual selected as a parent
	 */
	public Individual selectIndividual(Population population) {
		// Get individuals
		List<Individual> individuals = population.getPopulation();

		// Spin roulette wheel
		double populationFitness = population.getPopulationFitness();
		double rouletteWheelPosition = Math.random() * populationFitness;

		// Find parent
		double spinWheel = 0;
		for (Individual individual : individuals) {
			spinWheel += individual.getFitness();
			if (spinWheel >= rouletteWheelPosition) {
				return individual;
			}
		}
		return individuals.get(population.size() - 1);
	}

	/**
	 * Apply crossover to population
	 * 
	 * Crossover, more colloquially considered "mating", takes the population and
	 * blends individuals to create new offspring. It is hoped that when two
	 * individuals crossover that their offspring will have the strongest qualities
	 * of each of the parents. Of course, it's possible that an offspring will end
	 * up with the weakest qualities of each parent.
	 * 
	 * This method considers both the GeneticAlgorithm instance's crossoverRate and
	 * the elitismCount.
	 * 
	 * The type of crossover we perform depends on the problem domain. We don't want
	 * to create invalid solutions with crossover, so this method will need to be
	 * changed for different types of problems.
	 * 
	 * This particular crossover method selects random genes from each parent.
	 * 
	 * @param population
	 *            The population to apply crossover to
	 * @return The new population
	 */
	public Population crossoverPopulation(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		// Create new population
		List<Individual> newPopulation = new ArrayList<Individual>();

		newPopulation.clear();
		// Loop over current population by fitness
		for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
			Individual parent1 = population.getFittest(populationIndex);

			// Apply crossover to this individual?
			if (this.crossoverRate > Math.random()) {
				// Initialize offspring
				Individual offspring = new Individual(parent1.getChromosomeLength());

				// Find second parent
				Individual parent2 = selectIndividual(population);
				offspring = crossoverRandom(parent1, parent2);

				if (parent1.getFitness() <= calcFitness(offspring, fogDevices, cloudletList)
						&& !doesPopupationIncludeIndividual(population, offspring)) {
					newPopulation.add(offspring);
				} else {
					newPopulation.add(parent1);

				}
			} else {
				newPopulation.add(population.getFittest(populationIndex));
			}
		}
		population.getPopulation().clear();
		population.setPopulation(newPopulation);

		// System.out.println("--------AFTER CROSSOVER--------");
		// population.printPopulation();
		return population;
	}

	public Individual crossoverRandom(Individual parent1, Individual parent2) {
		Individual offspring = new Individual(parent1.getChromosomeLength());

		for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {

			if (Service.rand(0, 100) > 50) {
				offspring.setGene(geneIndex, parent1.getGene(geneIndex));
			} else {
				offspring.setGene(geneIndex, parent2.getGene(geneIndex));
			}
		}
		return offspring;
	}

	public Population searchFood(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		Individual individual;
		Individual newIndividual;

		// Number of scout bee
		int numScout;

		// Neighborhood search from elite bee (dynamic search)
		for (int populationIndex = 0; populationIndex < elitismCount; populationIndex++) {
			numScout = 100;
			individual = population.getFittest(populationIndex);

			for (int j = 0; j < numScout; j++) {

				// Create new individual identical to current one
				newIndividual = new Individual(individual.getChromosomeLength());
				for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++)
					newIndividual.setGene(geneIndex, individual.getGene(geneIndex));

				// Change random job to random node
				int job = Service.rand(0, cloudletList.size() - 1);
				int fogOld = newIndividual.getGene(job);
				int fogNew = Service.rand(0, fogDevices.size() - 1);
				newIndividual.setGene(job, fogNew);
				calcFitness(newIndividual, fogDevices, cloudletList);
				if (newIndividual.getFitness() > individual.getFitness()) {
					for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++)
						individual.setGene(geneIndex, newIndividual.getGene(geneIndex));
					calcFitness(individual, fogDevices, cloudletList);
					continue;
				}

				// Switch 2 random job from 2 random node
				int offset = Service.rand(0, fogDevices.size() - 1);
				for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
					if (newIndividual.getGene((geneIndex + offset) % fogDevices.size()) == fogNew) {
						newIndividual.setGene((geneIndex + offset) % fogDevices.size(), fogOld);
						break;
					}

				}
				calcFitness(newIndividual, fogDevices, cloudletList);
				if (newIndividual.getFitness() > individual.getFitness()) {
					for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++)
						individual.setGene(geneIndex, newIndividual.getGene(geneIndex));
					calcFitness(individual, fogDevices, cloudletList);
				}
			}

		}

		// Re-initailize lowest fitness bee
		for (int populationIndex = (int) (0.85 * populationSize); populationIndex < population
				.size(); populationIndex++) {
			individual = population.getFittest(populationIndex);

			for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++)
				individual.setGene(geneIndex, Service.rand(0, fogDevices.size() - 1));
		}
		return population;
	}

	public Population mutatePopulation(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {

			if (this.mutationRate > Math.random() && populationIndex >= this.elitismCount) {
				Individual individual = population.getFittest(populationIndex);
				individual.setGene(Service.rand(0, cloudletList.size() - 1), Service.rand(0, fogDevices.size() - 1));
			}
		}
		return population;
	}

	public Population changeClone(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		int currentIndividual = 0;
		int nextIndividual = 1;
		Individual first;
		Individual second;
		boolean isClone = true;
		while (currentIndividual < population.size() && nextIndividual < population.size()) {
			first = population.getIndividual(currentIndividual);
			second = population.getIndividual(nextIndividual);

			if (first.getFitness() == second.getFitness()) {
				for (int geneIndex = 0; geneIndex < first.getChromosomeLength(); geneIndex++) {
					if (first.getGene(geneIndex) != second.getGene(geneIndex)) {
						isClone = false;
						break;
					}
				}
			} else {
				isClone = false;
			}
			if (isClone == true) {
				for (int geneIndex = 0; geneIndex < second.getChromosomeLength(); geneIndex++)
					second.setGene(geneIndex, Service.rand(0, fogDevices.size() - 1));
				nextIndividual++;

			} else {
				currentIndividual = nextIndividual;
				nextIndividual++;
				isClone = true;
			}
		}
		return population;
	}

	public boolean isCloneExist(Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		int currentIndividual = 0;
		int nextIndividual = 1;
		Individual first;
		Individual second;
		boolean isClone = true;
		while (currentIndividual < population.size() && nextIndividual < population.size()) {
			first = population.getIndividual(currentIndividual);
			second = population.getIndividual(nextIndividual);

			if (first.getFitness() == second.getFitness()) {
				for (int geneIndex = 0; geneIndex < first.getChromosomeLength(); geneIndex++) {
					if (first.getGene(geneIndex) != second.getGene(geneIndex)) {
						isClone = false;
						break;
					}
				}
			} else {
				isClone = false;
			}
			if (isClone == true) {
				return true;

			} else {
				currentIndividual++;
				nextIndividual++;
				isClone = true;
			}
		}
		return false;
	}

	public double getMinTime() {
		return this.minTime;
	}

	public double getMinCost() {
		return this.minCost;
	}

	public boolean doesPopupationIncludeIndividual(Population population, Individual individual) {
		boolean include = false;
		for (int index = 0; index < population.size(); index++) {
			boolean similar = true;
			if (individual.getFitness() == population.getIndividual(index).getFitness()) {
				for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
					if (individual.getGene(geneIndex) != population.getIndividual(index).getGene(geneIndex)) {
						similar = false;
					}
				}
				if (similar == true) {
					include = true;
					break;
				}
			}
			if (include == true)
				break;
		}
		return include;
	}

	public void incNumIteration() {
		numIteration++;
		return;

	}

	public Population randommizePopulation(int offset, Population population, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		List<Individual> newPopulation = new ArrayList<Individual>();

		newPopulation.clear();

		for (int populationIndex = 0; populationIndex < offset; populationIndex++) {

			newPopulation.add(population.getFittest(populationIndex));

		}

		while (newPopulation.size() < populationSize)
			newPopulation.add(new Individual(cloudletList.size(), fogDevices.size() - 1));
		population.getPopulation().clear();
		population.setPopulation(newPopulation);

		return population;

	}

	public static void main(String[] args) {

	}
}
