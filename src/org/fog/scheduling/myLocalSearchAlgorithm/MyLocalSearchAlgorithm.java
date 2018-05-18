package org.fog.scheduling.myLocalSearchAlgorithm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.MyStatistic;
import org.fog.scheduling.myGAEntities.MyIndividual;
import org.fog.scheduling.myGAEntities.MyService;

/**
 * @author DungBachViet
 *
 */
public class MyLocalSearchAlgorithm {
	private double minTime;
	private double minCost;

	public MyLocalSearchAlgorithm() {
	}
	
	
	/**
	 * Local Search using Sacrificed Hill Climbing Strategy
	 * 
	 * @param individual : The initial individual of Local Search
	 * @param INITIAL_SACRIFICE : The value of initial sacrifice
	 * @param DESCENDING_SPEED : The value of descending speed of the sacrifice
	 * @param numIterations : The max iterations of Local Search
	 * @param fogDevices : List of Fog Devices
	 * @param cloudletList : List of Tasks (need to be assigned on a specific Fog Device)
	 * 
	 * @return maxIndividual : The best individual in Local Search process
	 */
	public MyIndividual sacrificedHillClimbing(MyIndividual individual, double INITIAL_SACRIFICE,
			double DESCENDING_SPEED, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		
		// Sets up the initial sacrifice (it descends over time)
		double currentSacrifice = INITIAL_SACRIFICE; // 0.005
		
		// Calculate and save the current fitness of individual 
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);;
		
		// Saves the best individual and its best fitness over time
		double bestFitness = currentFitness;
		MyIndividual bestIndividual = MyService.clonedIndividual(individual);
		
		// Neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// Get the length of gene and its value domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
		ArrayList<Double> sacrificedHillClimbingFitness = new ArrayList<Double>();
		
		// Do iterations of local search and set up the initial time
		double start = System.currentTimeMillis();
		int iterationIndex = 0;
		do {
			
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Descending the degree of sacrifice over iterations
			if (iterationIndex % 300 == 0)
				currentSacrifice = currentSacrifice * DESCENDING_SPEED;
			
			// Clear all local moves at previous step
			listLocalMoves.clear();

			// Save the better assigning neighbors among assigning neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					
					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					// Checks condition of taking part in the better neighbors
					if (neighborFitness > (bestFitness - currentSacrifice)) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Move back to the state of current individual
					oldGeneValue = individual.getGene(geneIndex);
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}
			
			// Save the better swapping neighbors among swapping neighbors
			int geneIndex1, geneIndex2;
			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					// Checks condition of taking part in the better neighbors
					if (neighborFitness > (bestFitness - currentSacrifice)) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Move back to state of current individual
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());
				
			// Chooses randomly a better neighbor to become the next individual
			if (!listLocalMoves.isEmpty()) {

				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);

				if (abstractLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

			} else { // Restarts if list of local moves is empty
				System.out.println("\nRestart : Stop  Sacrificed Hill Climbing...\n");
				break; 
			}

			// Get fitness of the new individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Update the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
				MyService.copyAttributes(bestIndividual, individual);
			}
			
			sacrificedHillClimbingFitness.add(iterationIndex, currentFitness);
			
			// Create a new neighbor of this new individual
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}
			
			individual.printGene();
			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			iterationIndex++;
			
		} while (iterationIndex < numIterations);

		System.out.println("Running Time : " + (System.currentTimeMillis() - start));
		MyStatistic.sacrificedHillClimbingFitness = sacrificedHillClimbingFitness;
		
		return bestIndividual;
	}

	

	/**
	 * Local Search using Hill Climbing Strategy
	 * 
	 * @param individual : The initial individual of Local Search
	 * @param numIterations : The max iterations of Local Search
	 * @param fogDevices : List of Fog Devices
	 * @param cloudletList : List of Tasks (need to be assigned on a specific Fog Device)
	 *
	 * @return maxIndividual : The best individual in Local Search process
	 */
	public MyIndividual hillClimbing(MyIndividual individual, int numIterations, 
			List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		
		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// Calculate and save the current fitness of individual 
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);;
		
		// Saves the best individual and its best fitness over time
		double bestFitness = currentFitness;
		MyIndividual bestIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Create neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Get the length of gene and its value domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
		
		ArrayList<Double> hillClimbingFitness = new ArrayList<Double>();

		
		// Do iterations of local search and set up the initial time
		double start = System.currentTimeMillis();
		int iterationIndex = 0;
		do {
			
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();

			// Save the better assigning neighbors among assigning neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					
					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					// Checks condition of taking part in the better neighbors
					if (neighborFitness > currentFitness) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Move back to the state of current individual
					oldGeneValue = individual.getGene(geneIndex);
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			
			// Save the better swapping neighbors among swapping neighbors
			int geneIndex1, geneIndex2;
			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					// Checks condition of joining the better neighbors
					if (neighborFitness > currentFitness) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Move back to state of current individual
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			// Chooses randomly a better neighbor to become the next individual
			if (!listLocalMoves.isEmpty()) {

				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);

				if (abstractLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

			} else { // Restarts if list of local moves is empty
				System.out.println("\nRestart ...\n");
				break;
			}

			// Get fitness of the new individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Update the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
				MyService.copyAttributes(bestIndividual, individual);
			}
			
			hillClimbingFitness.add(iterationIndex, currentFitness);

			
			// Create a new neighbor of this new individual
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}
			
			individual.printGene();
			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			iterationIndex++;
			
		} while (iterationIndex < numIterations);

		System.out.println("Running Time : " + (System.currentTimeMillis() - start));
		MyStatistic.hillClimbingFitness = hillClimbingFitness;

		return bestIndividual;
	}

	
	
	
	/**
	 * Local Search using Simulated Annealing Strategy using many new records (sd nhiều kỷ lục mới so với kỷ lục trước)
	 * 
	 * @param individual : The initial individual of Local Search
	 * @param INITIAL_TEMPERATURE : The initial temperature of frozen schedule
	 * @param ENDING_TEMPERATURE : The ending temperature of frozen schedule
	 * @param numIterations : The max iterations of Local Search
	 * @param fogDevices : List of Fog Devices
	 * @param cloudletList : List of Tasks (need to be assigned on a specific Fog Device)
	 * 
	 * @return maxIndividual : The best individual in Local Search process
	 */
	public MyIndividual simulatedAnnealing(MyIndividual individual, double INITIAL_TEMPERATURE,
			double ENDING_TEMPERATURE, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// Set up parameters of Simulated Annealing
		double currentTemperature = INITIAL_TEMPERATURE;
		double descendingSpeed = 1 - (Math.log(INITIAL_TEMPERATURE) - Math.log(ENDING_TEMPERATURE)) / numIterations;
		
		// Get the length of ADN and its domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
		// The current fitness of the individual
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		
		// The best fitness (current best record - Kỷ lục tốt nhất hiện tại)
		double bestFitness = currentFitness;
		MyIndividual bestIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Initialize the neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Contains good local moves (good neighbors) at each iteration
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// The best neighbor (new records) surpassing the current best fitness
		List<MyAbstractLocalMove> localMovesOfManyNewRecords = new ArrayList<MyAbstractLocalMove>();
		
		ArrayList<Double> simulatedAnnealingFitness = new ArrayList<Double>();
		
		double random;
		int iterationIndex = 0;
		while (iterationIndex < numIterations) {
			
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();
			
			// Clear all previous records
			localMovesOfManyNewRecords.clear();
			
			// Choose the better neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					
					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					random = (Math.random());
					
					// Checks condition of joining in list of new records (kỷ lục mới)
					if (neighborFitness > bestFitness) {
						localMovesOfManyNewRecords.add(new MyAssignLocalMove(geneIndex, geneValue));
					} 
					// Simulated Annealing Mechanism - Accept the sacrifices
					else if (random <= Math.exp(-(currentFitness - neighborFitness) / currentTemperature)) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Moves the new neighbor back to current state
					oldGeneValue = individual.getGene(geneIndex);
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			int geneIndex1, geneIndex2;
			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					random = (Math.random());
					
					// Checks condition of joining in list of new records (kỷ lục mới)
					if (neighborFitness > bestFitness) {
						localMovesOfManyNewRecords.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}
					// Simulated Annealing Mechanism - Accept the sacrifices
					else if (random <= Math.exp(-(currentFitness - neighborFitness) / currentTemperature)) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Moves the new neighbor back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			// If exists any neighbor surpassing current best fitness, choose it immediately
			if (localMovesOfManyNewRecords.size() > 0) {
				
				int newRecordIndex = MyService.rand(0, localMovesOfManyNewRecords.size() - 1);
				MyAbstractLocalMove newRecord = localMovesOfManyNewRecords.get(newRecordIndex);
				
				if (newRecord instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) newRecord).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) newRecord).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (newRecord instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) newRecord).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) newRecord).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}
			}

			// Chooses randomly a better neighbor from listLocalMoves
			else if (!listLocalMoves.isEmpty()) {

				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);

				if (abstractLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

			} else { // Restarts if list of local moves is empty
				System.out.println("\nRestart ...\n");
				break;
			}

			// Descend current temperature according to gradually frozen schedule
//			currentTemperature *= descendingSpeed;

			// Get fitness of the individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Updates the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
				MyService.copyAttributes(bestIndividual, individual);
			}
				
			simulatedAnnealingFitness.add(iterationIndex, currentFitness);
			
			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}
			
			individual.printGene();
			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Current Temperature : " + currentTemperature);

			iterationIndex++;
		}
		
		MyStatistic.simulatedAnnealingFitness = simulatedAnnealingFitness;
		
		return bestIndividual;
	}

	
	
	
	

	public MyIndividual degradedCeiling(MyIndividual individual, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		
		// Set up the parameters of Degrated Ceiling
		double currentCeiling = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		double descendingSpeed = (0.9 - currentCeiling) / numIterations;
		
		// Get length of ADN and its domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
		// Calculate the current fitness of the individual
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		
		// Save the best individual and its best fitness
		double bestFitness = currentFitness;
		MyIndividual bestIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Create new neighbor of the individual 
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// Save the new record (if exists) and its fitness at each iteration
		ArrayList<MyAbstractLocalMove> newLocalMoveRecords = new ArrayList<MyAbstractLocalMove>();
//		double newLocalMoveRecordFitness = Double.MIN_VALUE;
		
		System.out.println("\nInitial Fitness : " + currentFitness);
		
		ArrayList<Double> degratedCeilingFitness = new ArrayList<Double>();
		
		int iterationIndex = 0;
		while (iterationIndex < numIterations) {
			
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();
			
			// Clear all records of previous iteration
			newLocalMoveRecords.clear();

			// Choose the better assigning neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					
					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					// Checks condition of taking part in the better neighbors
					if (neighborFitness > bestFitness) {
						newLocalMoveRecords.add(new MyAssignLocalMove(geneIndex, geneValue));

					} else if (neighborFitness >= currentCeiling) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Moves back to current state
					oldGeneValue = individual.getGene(geneIndex);
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			// Choose the better swapping neighbors
			int geneIndex1, geneIndex2;
			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					if (neighborFitness > bestFitness) {
						newLocalMoveRecords.add(new MySwapLocalMove(geneIndex1, geneIndex2));

					} else if (neighborFitness >= currentCeiling) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Moves back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			// If exists new record, choose the record immediately
			if (newLocalMoveRecords.size() > 0) {
				int recordIndex = MyService.rand(0, newLocalMoveRecords.size() - 1);
				MyAbstractLocalMove abstractLocalMoveRecord = newLocalMoveRecords.get(recordIndex);

				if (abstractLocalMoveRecord instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMoveRecord).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMoveRecord).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMoveRecord instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) abstractLocalMoveRecord).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) abstractLocalMoveRecord).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

			}

			// Chooses randomly a better neighbor from listLocalMoves
			else if (!listLocalMoves.isEmpty()) {

				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);

				if (abstractLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

			} else { // Restarts if list of local moves is empty
				System.out.println("\nRestart : Stop Degrated Ceiling\n");
				break;
			}

			currentCeiling += descendingSpeed;
			individual.printGene();

			// Calculates fitness again
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Updates the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
				MyService.copyAttributes(bestIndividual, individual);
			}
			
			degratedCeilingFitness.add(iterationIndex, currentFitness);

			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}

			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Current Ceiling : " + currentCeiling);
			
			iterationIndex++;

		}
		
		MyStatistic.degratedCeilingFitness = degratedCeilingFitness;
		
		return bestIndividual;
	}


	
	public MyIndividual tabuSearch(MyIndividual individual, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList, int maxStable, int maxInteration, int maxTime, int tabuLength) {

		// Get the length of ADN and its domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();

		// Initializes the Tabu list for assigning
		int[][] tabuAssignList = new int[chromosomeLength][maxGeneValue + 1];
		for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
			for (int geneValue = 0; geneValue <= individual.getMaxValue(); geneValue++) {
				tabuAssignList[geneIndex][geneValue] = -1;
			}
		}

		// Initializes the Tabu list for swapping
		int[][] tabuSwapList = new int[chromosomeLength][chromosomeLength];
		for (int geneIndex1 = 0; geneIndex1 < chromosomeLength; geneIndex1++) {
			for (int geneIndex2 = 0; geneIndex2 < chromosomeLength; geneIndex2++) {
				tabuSwapList[geneIndex1][geneIndex2] = -1;
			}
		}
		
		// Calculate the fitness of the individual
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		
		// Save the best individual and its best fitness
		double bestFitness = currentFitness;
		MyIndividual bestIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Create the new neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Contains good local moves for upcoming states
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		System.out.println("Initial Fitness : " + currentFitness);

		ArrayList<Double> tabuFitness = new ArrayList<Double>();
		
		// Do iterations for Tabu Search
		maxTime = maxTime * 1000;
		int notImprovingCount = 0;
		int iterationIndex = 0;
		double start = System.currentTimeMillis();
		
		while (iterationIndex < maxInteration) {

			// Clear all neighbors of previous iteration
			listLocalMoves.clear();

			int oldGeneValue;
			double neighborFitness;

			// Choose better assigning neighbors to list of local moves for coming iteration
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {

					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					if (tabuAssignList[geneIndex][geneValue] <= iterationIndex || neighborFitness > currentFitness) {
						if (neighborFitness >= currentFitness) {
							listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
						}
					}

					// Move back the current state
					oldGeneValue = individual.getGene(geneIndex);
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			// Choose better swapping neighbors to list of local moves for coming iteration
			for (int geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (int geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					if (tabuSwapList[geneIndex1][geneIndex2] <= iterationIndex || neighborFitness > currentFitness) {
						if (neighborFitness >= currentFitness) {
							listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
						}
					}

					// Moves back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			if (listLocalMoves.size() > 0) {

				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);

				if (abstractLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					tabuAssignList[geneIndex][geneValue] = iterationIndex + tabuLength;
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					int geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					int geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					tabuSwapList[geneIndex1][geneIndex2] = iterationIndex + tabuLength;
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

				// Calculates fitness again
				currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
				
				// Update best fitness (if exists new record)
				if (currentFitness > bestFitness) {
					notImprovingCount = 0;
					bestFitness = currentFitness;
					MyService.copyAttributes(bestIndividual, individual);
				} 
				// Update notImprovingCount (if not any improve)
				else if (currentFitness <= bestFitness) {
					notImprovingCount++;
					if (notImprovingCount > maxStable) {
						System.out.println("Tabu restart: Stop Tabu Search");
						break;
					}
				}
				
			// If list of local move is empty
			} else {
				System.out.println("Tabu restart: Stop Tabu Search");
				break;
			}

			
			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}

			System.out.println("\nIterations : " + iterationIndex);
			System.out.println("Fitness value: " + individual.getFitness());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Not Changed Succession : " + notImprovingCount);
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			
			tabuFitness.add(iterationIndex, currentFitness);
			
			iterationIndex++;
		}

		MyStatistic.tabuFitness = tabuFitness;
		return bestIndividual;

	}

	/**
	 * Calculates the lower boundary of Time and Cost
	 */
	public void calcMinTimeCost(List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		this.minTime = MyService.calcMinTime(fogDevices, cloudletList);
		this.minCost = MyService.calcMinCost(fogDevices, cloudletList);

	}

	public double getMinTime() {
		return minTime;
	}

	public void setMinTime(double minTime) {
		this.minTime = minTime;
	}

	public double getMinCost() {
		return minCost;
	}

	public void setMinCost(double minCost) {
		this.minCost = minCost;
	}

}
