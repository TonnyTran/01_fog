package org.fog.scheduling.myLocalSearchAlgorithm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.cloudbus.cloudsim.Cloudlet;
import org.fog.entities.FogDevice;
import org.fog.scheduling.myGAEntities.MyIndividual;
import org.fog.scheduling.myGAEntities.MyService;

public class MyLocalSearchAlgorithm {
	private double minTime;
	private double minCost;

	public MyLocalSearchAlgorithm() {
	}

	public MyIndividual sacrificedHillClimbing(MyIndividual individual, double INITIAL_SACRIFICE,
			double DESCENDING_SPEED, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {
		
		// Sets up the initial sacrifice (it descends over time)
		double currentSacrifice = INITIAL_SACRIFICE; // 0.005
		
		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// Initialize the list contain the solutions for restarting
		int numRestarts = 10;
		MyIndividual[] restartingList = new MyIndividual[numRestarts];
		for (int index = 0; index < numRestarts; index++)
			restartingList[index] = individual;
	
		// Calculate and save the current fitness of individual 
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);;
		
		// Saves the best fitness over time
		double bestFitness = currentFitness;
		
		// Neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Get the length of gene and its value domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
		// Do iterations of local search and set up the initial time
		double start = System.currentTimeMillis();
		int restartIteration = 0;
		int iterationIndex = 0;
		do {
			
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Descending the degree of sacrifice over iterations
			if (iterationIndex % 300 == 0)
			{
				currentSacrifice = currentSacrifice * DESCENDING_SPEED;
			}
				

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
					
					// 
					if (neighborFitness > (bestFitness - currentSacrifice)) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Move back to state of current individual
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());
			
//			// Update the restarting list
//			if (!listLocalMoves.isEmpty()) {
//
//				int localMoveIndex = MyService.rand(0, listLocalMoves.size() - 1);
//				MyAbstractLocalMove abstractLocalMove = listLocalMoves.get(localMoveIndex);
//
//				if (abstractLocalMove instanceof MyAssignLocalMove) {
//					int geneIndex = ((MyAssignLocalMove) abstractLocalMove).getGeneIndex();
//					int geneValue = ((MyAssignLocalMove) abstractLocalMove).getGeneValue();
//					MyIndividual myRestart = MyService.clonedIndividual(individual);
//					myRestart.setGene(geneIndex, geneValue);
//					restartingList[iterationIndex % numRestarts] = myRestart;
//
//				} else if (abstractLocalMove instanceof MySwapLocalMove) {
//					geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
//					geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
//					MyIndividual myRestart = MyService.clonedIndividual(individual);
//					MyService.swapIndividual(myRestart, geneIndex1, geneIndex2);
//					restartingList[iterationIndex % numRestarts] = myRestart;
//				}
//
//			}
//			
			
			
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
//				individual = restartingList[MyService.rand(0, numRestarts-1)];
				individual = new MyIndividual(chromosomeLength, maxGeneValue, true);
			}

			// Get fitness of the new individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Update the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
			}
				
				

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
		return individual;
	}


	public MyIndividual hillClimbing(MyIndividual individual, int numIterations, 
			List<FogDevice> fogDevices, List<? extends Cloudlet> cloudletList) {
		
		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		// Calculate and save the current fitness of individual 
		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);;
		
		// Saves the best fitness over time
		double bestFitness = currentFitness;
		
		// Neighbor of the individual
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		
		// Get the length of gene and its value domain
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		
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
					
					// 
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
//				System.out.println("\nRestart ...\n");
//				individual = new MyIndividual(chromosomeLength, maxGeneValue, true);
				break;
			}

			// Get fitness of the new individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Update the best fitness
			if (currentFitness > bestFitness) {
				bestFitness = currentFitness;
			}

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
		return individual;
	}

	
	public MyIndividual simulatedAnnealing(MyIndividual individual, double INITIAL_TEMPERATURE,
			double ENDING_TEMPERATURE, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// Contains good local moves for upcoming states
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		
		MyAbstractLocalMove bestLocalMove = null;
		

		double currentTemperature = INITIAL_TEMPERATURE;
		double descendingSpeed = 1 - (Math.log(INITIAL_TEMPERATURE) - Math.log(ENDING_TEMPERATURE)) / numIterations;
		double sacrifyingDegree;

		double currentFitness;
		double bestFitness = Double.MIN_VALUE;
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		int iterationIndex = 0;
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		do {
			iterationIndex++;
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();

			// Get fitness of the individual
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Updates the best fitness
			if (currentFitness > bestFitness)
				bestFitness = currentFitness;

			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}
			
			double bestLocalFitness = Double.MIN_VALUE;
			double random;
			
			// Choose the better neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					oldGeneValue = individual.getGene(geneIndex);
					
					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					random = (Math.random());
					
					// Checks condition of taking part in the better neighbors
					if (neighborFitness > bestFitness) {
						bestLocalFitness = neighborFitness;
						bestLocalMove = new MyAssignLocalMove(geneIndex, geneValue);
					}
					else if ((random) <= Math.exp(-(currentFitness - neighborFitness) / currentTemperature)) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
//						if (iterationIndex % 50 == 0)
//							System.out.println(" Math.random() : " + random + " Do Hy sinh : " +  Math.exp(-(currentFitness - neighborFitness) / currentTemperature));
					}

					// Moves back to current state
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			int geneIndex1, geneIndex2;
			int oldGeneValue1, oldGeneValue2;
			
			
			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					
					random = (Math.random());
					
					if (neighborFitness > bestFitness) {

						bestLocalFitness = neighborFitness;
						bestLocalMove = new MySwapLocalMove(geneIndex1, geneIndex2);
						

					} 			
					else if ((random) <= Math.exp(-(currentFitness - neighborFitness) / currentTemperature)) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
						
//						if (iterationIndex % 50 == 0)
//							System.out.println(" Math.random() : " + random + " Do Hy sinh : " + Math.exp(-(currentFitness - neighborFitness) / currentTemperature));
					}

					// Moves back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			if (bestLocalFitness > bestFitness) {
				bestFitness = bestLocalFitness;

				if (bestLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) bestLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) bestLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);
				} else if (bestLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) bestLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) bestLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}
			}

			// Chooses randomly a better neighbor from listLocalMoves
			else 
				
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
//				System.out.println("\nRestart ...\n");
//				individual = new MyIndividual(individual.getChromosomeLength(), individual.getMaxValue(), true);
			}

			
			
			
			
//			if (iterationIndex % 10 == 0) {
//				descendingSpeed = 1 - (Math.log(INITIAL_TEMPERATURE) - Math.log(ENDING_TEMPERATURE)) / (numIterations-iterationIndex);
//				currentTemperature *= descendingSpeed;
//			}
			
//			currentTemperature *= descendingSpeed;
			
			
			
			
			
			
			
			individual.printGene();
			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Current Temperature : " + currentTemperature);

		} while (iterationIndex < numIterations);
		return individual;
	}

	
	
	public MyIndividual simulatedAnnualing1(MyIndividual individual, double INITIAL_TEMPERATURE,
			double ENDING_TEMPERATURE, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// Contains good local moves for upcoming states
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		double bestFitness = Double.MIN_VALUE;
		MyAbstractLocalMove bestLocalMove = null;
		double bestLocalFitness = Double.MIN_VALUE;

		double currentTemperature = INITIAL_TEMPERATURE;
		double descendingSpeed = 1 - (Math.log(INITIAL_TEMPERATURE) - Math.log(ENDING_TEMPERATURE)) / numIterations;
		double sacrifyingDegree;

		double fitness;
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		int iterationIndex = 0;
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();
		do {
			iterationIndex++;
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();

			// Get fitness of the individual
			fitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			// Updates the best fitness
			if (fitness > bestFitness)
				bestFitness = fitness;

			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}
			// Choose the better neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					oldGeneValue = individual.getGene(geneIndex);

					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					// Checks condition of taking part in the better neighbors
					if (neighborFitness >= bestFitness) {

						bestLocalFitness = neighborFitness;
						bestLocalMove = new MyAssignLocalMove(geneIndex, geneValue);

					} else if (Math.random() <= Math.exp((neighborFitness - bestFitness) / currentTemperature)) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Moves back to current state
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			int geneIndex1, geneIndex2;
			int oldGeneValue1, oldGeneValue2;

			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					if (neighborFitness > bestFitness) {

						bestLocalFitness = neighborFitness;
						bestLocalMove = new MySwapLocalMove(geneIndex1, geneIndex2);

					} else if (Math.random() <= Math.exp((neighborFitness - bestFitness) / currentTemperature)) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Moves back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			if (bestLocalFitness > bestFitness) {
				if (bestLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) bestLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) bestLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);
				} else if (bestLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) bestLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) bestLocalMove).getGeneIndex2();
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
				individual = new MyIndividual(individual.getChromosomeLength(), individual.getMaxValue(), true);
			}

			currentTemperature *= descendingSpeed;
			individual.printGene();
			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Current Temperature : " + currentTemperature);

		} while (iterationIndex < numIterations);
		return individual;
	}
	
	
	
	public MyIndividual degradedCeiling(MyIndividual individual, int numIterations, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList) {

		// Contains local moves to good neighbors
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();
		double bestFitness = Double.MIN_VALUE;
		MyAbstractLocalMove bestLocalMove = null;
		double bestLocalFitness = Double.MIN_VALUE;

		double currentCeiling = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		double descendingSpeed = (0.9 - currentCeiling) / numIterations;
		double sacrifyingDegree;

		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);
		int iterationIndex = 0;
		int chromosomeLength = individual.getChromosomeLength();
		int maxGeneValue = individual.getMaxValue();

		double currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		System.out.println("\nInitial Fitness : " + currentFitness);

		do {
			iterationIndex++;
			System.out.println("\n--------------------------------------");
			System.out.println("Iteration : " + iterationIndex + " : ");

			// Clear all local moves at previous step
			listLocalMoves.clear();

			// Updates the best fitness
			if (currentFitness > bestFitness)
				bestFitness = currentFitness;

			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}

			// Choose the better neighbors
			int oldGeneValue;
			double neighborFitness;
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {
					oldGeneValue = individual.getGene(geneIndex);

					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					// Checks condition of taking part in the better neighbors
					if (neighborFitness >= bestFitness) {

						bestLocalFitness = neighborFitness;
						bestLocalMove = new MyAssignLocalMove(geneIndex, geneValue);

					} else if (neighborFitness >= currentCeiling) {
						listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
					}

					// Moves back to current state
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			int geneIndex1, geneIndex2;
			int oldGeneValue1, oldGeneValue2;

			for (geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);
					if (neighborFitness > bestFitness) {

						bestLocalFitness = neighborFitness;
						bestLocalMove = new MySwapLocalMove(geneIndex1, geneIndex2);

					} else if (neighborFitness >= currentCeiling) {
						listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
					}

					// Moves back to current state
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
				}
			}

			System.out.println("Number of changelist: " + listLocalMoves.size());

			if (bestLocalFitness > bestFitness) {
				if (bestLocalMove instanceof MyAssignLocalMove) {
					int geneIndex = ((MyAssignLocalMove) bestLocalMove).getGeneIndex();
					int geneValue = ((MyAssignLocalMove) bestLocalMove).getGeneValue();
					individual.setGene(geneIndex, geneValue);
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);
				} else if (bestLocalMove instanceof MySwapLocalMove) {
					geneIndex1 = ((MySwapLocalMove) bestLocalMove).getGeneIndex1();
					geneIndex2 = ((MySwapLocalMove) bestLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}
			}

			// Chooses randomly a better neighbor from listLocalMoves
			else 
				
				
				
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
				individual = new MyIndividual(individual.getChromosomeLength(), individual.getMaxValue(), true);
			}

			currentCeiling += descendingSpeed;
			individual.printGene();

			// Calculates fitness again
			currentFitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

			System.out.println("\nFitness value: " + individual.getFitness());
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Current Ceiling : " + currentCeiling);

		} while (iterationIndex < numIterations);
		return individual;
	}

	
	
	public MyIndividual tabuSearch(MyIndividual individual, List<FogDevice> fogDevices,
			List<? extends Cloudlet> cloudletList, int maxStable, int maxInteration, int maxTime, int tabuLength) {

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

		MyIndividual bestSolution = (MyIndividual) MyService.clonedIndividual(individual);
		double bestFitness = Double.MIN_VALUE;
		double fitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);
		System.out.println("Initial Fitness : " + fitness);

		// Contains good local moves for upcoming states
		List<MyAbstractLocalMove> listLocalMoves = new ArrayList<MyAbstractLocalMove>();

		double start = System.currentTimeMillis();
		int count = 0;
		maxTime = maxTime * 1000;
		Random R = new Random();
		int notChangedSuccession = 0;
		MyIndividual neighborIndividual = (MyIndividual) MyService.clonedIndividual(individual);

		while ((System.currentTimeMillis() - start) < maxTime && count < maxInteration) {

			listLocalMoves.clear();

			// Create a new neighbor
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				neighborIndividual.setGene(geneIndex, individual.getGene(geneIndex));
			}

			int oldGeneValue;
			double neighborFitness;
			double bestLocalFitness = Double.MIN_VALUE;

			// consider which gene changed makes individual better
			for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
				for (int geneValue = 0; geneValue <= maxGeneValue; geneValue++) {

					oldGeneValue = individual.getGene(geneIndex);

					// Moves to new neighbor and calculates its fitness
					neighborIndividual.setGene(geneIndex, geneValue);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					if (tabuAssignList[geneIndex][geneValue] <= count || neighborFitness >= bestLocalFitness) {
						if (neighborFitness > bestLocalFitness) {
							listLocalMoves.clear();
							listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
							bestLocalFitness = neighborFitness;
						} else if (neighborFitness == bestFitness) {
							listLocalMoves.add(new MyAssignLocalMove(geneIndex, geneValue));
						}
					}

					// Move back the current state
					neighborIndividual.setGene(geneIndex, oldGeneValue);
				}
			}

			for (int geneIndex1 = 0; geneIndex1 < chromosomeLength - 1; geneIndex1++) {
				for (int geneIndex2 = geneIndex1 + 1; geneIndex2 < chromosomeLength; geneIndex2++) {
					// Moves to new neighbor and calculates its fitness
					MyService.swapIndividual(neighborIndividual, geneIndex1, geneIndex2);
					neighborFitness = MyService.calcFitness(neighborIndividual, this.minTime, this.minCost, fogDevices,
							cloudletList);

					if (tabuSwapList[geneIndex1][geneIndex2] <= count || neighborFitness >= bestLocalFitness) {
						if (neighborFitness > bestLocalFitness) {
							listLocalMoves.clear();
							listLocalMoves.add(new MySwapLocalMove(geneIndex1, geneIndex2));
							bestLocalFitness = neighborFitness;
						} else if (neighborFitness == bestFitness) {
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
					tabuAssignList[geneIndex][geneValue] = count + tabuLength;
					System.out.println("Local Move : Assign Individual[" + geneIndex + "] = " + geneValue);

				} else if (abstractLocalMove instanceof MySwapLocalMove) {
					int geneIndex1 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex1();
					int geneIndex2 = ((MySwapLocalMove) abstractLocalMove).getGeneIndex2();
					MyService.swapIndividual(individual, geneIndex1, geneIndex2);
					tabuSwapList[geneIndex1][geneIndex2] = count + tabuLength;
					System.out.println("Local Move : Swap geneIndex1 = " + geneIndex1 + ", geneIndex2 = " + geneIndex2);
				}

				// Calculates fitness again
				fitness = MyService.calcFitness(individual, this.minTime, this.minCost, fogDevices, cloudletList);

				if (fitness > bestFitness) {
					bestFitness = fitness;
					for (int geneIndex = 0; geneIndex < chromosomeLength; geneIndex++) {
						bestSolution.setGene(geneIndex, individual.getGene(geneIndex));
					}
				}

				if (fitness <= bestFitness) {
					notChangedSuccession++;
					if (notChangedSuccession > maxStable) {
						notChangedSuccession = 0;
						System.out.println("Tabu restart:");
						individual = new MyIndividual(chromosomeLength, maxGeneValue, true);

						// Reset the assigning tabu
						for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
							for (int geneValue = 0; geneValue < individual.getMaxValue() + 1; geneValue++) {
								tabuAssignList[geneIndex][geneValue] = -1;
							}
						}

						// Reset the swapping Tabu
						for (int geneIndex1 = 0; geneIndex1 < chromosomeLength; geneIndex1++) {
							for (int geneIndex2 = 0; geneIndex2 < chromosomeLength; geneIndex2++) {
								tabuSwapList[geneIndex1][geneIndex2] = -1;
							}
						}
					}
				}

			} else {
				notChangedSuccession = 0;
				System.out.println("Tabu restart:");
				individual = new MyIndividual(chromosomeLength, maxGeneValue, true);
				for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
					for (int geneValue = 0; geneValue < individual.getMaxValue() + 1; geneValue++) {
						tabuAssignList[geneIndex][geneValue] = -1;
					}
				}

				// Reset the swapping Tabu
				for (int geneIndex1 = 0; geneIndex1 < chromosomeLength; geneIndex1++) {
					for (int geneIndex2 = 0; geneIndex2 < chromosomeLength; geneIndex2++) {
						tabuSwapList[geneIndex1][geneIndex2] = -1;
					}
				}
			}

			count++;

			System.out.println("\nIterations : " + count);
			System.out.println("Fitness value: " + individual.getFitness());
			System.out.println("Best Fitness : " + bestFitness);
			System.out.println("Not Changed Succession : " + notChangedSuccession);
			System.out.println("Min Time: " + this.getMinTime() + "/// Makespan: " + individual.getTime());
			System.out.println("Min Cost: " + this.getMinCost() + "/// TotalCost: " + individual.getCost());

		}

		return bestSolution;

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
