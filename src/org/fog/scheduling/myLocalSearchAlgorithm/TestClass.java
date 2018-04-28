package org.fog.scheduling.myLocalSearchAlgorithm;

import org.fog.scheduling.myGAEntities.MyIndividual;
import org.fog.scheduling.myGAEntities.MyService;

public class TestClass {

	
	public static void test(int numIteration) {
		
		MyIndividual individual = new MyIndividual(1000, 1000, true);
		double start = System.currentTimeMillis();
		int countIteration = 0;
		
		while (countIteration < numIteration) {
//			MyIndividual cloneIndividual = new MyIndividual(1000, 1000, false);
//			for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
//		
//				cloneIndividual.setGene(geneIndex, cloneIndividual.getGene(geneIndex));
//				
//			}
			
			MyIndividual cloneIndividual = (MyIndividual) MyService.deepCopy(individual);
			System.out.println("iteration : " + countIteration);
			countIteration++;
		}
		
		System.out.println("Time : " + (System.currentTimeMillis() - start)/1000);
		
		
	}
	
	public static void main(String[] args) {
		
		TestClass.test(500000);
		
	}
	
	
	
	
	
}
