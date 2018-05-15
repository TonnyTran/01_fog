package org.fog.entities;

import java.util.List;

import org.cloudbus.cloudsim.core.SimEvent;
import org.cloudbus.cloudsim.power.PowerDatacenterBroker;
import org.fog.scheduling.MySchedulingAlgorithm;
import org.fog.scheduling.SchedulingAlgorithm;

public class MyFogBroker extends PowerDatacenterBroker{

	private List<FogDevice> fogDevices;
	
	public MyFogBroker(String name) throws Exception {
		super(name);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void startEntity() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void processEvent(SimEvent ev) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void shutdownEntity() {
		// TODO Auto-generated method stub
		
	}

	public List<FogDevice> getFogDevices() {
		return fogDevices;
	}

	public void setFogDevices(List<FogDevice> fogDevices) {
		this.fogDevices = fogDevices;
	}
	
	public void assignCloudlet(String schedulingStrategy) {
		switch (schedulingStrategy) {
			case SchedulingAlgorithm.GA:
				SchedulingAlgorithm.runGeneticAlgorithm(fogDevices, cloudletList);
				break;
			case SchedulingAlgorithm.LOCAL_SEARCH:
				SchedulingAlgorithm.runLocalSearchAlgorithm(fogDevices, cloudletList);
				break;
			case SchedulingAlgorithm.TABU_SEARCH:
				SchedulingAlgorithm.runTabuSearchAlgorithm(fogDevices, cloudletList);
				break;
		}
	}
	
	public void myAssignCloudlet(String schedulingStrategy) {
		switch (schedulingStrategy) {
			case MySchedulingAlgorithm.GA:
				MySchedulingAlgorithm.runGeneticAlgorithm(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.HILL_CLIMBING:
				MySchedulingAlgorithm.runHillClimbing(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.TABU_SEARCH:
				MySchedulingAlgorithm.runTabuSearchAlgorithm(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.SACRIFICED_HILL_CLIMBING:
				MySchedulingAlgorithm.runSacrificedHillClimbing(fogDevices, cloudletList);
				break;
			
			case MySchedulingAlgorithm.SIMULATED_ANNEALING_BEST_NEW_RECORD:
				MySchedulingAlgorithm.runSimulatedAnnealing_BestNewRecord(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.SIMULATED_ANNEALING_MANY_NEW_RECORDS:
				MySchedulingAlgorithm.runSimulatedAnnealing_ManyNewRecords(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.DEGRATED_CEILING_BEST_NEW_RECORD:
				MySchedulingAlgorithm.runDegratedCeiling_BestNewRecord(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.DEGRATED_CEILING_MANY_NEW_RECORDS:
				MySchedulingAlgorithm.runDegratedCeiling_ManyNewRecords(fogDevices, cloudletList);
				break;
				
			case MySchedulingAlgorithm.COMBINATION_GA_LOCALSEARCH:
				MySchedulingAlgorithm.runGA_LocalSearch(fogDevices, cloudletList);
				break;
		}
	}

}
