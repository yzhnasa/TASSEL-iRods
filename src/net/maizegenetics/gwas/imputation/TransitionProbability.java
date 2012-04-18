package net.maizegenetics.gwas.imputation;

public class TransitionProbability {
	protected double[][] probabilityOfATransition;
	protected int[] positions;
	
	public void setTransitionProbability(double[][] probabilityMatrix) {
		probabilityOfATransition = probabilityMatrix;
	}
	
	public double getTransitionProbability(int state1, int state2) {
		return probabilityOfATransition[state1][state2];
	}
	
	public double getLnTransitionProbability(int state1, int state2) {
		return Math.log(getTransitionProbability(state1, state2));
	}
	
	public double getTransitionProbability(int state1, int state2, int node) {
		return probabilityOfATransition[state1][state2];
	}
	
	public double getLnTransitionProbability(int state1, int state2, int node) {
		return Math.log(getTransitionProbability(state1, state2, node));
	}
	
	public int getNumberOfStates() {
		return probabilityOfATransition.length;
	}
	
	public void setPositions(int[] positions) {
		this.positions = positions;
	}
}
