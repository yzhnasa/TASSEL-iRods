package net.maizegenetics.gwas.imputation;

import net.maizegenetics.pal.alignment.Alignment;

public class PhasedEmissionProbability extends EmissionProbability {
//	Alignment myGenotype;
	float[][][] parentAlleleProbability;  //1st dim is site (node), 2nd dim is haplotype, 3rd dim is nucleotide (A,C,G,T)
	float probMissing;
	float logProbMissing;
	byte[] knownHaplotype;
	float probHetObsAsHet = 0.2f;
	int myHaplotype;
	
	public PhasedEmissionProbability() {

		
	}
	
	@Override
	public double getProbObsGivenState(int state, int obs, int node) {
		//obs equals the nucleotide byte code for the observation
		//state in (0,1,2,3) indicating which parental haplotype
		//node is the position index in the alignment
		float prob = 0.0f;
		int knownhapValue = knownHaplotype[node];
		for (int i = 0; i < 4; i++) { //loop through possible state values
			float thisprob;
			if (obs > 3) { //obs is het
				if (knownhapValue == i) thisprob = .01f; //genotype state is homozygous, this is prob of observing a het when the state is homozygous
				else thisprob = probHetObsAsHet; //prob of observing a het if the genotype is a het
			} else { //obs is homozygote
				if (knownhapValue == i && knownhapValue == obs) thisprob = .998f; //genotype state is homozygous and equal to obs
				else if (knownhapValue == i) thisprob = .001f; //genotype state is homozygous and not equal to obs
				else thisprob = 1 - probHetObsAsHet; //genotype state is het
			}
			prob =+ thisprob * parentAlleleProbability[node][myHaplotype][i];
		}

		return prob;
	}

	@Override
	public double getLnProbObsGivenState(int state, int obs, int node) {
		//obs equals the nucleotide byte code for the observation
		//state in (0,1,2,3) indicating which parental haplotype
		//node is the position index in the alignment
		 
		
		return Math.log(getProbObsGivenState(state, obs, node));
	}

	public void setParentAlleleProbability(float[][][] probability) {
		
	}
	
	public void setProbabilityMissing(float pMissing) {
		probMissing = pMissing;
	}
	
	public void setKnownHaplotype(byte[] knownHaplotype) {
		this.knownHaplotype = knownHaplotype;
	}
	
	public void setProbabilityHetisObservedAsHet(float prob) {
		probHetObsAsHet = prob;
	}
	
	public void setKnownHaplotypeNumber(int whichIsKnown) {  // 0 or 1
		myHaplotype = 1 - whichIsKnown;
	}
	
}
