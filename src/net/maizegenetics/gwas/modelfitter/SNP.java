package net.maizegenetics.gwas.modelfitter;

import java.util.ArrayList;

import net.maizegenetics.pal.position.Chromosome;

public class SNP {
	public String name;
	public Chromosome locus;
	public int position;
	ArrayList<Object> alleles;
	int index;
	
	public SNP(String name, Chromosome locus, int positionInLocus, int index) {
		this.name = name;
		this.locus = locus;
		position = positionInLocus;
		this.index = index;
	}
	
	public SNP() {
		
	}

	@Override
	public String toString() {
		return name;
	}
	
	
}
