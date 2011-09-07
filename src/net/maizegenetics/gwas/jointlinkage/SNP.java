package net.maizegenetics.gwas.jointlinkage;

public class SNP {
	double[] score;
	int chr;
	double pos;
	String name;
	int index;
	
	public SNP(int chr, double pos, double[] score, int index) {
		this.chr = chr;
		this.pos = pos;
		this.score = score;
		this.index = index;
	}
	
	public SNP(int chr, double pos, double[] score, String name, int index) {
		this.chr = chr;
		this.pos = pos;
		this.score = score;
		this.name = name;
		this.index = index;
	}
}
