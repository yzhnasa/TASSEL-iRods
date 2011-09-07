package net.maizegenetics.gwas.jointlinkage;

public interface SNPdata {
	public boolean hasNext();
	public SNP next();
	public SNP getSnp(int index);
	public int getNumberOfSNPs();
	public void setStartingSNP(int start);
	public void setEndingSNP(int end);
	public int getChromosome(int index);
	public int getPosition(int index);
	public void resetSNPs();
	public double[] getPhenotype();
	public String[] getPopulations();
	public String getPhenotypeName();
}


