package net.maizegenetics.gwas.LeafTraits;

public class SnpDataImputed extends SnpData {
	
	public SnpDataImputed(int chromosome) {
		super(chromosome);
	}
	
	public SnpDataImputed(int chromosome, FileNames files) {
		super(chromosome, files);
	}
	
	public double[] getGenotype() {
		int count = 0;
		double[] geno = new double[26];
		for (int i = 12; i < 38; i++) {
			geno[count++] = Double.parseDouble(parsedLine[i]);
		}
		return geno;
	}

	public String createInputFileName() {
		if (files != null) return super.createInputFileName();
		StringBuilder sb = new StringBuilder();
		
		//fastphase
		sb.append("C:/Projects/NAM/association/Fastphase/fastphase_chr");
		sb.append(chromosome);
		sb.append(".txt");
//		sb.append("C:/Projects/NAM/association/Fastphase/fastphase_chr3_angletest.txt");
		
		//average
//		sb.append("C:/Projects/NAM/association/average/chr");
//		sb.append(chromosome);
//		sb.append("average.txt");
		
		//minor
//		sb.append("C:/Projects/NAM/association/minor freq/chr");
//		sb.append(chromosome);
//		sb.append("snp_minor.txt");
		
		
		return sb.toString();
	}

	public boolean IsB73Alt() {
		if (parsedLine[11].charAt(0) == '1') return true;
		return false;
	}

}
