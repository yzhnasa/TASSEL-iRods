package net.maizegenetics.gwas.jointlinkage;

public class FileNames {
	public static String snpFilename = "C:/Projects/NAM/data/ImputedMarkerGenotypes_flowering_traits_092909.txt";
	public static String modelFilename = "C:/users/peter/temp/modelout.txt";
	public static String stepFilename = "C:/users/peter/temp/stepout.txt";
	public static String scanFilename = "C:/users/peter/temp/scanout.txt";
	public static String permutationFilename = "C:/users/peter/temp/permutationout.txt";

	public static void setSnpFilename(String snpFileName) { snpFilename = snpFileName; }
	public static void setModelFilename(String modelFileName) { modelFilename = modelFileName; }
	public static void setStepFilename(String stepFileName) { stepFilename = stepFileName; }
	public static void setScanFilename(String scanFileName) { scanFilename = scanFileName; }
	public static void setPermutationFilename(String permutationFileName) { permutationFilename = permutationFileName; }
	
}
