package net.maizegenetics.gwas.NAMgwas;

import java.io.File;
import java.io.FilenameFilter;
import java.util.regex.Pattern;

import net.maizegenetics.gwas.jointlinkage.StepwiseModelFitterForNAMJointLinkage;

public class Analysis {
	public static final String ANALYSIS_STEPWISE_NOMISSING = "stepwise";
	public static final String ANALYSIS_BOOTSTRAP = "bootstrap";
	public static final String ANALYSIS_HAPLOTYPE = "haplotype";
	public static final String ANALYSIS_SINGLESNP = "singlesnp";
	public static final String ANALYSIS_PERMUTATIONS = "permutations";

	public static void main(String[] args) {
		if (args == null || args.length == 0) {
//			imputeParentsFromGBS();
//			callHets();
//			NamGBSData.compareToNAMarray("C:/Users/Peter/temp/gbs/cov20_chr10_Z001.hets.txt", 10);
			
//			NamGBSData.createImputedSNPDatasets();
//			NamGBSData.concatenateSomeFiles();
			StepwiseModelFitterForNAMJointLinkage fitter = new StepwiseModelFitterForNAMJointLinkage();
		}
		fitModels(args);
	}
	
	public static void imputeParentsFromGBS() {
		final File dir = new File("C:/Projects/NAM/namgbs/alignment071811cov20LDE");
//		final File dir = new File("C:/Projects/NAM/namgbs/allfusion_110401.LD/allLD");
//		final File dir = new File("C:/Projects/NAM/namgbs/Archive");
		final File[] gbsfiles = dir.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File dr, String fname) {
				if (dr.equals(dir) && fname.endsWith(".hmp.txt")) return true;
				return false;
			}
		});
		
		for (int ndx = 0; ndx < 10; ndx++) {
			System.out.println("Imputing parents for " + gbsfiles[ndx].getPath());
			NamGBSData namgbs = new NamGBSData(gbsfiles[ndx].getPath());
		}

		//		NamGBSData namgbs = new NamGBSData(10, "C:/Projects/NAM/namgbs/maize071811.cov10.fT1hLDE1.mgNoHet.imp.c10.hmp.txt/maize071811.cov10.fT1hLDE1.mgNoHet.imp.c10.hmp.txt");
	}
	
	public static void callHets() {
		final File dir = new File("C:/Projects/NAM/namgbs/datasets/cov20");
		final File[] gbsfiles = dir.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File dr, String fname) {
				if (dr.equals(dir) && fname.startsWith("cov20_chr")) return true;
				return false;
			}
		});
		
		//one of the files might be cov20_chr1_Z001.txt
		for (File file:gbsfiles) {
			String infileName = file.getPath();
			String basename = infileName.substring(0, infileName.lastIndexOf('.'));
			String outfileName = basename + ".hets.txt";
			NamGBSData.callHetsFromABHFile(infileName, outfileName);
			System.out.println("calling hets for " + infileName);
		}
		
		System.out.println("Finished calling hets.");
	}
	
	public static void fitModels(String[] args) {
		String filename = null;
		Pattern equal = Pattern.compile("=");

		if (args.length == 0) {
			return;
		}
		else {
			for (String arg : args) {
				String[] parsed = equal.split(arg);
				if (parsed[0].startsWith("file")) filename = parsed[1];
			}

			if (filename == null) {
				System.err.println("No file name specified. Exiting program.");
				System.exit(-1);
			}
			else {
				FileNames files = new FileNames(filename);
				for (int chr = 1; chr <=10; chr++) {
					if (files.analysis.equalsIgnoreCase("project")) {
						if (files.projectedFile[chr - 1] != null) {
							ModelFitter fitter = new ModelFitter(chr, files);
							fitter.init();
							fitter.createProjectionData();
						}
					} else if (files.chrmodel[chr - 1] != null) {
						if (files.threaded && !files.analysis.equalsIgnoreCase(ANALYSIS_BOOTSTRAP)) {
							createModelFitterInstance(chr, files).start();
						} else {
							//debug
							System.out.println("testing chromosome " + chr);
							
							ModelFitter theFitter = createModelFitterInstance(chr, files);
							theFitter.init();
							theFitter.testSnps();
						}
					}
				}
			}
		}

	}

	public static ModelFitter createModelFitterInstance(int chr, FileNames files) {
		if (files.analysis.length() == 0) return  null;
		String name = files.analysis.substring(0, 4);
		if (name.equalsIgnoreCase("step")) {
			files.analysis = ANALYSIS_STEPWISE_NOMISSING;
			return new ModelFitterForwardRegression(chr, files);
		}
		if (name.equalsIgnoreCase("hapl")) {
			//not implemented
			System.out.println("Haplotype analysis is not implemented.");
			System.exit(0);
//			files.analysis = ANALYSIS_HAPLOTYPE;
		}
		if (name.equalsIgnoreCase("sing")) {
			files.analysis = ANALYSIS_SINGLESNP;
			return new ModelFitterNoMissingData(chr, files);
		}
		if (name.equalsIgnoreCase("perm")) {
//			files.analysis = ANALYSIS_PERMUTATIONS;
			System.out.println("Permutation testing is not implemented.");
			System.exit(0);
		}
		if (name.equalsIgnoreCase("boot")) {
			return new ModelFitterBootstrapForward(chr, files);
		}
		return null;
	}

}
