package net.maizegenetics.gwas.LeafTraits;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.regex.Pattern;

import javax.swing.Timer;

public class Analysis {
	public static final String ANALYSIS_STEPWISE_NOMISSING = "stepwise";
	public static final String ANALYSIS_BOOTSTRAP = "bootstrap";
	public static final String ANALYSIS_HAPLOTYPE = "haplotype";
	public static final String ANALYSIS_SINGLESNP = "singlesnp";
	public static final String ANALYSIS_PERMUTATIONS = "permutations";
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		fitModels(args);
//		Utilities.mergeFiles();
//		testChr3();
//		LD(10);
//		testlg1();
	}
	
	public static void fitAModel() {
		
		final ModelFitter fitter = new ModelFitterNoMissingData(8,0);
		fitter.start();
		
		Timer aTimer = new Timer(5000, new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
				System.out.println("progress = " + fitter.getProgress() + ", state is " + fitter.getState());
			} 
			
		});
		aTimer.start();
	}

	
	public static void fitChromosomeModels() {
		int traitno = 0;
		final ModelFitterNoMissingData[] fitters = new ModelFitterNoMissingData[10];
		for (int i = 0; i < 10; i++) {
			int chr = i + 1;
			fitters[i] = new ModelFitterNoMissingData(chr, traitno);
			fitters[i].start();
		}
		System.out.println(System.currentTimeMillis());
		Timer aTimer = new Timer(5000, new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
				for (int i = 0; i < 10; i++) {
					int chr = i + 1;
					if (fitters[i].isAlive()) System.out.println("Chromosome " + chr + " progress = " + fitters[i].getProgress());
				}
				System.out.println(System.currentTimeMillis());

			} 
			
		});
		aTimer.start();
	}
	
	public static void calculatePermutations() {
		int traitno = 0;
		final ModelFitterWithPermutations[] fitters = new ModelFitterWithPermutations[10];
		for (int i = 0; i < 10; i++) {
			int chr = i + 1;
			fitters[i] = new ModelFitterWithPermutations(chr, traitno);
			fitters[i].setNumberOfPermutations(1000);
			fitters[i].start();
		}
		System.out.println("starting permutations");
		System.out.println(System.currentTimeMillis());
		Timer aTimer = new Timer(5000, new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
				for (int i = 0; i < 10; i++) {
					int chr = i + 1;
					if (fitters[i].isAlive()) System.out.println("Chromosome " + chr + " progress = " + fitters[i].getProgress());
				}
				System.out.println(System.currentTimeMillis());
			} 
		});
		aTimer.start();
	}
	
	public static void fitChromosomeTenModel() {
		
		final ModelFitterStepwiseNoMissing fitter = new ModelFitterStepwiseNoMissing(10, 0);
		fitter.setEnterLimit(getLimit(0, 10));
		long start = System.currentTimeMillis();
		fitter.start();
		System.out.println("Finished fitting chromosome 10 using old gwas method, elapsed time = " + (System.currentTimeMillis() - start));
	}
	
	public static void fitStepwiseModels() {
		int traitnum = 1;
		for (int chr = 1; chr <=10; chr++) {
			ModelFitterStepwiseNoMissing fitter = new ModelFitterStepwiseNoMissing(chr, traitnum);
			fitter.setEnterLimit(getLimit(traitnum, chr));
			fitter.start();
		}
	}

	public static void fitStepwiseModels(String[] args) {
		String filename = null;
		Pattern equal = Pattern.compile("=");

		if (args.length == 0) {
			fitStepwiseModels();
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
					if (files.chrmodel[chr - 1] != null) {
						ModelFitterStepwiseNoMissing fitter = new ModelFitterStepwiseNoMissing(chr, files);
						fitter.start();
					}
				}
			}
		}

	}

	public static void fitModels(String[] args) {
		String filename = null;
		Pattern equal = Pattern.compile("=");

		if (args.length == 0) {
			fitStepwiseModels();
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
					if (files.chrmodel[chr - 1] != null) {
						if (files.threaded) createModelFitterInstance(chr, files).start();
						else {
							ModelFitter theFitter = createModelFitterInstance(chr, files);
							theFitter.init();
							theFitter.testSnps();
						}
					}
				}
			}
		}

	}
	
	public static void fitStepwiseCofactorModel() {
		int traitnum = 0;
		int chr = 8;
		final ModelFitterStepwiseWithCofactors fitter = new ModelFitterStepwiseWithCofactors(chr, traitnum);
		fitter.setEnterLimit(getLimit(traitnum, chr));
		fitter.start();
		
		Timer aTimer = new Timer(20000, new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
					if (fitter.isAlive()) System.out.println(fitter.showProgress());
			} 
		});
		aTimer.start();

	}
	
	public static void fitSomeStepwiseCofactorModels() {
		int traitnum = 0;
		for (int chr = 1; chr <= 10; chr++) {
			ModelFitterStepwiseWithCofactors fitter = new ModelFitterStepwiseWithCofactors(chr, traitnum);
			fitter.setEnterLimit(getLimit(traitnum, chr));
			fitter.start();
		}
	}

	public static void fitGenomewideSnpModel() {
		StepwiseGenomeWideSnpFinder snpfinder = new StepwiseGenomeWideSnpFinder();
		snpfinder.findSnps("C:/Projects/NAM/leaf traits/filenames_whole_genome.txt");
	}
	
	public static ModelFitter createModelFitterInstance(int chr, FileNames files) {
		if (files.analysis.length() == 0) return  new ModelFitterStepwiseNoMissing(chr, files);
		String name = files.analysis.substring(0, 4);
		if (name.equalsIgnoreCase("step")) {
			files.analysis = ANALYSIS_STEPWISE_NOMISSING;
			return new ModelFitterStepwiseNoMissing(chr, files);
		}
		if (name.equalsIgnoreCase("hapl")) {
			files.analysis = ANALYSIS_HAPLOTYPE;
			return new ModelFitterStepwiseNoMissing(chr, files);
		}
		if (name.equalsIgnoreCase("sing")) {
			files.analysis = ANALYSIS_SINGLESNP;
			return new ModelFitterNoMissingData(chr, files);
		}
		if (name.equalsIgnoreCase("perm")) {
			files.analysis = ANALYSIS_PERMUTATIONS;
			return new ModelFitterWithPermutations(chr, files);
		}
		if (name.equalsIgnoreCase("boot")) {
			return new ModelFitterBootstrapStepwise(chr, files);
		}
		return null;
	}
	
	public static double getLimit(int trait, int chromosome) {
		switch (trait) {
		case 0:
			switch (chromosome) {
			case 1: return 5.94E-06;
			case 2: return 5.22E-06;
			case 3: return 5.51E-06;
			case 4: return 6.04E-06;
			case 5: return 6.69E-06;
			case 6: return 7.78E-06;
			case 7: return 7.66E-06;
			case 8: return 8.52E-06;
			case 9: return 1.13E-05;
			case 10: return 1.18E-05;
			default: return 0;
			}
		case 1:
			switch (chromosome) {
			case 1: return 4.61E-06;
			case 2: return 8.11E-06;
			case 3: return 6.38E-06;
			case 4: return 7.47E-06;
			case 5: return 5.73E-06;
			case 6: return 8.40E-06;
			case 7: return 9.64E-06;
			case 8: return 1.26E-05;
			case 9: return 8.58E-06;
			case 10: return 1.24E-05;
			default: return 0;
			}
		case 2: //									
			switch (chromosome) {
			case 1: return 6.44E-06;
			case 2: return 6.54E-06;
			case 3: return 6.38E-06;
			case 4: return 7.60E-06;
			case 5: return 6.46E-06;
			case 6: return 6.55E-06;
			case 7: return 7.41E-06;
			case 8: return 9.75E-06;
			case 9: return 8.92E-06;
			case 10: return 1.08E-05;
			default: return 0;
			}
		default:
			return 0;
		}
	}
	
	public static void testChr3() {
		ModelFitterNoMissingData mfit = new ModelFitterNoMissingData(3,0);
		mfit.run();
	}
	
	public static void LD(int chromosome) {
		System.out.println("Calculating r2 for subsample.");
		LinkageDisequilibrium ld = new LinkageDisequilibrium(chromosome);
		System.out.println("Finished.");
		
	}
	
	public static void testlg1() {
		FileNames fn = new FileNames("C:/Projects/NAM/leaf traits/filenames_leaflength_bychr.txt");
//		ModelFitterMultiSNP mf = new ModelFitterMultiSNP(2, fn, 3209000, 5212000);
		long start = System.currentTimeMillis();
		ModelFitterMultiSNP mf = new ModelFitterMultiSNP(2, fn, 4000000, 4400000);
		mf.run();
		System.out.println("Analysis finished. Elapsed time = " +(System.currentTimeMillis() - start));
	}
}
