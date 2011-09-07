package net.maizegenetics.gwas.jointlinkage;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.NestedCovariateModelEffect;
import net.maizegenetics.jGLiM.dm.PartitionedLinearModel;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;

public class StepwiseModelFitterForNAMJointLinkage {
	double[] data;
	String[] pops;
	SweepFastLinearModel sflm;
	ArrayList<ModelEffect> model = new ArrayList<ModelEffect>();
	FactorModelEffect popEffect;
	int numberOfFiles = 10;
	double enterlimit = 1e-6;
	double exitlimit = 2e-6;
	boolean multithreaded = false;
	int numberOfProcessors;
	SNPdata snpdata;
	int[][] snpIndices;
	BufferedWriter bwstep = null;
	BufferedWriter bwmodel = null;
	ExecutorService exec;
	NamMap nammap;

	public StepwiseModelFitterForNAMJointLinkage() {
		loadData();
		numberOfProcessors = Runtime.getRuntime().availableProcessors();
		if (numberOfProcessors == 1) multithreaded = false;
		//		multithreaded = false;
		if (multithreaded) {
			exec = Executors.newFixedThreadPool(numberOfProcessors); 
			subdivideSnps();
			long start = System.currentTimeMillis();
			fitModel();
			System.out.println("time to fit model: " + (System.currentTimeMillis() - start));
			exec.shutdown();
			scanAndFindCI();
		} else {
			fitModel();
			scanAndFindCI();
		}
	}

	private void loadData() {
//		snpdata = new SNPdataFromArray();
		snpdata = new SNPdataNAMGBS();
		data = snpdata.getPhenotype();
		pops = snpdata.getPopulations();
	}

	private void subdivideSnps() {
		int totalNumberOfSnps = snpdata.getNumberOfSNPs();
		int snpsPerSet = totalNumberOfSnps / numberOfProcessors;
		int start = 0;
		int end = snpsPerSet - 1;
		snpIndices = new int[numberOfProcessors][2];

		for (int i = 0; i < numberOfProcessors; i++) {
			snpIndices[i][0] = start;
			snpIndices[i][1] = end;
			start += snpsPerSet;
			end += snpsPerSet;
		}

		snpIndices[numberOfProcessors - 1][1] = totalNumberOfSnps - 1;
	}

	public void fitModel() {
		int numberOfTaxa = data.length;
		int[] mean = new int[numberOfTaxa];
		ArrayList<Object> popLevelIds = new ArrayList<Object>();
		int[] poplevels = ModelEffectUtils.getIntegerLevels(pops, popLevelIds);
		FactorModelEffect meanEffect = new FactorModelEffect(mean, false);
		meanEffect.setID("mean");
		popEffect = new FactorModelEffect(poplevels, true);
		popEffect.setID("Population");

		model.add(meanEffect);
		model.add(popEffect);

		sflm = new SweepFastLinearModel(model, data);

		//initialize file for step output
		openStepOutputFile(true);
		closeStepOutputFile();

		//fit a model
		if (multithreaded) {
			while(forwardStepMT()) {
				while(backwardStep());
			}
		} else {
			while(forwardStep()) {
				while(backwardStep());
			}
		}
		writeModelOutputToFile();
	}

	public void writeModelOutputToFile() {
		openModelOutputFile(true, false);
		double[] errorssdf = sflm.getResidualSSdf();
		double totalSS = sflm.getModelcfmSSdf()[0];
		for (int i = 0; i < model.size(); i++) {
			String chr, pos;
			String name = " ";
			if (i == 0) {
				chr = "mean";
				pos = " ";
			} else if (i == 1) {
				chr = "population";
				pos = " ";
			} else {
				SNP snp = (SNP) model.get(i).getID();
				chr = Integer.toString(snp.chr);
				pos = Double.toString(snp.pos);
				name = snp.name;
			}

			double[] ssdf = sflm.getMarginalSSdf(i);
			double SS = ssdf[0];
			int df = (int) ssdf[1];
			double MS = ssdf[0] / ssdf[1];
			double F = MS / errorssdf[0] * errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, ssdf[1], errorssdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}
			double rsq = SS / totalSS;
			outputModel(chr, pos, name, SS, df, MS, F, p, rsq);
		}
		closeModelOutputFile();
	}

	//returns true if a new term is added to the model
	public boolean forwardStep() {
		double bestss = 0;
		NestedCovariateModelEffect besteffect = null;

		//test all the snps. Select the one that gives the highest model SS
		PartitionedLinearModel plm = new PartitionedLinearModel(model, sflm);
		int numberOfSites = snpdata.getNumberOfSNPs(); //subtract one for the header row
		for (int s = 0; s < numberOfSites; s++) {
			SNP snp = snpdata.getSnp(s);
			if (snp != null) {
				NestedCovariateModelEffect nested = new NestedCovariateModelEffect(snp.score, popEffect);
				nested.setID(snp);
				plm.testNewModelEffect(nested);
				double modelss = plm.getModelSS();
				if (modelss > bestss) {
					bestss = modelss;
					besteffect = nested;
				}
			}
		}

		//if the p-value for the select SNP is less than the enter limit, add it to the model and recalculate the model solution
		plm.testNewModelEffect(besteffect);
		double[] Fp = plm.getFp();
		if (Fp[1] < enterlimit) {
			SNP snp = (SNP) besteffect.getID();
			outputStep(true, snp, Fp[0], Fp[1]);
			model.add(besteffect);
			sflm = new SweepFastLinearModel(model, data);
			return true;
		} else {
			return false;
		}
	}

	//returns true if a new term is added to the model
	//multithreaded version of forward step
	public boolean forwardStepMT() {

		//submit each subset of snps to a different thread
		ArrayList<Future<ForwardStepResult>> theFutures = new ArrayList<Future<ForwardStepResult>>();
		for (int i = 0; i < numberOfProcessors; i++) {
			//copy model
			ArrayList<ModelEffect> modelcopy = new ArrayList<ModelEffect>();
			for (ModelEffect effect : model) modelcopy.add(effect.getCopy());
			FactorModelEffect popeffectcopy = (FactorModelEffect) modelcopy.get(1);
			PartitionedLinearModel plm = new PartitionedLinearModel(modelcopy, sflm);
			theFutures.add(exec.submit(new ForwardStepJointLinkage(plm, snpdata, popeffectcopy, snpIndices[i][0], snpIndices[i][1])));
		}

		ForwardStepResult bestResult = null;
		for (Future<ForwardStepResult> f : theFutures) {
			ForwardStepResult thisResult;
			try {
				thisResult = f.get();
			} catch (Exception e) {
				thisResult = null;
			}

			if (bestResult == null) {
				bestResult = thisResult;
			} else if (thisResult != null) {
				if (thisResult.bestF > bestResult.bestF) bestResult = thisResult;
			}
		}

		//if the p-value for the select SNP is less than the enter limit, add it to the model and recalculate the model solution
		if (bestResult.bestp < enterlimit) {
			SNP snp = (SNP) bestResult.bestEffect.getID();
			outputStep(true, snp, bestResult.bestF, bestResult.bestp);

			model.add(bestResult.bestEffect);
			sflm = new SweepFastLinearModel(model, data);
			return true;
		} else {
			return false;
		}
	}

	//returns true if a term is removed from the model
	public boolean backwardStep() {
		int numberOfTerms = model.size();
		if (numberOfTerms <= 2) return false;

		//find the model term (snps only) with the largest p-value
		double maxp = 0;
		double minF= -1;
		int maxterm = 0;
		double[] errorssdf = sflm.getResidualSSdf();

		for (int t = 2; t < numberOfTerms; t++) {
			double[] termssdf = sflm.getIncrementalSSdf(t);
			double F = termssdf[0]/termssdf[1]/errorssdf[0]*errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, termssdf[1], errorssdf[1]);
				if (p > maxp) {
					maxterm = t;
					maxp = p;
					minF = F; 
				}
			} catch(Exception e){
				p = Double.NaN;
			}
		}

		//if that maxp is >= exitlimit, then remove maxterm from the model, recalculate the model, and return true;
		if (maxp >= exitlimit) {
			ModelEffect me = model.remove(maxterm);
			SNP snp = (SNP) me.getID();
			outputStep(false, snp, minF, maxp);
			sflm = new SweepFastLinearModel(model, data);
			return true;
		}

		return false;
	}

	public void scanAndFindCI() {
		//rescan the model to refine position and find CI's
		//rescan requires a map

		//for each snp in the model:
		//1. add an adjacent snp to left of the original
		//2. solve the new model
		//3. if the marginal p of the model snp is sig stop and repeat on right of snp.
		//4. find the max ss of any point with the CI interval (not including the original marker)
		//5. if the max ss is greater than the original marker ss, replace the original and go to 1.
		//report steps as process proceeds

		int numberOfTerms = model.size();
		int[] upperbound = new int[numberOfTerms - 2];
		int[] lowerbound = new int[numberOfTerms - 2];
		for (int t = 2; t < numberOfTerms; t++){ //0 is the mean, 1 is populations
			//find the CI bounds
			lowerbound[t - 2] = scanASide(true, t);
			upperbound[t - 2] = scanASide(false, t);

			//find the marker with highest ss in the CI
			SNP bestsnp = null;
			double bestss = 0;
			ModelEffect besteffect = null;
			ModelEffect currentme = model.remove(t);
			sflm = new SweepFastLinearModel(model, data);
			PartitionedLinearModel plm = new PartitionedLinearModel(model, sflm);
			for (int m = lowerbound[t - 2]; m <= upperbound[t - 2]; m++) {
				SNP testsnp = snpdata.getSnp(m);
				NestedCovariateModelEffect snpeffect = new NestedCovariateModelEffect(testsnp.score, popEffect);
				snpeffect.setID(testsnp);
				plm.testNewModelEffect(snpeffect);
				double ss = plm.getModelSS();
				if (ss > bestss) {
					bestss = ss;
					bestsnp = testsnp;
					besteffect = snpeffect;
				}
			}
			model.add(t, besteffect);

			//did the best marker change?
			boolean markerchanged = false;
			SNP currentSnp = (SNP) currentme.getID();
			if (currentSnp.index != bestsnp.index) markerchanged = true;

			//if this marker is different than the one tested, reset the CI
			if (markerchanged) {
				lowerbound[t - 2] = scanASide(true, t);
				upperbound[t - 2] = scanASide(false, t);
			}
		}

		//report the results including CI's
		sflm = new SweepFastLinearModel(model, data);
		writeModelWithCI2File(lowerbound, upperbound);
	}

	private void writeModelWithCI2File(int[] lowerbound, int[] upperbound) {
		openModelOutputFile(true, true);
		double[] errorssdf = sflm.getResidualSSdf();
		double totalSS = sflm.getModelcfmSSdf()[0];
		for (int i = 0; i < model.size(); i++) {
			String chr, pos;
			String name = " ";
			double lower = 0;
			double upper = 0;
			if (i == 0) {
				chr = "mean";
				pos = " ";
			} else if (i == 1) {
				chr = "population";
				pos = " ";
			} else {
				SNP snp = (SNP) model.get(i).getID();
				chr = Integer.toString(snp.chr);
				pos = Double.toString(snp.pos);
				name = snp.name;
				snp = snpdata.getSnp(lowerbound[i - 2]);
				lower = snp.pos;
				snp = snpdata.getSnp(upperbound[i - 2]);
				upper = snp.pos;
			}

			double[] ssdf = sflm.getMarginalSSdf(i);
			double SS = ssdf[0];
			int df = (int) ssdf[1];
			double MS = ssdf[0] / ssdf[1];
			double F = MS / errorssdf[0] * errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, ssdf[1], errorssdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}
			double rsq = SS / totalSS;
			outputModel(chr, pos, name, SS, df, MS, F, p, rsq, lower, upper);
		}
		closeModelOutputFile();
	}
	
	private int scanASide(boolean left, int whichModelTerm) {
		BufferedWriter bwlog = null;
		try {
			bwlog = new BufferedWriter(new FileWriter(FileNames.scanFilename, true));
		} catch (IOException e) {
			System.out.println("unable to open " + FileNames.scanFilename);
		}
		
		double alpha = 0.05;
		int minIndex = 0;
		int maxIndex = snpdata.getNumberOfSNPs() - 1;
		int incr;
		if (left) {
			incr = -1;
		} else {
			incr = 1;
		}

		SNP modelsnp = (SNP) model.get(whichModelTerm).getID();
		int markerIndex = modelsnp.index;
		int chr = modelsnp.chr;
		boolean boundfound = false;
		int testIndex = markerIndex;
		int lastterm = model.size();
		do {
			testIndex += incr;

			SNP testsnp = snpdata.getSnp(testIndex);
			if (testsnp.chr != chr) {
				testIndex -= incr;
				boundfound = true;
			} else {
				NestedCovariateModelEffect snpeffect = new NestedCovariateModelEffect(testsnp.score, popEffect);
				snpeffect.setID(testsnp);
				model.add(snpeffect);
				SweepFastLinearModel sflm = new SweepFastLinearModel(model, data);
				double[] snpssdf = sflm.getMarginalSSdf(whichModelTerm);
				double[] errorssdf = sflm.getResidualSSdf();
				double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
				double p;

				try {
					p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
				} catch(Exception e) {
					p = 1;
				}

				//write result to log
				StringBuffer sb = new StringBuffer();
				sb.append(modelsnp.chr);
				sb.append("\t").append(modelsnp.pos);
				sb.append("\t").append(testsnp.pos);
				sb.append("\t").append(testsnp.name);
				sb.append("\t").append(F);
				sb.append("\t").append(p);
				System.out.println(sb);
				try {
					bwlog.write(sb.toString());
					bwlog.newLine();
				} catch(IOException e) {}

				if (p < alpha) {
					boundfound = true;
				}
			}
			model.remove(lastterm);
		} while (!boundfound && testIndex > minIndex && testIndex < maxIndex);

		try {
			bwlog.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return testIndex;
	}

	private void outputStep(boolean added, SNP snp, double F, double p){
		StringBuilder sb = new StringBuilder(snpdata.getPhenotypeName());

		if (added) {
			sb.append("\t").append("added");
		} else {
			sb.append("\t").append("removed");
		}

		sb.append("\t").append(snp.chr);
		sb.append("\t").append(snp.pos);
		sb.append("\t").append(snp.name);
		sb.append("\t").append(F);
		sb.append("\t").append(p);
		System.out.println(sb);

		openStepOutputFile(false);
		try {
			bwstep.write(sb.toString());
			bwstep.newLine();
		} catch (IOException e) {
			System.err.println("Failed writing to step output: " + sb);
			e.printStackTrace();
		}
		closeStepOutputFile();
	}

	private void outputModel(String chr, String pos, String name, double SS, int df, double MS, double F, double p, double rsq) {
		StringBuilder sb = new StringBuilder(snpdata.getPhenotypeName());
		sb.append("\t").append(chr);
		sb.append("\t").append(pos);
		sb.append("\t").append(name);
		sb.append("\t").append(SS);
		sb.append("\t").append(df);
		sb.append("\t").append(MS);
		sb.append("\t").append(F);
		sb.append("\t").append(p);
		sb.append("\t").append(rsq);
		try {
			bwmodel.write(sb.toString());
			bwmodel.newLine();
		} catch (IOException e) {
			System.err.println("Failed writing to model output: " + sb);
			e.printStackTrace();
		}
	}

	private void outputModel(String chr, String pos, String name, double SS, int df, double MS, double F, double p, double rsq, double lowerCI, double upperCI) {
		StringBuilder sb = new StringBuilder(snpdata.getPhenotypeName());
		sb.append("\t").append(chr);
		sb.append("\t").append(pos);
		sb.append("\t").append(name);
		sb.append("\t").append(SS);
		sb.append("\t").append(df);
		sb.append("\t").append(MS);
		sb.append("\t").append(F);
		sb.append("\t").append(p);
		sb.append("\t").append(rsq);
		sb.append("\t").append(lowerCI);
		sb.append("\t").append(upperCI);
		try {
			bwmodel.write(sb.toString());
			bwmodel.newLine();
		} catch (IOException e) {
			System.err.println("Failed writing to model output: " + sb);
			e.printStackTrace();
		}
	}

	private void openStepOutputFile(boolean writeHeader) {
		try {
			bwstep = new BufferedWriter(new FileWriter(FileNames.stepFilename, true));
			if (writeHeader) {
				bwstep.write("Trait\tAction\tChromosome\tPosition\tName\tF\tpvalue");
				bwstep.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void openModelOutputFile(boolean writeHeader, boolean withCI) {
		try {
			bwmodel = new BufferedWriter(new FileWriter(FileNames.modelFilename, true));
			if (writeHeader && withCI) {
				bwmodel.write("Trait\tChromosome\tPosition\tName\tSS\tdf\tMS\tF\tpvalue\tRsq\tlowerCI\tupperCI");
				bwmodel.newLine();
			} else if (writeHeader) {
				bwmodel.write("Trait\tChromosome\tPosition\tName\tSS\tdf\tMS\tF\tpvalue\tRsq");
				bwmodel.newLine();
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void closeStepOutputFile() {
		try {
			bwstep.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void closeModelOutputFile() {
		try {
			bwmodel.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void setEnterLimit(double limit) { enterlimit = limit; }
	public void setExitLimit(double limit) { exitlimit = limit; }
}
