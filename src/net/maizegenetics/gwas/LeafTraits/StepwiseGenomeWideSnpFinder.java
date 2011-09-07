package net.maizegenetics.gwas.LeafTraits;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EventListener;

import javax.swing.Timer;

import cern.colt.matrix.DoubleMatrix1D;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.LinearModelWithSweep;
import net.maizegenetics.jGLiM.ModelEffect;

public class StepwiseGenomeWideSnpFinder implements EventListener {
	BestSnpFinder[] finders;
	int chromosomesRunning;
	BufferedWriter bwsteps;
	BufferedWriter bwmodel;
	ArrayList<ModelEffect> modelEffects;
	double[] phenotypes;
	FileNames files;
	String[] samples;
	int[] popIndex;
	
	double entryLimit = 1e-4;
	int maxsnps = 500;
	
	public void findSnps(String fileList) {
		
		files = new FileNames(fileList);
		
		//import phenotypes, sample names
		ArrayList<String> sampleList = new ArrayList<String>();
		ArrayList<Double> phenoList = new ArrayList<Double>();
		ArrayList<Integer> popList = new ArrayList<Integer>();
		
		try {
			//import phenotypes
			BufferedReader br = new BufferedReader(new FileReader(files.phenotypes));
			br.readLine();
			String input = br.readLine();
			while (input != null) {
				String[] parsedLine = input.split("\t");
				try {
					String sample = parsedLine[0];
					Integer pop = Integer.decode(parsedLine[1]);
					Double trait = Double.valueOf(parsedLine[2]);
					sampleList.add(sample);
					phenoList.add(trait);
					popList.add(pop);
				}
				catch(Exception e) {}
				input = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		samples = new String[sampleList.size()];
		sampleList.toArray(samples);
		int nSamples = samples.length;
		popIndex = new int[nSamples];
		phenotypes = new double[nSamples];
		
		for (int i = 0; i < nSamples; i++) {
			popIndex[i] = popList.get(i) - 1;
			phenotypes[i] = phenoList.get(i);
		}
		
		//build the model
		int[] mean = new int[nSamples];
		ModelEffect meanEffect = new ModelEffect(mean);
		ModelEffect popEffect = new ModelEffect(ModelEffect.getIntegerLevels(popList));
		modelEffects = new ArrayList<ModelEffect>();
		modelEffects.add(meanEffect);
		modelEffects.add(popEffect);
		
		openFilesForOutput();
		fitNextSnp();
		
	}
	
	public synchronized void thisChromosomeIsFinished(int chromosome) {
		chromosomesRunning--;
		System.out.println("Chromosome " + chromosome + " is finished. chromosomesRunning = " + chromosomesRunning);
		if (chromosomesRunning == 0) pickBestSnp();
	}
	
	public void pickBestSnp() {
		System.out.println("picking best snp...");
		int bestchr = 0;
		double minp = finders[0].minp;
		double maxF = finders[0].maxF;
		
		for (int i = 1; i < 10; i++) {
			if (finders[i].minp < minp || (finders[i].minp == minp && finders[i].maxF > maxF)) {
				minp = finders[i].minp;
				maxF = finders[i].maxF;
				bestchr = i;
			}
		}
		
		if (minp < entryLimit && modelEffects.size() < 2 + maxsnps) {
			
			try {
				StringBuilder sb = new StringBuilder();
				sb.append(finders[bestchr].chromosome);
				int pos = finders[bestchr].getBestPosition();
				sb.append("\t").append(finders[bestchr].getBestPosition());
				sb.append("\t").append(finders[bestchr].getBestCm());
				sb.append("\t").append(maxF);
				sb.append("\t").append(minp);
				System.out.println(sb.toString());
				bwsteps.write(sb.toString());
				bwsteps.newLine();
				bwsteps.flush();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			
			modelEffects.add(finders[bestchr].bestEffect);
			System.out.println("fitting next snp...");
			fitNextSnp();
		}
		else {
			System.out.println("Finished fitting snps, summarizing model...");
			LinearModelWithSweep lmws = new LinearModelWithSweep(modelEffects, phenotypes);
			double errss = lmws.getErrorSS();
			double errdf = lmws.getErrordf();
			
			int neffects = modelEffects.size();
			for (int i = 2; i < neffects; i++) {
				ModelEffect me = modelEffects.get(i);
				StringBuilder sb = new StringBuilder();
				int[] id = (int[]) me.getId();
				sb.append(id[0]).append("\t").append(id[1]);
				double cmpos = finders[0].theAGPMap.getCmFromPosition(id[0],id[1]);
				sb.append("\t").append(cmpos);
				
				//append effect
				DoubleMatrix1D beta = lmws.getBeta();
				int nparam = beta.size();
				int whichParam = nparam - neffects + i;
				sb.append("\t").append(beta.getQuick(whichParam));
				
				//append F, p
				double[] ssdf = lmws.marginalEffectSSdf(i);
				lmws.incrementalEffectSS();
				double F = ssdf[0] / ssdf[1] / errss * errdf;
				double p;
				try {p = AbstractLinearModel.Ftest(F, ssdf[1], errdf);}
				catch (Exception e) {p = Double.NaN;}
				sb.append("\t").append(F);
				sb.append("\t").append(p);
				
				try {
					bwmodel.write(sb.toString());
					bwmodel.newLine();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
			
			try {
				bwmodel.close();
				bwsteps.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
	public void fitNextSnp() {
		//instantiate the chromosome threads
		chromosomesRunning = 0;
		finders = new BestSnpFinder[10];
		for (int i = 0; i < 10; i++) {
			finders[i] = new BestSnpFinder(this, files, i + 1, phenotypes, samples, popIndex);
			finders[i].setBaseEffects(modelEffects);
			finders[i].start();
			chromosomesRunning++;
		}
	}
	
	public void openFilesForOutput() {
		try {
			bwsteps = new BufferedWriter(new FileWriter(files.steps));
			bwmodel = new BufferedWriter(new FileWriter(files.model));
			bwsteps.write("Chromosome\tposition\tcM\tF\tpvalue");
			bwsteps.flush();
			bwsteps.newLine();
			bwmodel.write("Chromosome\tposition\tcM\teffect\tF\tpvalue");
			bwsteps.flush();
			bwmodel.newLine();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void closeOutput() {
		try {
			bwsteps.close();
			bwmodel.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
}
