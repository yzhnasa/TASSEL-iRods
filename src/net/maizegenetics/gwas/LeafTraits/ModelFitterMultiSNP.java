package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedWriter;
import java.util.ArrayList;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.SweepFast;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class ModelFitterMultiSNP extends ModelFitter {
	
	int[] sampleIndex;
	ArrayList<String> popLabels;
	int maxmarker;
	ArrayList<String> pops;
	int[] popIndex;
	ArrayList<Integer> snpList = null;
	ArrayList<double[]> parentSNPGenotypes = new ArrayList<double[]>();
	ArrayList<Integer> snpPosition = new ArrayList<Integer>();
	ArrayList<String> snpAllele = new ArrayList<String>();
	int startPosition, endPosition;

	public ModelFitterMultiSNP(int chromosome, FileNames files, int startPosition, int endPosition) {
		super(chromosome, files);
		this.startPosition = startPosition;
		this.endPosition = endPosition;
	}
	
	public void testSnps() {
		testSnpsNoInteraction("twosnptest.txt", true);
	}
	
	public void testSnpsNoInteraction(String output, boolean pairs) {

		// load the desired snps
		SnpData snpdata = new SnpDataImputed(chromosome, files);
		while (snpdata.next()) {
			int pos = snpdata.getPosition();
			if (pos >= startPosition && pos <= endPosition) {
				parentSNPGenotypes.add(snpdata.getGenotype());
				snpPosition.add(new Integer(pos));
				snpAllele.add(snpdata.getAllele());
			}
		}

		totalSnps = parentSNPGenotypes.size();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		popIndex = new int[nSamples];
		pops = new ArrayList<String>();

		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);

		//get the data
		double[] y = new double[nSamples];
		for (int i= 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];

		//setup the output
		BufferedWriter bw;
		bw = openOutputFile("C:/Projects/NAM/leaf traits/liguless/" + output);
		StringBuilder sb = new StringBuilder("SNP1_pos\tSNP1_allele\tSNP2_pos\tSNP2_allele\tF\tp\tlog(1/p)");
		writeToOutput(sb.toString(), bw);

		snpcount = 0;


		//set up the reduced model with no snps
		FactorModelEffect meanModelEffect = new FactorModelEffect(new int[nSamples], false);
		int[] poplevels = ModelEffectUtils.getIntegerLevels(pops);
		FactorModelEffect popModelEffect = new FactorModelEffect(poplevels, true);

		DoubleMatrix[][] rxtx = new DoubleMatrix[2][2];
		rxtx[0][0] = meanModelEffect.getXtX();
		rxtx[0][1] = meanModelEffect.getXtX2(popModelEffect);
		rxtx[1][1] = popModelEffect.getXtX();

		DoubleMatrix[] rxty = new DoubleMatrix[2];
		rxty[0] = meanModelEffect.getXty(y);
		rxty[1] = popModelEffect.getXty(y);

		int numberOfxtxColumns = rxtx[0][0].numberOfColumns() + rxtx[0][1].numberOfColumns();
		double yty = 0;
		for (double d:y) yty += d*d;

		SweepFast sf = new SweepFast(rxtx, rxty, yty);
		double reducedModelDF = 0;
		for (int i = 0; i < numberOfxtxColumns; i++) {
			if (sf.revg2sweep(i)) reducedModelDF++;
		}
		double reducedModelResidualSS = sf.getResidualSS();

		System.out.println("Total snps to be tested = " + totalSnps);

		//test all individual snps
		if (!pairs) {
			DoubleMatrix[][] xtx = new DoubleMatrix[3][3];
			xtx[0][0] = rxtx[0][0];
			xtx[0][1] = rxtx[0][1];
			xtx[1][1] = rxtx[1][1];
			DoubleMatrix[] xty = new DoubleMatrix[3];
			xty[0] = rxty[0];
			xty[1] = rxty[1];

			numberOfxtxColumns++;
			for (int s = 0; s < totalSnps; s++) {
				if (s%1000 == 0 ) System.out.println("Testing snp " + s);
				double[] snp = projectSnp(parentSNPGenotypes.get(s), snpPosition.get(s), popIndex);
				//calculate sum of squares
				double sumsq = 0;
				for (double d:snp) sumsq += d*d;

				xtx[0][2] = meanModelEffect.getXty(snp);
				xtx[1][2] = popModelEffect.getXty(snp);
				xtx[2][2] = DoubleMatrixFactory.DEFAULT.make(1, 1, sumsq);

				double sumprod = 0;
				for (int i = 0; i < nSamples; i++) sumprod += y[i] * snp[i];
				xty[2] = DoubleMatrixFactory.DEFAULT.make(1,1,sumprod);

				sf = new SweepFast(xtx, xty, yty);
				double fullModelDF = 0;
				for (int i = 0; i < numberOfxtxColumns; i++) {
					if (sf.revg2sweep(i)) fullModelDF++;
				}
				double fullModelResidualSS = sf.getResidualSS();

				double F,p;
				double snpDF = fullModelDF - reducedModelDF;
				double errorDF = nSamples - fullModelDF;
				double snpMS = (reducedModelResidualSS - fullModelResidualSS) / snpDF;
				F =  snpMS / fullModelResidualSS * errorDF;
				try {
					p = AbstractLinearModel.Ftest(F, snpDF, errorDF);
				} catch (Exception e) {
					p = Double.NaN;
				}

				sb = new StringBuilder();
				sb.append(snpPosition.get(s)).append("\t");
				sb.append(snpAllele.get(s)).append("\t");
				sb.append("none").append("\t");
				sb.append("none").append("\t");
				sb.append(F).append("\t");
				sb.append(p).append("\t");
				sb.append(-Math.log10(p));
				writeToOutput(sb.toString(), bw);

			}
		} else {
			
			//second, test all snp pairs
			DoubleMatrix[][] xtx = new DoubleMatrix[4][4];
			xtx[0][0] = rxtx[0][0];
			xtx[0][1] = rxtx[0][1];
			xtx[1][1] = rxtx[1][1];
			DoubleMatrix[] xty = new DoubleMatrix[4];
			xty[0] = rxty[0];
			xty[1] = rxty[1];
			numberOfxtxColumns += 2;

			for (int s1 = 0; s1 < totalSnps; s1++) {
				if (s1%1000 == 0) System.out.println("snp1 = " + s1);
				double[] snp1 = projectSnp(parentSNPGenotypes.get(s1), snpPosition.get(s1), popIndex);
				double sumsq = 0;
				for (double d:snp1) sumsq += d*d;

				xtx[0][2] = meanModelEffect.getXty(snp1);
				xtx[1][2] = popModelEffect.getXty(snp1);
				xtx[2][2] = DoubleMatrixFactory.DEFAULT.make(1, 1, sumsq);

				double sumprod = 0;
				for (int i = 0; i < nSamples; i++) sumprod += y[i] * snp1[i];
				xty[2] = DoubleMatrixFactory.DEFAULT.make(1,1,sumprod);

				int pos1 = snpPosition.get(s1);
				String allele1 = snpAllele.get(s1);

				for (int s2 = s1 + 1; s2 < totalSnps; s2++) {
					double[] snp2 = projectSnp(parentSNPGenotypes.get(s2), snpPosition.get(s2), popIndex);
					sumsq = 0;
					sumprod = 0;
					double sumprod12 = 0;
					for (int i = 0; i < nSamples; i++) {
						double val = snp2[i];
						sumsq += val * val;
						sumprod += val * y[i];
						sumprod12 += val * snp1[i];
					}

					xtx[0][3] = meanModelEffect.getXty(snp2);
					xtx[1][3] = popModelEffect.getXty(snp2);
					xtx[2][3] = DoubleMatrixFactory.DEFAULT.make(1, 1, sumprod12);
					xtx[3][3] = DoubleMatrixFactory.DEFAULT.make(1, 1, sumsq);
					xty[3] = DoubleMatrixFactory.DEFAULT.make(1, 1, sumprod);

					sf = new SweepFast(xtx, xty, yty);
					double fullModelDF = 0;
					for (int i = 0; i < numberOfxtxColumns; i++) {
						if (sf.revg2sweep(i)) fullModelDF++;
					}
					double fullModelResidualSS = sf.getResidualSS();

					double F,p;
					double snpDF = fullModelDF - reducedModelDF;
					double errorDF = nSamples - fullModelDF;
					double snpMS = (reducedModelResidualSS - fullModelResidualSS) / snpDF;
					F =  snpMS / fullModelResidualSS * errorDF;
					try {
						p = AbstractLinearModel.Ftest(F, snpDF, errorDF);
					} catch (Exception e) {
						p = Double.NaN;
					}

					sb = new StringBuilder();
					sb.append(pos1).append("\t");
					sb.append(allele1).append("\t");
					sb.append(snpPosition.get(s2)).append("\t");
					sb.append(snpAllele.get(s2)).append("\t");
					sb.append(F).append("\t");
					sb.append(p).append("\t");
					sb.append(-Math.log10(p));
					writeToOutput(sb.toString(), bw);

				}
			}
		}

		closeOutput(bw);

	}

	public double[] projectSnp(double[] parents, int pos, int[] popIndex) {
		Object[][] markers = theAGPMap.getInterval(chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = chromosomeLength[chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxmarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;
		
		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		for (int i = 0; i < nSamples; i++) {
			snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
			snpvalues[i] = snpvalues[i] * parents[popIndex[i]];
		}
		return snpvalues;
	}

}
