package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.NestedCovariateModelEffect;
import net.maizegenetics.stats.linearmodels.PartitionedLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.trait.MarkerPhenotypeAdapter;
import net.maizegenetics.trait.MarkerPhenotypeAdapterUtils;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.LinkedList;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.BasicShuffler;

public class StepwiseOLSModelFitter {
	private MarkerPhenotypeAdapter myData;

	//settable parameters
	private double[] enterlimits = null;
	private double[] exitlimits = null;
	private double enterlimit = 1e-5;
	private double exitlimit = 2e-5;
	private int maxNumberOfMarkers = 1000;
	private FactorModelEffect nestingEffect;
	private int nestingFactorIndex;
	private boolean isNested;
	private ArrayList<String> nestingFactorNames;

	//global variables used by the analysis
	private ArrayList<ModelEffect> currentModel;
	private int currentPhenotypeIndex;
	private int numberOfBaseModelEffects;
	private double[] y;
	private boolean[] missing;
	int numberNotMissing;
	private ArrayList<String[]> factorList;
	private ArrayList<double[]> covariateList;
	private LinkedList<Object[]> resultRowsAnova = new LinkedList<Object[]>();
	private LinkedList<Object[]> resultRowsAnovaWithCI = new LinkedList<Object[]>();        
	private LinkedList<Object[]> rowsSiteEffectTable = new LinkedList<Object[]>();
	private LinkedList<Object[]> rowsSiteEffectTableWithCI = new LinkedList<Object[]>();
	private MODEL_TYPE modelType = MODEL_TYPE.mbic;
	private String datasetName;
	private double globalbestbic = Double.MAX_VALUE; //Rationale: the optimal model has minimum BIC. Thus, initialize bestbic with the largest possible double java can handle.
	private double globalbestmbic = Double.MAX_VALUE;
	private double globalbestaic = Double.MAX_VALUE;
	public enum MODEL_TYPE {pvalue, bic, mbic, aic};
	private boolean test;
    private int[][] theUpperAndLowerBound;

    private int numberOfPermutations = 0;
    private double alpha = 0.05;
    private double[] minPvalues = null;

	private final String[] anovaReportHeader = new String[]{"Trait", "Name","Locus","Position","df","SS","MS", "F", "pr>F", "BIC", "mBIC", "AIC"};
    private final String[] anovaReportWithCIHeader = new String[]{"Trait", "Name","Locus","Position","df","SS","MS", "F", "pr>F", "BIC", "mBIC", "AIC", "SuppLeft", "SuppRight"};
	public void setModelType(MODEL_TYPE modelType){
		this.modelType = modelType;
	}

	public StepwiseOLSModelFitter(MarkerPhenotypeAdapter anAdapter, String datasetName) {
		myData = anAdapter;
		this.datasetName = datasetName;
	}

	public DataSet runAnalysis() {
		//numbers of various model components
		int numberOfPhenotypes = myData.getNumberOfPhenotypes();

		//cycle through the phenotypes
		//notation: 
		//X is the design matrix without the markers, rows of X will be deleted if marker data is missing
		//Xm is the design matrix with markers
		//y is the data

		for (int ph = 0; ph < numberOfPhenotypes; ph++) {
			currentPhenotypeIndex = ph;
			if (enterlimits != null) enterlimit = enterlimits[ph];
			if (exitlimits != null) exitlimit = exitlimits[ph];

			//get phenotype data
			double[] phenotypeData = myData.getPhenotypeValues(ph);

			//keep track of missing rows
			missing = myData.getMissingPhenotypes(ph);

			//get factors
			factorList = MarkerPhenotypeAdapterUtils.getFactorList(myData, ph, missing);

			//get covariates
			covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(myData, ph, missing);

			//remove missing values from the arrays
			numberNotMissing = 0;
			int totalNumber = missing.length;
			for (boolean b:missing) if (!b) numberNotMissing++;

			int ptr = 0;
			y = new double[numberNotMissing];
			for (int i = 0; i < totalNumber; i++) {
				if (!missing[i]) y[ptr++] = phenotypeData[i];
			}

			if (factorList != null) {
				int n = factorList.size();
				for (int f = 0 ; f < n; f++) {
					String[] newfactor = new String[numberNotMissing];
					String[] oldfactor = factorList.get(f);
					ptr = 0;
					for (int i = 0; i < totalNumber; i++) {
						if (!missing[i]) newfactor[ptr++] = oldfactor[i];
					}
					factorList.set(f, newfactor);
				}
			}

			if (covariateList != null) {
				int n = covariateList.size();
				for (int f = 0 ; f < n; f++) {
					double[] newcov = new double[numberNotMissing];
					double[] oldcov = covariateList.get(f);
					ptr = 0;
					for (int i = 0; i < totalNumber; i++) {
						if (!missing[i]) newcov[ptr++] = oldcov[i];
					}
					covariateList.set(f, newcov);
				}
			}
			fitModel();
			scanAndFindCI();
		}

		return null;
	}

	public void fitModel() {
		//build the base model
		currentModel = new ArrayList<ModelEffect>();
		int numberOfTaxa = y.length;
		int[] mean = new int[numberOfTaxa];

		FactorModelEffect meanEffect = new FactorModelEffect(mean, false);
		meanEffect.setID("mean");
		currentModel.add(meanEffect);

		//add the factor effects
		if (factorList != null) {
			for (int f = 0; f < factorList.size(); f++) {
				ArrayList<String> ids = new ArrayList<String>();
				int[] levels = ModelEffectUtils.getIntegerLevels(factorList.get(f), ids);
				FactorModelEffect fme = new FactorModelEffect(levels, true, new Object[]{myData.getFactorName(f), ids});
				currentModel.add(fme);
				if (isNested && f == nestingFactorIndex) {
					nestingEffect = fme;
					nestingFactorNames = ids; 
				}
			}
		}

		//add the covariate effects
		if (covariateList != null) {
			for (int c = 0; c < covariateList.size(); c++) {
				CovariateModelEffect cme = new CovariateModelEffect(covariateList.get(c), myData.getCovariateName(c));
				currentModel.add(cme);
			}
		}
		numberOfBaseModelEffects = currentModel.size();
		//permutations should be run right here
		System.out.println("-----Number of Permutations = " + numberOfPermutations + "------------");
		if(numberOfPermutations > 0) runPermutationTestNoMissingData(); 
		System.out.println("--------------the Enter limit is " + enterlimit);
		System.out.println("--------------the Exit limit is " + exitlimit);  			

		while(forwardStep()){
			while(backwardStep());
		}

		SweepFastLinearModel theLinearModel = new SweepFastLinearModel(currentModel, y);
		appendAnovaResults(theLinearModel);
		appendSiteEffectEstimates(theLinearModel);

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
        SweepFastLinearModel sflm;
        int numberOfTerms = currentModel.size();
        int firstMarkerIndex = 1;
        if (factorList != null){
            firstMarkerIndex = firstMarkerIndex+factorList.size();
        }
        if (covariateList != null){
            firstMarkerIndex = firstMarkerIndex+covariateList.size();
        }        

        System.out.println("firstMarkerIndex " + firstMarkerIndex);
        int[] upperbound = new int[numberOfTerms - firstMarkerIndex];
        int[] lowerbound = new int[numberOfTerms - firstMarkerIndex];
        theUpperAndLowerBound = new int[numberOfTerms - firstMarkerIndex][2];
        for (int t = firstMarkerIndex; t < numberOfTerms; t++){ //0 is the mean, 1 is populations
            //find the CI bounds
            lowerbound[t - firstMarkerIndex] = scanASide(true, t);
            System.out.println("lowerbound[t - firstMarkerIndex]");
            System.out.println(lowerbound[t - firstMarkerIndex]);
            upperbound[t - firstMarkerIndex] = scanASide(false, t);
            System.out.println("upperbound[t - firstMarkerIndex]");
            System.out.println(upperbound[t - firstMarkerIndex]);
            theUpperAndLowerBound[t-firstMarkerIndex][0] = lowerbound[t-firstMarkerIndex];
            theUpperAndLowerBound[t-firstMarkerIndex][1] = upperbound[t-firstMarkerIndex];
             
            //find the marker with highest ss in the CI
            SNP bestsnp = null;
            double bestss = 0;
            ModelEffect besteffect = null;
            ModelEffect currentme = currentModel.remove(t);
            sflm = new SweepFastLinearModel(currentModel, y);
            PartitionedLinearModel plm = new PartitionedLinearModel(currentModel, sflm);
            
            			//NEW CODE: create the appropriate marker effect
            for (int m = lowerbound[t - firstMarkerIndex]; m <= upperbound[t - firstMarkerIndex]; m++) {
                ModelEffect markerEffect = null;
                SNP testsnp = new SNP(myData.getMarkerName(m), new Chromosome(myData.getLocusName(m)), (int) myData.getMarkerChromosomePosition(m), m);
                Object[] markerValues = myData.getMarkerValue(currentPhenotypeIndex, m);
                Object[] markersForNonMissing = new Object[numberNotMissing];
                int allCount = 0;
                int nonmissingCount = 0;
                for (Object marker:markerValues) {
                        if (!missing[allCount++]) markersForNonMissing[nonmissingCount++] = marker;
                }

                if (myData.isMarkerDiscrete(t)) {
                        ArrayList<Object> markerIds = new ArrayList<Object>();
                        int[] levels = ModelEffectUtils.getIntegerLevels(markersForNonMissing, markerIds);
                        testsnp.alleles = markerIds;

                        if (isNested) {
                                //not implemented yet
                                markerEffect = new FactorModelEffect(levels, true, testsnp);
                        } else {
                                markerEffect = new FactorModelEffect(levels, true, testsnp);
                        }
                } else {
                        //int n = markerValues.length;
                        double[] cov = new double[numberNotMissing];
                        for (int i = 0; i < numberNotMissing; i++) cov[i] = ((Double) markersForNonMissing[i]).doubleValue();

                        if (isNested) {
                                CovariateModelEffect cme = new CovariateModelEffect(cov);
                                markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
                                markerEffect.setID(testsnp);
                        } else {
                        markerEffect = new CovariateModelEffect(cov, testsnp);
                        }
                       }

                plm.testNewModelEffect(markerEffect);
                double modelss = plm.getModelSS();
                if (modelss > bestss) {
                    bestss = modelss;
                    bestsnp = testsnp;
                    besteffect = markerEffect;
                }
            }
            
            currentModel.add(t, besteffect);
            ////
            //did the best marker change?
            boolean markerchanged = false;
            SNP currentSnp = (SNP) currentme.getID();
            if (currentSnp.index != bestsnp.index) markerchanged = true;

            //if this marker is different than the one tested, reset the CI
            if (markerchanged) {
                lowerbound[t - firstMarkerIndex] = scanASide(true, t);
                upperbound[t - firstMarkerIndex] = scanASide(false, t);
                theUpperAndLowerBound[t-firstMarkerIndex][0] = lowerbound[t-firstMarkerIndex];
                theUpperAndLowerBound[t-firstMarkerIndex][1] = upperbound[t-firstMarkerIndex];
            }
            System.out.println("upperAndLowerBound[t-firstMarkerIndex][0]");
            System.out.println(theUpperAndLowerBound[t-firstMarkerIndex][0]); 
            System.out.println("upperAndLowerBound[t-firstMarkerIndex][1]");
            System.out.println(theUpperAndLowerBound[t-firstMarkerIndex][1]);
        }
        
        //report the results including CI's
        //writeModelWithCI2File(lowerbound, upperbound); Turn this into a global matrix of upper and lower bounds
        SweepFastLinearModel theLinearModel = new SweepFastLinearModel(currentModel, y);
        appendAnovaResultsWithCI(theLinearModel);
        appendSiteEffectEstimatesWithCI(theLinearModel);

 
    }        
        
  
    
    private int scanASide(boolean left, int whichModelTerm) {
        /*BufferedWriter bwlog = null;
        try {
            bwlog = new BufferedWriter(new FileWriter(scanFilename, true));
        } catch (IOException e) {
            System.out.println("unable to open " + scanFilename);
        }*/
       
        double alpha = 0.05;
        int minIndex = 0;
        int maxIndex = myData.getNumberOfMarkers() - 1;
        int incr;
        if (left) {
            incr = -1;
        } else {
            incr = 1;
        }

        SNP modelsnp = (SNP) currentModel.get(whichModelTerm).getID();
        int markerIndex = modelsnp.index;
        Chromosome chr = modelsnp.locus;
        int chrAsInt = chr.getChromosomeNumber();
        boolean boundfound = false;
        int testIndex = markerIndex;
        int lastterm = currentModel.size();
        
        ///
        do{
            testIndex += incr;
            ModelEffect markerEffect = null;
            SNP snp = new SNP(myData.getMarkerName(testIndex), new Chromosome(myData.getLocusName(testIndex)), (int) myData.getMarkerChromosomePosition(testIndex), testIndex);
            Chromosome chrOfTestedSnp = snp.locus;
            int chrOfTestedSnpAsInt = chrOfTestedSnp.getChromosomeNumber();
            if (snp == null || chrOfTestedSnpAsInt != chrAsInt) {
                testIndex -= incr;
                boundfound = true;
            } else {
 		//create the appropriate marker effect

		Object[] markerValues = myData.getMarkerValue(currentPhenotypeIndex, testIndex);
		Object[] markersForNonMissing = new Object[numberNotMissing];
		int allCount = 0;
		int nonmissingCount = 0;
		for (Object marker:markerValues) {
			if (!missing[allCount++]) markersForNonMissing[nonmissingCount++] = marker;
		}
		
		if (myData.isMarkerDiscrete(testIndex)) {
			ArrayList<Object> markerIds = new ArrayList<Object>();
			int[] levels = ModelEffectUtils.getIntegerLevels(markersForNonMissing, markerIds);
			snp.alleles = markerIds;
			
			if (isNested) {
			//not implemented yet
                		markerEffect = new FactorModelEffect(levels, true, snp);
			} else {
				markerEffect = new FactorModelEffect(levels, true, snp);
			}
		} else {
			//int n = markerValues.length;
			double[] cov = new double[numberNotMissing];
			for (int i = 0; i < numberNotMissing; i++) cov[i] = ((Double) markersForNonMissing[i]).doubleValue();

			if (isNested) {
				CovariateModelEffect cme = new CovariateModelEffect(cov);
				markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
				markerEffect.setID(snp);
			} else {
				markerEffect = new CovariateModelEffect(cov, snp);
			}
		}
     
               ///
                currentModel.add(markerEffect);
                SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
                double[] snpssdf = sflm.getMarginalSSdf(whichModelTerm);
                double[] errorssdf = sflm.getResidualSSdf();
                double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
                double p;
                try {
                    p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
                } catch(Exception e) {
                    p = 1;
                }     
           
                if(p < alpha){
                    boundfound = true;
                }
                currentModel.remove(lastterm);
            }    
       }  while (!boundfound && testIndex > minIndex && testIndex < maxIndex);      
        ////      
        return testIndex;
    }      
        
	public boolean forwardStep() {
		double bestss = 0;
		double bestbic = Double.MAX_VALUE;
		double bestmbic = Double.MAX_VALUE;
		double bestaic = Double.MAX_VALUE;
		ModelEffect besteffect = null;

		SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
		PartitionedLinearModel plm = new PartitionedLinearModel(currentModel, sflm);
		int numberOfSites = myData.getNumberOfMarkers();

		System.out.println("We are in forwardStep()");

		double[] temperrorssdf = sflm.getResidualSSdf();
		System.out.println("The value of errorss before the loop in forwardStep() is: " + temperrorssdf[0]);
		for (int s = 0; s < numberOfSites; s++) {
			//create the appropriate marker effect

			ModelEffect markerEffect = null;
			SNP snp = new SNP(myData.getMarkerName(s), new Chromosome(myData.getLocusName(s)), (int) myData.getMarkerChromosomePosition(s), s);
			Object[] markerValues = myData.getMarkerValue(currentPhenotypeIndex, s);
			Object[] markersForNonMissing = new Object[numberNotMissing];
			int allCount = 0;
			int nonmissingCount = 0;
			for (Object marker:markerValues) {
				if (!missing[allCount++]) markersForNonMissing[nonmissingCount++] = marker;
			}
			
			if (myData.isMarkerDiscrete(s)) {
				ArrayList<Object> markerIds = new ArrayList<Object>();
				int[] levels = ModelEffectUtils.getIntegerLevels(markersForNonMissing, markerIds);
				snp.alleles = markerIds;
				
				if (isNested) {
					//not implemented yet
					markerEffect = new FactorModelEffect(levels, true, snp);
				} else {
					markerEffect = new FactorModelEffect(levels, true, snp);
				}
			} else {
				//int n = markerValues.length;
				double[] cov = new double[numberNotMissing];
				for (int i = 0; i < numberNotMissing; i++) cov[i] = ((Double) markersForNonMissing[i]).doubleValue();

				if (isNested) {
					CovariateModelEffect cme = new CovariateModelEffect(cov);
					markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
					markerEffect.setID(snp);
				} else {
					markerEffect = new CovariateModelEffect(cov, snp);
				}
			}


			plm.testNewModelEffect(markerEffect);
			double modelss = plm.getModelSS();

			currentModel.add(markerEffect); //Temporary; Remove if wrong
			SweepFastLinearModel sflm2 = new SweepFastLinearModel(currentModel, y);
			
			int n = numberNotMissing;
			//Calculate the BIC
			double [] errorss = sflm2.getResidualSSdf();
			double [] modeldf = sflm2.getFullModelSSdf();
			double pForBIC  = modeldf[1]; 
			double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ;
			double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;
			if(s==2) System.out.println("The value of errorss in forwardStep() is: " + errorss[0]);

			//Calculate the mBIC
			int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
			double pForMbic = modeldf[1];
			double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
			double mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
					(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) ));

			switch(modelType){
			case pvalue:
				test = modelss > bestss;
				break;
			case bic:
				test = bic < bestbic;
				break;
			case mbic:
				test = mbic < bestmbic;
				break;
			case aic:
				test = aic < bestaic;
				break;
			}
			if (test) {
				bestmbic = mbic;
				bestbic = bic;
				bestaic = aic;
				bestss = modelss;
				besteffect = markerEffect;
			}
			currentModel.remove(markerEffect); //Temporary; Remove if wrong
		}

		//if the p-value for the select SNP is less than the enter limit, add it to the model and recalculate the model solution
		plm.testNewModelEffect(besteffect);
		double[] Fp = plm.getFp();
		if(modelType == MODEL_TYPE.pvalue){
			if ( Fp[1] < enterlimit) {
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else { 
				return false; 
			}     
		} else if(modelType == MODEL_TYPE.bic){
			if(bestbic < globalbestbic){
				globalbestbic = bestbic;
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			} 
		} else if(modelType == MODEL_TYPE.mbic){
			if(bestmbic < globalbestmbic){
				globalbestmbic = bestmbic;
				//System.out.println("***********The value of globalbestmbic at the end of forwardStep() is " + globalbestmbic);
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			}
		} else if(modelType == MODEL_TYPE.aic){
			if(bestaic < globalbestaic){
				globalbestaic = bestaic;
				//System.out.println("***********The value of globalbestmbic at the end of forwardStep() is " + globalbestmbic);
				currentModel.add(besteffect);
				if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
				return true;
			} else{
				return false;
			} 
		}  else{
			return false;
		}  


	}

	public boolean backwardStep() {
		int numberOfTerms = currentModel.size();
		if (numberOfTerms <= numberOfBaseModelEffects) return false;
		double bestbic = Double.MAX_VALUE;
		double bestmbic = Double.MAX_VALUE;
		double bestaic = Double.MAX_VALUE;

		int n = y.length;
		int numberOfSites = myData.getNumberOfMarkers();

		SweepFastLinearModel sflm0 = new SweepFastLinearModel(currentModel, y);

		//find the model term (snps only) with the largest p-value
		double maxp = 0;
		double minF= -1;
		int maxterm = 0;
		ModelEffect worstMarkerEffect = null;
		double[] errorssdf = sflm0.getResidualSSdf();

		for (int t = numberOfBaseModelEffects; t < numberOfTerms; t++) {
			double bic;
			double mbic;
			double aic;
			ModelEffect markerEffect = null;
			double[] termssdf = sflm0.getIncrementalSSdf(t);
			double F = termssdf[0]/termssdf[1]/errorssdf[0]*errorssdf[1];
			double p;
			if(modelType != MODEL_TYPE.pvalue){

				ModelEffect meTest1 = currentModel.remove(t);
				SNP snpRemoved = (SNP) meTest1.getID(); 
				System.out.println("CORRECT The SNP just removed is: " + snpRemoved);

				SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y); 
				//Calculate the BIC
				double [] errorss = sflm.getResidualSSdf();
				double [] modeldf = sflm.getFullModelSSdf();
				double pForBIC  = modeldf[1]; 
				bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
				aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

				//Calculate the mBIC
				int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
				double pForMbic = modeldf[1];
				double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
				mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
						(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) ));


				currentModel.add(t, meTest1);

				sflm = new SweepFastLinearModel(currentModel, y); 


			} else{
				bic = Double.MAX_VALUE;
				mbic = Double.MAX_VALUE;
				aic = Double.MAX_VALUE;
			}

			try{
				p = LinearModelUtils.Ftest(F, termssdf[1], errorssdf[1]);  
			} catch (Exception e){
				p = Double.NaN;
			}


			switch(modelType){
			case pvalue:
				try {
					if (p > maxp) {
						maxterm = t;
						maxp = p;
						minF = F; 
					}
				} catch(Exception e){

				}
				break;
			case bic:
				if(bic < bestbic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;  //Actually, this should be "minterm" becasue we are finding the model with the minimum BIC, but I am keeping this as "maxterm" for simpilicity of the calculations
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F; 
				}
				break;
			case mbic:
				if(mbic < bestmbic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F;
				}
				break;
			case aic:
				if(aic < bestaic){
					bestbic = bic;
					bestaic = aic;
					bestmbic = mbic;
					maxterm = t;
					worstMarkerEffect = markerEffect;
					maxp = p;
					minF = F;
				}
				break;
			} 

		}


		//if that maxp is >= exitlimit, then remove maxterm from the model, recalculate the model, and return true;
		switch(modelType){
		case pvalue:
			test = maxp >= exitlimit;
			break;
		case bic:
			test = bestbic < globalbestbic;
			break;
		case mbic:
			test = bestmbic < globalbestmbic;
			break;
		case aic:
			test = bestaic < globalbestaic;
			break;
		}


		if ((test)&&(maxterm != 0)) {
			ModelEffect me = currentModel.remove(maxterm);
			globalbestbic = bestbic;
			globalbestmbic = bestmbic;
			globalbestaic = bestaic;
			return true;
		}

		return false;
	}

        public void runPermutationTestNoMissingData() {
                ArrayList<double[]> permutedData = new ArrayList<double[]>();
                DoubleMatrix PvalueVectorAcrossMarkers = null;
                
                int indexOfThreshold = (int) (alpha*numberOfPermutations);
                        
                int numberOfSites = myData.getNumberOfMarkers();
                DoubleMatrix[][] Xmatrices;
		System.out.println("-----------------Running permutations...----------------");

		int numberOfObs = y.length;
		double totalSS = 0;
		for (int i = 0; i < numberOfObs; i++) totalSS += y[i] * y[i];
		
		ArrayList<double[]> baseModelSSdf = new ArrayList<double[]>();
                double[] permTotalSS = new double[numberOfPermutations];
		
		/*
                 * I commented out this code because this step is already done in StepwiseOLSModelFitter
                int[] mean = new int[numberOfObs];
		FactorModelEffect meanME = new FactorModelEffect(mean, false, "mean");
		FactorModelEffect popME = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops), true);
		ArrayList<ModelEffect> effectList = new ArrayList<ModelEffect>();
		effectList.add(meanME);
		effectList.add(popME);
		DoubleMatrix[][] Xmatrices = new DoubleMatrix[1][3];
		Xmatrices[0][0] = meanME.getX();
		Xmatrices[0][1] = popME.getX();
                */
                int[] mean = new int[numberOfObs];
		FactorModelEffect meanME = new FactorModelEffect(mean, false, "mean");
                int columnNumber = 2;
                if(covariateList != null){
                    columnNumber = columnNumber + covariateList.size();
                    //(any other covariates)+one marker = currentModel.size()+1
                }if(factorList != null){
                    columnNumber = columnNumber + factorList.size();                
                }
                
                Xmatrices = new DoubleMatrix[1][columnNumber];//Intercept+family term+                
	
		Xmatrices[0][0] = meanME.getX();
                int indexForFactorLevelinXmatrices = 0;
                if(covariateList != null){
                    indexForFactorLevelinXmatrices = covariateList.size(); //This is not equal to
                    // covariateList.size()+1 because java starts indices at 0; not 1
                    for(int i=0;  i < covariateList.size(); i++){
                        CovariateModelEffect cme = new CovariateModelEffect(covariateList.get(i), 
                                                  myData.getCovariateName(i));
                        Xmatrices[0][i+1] = cme.getX();
                    }
                    //(any other covariates)+one marker = currentModel.size()+1
                }if(factorList != null){
			for (int f = 0; f < factorList.size(); f++) {
				ArrayList<String> ids = new ArrayList<String>();
				int[] levels = ModelEffectUtils.getIntegerLevels(factorList.get(f), ids);
				FactorModelEffect fme = new FactorModelEffect(levels, true, new Object[]{myData.getFactorName(f), ids});
				if (isNested && f == nestingFactorIndex) {
					nestingEffect = fme;
					nestingFactorNames = ids; 
				}
                                Xmatrices[0][f+indexForFactorLevelinXmatrices+1] = fme.getX();
			}                                   
                }            
                
                SweepFastLinearModel baseSflm = new SweepFastLinearModel(currentModel, y);

                        
                DoubleMatrix resAsDoubleMatrix = baseSflm.getResiduals();
                DoubleMatrix predAsDoubleMatrix = baseSflm.getPredictedValues();
 		double[] res = new double[numberOfObs];
                double[] pred = new double[numberOfObs];
                for(int i = 0; i < numberOfObs; i++){
                    res[i] = resAsDoubleMatrix.get(i,0);
                    pred[i] = predAsDoubleMatrix.get(i,0);
                }
//                System.out.println("resAsDoubleMatrix are: " + resAsDoubleMatrix.toString());
//                System.out.println("predAsDoubleMatrix are: " + predAsDoubleMatrix.toString());
  
//                System.out.println("residuals[0] are: " + res[0]);
//                System.out.println("predicted[0] values are: " + pred[0]);
//TEMPORARY CODE
/*		double[][] permsFromSAS = ImportAsTextFile.ImportTwoDimensionalArrayAsText(numberOfObs,numberOfPermutations, 
                                           "G:\\Lipka_Hal\\Java_Coding_Project\\Permutation_Code\\Validatation\\Input_Files_for_R\\Permuted.Kernel.Color.Phenotypes.from.R.Set5.txt");
                
                System.out.println("permsFromSAS[0][0]: " + permsFromSAS[0][0]);
                System.out.println("permsFromSAS[0][999]: " + permsFromSAS[0][999]);
                System.out.println("permsFromSAS[1683][0]: " + permsFromSAS[1683][0]);
                System.out.println("permsFromSAS[1683][999]: " + permsFromSAS[1683][999]);
                
                
                double[] minP = new double[numberOfPermutations];
                DoubleMatrix minPAsDoubleMatrix = null;
		for (int p = 0; p < numberOfPermutations; p++) {
			minP[p] = 1;
			double[] pdata = new double[numberOfObs];
			//System.arraycopy(res, 0, pdata, 0, numberOfObs);
			//BasicShuffler.shuffleEd(pdata);
                        //System.out.println("pdata");
                        for(int i = 0; i < numberOfObs; i++){
                            pdata[i] = permsFromSAS[i][p];
                            //System.out.println(pdata[i]);
                        } //Add the predicted value back to the residual
                        totalSS = 0;
                        for (int i = 0; i < numberOfObs; i++) totalSS += pdata[i] * pdata[i];
                        permTotalSS[p] = totalSS;
                        //System.out.println("pdata[0] after pred added back are: " + pdata[0]);
			permutedData.add(pdata); 
			
			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, pdata);
			baseModelSSdf.add(sflm.getFullModelSSdf());
		}    */          
                
 //END TEMPORARY CODE            
                
                
		//permute data
                //This for loop below is temporarily commented out. UNCOMMENT THIS ONCE
                //YOU FIGURED OUT WHY SAS AND TASSEL ARE GIVING DIFFERENT ANSWERS
		double[] minP = new double[numberOfPermutations];
                
                DoubleMatrix minPAsDoubleMatrix = null;
		for (int p = 0; p < numberOfPermutations; p++) {
			minP[p] = 1;
			double[] pdata = new double[numberOfObs];
			System.arraycopy(res, 0, pdata, 0, numberOfObs);
			BasicShuffler.shuffle(pdata);

                        //System.out.println("pdata");
                        for(int i = 0; i < numberOfObs; i++){
                            pdata[i] = pdata[i] + pred[i];
                            //System.out.println(pdata[i]);
                        } //Add the predicted value back to the residual
                        totalSS = 0;
                        for (int i = 0; i < numberOfObs; i++) totalSS += pdata[i] * pdata[i];
                        permTotalSS[p] = totalSS;
                        //System.out.println("pdata[0] after pred added back are: " + pdata[0]);
			permutedData.add(pdata); 
			
			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, pdata);
			baseModelSSdf.add(sflm.getFullModelSSdf());
		}
		
                
		for (int m = 0; m < numberOfSites; m++) {
			if (m % 50 == 0) System.out.println("Testing marker " + m);                   
  			ModelEffect markerEffect = null;                      
			SNP snp = new SNP(myData.getMarkerName(m), new Chromosome(myData.getLocusName(m)), (int) myData.getMarkerChromosomePosition(m), m);;
			Object[] markerValues = myData.getMarkerValue(currentPhenotypeIndex, m);
			Object[] markersForNonMissing = new Object[numberNotMissing];
			int allCount = 0;
			int nonmissingCount = 0;
			for (Object marker:markerValues) {
				if (!missing[allCount++]) markersForNonMissing[nonmissingCount++] = marker;
			}
                        
        			if (myData.isMarkerDiscrete(m)) {
				ArrayList<Object> markerIds = new ArrayList<Object>();
				int[] levels = ModelEffectUtils.getIntegerLevels(markersForNonMissing, markerIds);
				snp.alleles = markerIds;
				
				if (isNested) {
					//not implemented yet
					markerEffect = new FactorModelEffect(levels, true, snp);
				} else {
					markerEffect = new FactorModelEffect(levels, true, snp);
				}
			} else {
				//int n = markerValues.length;
				double[] cov = new double[numberNotMissing];
				for (int i = 0; i < numberNotMissing; i++) cov[i] = ((Double) markersForNonMissing[i]).doubleValue();

				if (isNested) {
					CovariateModelEffect cme = new CovariateModelEffect(cov);
					markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
					markerEffect.setID(snp);
				} else {
					markerEffect = new CovariateModelEffect(cov, snp);
				}
			}
 
                        /*Begin old code that may be commented out
                        float[] snpscore = snp.score;
			double[] dscore = new double[numberOfObs];
			for (int i = 0; i < numberOfObs; i++) {
				dscore[i] = snpscore[i];
			}
			
			NestedCovariateModelEffect ncme = new NestedCovariateModelEffect(dscore, popME);
                        */
                        
                        //System.out.println("columnNumber - 1: " + (columnNumber-1));
			Xmatrices[0][columnNumber-1] = markerEffect.getX(); //columnNumber-1 because java
                        //starts counting things from 0
			DoubleMatrix X = DoubleMatrixFactory.DEFAULT.compose(Xmatrices);

                        //Temporary Code
                      /*  ArrayList<double[]> baseModelSSdfTest = new ArrayList<double[]>();
                        SweepFastLinearModel sflmTest = new SweepFastLinearModel(currentModel, y);
                        baseModelSSdfTest.add(sflmTest.getFullModelSSdf());*/
                        //End temporary code 
                        
			currentModel.add(markerEffect);
			SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);         
			double[] modelSSdf = sflm.getFullModelSSdf();
			DoubleMatrix G = sflm.getInverseOfXtX();
                        
                        //Temporary code
 
             /*           DoubleMatrix ypermTemp = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, y);
                        DoubleMatrix XtyTemp = X.crossproduct(ypermTemp);
                        double[] reducedSSdfTemp = sflmTest.getFullModelSSdf();
                        double fullSSTemp = XtyTemp.crossproduct(G.mult(XtyTemp)).get(0, 0);
                        System.out.println("fullSSTemp: " + fullSSTemp);
                        double fulldfTemp = modelSSdf[1];
                        System.out.println("modelSSdf[1] " + modelSSdf[1]);
                        double markerSSTemp = fullSSTemp - reducedSSdfTemp[0];
                        System.out.println("markerSS " + markerSSTemp);
                        double markerdfTemp = fulldfTemp - reducedSSdfTemp[1];
                        System.out.println("markerdf " + markerdfTemp);
                        double errorSSTemp = totalSS - fullSSTemp;
                        System.out.println("errorSS " + errorSSTemp);
                        double errordfTemp = numberOfObs - fulldfTemp;
                        System.out.println("errordf " + errordfTemp);
                        double Ftemp = markerSSTemp / markerdfTemp / errorSSTemp * errordfTemp;
                        System.out.println("Ftemp " + Ftemp);*/

                        //Resuts match PROC GLM exaclty
                        //End temporary code
			
			for (int p = 0; p < numberOfPermutations; p++) {
				double[] pdata = permutedData.get(p);
                                DoubleMatrix yperm = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, pdata);
                                totalSS = permTotalSS[p];
                                
				DoubleMatrix Xty = X.crossproduct(yperm);
				double[] reducedSSdf = baseModelSSdf.get(p);
                                //System.out.println("reducedSSdf[0]: " + reducedSSdf[0]);
                                //System.out.println("reducedSSdf[1]: " + reducedSSdf[1]);
				double fullSS = Xty.crossproduct(G.mult(Xty)).get(0, 0);
                                //System.out.println("fullSS: " + fullSS);
				double fulldf = modelSSdf[1];
                                //System.out.println("fulldf: " + fulldf);
				double markerSS = fullSS - reducedSSdf[0];
                                //System.out.println("markerSS: " + markerSS);
				double markerdf = fulldf - reducedSSdf[1];
                                //System.out.println("markerdf: " + markerdf);
				double errorSS = totalSS - fullSS;
                                //System.out.println("errorSS: " + errorSS);
				double errordf = numberOfObs - fulldf;
                                //System.out.println("errordf: " + errordf);
				double F = markerSS / markerdf / errorSS * errordf;
                                //System.out.println("F: " + F);
				double probF;
				try {
					probF = LinearModelUtils.Ftest(F, markerdf, errordf);
				} catch(Exception e) {
					probF = 1;
				}
				//System.out.println("probF: "+probF);
				//minP[p] = Math.min(minP[p], probF);
                                minP[p] = probF;
                                //System.out.println("minP[p]:"+ minP[p]);
				minPAsDoubleMatrix = DoubleMatrixFactory.DEFAULT.make(minP.length, 1, minP);
			}
                        if(m == 0){
                            PvalueVectorAcrossMarkers = minPAsDoubleMatrix;
                        }else{
                            PvalueVectorAcrossMarkers = PvalueVectorAcrossMarkers.concatenate(minPAsDoubleMatrix, false);//Change this to a DoubleMatrix; transpose it
                            //for the next step                            
                        }

			currentModel.remove(markerEffect);
		}
		
                System.out.println(PvalueVectorAcrossMarkers.toString());
                for(int p=0; p<numberOfPermutations; p++){
                    double minFromOnePermutation = 1;
                    for(int m=0; m < numberOfSites; m++){
                       //System.out.println("Marker:" + m + ":Permutation:"+ p +":"+ PvalueVectorAcrossMarkers.get(p,m));
                       minFromOnePermutation = Math.min(PvalueVectorAcrossMarkers.get(p,m), minFromOnePermutation);
                    }
                   //obtain the minimum P-value across each permutation
                    //System.out.println("--------------p is " + p);
                    minPvalues[p] = minFromOnePermutation;
                    //System.out.println(" minPvalues[p] " +  minPvalues[p]);
                  //PvalueVectorAcrossMarkers.
                }
		
                //System.out.println("minPvalues[0]: " + minPvalues[0]);
                //System.out.println("minPvalues[1]: " + minPvalues[1]);                
                
                //sort the P-values from smallest to largest
                Arrays.sort(minPvalues);

                //System.out.println("minPvalues[0] after permuting: " + minPvalues[0]);
                //System.out.println("minPvalues[1] after permuting: " + minPvalues[1]); 
                
                //select the (alpha*nperm)th smallest P-value
                enterlimit = minPvalues[indexOfThreshold];
                exitlimit = 2*enterlimit;
                
                //System.out.println("indexOfThreshold: " + indexOfThreshold);
                //System.out.println("enterlimit: " + enterlimit);
                //System.out.println("exitlimit: " + exitlimit); 
		//String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		
                //Turn this into a new method. Put it into the plugger. Have a logical statement to run this 
                //only if permutations are selected.
                //Create a table that has the permuted P-values

                

                
                System.out.println("--------------the Enter limit is " + enterlimit);
                System.out.println("--------------the Exit limit is " + exitlimit); 
                System.out.println("--------------indexOfThreshold is " + indexOfThreshold); 
	}

        
 	public TableReport getPermutationReport() {
		if (numberOfPermutations == 0) return null;
                
                LinkedList<Object[]> pValueReportTable = new LinkedList<Object[]>();;
                //Object[] pValueReportRow = new Object[numberOfPermutations];
                for(int i=0; i < numberOfPermutations; i++){
                    Object[] pValueReportRow = new Object[1];
                    pValueReportRow[0] = new Double(minPvalues[i]);          
                    //System.out.println("minPvalues[i]" + minPvalues[i]);
                    pValueReportTable.add(pValueReportRow);
                }
                
                //pValueReportTable.add(pValueReportRow);
               	Object[][] table = new Object[pValueReportTable.size()][];
		pValueReportTable.toArray(table);                
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = numberOfPermutations + " Permutations for "  + datasetName;
		String[] reportHeader = new String[]{"P-value"};

		return new SimpleTableReport(reportName, reportHeader, table);
	}       
        
        
	public LinkedList<Object[]> createReportRowsFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		int ncol = anovaReportHeader.length;
		LinkedList<Object[]> reportTable = new LinkedList<Object[]>();
		double[] residualSSdf = sflm.getResidualSSdf();
		int n = y.length;
		int numberOfSites = myData.getNumberOfMarkers();


		//Calcualte the BIC
		double [] errorss = sflm.getResidualSSdf();
		double [] modeldf = sflm.getFullModelSSdf();
		double pForBIC  = modeldf[1];
		double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
		double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

		System.out.println("error ss is " + errorss[0]);
		System.out.println("pForBIC is " + pForBIC);
		System.out.println("n is " + n);


		//Calculate the mBIC
		int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
		double pForMbic = modeldf[1];
		double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
		double mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
				(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) ));

		int effectPtr = 0;
		for (ModelEffect me : currentModel) {
			Object[] reportRow = new Object[ncol];
			int ptr = 0;
			reportRow[ptr++] = traitname;

			if (me.getID() instanceof SNP) {
				SNP snp = (SNP) me.getID();
				reportRow[ptr++] = snp.name;
				reportRow[ptr++] = snp.locus.getName();
				reportRow[ptr++] = Integer.toString(snp.position);
			} else {
				if (me.getID() instanceof String) reportRow[ptr++] = me.getID().toString();
				else if (me instanceof FactorModelEffect) reportRow[ptr++] = ((Object[]) me.getID())[0].toString();
				else reportRow[ptr++] = me.getID().toString();
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			double[] effectSSdf = sflm.getMarginalSSdf(effectPtr);
			double ms = effectSSdf[0] / effectSSdf[1];
			double Fval = ms / residualSSdf[0] * residualSSdf[1];
			double pval;
			try {
				pval = LinearModelUtils.Ftest(Fval, effectSSdf[1], residualSSdf[1]);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			reportRow[ptr++] = new Integer((int) effectSSdf[1]);
			reportRow[ptr++] = new Double(effectSSdf[0]);
			reportRow[ptr++] = new Double(ms);
			reportRow[ptr++] = new Double(Fval);
			reportRow[ptr++] = new Double(pval);
			reportRow[ptr++] = new Double(bic);
			reportRow[ptr++] = new Double(mbic);
			reportRow[ptr++] = new Double(aic);
			reportTable.add(reportRow);
			effectPtr++;
		}
		int ptr = 0;
		Object[] reportRow = new Object[ncol];
		reportRow[ptr++] = traitname;
		reportRow[ptr++] = "Error";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = new Integer((int) residualSSdf[1]);
		reportRow[ptr++] = new Double(residualSSdf[0]);
		reportRow[ptr++] = new Double(residualSSdf[0]/residualSSdf[1]);
		reportRow[ptr++] = new Double(Double.NaN);
		reportRow[ptr++] = new Double(Double.NaN);
		reportTable.add(reportRow);

		return reportTable;
	}

        public LinkedList<Object[]> createReportRowsFromCurrentModelAfterScanCI(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		int ncol = anovaReportWithCIHeader.length;
		LinkedList<Object[]> reportTable = new LinkedList<Object[]>();
		double[] residualSSdf = sflm.getResidualSSdf();
		int n = y.length;
		int numberOfSites = myData.getNumberOfMarkers();
                int firstMarkerIndex = 1;
                if (factorList != null){
                    firstMarkerIndex = firstMarkerIndex+factorList.size();
                }
                if (covariateList != null){
                    firstMarkerIndex = firstMarkerIndex+covariateList.size();
                }  

		//Calcualte the BIC
		double [] errorss = sflm.getResidualSSdf();
		double [] modeldf = sflm.getFullModelSSdf();
		double pForBIC  = modeldf[1];
		double bic = (n*Math.log(errorss[0]/n)) + (pForBIC*Math.log(n)) ; 
		double aic = (n*Math.log(errorss[0]/n)) + (2*pForBIC) ;

		System.out.println("error ss is " + errorss[0]);
		System.out.println("pForBIC is " + pForBIC);
		System.out.println("n is " + n);


		//Calculate the mBIC
		int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
		double pForMbic = modeldf[1];
		double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
		double mbic = (n*Math.log(errorss[0])) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
				(2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) ));

		int effectPtr = 0;
		for (ModelEffect me : currentModel) {
			Object[] reportRow = new Object[ncol];
			int ptr = 0;
			reportRow[ptr++] = traitname;

			if (me.getID() instanceof SNP) {
				SNP snp = (SNP) me.getID();
				reportRow[ptr++] = snp.name;
				reportRow[ptr++] = snp.locus.getName();
				reportRow[ptr++] = Integer.toString(snp.position);
			} else {
				if (me.getID() instanceof String) reportRow[ptr++] = me.getID().toString();
				else if (me instanceof FactorModelEffect) reportRow[ptr++] = ((Object[]) me.getID())[0].toString();
				else reportRow[ptr++] = me.getID().toString();
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			double[] effectSSdf = sflm.getMarginalSSdf(effectPtr);
			double ms = effectSSdf[0] / effectSSdf[1];
			double Fval = ms / residualSSdf[0] * residualSSdf[1];
			double pval;
			try {
				pval = LinearModelUtils.Ftest(Fval, effectSSdf[1], residualSSdf[1]);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			reportRow[ptr++] = new Integer((int) effectSSdf[1]);
			reportRow[ptr++] = new Double(effectSSdf[0]);
			reportRow[ptr++] = new Double(ms);
			reportRow[ptr++] = new Double(Fval);
			reportRow[ptr++] = new Double(pval);
			reportRow[ptr++] = new Double(bic);
			reportRow[ptr++] = new Double(mbic);
			reportRow[ptr++] = new Double(aic);
                        //System.out.println("The value of effectPtr is "+effectPtr);

                        if (me.getID() instanceof SNP) {
				reportRow[ptr++] = new String(myData.getMarkerName(theUpperAndLowerBound[effectPtr-firstMarkerIndex][0]));
				reportRow[ptr++] = new String(myData.getMarkerName(theUpperAndLowerBound[effectPtr-firstMarkerIndex][1]));
				
			} else {
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			reportTable.add(reportRow);
			effectPtr++;
		}
		int ptr = 0;
		Object[] reportRow = new Object[ncol];
		reportRow[ptr++] = traitname;
		reportRow[ptr++] = "Error";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = new Integer((int) residualSSdf[1]);
		reportRow[ptr++] = new Double(residualSSdf[0]);
		reportRow[ptr++] = new Double(residualSSdf[0]/residualSSdf[1]);
		reportRow[ptr++] = new Double(Double.NaN);
		reportRow[ptr++] = new Double(Double.NaN);
		reportTable.add(reportRow);

		return reportTable;
	}

	public Datum createReportFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		LinkedList<Object[]> reportTable = createReportRowsFromCurrentModel(sflm);
		Object[][] table = new Object[reportTable.size()][];
		reportTable.toArray(table);
		String reportName = "ANOVA for " + traitname + ", " + datasetName;
		TableReport tr = new SimpleTableReport(reportName, anovaReportHeader, table);
		return new Datum(reportName, tr, "");
	}

	public Datum createReportFromCurrentModelWithCI(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		LinkedList<Object[]> reportTable = createReportRowsFromCurrentModelAfterScanCI(sflm);
		Object[][] table = new Object[reportTable.size()][];
		reportTable.toArray(table);
		String reportName = "ANOVA With CI for " + traitname + ", " + datasetName;
		TableReport tr = new SimpleTableReport(reportName, anovaReportWithCIHeader, table);
		return new Datum(reportName, tr, "");
	}
        
	public void appendAnovaResults(SweepFastLinearModel sflm) {
		resultRowsAnova.addAll(createReportRowsFromCurrentModel(sflm));
	}

        public void appendAnovaResultsWithCI(SweepFastLinearModel sflm) {
            resultRowsAnovaWithCI.addAll(createReportRowsFromCurrentModelAfterScanCI(sflm));
	}

	public TableReport getAnovaReport() {
		String reportName = "ANOVA table for " + datasetName;
		Object[][] table = new Object[resultRowsAnova.size()][];
		resultRowsAnova.toArray(table);
		return new SimpleTableReport(reportName, anovaReportHeader, table);
	}

	public TableReport getAnovaReportWithCI() {
		String reportName = "ANOVA table with CI scan for " + datasetName;
		Object[][] table = new Object[resultRowsAnovaWithCI.size()][];
		resultRowsAnovaWithCI.toArray(table);
		return new SimpleTableReport(reportName, anovaReportWithCIHeader, table);
	}
        
	public TableReport getMarkerEffectReport() {
		if (rowsSiteEffectTable.size() == 0) return null;
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = "Marker effects for " + datasetName;
		String[] reportHeader = new String[]{"Trait","Snp","Locus","Position","Within","Estimate"};
		Object[][] table = new Object[rowsSiteEffectTable.size()][];
		rowsSiteEffectTable.toArray(table);
		return new SimpleTableReport(reportName, reportHeader, table);
	}

	public TableReport getMarkerEffectReportWithCI() {
		if (rowsSiteEffectTableWithCI.size() == 0) return null;
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = "Marker effects for " + datasetName;
		String[] reportHeader = new String[]{"Trait","Snp","Locus","Position","Within","Estimate"};
		Object[][] table = new Object[rowsSiteEffectTableWithCI.size()][];
		rowsSiteEffectTableWithCI.toArray(table);
		return new SimpleTableReport(reportName, reportHeader, table);
	}

	public void appendSiteEffectEstimates(SweepFastLinearModel sflm) {

		//columns are trait, snp name, locus, position, factor value, estimate
		int nBaseModelParameters = 0;
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		for (int i = 0; i < numberOfBaseModelEffects; i++) {
			nBaseModelParameters += currentModel.get(i).getEffectSize();
		}

		double[] beta = sflm.getBeta();
		int parameterIndex = nBaseModelParameters;
		for (int s = numberOfBaseModelEffects; s < currentModel.size(); s++) {
			ModelEffect me = currentModel.get(s);
			if (me instanceof NestedCovariateModelEffect) {
				NestedCovariateModelEffect ncme = (NestedCovariateModelEffect) me;
				FactorModelEffect fme = ncme.getFactorModelEffect();
				SNP snp = (SNP) ncme.getID();
				int n = nestingFactorNames.size();
				for (int i = 0; i < n; i++) {
					Object[] rowValues = new Object[6];
					int ptr = 0;
					rowValues[ptr++] = traitname;
					rowValues[ptr++] = snp.name;
					rowValues[ptr++] = snp.locus.getName();
					rowValues[ptr++] = new Integer(snp.position);
					rowValues[ptr++] = nestingFactorNames.get(i);
					rowValues[ptr++] = new Double(beta[parameterIndex++]);
					rowsSiteEffectTable.add(rowValues);
				}
			} else if (me instanceof CovariateModelEffect) {
				SNP snp = (SNP) me.getID();
				Object[] rowValues = new Object[6];
				int ptr = 0;
				rowValues[ptr++] = traitname;
				rowValues[ptr++] = snp.name;
				rowValues[ptr++] = snp.locus.getName();
				rowValues[ptr++] = new Integer(snp.position);
				rowValues[ptr++] = "NA";
				rowValues[ptr++] = new Double(beta[parameterIndex++]);
				rowsSiteEffectTable.add(rowValues);
			} else if (me instanceof FactorModelEffect) {

			}
		}
	}

        public void appendSiteEffectEstimatesWithCI(SweepFastLinearModel sflm) {

		//columns are trait, snp name, locus, position, factor value, estimate
		int nBaseModelParameters = 0;
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		for (int i = 0; i < numberOfBaseModelEffects; i++) {
			nBaseModelParameters += currentModel.get(i).getEffectSize();
		}

		double[] beta = sflm.getBeta();
		int parameterIndex = nBaseModelParameters;
		for (int s = numberOfBaseModelEffects; s < currentModel.size(); s++) {
			ModelEffect me = currentModel.get(s);
			if (me instanceof NestedCovariateModelEffect) {
				NestedCovariateModelEffect ncme = (NestedCovariateModelEffect) me;
				FactorModelEffect fme = ncme.getFactorModelEffect();
				SNP snp = (SNP) ncme.getID();
				int n = nestingFactorNames.size();
				for (int i = 0; i < n; i++) {
					Object[] rowValues = new Object[6];
					int ptr = 0;
					rowValues[ptr++] = traitname;
					rowValues[ptr++] = snp.name;
					rowValues[ptr++] = snp.locus.getName();
					rowValues[ptr++] = new Integer(snp.position);
					rowValues[ptr++] = nestingFactorNames.get(i);
					rowValues[ptr++] = new Double(beta[parameterIndex++]);
					rowsSiteEffectTableWithCI.add(rowValues);
				}
			} else if (me instanceof CovariateModelEffect) {
				SNP snp = (SNP) me.getID();
				Object[] rowValues = new Object[6];
				int ptr = 0;
				rowValues[ptr++] = traitname;
				rowValues[ptr++] = snp.name;
				rowValues[ptr++] = snp.locus.getName();
				rowValues[ptr++] = new Integer(snp.position);
				rowValues[ptr++] = "NA";
				rowValues[ptr++] = new Double(beta[parameterIndex++]);
				rowsSiteEffectTableWithCI.add(rowValues);
			} else if (me instanceof FactorModelEffect) {

			}
		}
	}

	public void setEnterlimits(double[] enterlimits) {
		this.enterlimits = enterlimits;
	}

	public void setExitlimits(double[] exitlimits) {
		this.exitlimits = exitlimits;
	}

	public void setEnterlimit(double enterlimit) {
		this.enterlimit = enterlimit;
	}

	public void setExitlimit(double exitlimit) {
		this.exitlimit = exitlimit;
	}

	public void setMaxNumberOfMarkers(int maxNumberOfMarkers) {
		this.maxNumberOfMarkers = maxNumberOfMarkers;
	}

	public void setNestingEffectIndex(int nestingFactorIndex) {
		this.nestingFactorIndex = nestingFactorIndex;
	}

	public void setNested(boolean isNested) {
		this.isNested = isNested;
	}
        
	public void setNumberOfPermutations(int numberOfPermutations) {
		this.numberOfPermutations = numberOfPermutations;
                minPvalues = new double[this.numberOfPermutations];
}
 
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}
}
