
package net.maizegenetics.pd;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.pal.alignment.BitNucleotideAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


public class PDAnnotation {

private static final String PHYSICAL_POSITIONS = "PhysPos";
private static final String CHROMOSOME = "_C";
private static final String GWAS_TRAIT = "_GT";
private static final String MINOR_ALLELE_FREQUENCY = "MAF";
private static final String HAS_DATA = "HasData"; // summary index where if any trait has a value at that location, value is set to 1

    private boolean myIsSBit = true;

    private String[] chromosomes = new String[] { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10" };
    private String hapMapFile_prefix = "/maizeHapMapV2_B73RefGenV2_201203028_";
    private String hapMapFile_suffix = ".hmp.txt";
    
    private static final int PHYSICAL_POSITION_COLUMN = 0;
    private static final int MINOR_ALLELE_FREQUENCY_COLUMN = 1;
    private static final int COLUMN_OFFSET = 1; // one for physical position column and another for minor allele frequency column



    public void init(){
        String hapMapPath = "/local/workdir/dek29/pd/HapMapV2RefGenV2/";
        File aHapMapDir = new File(hapMapPath);
        String pathGwas = "/local/workdir/dek29/pd/gwas_results/";
        
        File aGwasDir = new File(pathGwas);

        PDAnnotation pd = new PDAnnotation();

        pd.loadAllChromosomes(aGwasDir, aHapMapDir);
    }


    public void loadAllChromosomes(File gwasDirIn, File hapMapDir){
        IHDF5Writer writer = HDF5Factory.open("chromosomes10.h5");
        // Define the block size as 10 x 10.
// writer.createIntMatrix("results", 10, 10);

        int[] hasData = null;
        for(int i =0; i<chromosomes.length; i++){
            String chromosomeFile = hapMapDir + hapMapFile_prefix + chromosomes[i] + hapMapFile_suffix;
            

            BitNucleotideAlignment bna = readFile(chromosomeFile);
            
            int[] alignmentPhysPos = bna.getPhysicalPositions();
            hasData = new int[alignmentPhysPos.length];
            float[] maf = new float[alignmentPhysPos.length];
            for(int j=0; j< maf.length; j++){
             maf[j] = (float)bna.getMinorAlleleFrequency(j);
            }
           
           // todo: talk to Ed about what he meant here
            // get the alleles

          //delete bna
            bna = null;
            
          //write positions to hdf "pos"+chromosome
           writer.createIntArray(PDAnnotation.PHYSICAL_POSITIONS + PDAnnotation.CHROMOSOME + (i + 1), alignmentPhysPos.length);
           writer.writeIntArray(PDAnnotation.PHYSICAL_POSITIONS, alignmentPhysPos);
           
              
           //write alleles to hdf "allele"+chromosome
           // which version? String[][] ?
           
           // write minor allele frequencies
           writer.createFloatArray(PDAnnotation.MINOR_ALLELE_FREQUENCY, maf.length);
           writer.writeFloatArray(PDAnnotation.MINOR_ALLELE_FREQUENCY, maf);
            
           FolderParser fp = new FolderParser(gwasDirIn);
           String[] traits = fp.getAllTraits();

           int physPosIndex = 1;
           int pValueIndex = 6;
           for(int j = 0; j < traits.length; j++){
         File gwasFile = fp.getFile(chromosomes[i], traits[j]);
         System.out.println(gwasFile.toString());
         // pull out the physical location and p-value
         float[][] fVals = parseGWASFloatValues( gwasFile, physPosIndex, pValueIndex, true);
        
         if(fVals == null) continue;
        
         // find the matches
         float[] pVals = findMatches(alignmentPhysPos, fVals, hasData);

               // if the file is empty, move on to next trait
               if(fVals == null) continue;

               String dataSetName = "GW_c"+(i+1)+"_T"+traits[j];
               writer.createFloatArray(dataSetName, pVals.length);
               writer.writeFloatArray(dataSetName, pVals);

           } // end of traits loop
           
           // write
           writer.createIntArray(PDAnnotation.HAS_DATA + PDAnnotation.CHROMOSOME + (i + 1), hasData.length);
           writer.writeIntArray(PDAnnotation.PHYSICAL_POSITIONS, hasData);
           
        } // end of chromosome loop
        writer.close();
    }

    
    private float[] findMatches(int[] alignmentPhysPositions, float[][] gwasResults, int[] hasData){
     float[] results = new float[alignmentPhysPositions.length];
    
     int lastIndex = 0;
     for(int i = 0; i < alignmentPhysPositions.length; i++){
     int anIndex = lastIndex;
     while(anIndex<gwasResults[0].length ){
     if(((int)gwasResults[0][anIndex]) == alignmentPhysPositions[i])
     {
     results[i] = gwasResults[1][anIndex];
     hasData[i] = 1;
     System.out.println("\t" + alignmentPhysPositions[i] + "\t" + (int)gwasResults[0][anIndex] + "\t" + gwasResults[1][anIndex]);
     lastIndex = anIndex;
     }
     anIndex++;
     }
     }
     return results;
    }
    
    /**
* Load all results gwas results for a given chromosome
* @param bna
* @param gwasDirIn
* @param chromosomeIn
* @return
*/
    private int[][] loadGWAS(BitNucleotideAlignment bna, File gwasDirIn, String chromosomeIn){

        FolderParser fp = new FolderParser(gwasDirIn);
        String[] traits = fp.getAllTraits();

        int[] alignmentPhysPos = bna.getPhysicalPositions();

        int chrBaseCount = alignmentPhysPos.length;

// float[] maf = new float[chrBaseCount];
// for(int i=0; i< maf.length; i++){
// maf[i] = (float)bna.getMinorAlleleFrequency(i);
// }

        int[][] result = new int[chrBaseCount][traits.length + PDAnnotation.COLUMN_OFFSET];

        // add the alignment's physical positions to the (float) data array
        try{
            addNewColumns(alignmentPhysPos, result, PDAnnotation.PHYSICAL_POSITION_COLUMN);
        }catch (Exception e ) { /* ignore */ }

        // when starting a new chromosome, write out the results of the last chromosome
        // write 2^16 bases as a block in HDF5
        for(int j=0; j<traits.length; j++){

            File gwasFile = fp.getFile(chromosomeIn, traits[j]);
            System.out.println("Loading trait: " + j + " " + traits[j]);
            int fieldOfInterest = 1;
            int[]intVals = parseGWASIntValues(gwasFile, fieldOfInterest, true );

            // if the file is empty, move on to next trait
            if(intVals == null) continue;

            // get the gwas value(s) of interest in the appropriate column
            int gwasPos = 6;

            int[] vals = new int[alignmentPhysPos.length];
            //vals = parseGWASFloatValues( gwasFile, gwasPos, true);

            // add the gwas values to the final (float) data array
            try{
                addNewColumns(vals, result, j+ COLUMN_OFFSET, intVals);;
            }catch (Exception e) {
                System.err.println("Problem adding values to results. File: " + gwasFile.toString() + " trait: " + traits[j]);
            }
        }

        bna = null;
        fp = null;
        return result;// all gwas trait results for that chromosome
    }


    /**
*
* @param aFile
* @param fieldOfInterest
* @param hasHeader If true, skips the first line in the file
* @return
*/
    private int[] parseGWASIntValues(File aFile, int fieldOfInterest, boolean hasHeader){

        int[] fieldValue = null;
        StringBuffer sb = null;
        try {
            BufferedReader in = new BufferedReader(new FileReader(aFile));
            String str;
            while ((str = in.readLine()) != null) {
                if (sb == null) {
                    sb = new StringBuffer(str + "\n");
                } else {
                    sb.append(str.trim() + "\n");
                }
            }
            if (sb == null) { // if the StringBuffer was never initialized
                // then there are no records in a file
                String msg = "File " + aFile.getName() + " could not be loaded. Could it be empty?";
                return null;
            } else {

                String[] line = sb.toString().split("\n");
                int lineCount = line.length;

                int lineIndex = 0;
                fieldValue = new int[lineCount];
                if(hasHeader){
                    lineIndex = 1;
                    fieldValue = new int[lineCount-1];
                }

                int range = 0;

                while (lineIndex < lineCount) {
                    String[] val = line[lineIndex].split("\t");
                    int value = -1;
                    String toBeParsed = null;
                    try{
                        toBeParsed = val[fieldOfInterest];
                        if(hasHeader){
                            fieldValue[lineIndex-1] = Integer.parseInt(toBeParsed);
                        }else{
                            fieldValue[lineIndex] = Integer.parseInt(toBeParsed);
                        }
                    }catch(NumberFormatException nfe){
                        System.err.println("File " + aFile.getName() + " int parsing failed: " + toBeParsed + "lineIndex: " + lineIndex);
                    }
                    lineIndex++;
                }
            }
            in.close();
        } catch (IOException ioe) {
            String msg = "Failed while reading in file: " + aFile.getAbsolutePath() + "\n";
            System.err.println(msg);
            ioe.printStackTrace();
        }

        return fieldValue;
    }

    /**
*
* @param input float[] data to be added to the dataFinal
* @param dataFinal float[][] which is the recipient of added data
* @param column column in dataFinal to which the input is added
* @param inputPosition int[] physical positions of the input data.
* @throws Exception input[] must be of the same length as inputPosition
*/
    private void addNewColumns(int[] input, int[][] dataFinal, int column, int[] inputPosition) throws Exception{
        if(input.length != inputPosition.length){ throw new Exception("input[] must be the same length as the inputPostions[]"); }

        if(dataFinal[0].length < column) {
            throw new ArrayIndexOutOfBoundsException("Column does not exist: " + column);
        }

        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < dataFinal.length; i++){
            float alignmentPos = dataFinal[i][0];
            for(int j = 0; j < inputPosition.length; j++){
// if((alignmentPos -50 < (float)inputPosition[j]) && (alignmentPos +50 > (float)inputPosition[j])){ // for range limiting
                if((alignmentPos == (float)inputPosition[j])){
                    dataFinal[i][column] = input[j];
                }
            }
        }
    }

    
    /**
*
* @param input
* @param dataFinal
* @param column
*/
    private void addNewColumns(float[] input, float[][] dataFinal, int column) throws Exception{
        if(input.length != dataFinal.length) throw new Exception("The float[][] dataFinal must be the same length as the input");
        
        if(dataFinal[0].length < column) {
            throw new ArrayIndexOutOfBoundsException("Column does not exist: " + column);
        }
        int lastIndex = 0;
        for(int i = 0; i < dataFinal.length; i++){
            dataFinal[i][column] = input[i];
        }
    }

    /**
*
* @param input
* @param dataFinal
* @param column
*/
    private void addNewColumns(int[] input, int[][] dataFinal, int column) throws Exception{
        if(input.length != dataFinal.length) throw new Exception("The float[][] dataFinal must be the same length as the input");
        
        if(dataFinal[0].length < column) {
            throw new ArrayIndexOutOfBoundsException("Column does not exist: " + column);
        }
        int lastIndex = 0;
        for(int i = 0; i < dataFinal.length; i++){
            dataFinal[i][column] = input[i];
        }
    }

    /**
*
* @param aFile
* @param fieldOfInterest
* @param hasHeader If true, skips the first line in the file
* @return
*/
    private float[][] parseGWASFloatValues(File aFile, int physPositionIndex, int pValIndex, boolean hasHeader){

        float[][] fieldValue = null;
        StringBuffer sb = null;
        try {
            BufferedReader in = new BufferedReader(new FileReader(aFile));
            String str;
            while ((str = in.readLine()) != null) {
                if (sb == null) {
                    sb = new StringBuffer(str + "\n");
                } else {
                    sb.append(str.trim() + "\n");
                }
            }
            if (sb == null) { // if the StringBuffer was never initialized
                // then there are no records in a file
                String msg = "File " + aFile.getName() + " could not be loaded. Could it be empty?";
                return null;
            } else {

                String[] line = sb.toString().split("\n");
                int lineCount = line.length;
                int lineIndex = 0;
                fieldValue = new float[2][lineCount];
                if(hasHeader) {
                    lineIndex=1;
                    fieldValue = new float[2][lineCount-1];
                    
                }
                if(lineCount<2) return null;

                int field = 1; // field containing information of interest
                int range = 0;
                while (lineIndex < lineCount) {
                    String[] val = line[lineIndex].split("\t");
                    int value = -1;
                    try {

                        if (hasHeader) {
                         fieldValue[0][lineIndex - 1] = Float.parseFloat(val[physPositionIndex]);
                            fieldValue[1][lineIndex - 1] = Float.parseFloat(val[pValIndex]);
                            
                        } else {
                            fieldValue[0][lineIndex ] = Float.parseFloat(val[physPositionIndex]); //physical position
                            fieldValue[1][lineIndex] = Float.parseFloat(val[pValIndex]); //

                        }
                        
                    } catch (NumberFormatException nfe) {
                        System.err.println("File " + aFile.getName() + " float parsing failed: " + val[field] + "lineIndex: " + lineIndex);
                    }
                    lineIndex++;
                }
            }
            in.close();
        } catch (IOException ioe) {
            String msg = "Failed while reading in file: " + aFile.getAbsolutePath() + "\n";
            System.err.println(msg);
            ioe.printStackTrace();
        }

        return fieldValue;
    }

    public BitNucleotideAlignment readFile(String inFile){
        return (BitNucleotideAlignment) ImportUtils.readFromHapmap(inFile, myIsSBit, null /*progressListener*/);
    }

    public static void main(String[] args) {
        PDAnnotation p = new PDAnnotation();
        p.init();
    }
}