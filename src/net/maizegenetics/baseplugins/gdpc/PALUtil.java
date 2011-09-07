/*
 * PALUtil.java
 *
 * Created on July 29, 2003, 5:57 PM
 */
package net.maizegenetics.baseplugins.gdpc;

import gov.usda.gdpc.*;

import java.math.BigDecimal;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.StringTokenizer;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AllelePositionBLOBUtils;
import net.maizegenetics.pal.alignment.GdpdmBLOBUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import net.maizegenetics.pal.alignment.SimpleAlignment;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 * @author  terryc
 */
public class PALUtil {

    /** Creates a new instance of PALUtil */
    private PALUtil() {
    }

    /**
     * Converts the given GDPC GenotypeTable to a PAL SimpleAnnotatedAlignment.
     * This defaults to using the first column (index 0) in the genotype table.
     *
     * @param table the genotype table
     *
     * @return the simple annotated alignment
     */
    public static SimpleAlignment createAnnotationAlignment(GenotypeTable table) {
        return createAnnotationAlignment(table, 0);
    }

    /**
     * Converts the given GDPC TaxonParentGroup to a PAL KinshipMatrix.
     *
     * @param group the taxon parent group
     *
     * @return kinship matrix
     */
    public static DistanceMatrix createKinshipMatrix(TaxonParentGroup group) {

        String taxa[] = Convert.convertTaxonParentsToTableOrder(group);
        IdGroup idg = new SimpleIdGroup(taxa);
        int[][] nidTable = Import.threeColumnToNumerical(taxa);
        int[][] sortId = Import.sortID(nidTable);
        int[][] completedTable = Import.tableFiller(sortId, nidTable);
        double[][] kinshipMatrix = Import.kinshipRelation(completedTable);
        DistanceMatrix distanceMatrix = new DistanceMatrix(kinshipMatrix, idg);
        return distanceMatrix;

    }

    /**
     * Converts the given GDPC GenotypeTable to a PAL SimpleAnnotatedAlignment.
     *
     * @param table the genotype table
     * @param genotypeExpIndex Genotype Experiment index to use.
     *
     * @return the simple annotated alignment
     */
    public static SimpleAlignment createAnnotationAlignment(GenotypeTable table, int genotypeExpIndex) {

        SimpleAlignment theNewAlignment = null;
        // int totalGenes = dataTable.numColumns();
        int totalTaxa = table.numRows();
        int nonNullTaxa = 0;
        for (int j = 0; j < totalTaxa; j++) {
            if (table.get(genotypeExpIndex, j) != null) {
                nonNullTaxa++;
            }
        }

        String[] seq = new String[nonNullTaxa];
        String[] qualScore = new String[nonNullTaxa];
        net.maizegenetics.pal.ids.Identifier[] taxaName = new net.maizegenetics.pal.ids.Identifier[nonNullTaxa];
        //System.out.println(" Genes = "+totalGenes+" Taxa = "+totalTaxa);

        int currentNullTaxa = 0;
        for (int j = 0; j < totalTaxa; j++) {
            if (table.get(genotypeExpIndex, j) != null) {
                Taxon t = (Taxon) table.getRowHeading(j);
                // taxaName[currentNullTaxa] = new net.maizegenetics.pal.ids.Identifier(t.getProperty(TaxonProperty.ACCESSION).toString() +  "_" +
                // t.getProperty(TaxonProperty.SOURCE).toString());
                taxaName[currentNullTaxa] = new net.maizegenetics.pal.ids.Identifier(t.toString());
                AlleleList alleleList = (AlleleList) table.get(genotypeExpIndex, j);
                ListIterator itr2 = alleleList.listIterator();
                while (itr2.hasNext()) {
                    Allele allele = (Allele) itr2.next();

                    seq[currentNullTaxa] = (String) allele.getProperty(AlleleProperty.VALUE);
                    // System.out.println("seq[" + j + "] = " + seq[j]);

                    // Quality Score
                    qualScore[currentNullTaxa] = (String) allele.getProperty(AlleleProperty.QUALITY);
                // System.out.println("qualScore[" + j + "] = " + qualScore[j]);
                }
                currentNullTaxa++;
            }
        }

        GenotypeExperiment theGenotypeExperiment = (GenotypeExperiment) table.getColumnHeading(genotypeExpIndex);
        String locusName = createLocusName(theGenotypeExperiment);

        //theNewAlignment = new SimpleAnnotatedAlignment(taxaName, seq, qualScore, "-", new Nucleotides());
        theNewAlignment = new SimpleAlignment(new SimpleIdGroup(taxaName), seq, new IUPACNucleotides(), null, null, null, locusName, null, null);

        // Old method call.  No longer used except for reference.
        // createMetaData(dataTable, genotypeExpIndex, theNewAlignment.getDataType(0).getDescription(), theNewAlignment.getSiteCount(), nonNullTaxa);

        // createComments(theLocus, locusName, theNewAlignment.getDataType(0).getDescription(), theNewAlignment.getSiteCount(), nonNullTaxa);

        return theNewAlignment;

    }

    public static String createLocusName(GenotypeExperiment theGenotypeExperiment) {

        //What Tassel call a locus is really a marker
        return theGenotypeExperiment.getProperty(GenotypeExperimentProperty.NAME) + ".ID" +
                theGenotypeExperiment.getProperty(GenotypeExperimentProperty.ID);

    }

    public static String createComments(Locus theLocus, String locusName, String dataTypeDesc, int siteNumber, int totalTaxa) {

        return "Raw Sequence \nNumber of Sequences: " + totalTaxa +
                "\nNumber of Sites: " + siteNumber + "\nDataType: " + dataTypeDesc +
                "\nGenotypeExperimentName: " + locusName +
                "\nLocusName: " + theLocus.getProperty(LocusProperty.NAME);

    }

    public static SimplePhenotype createCharacterAlignment(PhenotypeTable phenoTable) {

        SimplePhenotype theNewAlignment = null;
        int totalTraits = phenoTable.numColumns();
        int totalTaxa = phenoTable.numRows();
        String[] taxaName = new String[totalTaxa];
        //   String[] envID = new String[totalTraits];
        //   String[] traitNames = new String[totalTraits];
        String[] headerNames = new String[]{"Trait", "Env"};  //todo replace with reading more headers in future, rep, block, etc.
        String[][] headers = new String[2][totalTraits];
        ArrayList<Trait> traits = new ArrayList();
        double[][] data = new double[totalTaxa][totalTraits];
        SimpleIdGroup theSimpleIdGroup = null;

        for (int i = 0; i < totalTaxa; i++) {
            Taxon t = (Taxon) phenoTable.getRowHeading(i);
            taxaName[i] = t.toString();
            // taxaName[i]= t.getProperty(TaxonProperty.ACCESSION).toString();
            //  if(t.getProperty(TaxonProperty.SOURCE)!=null)
            //           taxaName[i]+="_"+ t.getProperty(TaxonProperty.SOURCE).toString();
            // System.out.println(taxaName[i]+"  "+i);
            /* taxaName[i] = new net.maizegenetics.pal.ids.Identifier(t.getProperty(TaxonProperty.NAME).
            toString()
            + ".ID" +
            t.getProperty(TaxonProperty.ID).
            toString()).toString();
             */
            for (int j = 0; j < totalTraits; j++) {
                headers[1][j] = phenoTable.getColumnHeading(j).toString();
                headers[0][j] = phenoTable.getPhenotypeOntology().toString();
                traits.add(new Trait(headers[0][j], false, "data"));
                if (phenoTable.get(j, i) == null) {
                    data[i][j] = Double.NaN;
                } else {
                    String valueString = phenoTable.get(j, i).toString();

                    if (valueString.indexOf(':') != -1) {
                        int count = 0;
                        double total = 0;
                        StringTokenizer st = new StringTokenizer(valueString, ":");
                        while (st.hasMoreElements()) {
                            double value = 0;
                            String element = (String) st.nextElement();
                            try {
                                value = Double.parseDouble(element);
                                count++;
                            } catch (NumberFormatException nfe) {
                                System.err.println("\n*****GDPCDataLoader:\n" + nfe);
                            }
                            total += value;
                        }
                        BigDecimal bd = new BigDecimal((double) total / count);
                        bd = bd.setScale(4, BigDecimal.ROUND_HALF_UP);
                        data[i][j] = bd.doubleValue();
                    } else {
                        data[i][j] = Double.parseDouble(valueString);
                    }
                }
            }
        }
        theSimpleIdGroup = new SimpleIdGroup(taxaName);
        //theNewAlignment = new SimplePhenotype(theSimpleIdGroup, data, headerNames, headers);
        theNewAlignment = new SimplePhenotype(theSimpleIdGroup, traits, data);
        //String Comments = "Number of Sequences: " + totalTaxa + "\nNumber of Traits: " + totalTraits;

        return theNewAlignment;

    }

    public static Alignment[] createMultiLocusAnnotatedAlignment(GenotypeTable dataTable, ArrayList polyList) {

        int totalTaxa = dataTable.numRows();
        int totalSites = dataTable.numColumns();

        int packSiteSize = (totalSites % 2 == 1) ? ((totalSites + 1) / 2) : (totalSites / 2); //number of bytes to pack bases, add an extra site if odd
        byte[][] alleleBLOB = new byte[totalTaxa][GdpdmBLOBUtils.totalHeaderPadding + packSiteSize];

        for (int j = 0; j < totalTaxa; j++) {
            Taxon t = (Taxon) dataTable.getRowHeading(j);
            GdpdmBLOBUtils.setHeaderOfBLOB(alleleBLOB[j], GdpdmBLOBUtils.compressionAlgorithm, GdpdmBLOBUtils.alleleBLOBtype,
                    totalSites, "ERROR",
                    "UNKNOWN", 0, totalSites - 1, t.toString().toUpperCase());
        }

        for (int i = 0; i < totalSites; i++) {

            GenotypeExperiment theGenotypeExperiment = (GenotypeExperiment) dataTable.getColumnHeading(i);
            Object theDT = theGenotypeExperiment.getProperty(GenotypeExperimentProperty.POLY_TYPE);

            if (!theDT.equals(GenotypeExperimentProperty.POLY_TYPE_SNP)) {
                throw new IllegalArgumentException("PALUtil: createMultiLocusAnnotatedAlignment: Must be SNP values.");
            }

            for (int j = 0; j < totalTaxa; j++) {
                AlleleList alleleList = (AlleleList) dataTable.get(i, j);

                if (alleleList != null) {
                    List uniqueAlleles = new ArrayList();
                    for (int k = 0; k < alleleList.size(); k++) {
                        Allele allele = (Allele) alleleList.get(k);
                        char snp = ((String) allele.getProperty(AlleleProperty.VALUE)).charAt(0);
                        if (!uniqueAlleles.contains(snp)) {
                            uniqueAlleles.add(snp);
                        }
                    }
                    if (uniqueAlleles.size() == 0) {
                        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(alleleBLOB[j], i, DataType.UNKNOWN_CHARACTER);
                    } else if (uniqueAlleles.size() == 1) {
                        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(alleleBLOB[j], i, ((Character) uniqueAlleles.get(0)).charValue());
                    } else {
                        byte[] snps = new byte[uniqueAlleles.size()];
                        for (int b = 0; b < uniqueAlleles.size(); b++) {
                            snps[b] = (byte) ((Character) uniqueAlleles.get(b)).charValue();
                        }
                        AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(alleleBLOB[j], i, (char) AllelePositionBLOBUtils.getBaseFromSNPValue(snps));
                    }
                } else {
                    AllelePositionBLOBUtils.setHalfByteInAlleleBLOB(alleleBLOB[j], i, DataType.UNKNOWN_CHARACTER);
                }
            }

        }

        Alignment[] theNewAlignment = new Alignment[1];
        theNewAlignment[0] = new Pack1Alignment(alleleBLOB);

        return theNewAlignment;

    }

}
