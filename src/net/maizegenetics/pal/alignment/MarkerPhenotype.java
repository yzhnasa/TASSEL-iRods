package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;

import java.io.PrintWriter;
import java.io.StringWriter;

import java.util.List;


/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Dec 5, 2006
 */
public class MarkerPhenotype extends AbstractAlignment implements Phenotype {

    //
    // Private stuff
    //
    private Alignment aAlign;  //alignment to hold markers
    private Phenotype cAlign;  //alignment to hold characters

    /**
     * concatenate alignments
     * @param aa
     * @param ca
     * @param union
     * @throws IllegalArgumentException
     */
    private MarkerPhenotype(Alignment aa, Phenotype ca)
            throws IllegalArgumentException {
        super(aa);
        aAlign = aa;
        cAlign = ca;
    }

    public static MarkerPhenotype getInstance(Alignment aa, Phenotype ca, boolean union) {
        IdGroup idGroup = getIdGroup(aa.getIdGroup(), ca.getTaxa(), union);
        Alignment align = FilterAlignment.getInstance(aa, idGroup);
        Phenotype phenotype = FilterPhenotype.getInstance(ca, idGroup, null);
        return new MarkerPhenotype(align, phenotype);
    }

    public static MarkerPhenotype getInstance(MarkerPhenotype aac, IdGroup group) {
        Alignment aa = FilterAlignment.getInstance(aac, group);
        Phenotype ca = FilterPhenotype.getInstance(aac, group, null);
        return new MarkerPhenotype(aa,ca);
    }

    private static IdGroup getIdGroup(IdGroup group1, IdGroup group2, boolean union) {
    	if (union) {
    		return IdGroupUtils.getAllIds(group1, group2);
    	} else {
    		return IdGroupUtils.getCommonIds(group1, group2);
    	}
    }

    /** sequence alignment at (sequence, site) */
    public char getBaseChar(int seq, int site) {
        return aAlign.getBaseChar(seq, site);
    }

    /** Return the position along the locus (ignores gaps) */
    public int getPositionInLocus(int site) {
        return aAlign.getPositionInLocus(site);
    }

    /** Returns position type (eg.  I=intron, E-exon, P=promoter, 1=first, 2=second, 3=third, etc.*/
    public byte getPositionType(int site) {
        return 0;
    }

    //implement Phenotype
    public double getData(int taxon, int trait) {
        return cAlign.getData(taxon, trait);
    }

    /*
     * Return number of sites for each sequence in this alignment
     */
    public int getNumberOfTraits() {
        return cAlign.getNumberOfTraits();
    }

    /**
     * @param trait	an integer j, representing the jth column in the data set
     * @return		the Trait for this column
     */
    public Trait getTrait(int trait) {
        return cAlign.getTrait(trait);
    }

    /**
     * @return the traits or columns of this dataset, in order by column number
     */
    public List<Trait> getTraits() {
        return cAlign.getTraits();
    }

    public float getSiteScore(int seq, int site) {
        return aAlign.getSiteScore(seq, site);
    }

    public boolean hasSiteScores() {
        return aAlign.hasSiteScores();
    }

    /**
     * @param factor
     * @return the name of this factor
     */
    public String getFactorName(int factor) {
        return cAlign.getFactorName(factor);
    }

    /**
     * @return
     */
    public int getNumberOfFactors() {
        return cAlign.getNumberOfFactors();
    }

    /**
     * Use a copy of the factor names to create a new Phenotype to avoid making unintended changes to the original
     *
     * @return
     */
    public String[] getFactorNameCopy() {
        return cAlign.getFactorNameCopy();
    }

    // interface Report
    //public int[] getActiveFactorCopy() {
    //    return cAlign.getActiveFactorCopy();
    //}

    //public String getActiveFactorName(int factor) {
    //    return cAlign.getActiveFactorName(factor);
    //}

    //public String[] getActiveFactorNames() {
    //    return cAlign.getActiveFactorNames();
    //}

    //public int getNumberOfActiveFactors() {
    //    return cAlign.getNumberOfActiveFactors();
    //}

    //public boolean isFactorActive(int factor) {
    //    return cAlign.isFactorActive(factor);
    //}

    //public void setActiveFactor(int factor, boolean active) {
    //    cAlign.setActiveFactor(factor, active);
    //}
    /**
     * sort the sites by chromosome, then by chromosomal location, and final locusPosition
     */
    //public void report(PrintWriter out) {
    //    AlignmentUtils.report(this, out);
    //    out.println();
    //    aAlign.report(out);
    //    out.println();
    //    cAlign.report(out);
    //    out.println();
    //}
    /** returns representation of this alignment as a string */
    public String toString() {

        StringWriter sw = new StringWriter();

        PrintWriter out = new PrintWriter(sw);
        out.println("  " + getSequenceCount() + " " + getNumberOfTraits());
        out.print("Taxa/Trait\t");
        for (int j = 0; j < getNumberOfTraits(); j++) {
            out.print(getTrait(j).getName() + "\t");
        }
        out.println();
        //if (getEnvironmentName(0) != null) {
        //    out.print("Taxa/Environ\t");
        //    for (int j = 0; j < getNumberOfTraits(); j++) {
        //        out.print(getEnvironmentName(j) + "\t");
        //    }
        //    out.println();
        //}
        for (int i = 0; i < getSequenceCount(); i++) {
            out.print(getIdGroup().getIdentifier(i).getName() + "\t");
            for (int j = 0; j < getNumberOfTraits(); j++) {
                out.print(getData(i, j) + "\t");
            }
            out.println(getAlignedSequenceString(i));
        }
        return sw.toString();
    }

    //Implementation of TableReport Interface
    /**
     * Return column names for the table
     */
    public Object[] getTableColumnNames() {
        String[] basicLabels = new String[getNumberOfTraits() + 2];
        basicLabels[0] = "Taxa";
        for (int c = 0; c < getNumberOfTraits(); c++) {
            //basicLabels[c + 1] = getTrait(c).getName() + "." + getEnvironmentName(c);
            basicLabels[c + 1] = getTrait(c).getName();
        }
        basicLabels[getNumberOfTraits() + 1] = "Haplotype";
        return basicLabels;
    }

    /**
     * get the data elements
     *
     * @return the data elements
     */
    public Object[][] getTableData() {
        return getTableData(0, getRowCount() - 1);
    }

    /**
     * Get Table Data in specified range inclusive.
     * Indices start at 0.
     *
     * @param start start position
     * @param end end position
     * @return
     */
    public Object[][] getTableData(int start, int end) {

        if ((start < 0) || (end >= getRowCount())) {
            throw new IndexOutOfBoundsException("getTableData: start: " + start + "  end: " + end);
        }

        if (end < start) {
            return null;
        }

        Object[][] temp = new Object[end - start + 1][];

        for (int i = start; i <= end; i++) {
            temp[i] = getRow(i);
        }

        return temp;

    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        Object[] data;
        int haplotypeColumn = getNumberOfTraits() + 1;
        data = new String[getNumberOfTraits() + 2];
        data[0] = getIdGroup().getIdentifier(row).getName();
        for (int c = 0; c < getNumberOfTraits(); c++) {
            data[c + 1] = "" + getData(row, c);
        }
        int siteCount = Math.min(getSiteCount(), 10);
        StringBuilder builder = new StringBuilder();
        builder.append(getAlignedSequenceString(row, 0, siteCount));
        if (getSiteCount() > 10) {
            builder.append("...");
        }
        data[haplotypeColumn] = builder.toString();
        //data[haplotypeColumn] = getFormattedSequenceString(row);
        return data;

    }

    /**
     * Return the name for the title of the ANOVA
     */
    public String getTableTitle() {
        return "Phenotypes and Genotypes";
    }

    /**
     *
     * @return  the maximum number of factors or header rows in any of the character Alignments
     * that make up this concantenated alignment.
     */

    // public int getMaxFactors() {
    //  return cAlign.get
    // }
    /**
     * Provides a method for getting a copy all the trait values
     *
     * @return
     */
    //public double[][] getTraitValues() {
    //    return cAlign.getTraitValues();
    //}
    public int getRowCount() {
        return getSequenceCount();
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        return getNumberOfTraits() + 2;
    }

    public Alignment getAlignment() {
        return aAlign;
    }

    public Phenotype getPhenotype() {
        return cAlign;
    }

    public byte getBase(int taxon, int site) {
        return aAlign.getBase(taxon, site);
    }

    public byte getBase(int taxon, int site, int allele) {
        return aAlign.getBase(taxon, site, allele);
    }

    public byte getBase(int taxon, Locus locus, int physicalPosition, int allele) {
        return aAlign.getBase(taxon, locus, physicalPosition, allele);
    }

    public String[] getSNPIDs() {
        return aAlign.getSNPIDs();
    }

    public String getSNPID(int site) {
        return aAlign.getSNPID(site);
    }

    public int getSiteCount() {
        return aAlign.getSiteCount();
    }

    public int getLocusSiteCount(Locus locus) {
        return aAlign.getLocusSiteCount(locus);
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        return aAlign.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    public byte[] getPositionTypes() {
        return aAlign.getPositionTypes();
    }

    public Locus getLocus(int site) {
        return aAlign.getLocus(site);
    }

    public Locus[] getLoci() {
        return aAlign.getLoci();
    }

    public int getNumLoci() {
        return aAlign.getNumLoci();
    }

    public double getData(Identifier taxon, Trait trait) {
        return cAlign.getData(taxon, trait);
    }

    public void setData(int taxon, int trait, double value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setData(Identifier taxon, Trait trait, double value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int whichTrait(Trait trait) {
        return whichTrait(trait);
    }

    public int whichTaxon(Identifier taxon) {
        return cAlign.whichTaxon(taxon);
    }

    public Identifier getTaxon(int taxon) {
        return cAlign.getTaxon(taxon);
    }

    public int getNumberOfTaxa() {
        return cAlign.getNumberOfTaxa();
    }

    public IdGroup getTaxa() {
    	if (aAlign != null) return aAlign.getIdGroup();
    	else return cAlign.getTaxa();
    }

    public double[][] getData() {
        return cAlign.getData();
    }

    public void setActiveFactor(int factor, boolean active) {
        cAlign.setActiveFactor(factor, active);
    }

    public boolean isFactorActive(int factor) {
        return cAlign.isFactorActive(factor);
    }

    public String getActiveFactorName(int factor) {
        return cAlign.getActiveFactorName(factor);
    }

    public String[] getActiveFactorNames() {
        return cAlign.getActiveFactorNames();
    }

    public int getNumberOfActiveFactors() {
        return cAlign.getNumberOfActiveFactors();
    }

	@Override
	public String[] getRowNames() {
		return cAlign.getRowNames();
	}

	@Override
	public Object getValueAt(int row, int col) {
        int haplotypeColumn = cAlign.getColumnCount();
        if (col == haplotypeColumn) {
            int siteCount = Math.min(getSiteCount(), 10);
            StringBuilder builder = new StringBuilder();
            builder.append(getAlignedSequenceString(row, 0, siteCount));
            if (getSiteCount() > 10) {
                builder.append("...");
            }
            return builder.toString();
        }
        return cAlign.getValueAt(row, col) ;
	}

}
