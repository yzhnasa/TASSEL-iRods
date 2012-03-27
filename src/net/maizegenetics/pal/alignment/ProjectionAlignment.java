/*
 * SBitNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.OpenBitSet;

/**
 * This class projects high Density markers on large group of taxa through 
 * a look up table system.  The lookup table generally needs be built through
 * some imputation approach.
 * @author ed
 */
public class ProjectionAlignment extends AbstractAlignment {
    private int[][]  siteBreaks;  //temporary - saving not needed
    private int[][]  hdTaxa;  //taxa ids should be saved
    private int[][] posBreaks;  //positions should be saved
    private SBitAlignment hdAlign;  //high density marker alignment that is being projected.
    private OpenBitSet[] currData;  //this should be a buffer 
    private IdGroup myIdGroup;  //taxa ids should saved.

    public ProjectionAlignment(Alignment hdAlign, IdGroup ldIDGroup, int maxNumAlleles, boolean retainRareAlleles) {
        super(hdAlign, maxNumAlleles, retainRareAlleles);
        if(hdAlign instanceof SBitAlignment) {this.hdAlign=(SBitAlignment)hdAlign;}  //this is taking three minutes to convert.
        else {this.hdAlign=new SBitAlignment(hdAlign,2,false);}
        myIdGroup=ldIDGroup;  //this is messy there is one myIdGroup defined in the super and one here.
        this.
        siteBreaks=new int[getSequenceCount()][];
        posBreaks=new int[getSequenceCount()][];
        hdTaxa=new int[getSequenceCount()][];
        long currentTime = System.currentTimeMillis();
        
        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        System.out.println("Time to load alleles: " + ((currentTime - prevTime) / 1000));
    }
    
    public void setCompositionOfTaxon(int taxon, int[]  posBreaks, int[]  hdTaxa) {
        this.posBreaks[taxon]=posBreaks;
        this.hdTaxa[taxon]=hdTaxa;
        this.siteBreaks[taxon]=new int[this.posBreaks[taxon].length];
        for (int i = 0; i < this.posBreaks[taxon].length; i++) {
            int site=hdAlign.getSiteOfPhysicalPosition(posBreaks[i], null);
            if(site<0) site=-(site+1);  
            this.siteBreaks[taxon][i]=site;
        }
    }
    
    public void setCompositionOfTaxon(String taxa, int[]  posBreaks, int[]  hdTaxa) {
        int taxon=myIdGroup.whichIdNumber(taxa);
        setCompositionOfTaxon(taxon, posBreaks, hdTaxa);
    }
    
    @Override
    public int getSequenceCount() {
        return myIdGroup.getIdCount();
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }
    
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    /**
     * This is the slow implementation of this.  Most of these should be buffered bit
     * sets.  Note the imputation of the taxon is likely to be the the same for 64 or more sites
     * in a row (potentially, 10,000s of sites in many cases).
     */
    public byte getBase(int taxon, int site) {
        int b=Arrays.binarySearch(siteBreaks[taxon], site);
        if(b<0) b=-(b+1);
        return hdAlign.getBase(hdTaxa[taxon][b], site);
    }

    @Override
    public boolean isSBitFriendly() {
        return true;
    }

    @Override
    public boolean isTBitFriendly() {
        return false;
    }

}
