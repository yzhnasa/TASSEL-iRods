/*
 * MutableNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;


import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class MutableNucleotideAlignment extends MutableSingleEncodeAlignment implements MutableAlignment {

    private MutableNucleotideAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return new MutableNucleotideAlignment(a, maxTaxa, maxNumSites);
        } else {
            throw new IllegalArgumentException("MutableNucleotideAlignment: getInstance: alignment must be nucleotide data.");
        }
    }

    public static MutableNucleotideAlignment getInstance(Alignment a) {
        return getInstance(a, a.getSequenceCount(), a.getSiteCount());
    }

    private MutableNucleotideAlignment(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableNucleotideAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignment getInstance(IdGroup idGroup, int initNumSites) {
        return new MutableNucleotideAlignment(idGroup, initNumSites, idGroup.getIdCount(), initNumSites);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

}
