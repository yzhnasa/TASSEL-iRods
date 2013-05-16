/*
 * TOPMInterface
 */
package net.maizegenetics.gbs.maps;

import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.pal.alignment.Locus;

/**
 *
 * @author edbuckler
 */
public interface TOPMInterface extends Tags {

    public int addVariant(int tagIndex, byte offset, byte base);

    public int compare(int index1, int index2);

    public int getChromosome(int index);

    /**
     * Returns an array whose <i>values</i> are the distinct chromosomes in this
     * file, as stored in the chromosome[] array. The indices are arbitrary.
     */
    public int[] getChromosomes();

    public byte getDcoP(int index);

    public byte getDivergence(int index);

    public int getEndPosition(int index);

    public Locus[] getLoci();

    public Locus getLocus(int tagIndex);

    public byte getMapP(int index);

    public byte getMultiMaps(int index);

    public int[] getPositionArray(int index);

    public int getReadIndexForPositionIndex(int posIndex);

    public int getSize();

    public int getStartPosition(int index);

    public byte getStrand(int tagIndex);

    public int getMaxNumVariants();

    public byte getVariantDef(int tagIndex, int variantIndex);

    /**
     * Returns an array containing all variant definitions for the tag at the
     * supplied index.
     */
    public byte[] getVariantDefArray(int tagIndex);

    public byte getVariantPosOff(int tagIndex, int variantIndex);

    /**
     * Returns an array containing all variant position offsets for the tag at
     * the supplied index.
     */
    public byte[] getVariantPosOffArray(int tagIndex);

    public byte[][] getVariantOff();

    public byte[][] getVariantDef();

    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax);

    public void setDivergence(int index, byte divergence);

    public void setMapP(int index, byte mapP);

    public void setMapP(int index, double mapP);

    public void setVariantDef(int tagIndex, int variantIndex, byte def);

    public void setVariantPosOff(int tagIndex, int variantIndex, byte offset);
    
    public void clearVariants();
}
