/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import net.maizegenetics.pal.alignment.Locus;

/**
 *
 * @author edbuckler
 */
public interface TOPMInterface {

    int addVariant(int tagIndex, byte offset, byte base);

    int compare(int index1, int index2);

    int getChromosome(int index);

    /**
     * Returns an array whose <i>values</i> are the distinct chromosomes
     * in this file, as stored in the chromosome[] array.  The indices are arbitrary.
     */
    int[] getChromosomes();

    byte getDcoP(int index);

    byte getDivergence(int index);

    int getEndPosition(int index);

    Locus[] getLoci();

    Locus getLocus(int tagIndex);

    byte getMapP(int index);

    byte getMultiMaps(int index);

    int[] getPositionArray(int index);

    int getReadIndexForPositionIndex(int posIndex);

    int getSize();

    int getStartPosition(int index);

    byte getStrand(int tagIndex);

    byte getVariantDef(int tagIndex, int variantIndex);

    /**
     * Returns an array containing all variant definitions for the
     * tag at the supplied index.
     */
    byte[] getVariantDefArray(int tagIndex);

    byte getVariantPosOff(int tagIndex, int variantIndex);

    /**
     * Returns an array containing all variant position offsets for the
     * tag at the supplied index.
     */
    byte[] getVariantPosOffArray(int tagIndex);

    void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax);

    void setDivergence(int index, byte divergence);

    void setMapP(int index, byte mapP);

    void setMapP(int index, double mapP);

    void setVariantDef(int tagIndex, int variantIndex, byte def);

    void setVariantPosOff(int tagIndex, int variantIndex, byte offset);
    
}
