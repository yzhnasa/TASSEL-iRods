/*
 *  TaxaList
 */
package net.maizegenetics.pal.taxa;

import java.util.List;

/**
 *
 * @author terry
 */
public interface TaxaList extends List<AnnotatedTaxon> {

    /**
     * Returns number of sequences (taxa).
     *
     * @return number of sequences
     */
    public int getSequenceCount();

    /**
     * Returns number of taxa (same as getSequenceCount()
     *
     * @return number of taxa
     */
    public int getTaxaCount();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String getTaxaName(int index);

    /**
     * Return full taxa name at given index.
     *
     * @param index
     * @return full taxa name
     */
    public String getFullTaxaName(int index);
}
