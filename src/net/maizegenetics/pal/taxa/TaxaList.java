/*
 *  TaxaList
 */
package net.maizegenetics.pal.taxa;

import java.util.List;
import java.util.Set;

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

    /**
     * Return a list of all matching taxa indices for a given name.  Matches will be done with both fullName and regular name.
     * @param name
     * @return set of indices for matching taxa (Set is empty if no match).
     */
    public Set<Integer> getTaxaMatchingIndices(String name);
}
