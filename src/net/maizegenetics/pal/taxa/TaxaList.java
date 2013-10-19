/*
 *  TaxaList
 */
package net.maizegenetics.pal.taxa;

import java.util.List;

/**
 *
 * @author terry
 */
public interface TaxaList extends List<Taxon> {

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
     * Return a list of all matching taxa indices for a given name. Matches will
     * depend on the Tassel Preference ID Join Strict.
     *
     * @param name name
     *
     * @return Indices for matching taxa (Empty if no match).
     */
    public List<Integer> getIndicesMatchingTaxon(String name);
    
    /**
     * Return a list of all matching taxa indices for a given name. Matches will
     * depend on the Tassel Preference ID Join Strict.
     *
     * @param taxon taxon
     *
     * @return Indices for matching taxa (Empty if no match).
     */
    public List<Integer> getIndicesMatchingTaxon(Taxon taxon);
}
