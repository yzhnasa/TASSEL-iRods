/*
 *  TaxaList
 */
package net.maizegenetics.taxa;

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
    public int numberOfTaxa();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String taxaName(int index);

    /**
     * Return a list of all matching taxa indices for a given name. Matches will
     * depend on the Tassel Preference ID Join Strict.
     *
     * @param name name
     *
     * @return Indices for matching taxa (Empty if no match).
     */
    public List<Integer> indicesMatchingTaxon(String name);
    
    /**
     * Return a list of all matching taxa indices for a given name. Matches will
     * depend on the Tassel Preference ID Join Strict.
     *
     * @param taxon taxon
     *
     * @return Indices for matching taxa (Empty if no match).
     */
    public List<Integer> indicesMatchingTaxon(Taxon taxon);
}
