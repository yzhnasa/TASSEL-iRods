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
     * Returns number of taxa
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
     * Return a list of all matching taxa indices for a given name.
     *
     * @param name name
     *
     * @return Indices for matching taxa (-1 if no match).
     */
    public int indexOf(String name);
    
    /**
     * Return a list of all matching taxa indices for a given name.
     *
     * @param taxon taxon
     *
     * @return Indices for matching taxa (-1 if no match).
     */
    public int indexOf(Taxon taxon);
}
