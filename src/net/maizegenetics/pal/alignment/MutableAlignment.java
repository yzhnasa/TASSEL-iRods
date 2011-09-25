/*
 * MutableAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public interface MutableAlignment extends Alignment {

    /**
     * Sets base at given taxon and site.  newBase should
     * contain both diploid allele values (one in first four bits
     * and second in last four bits).  This method may not be
     * supported by every implementation.
     *
     * @param taxon taxon
     * @param site site
     * @param newBase new diploid allele values
     */
    public void setBase(int taxon, int site, byte newBase);

    /**
     * Sets bases starting at specified site for given taxon.
     * Each base should contain both diploid allele values
     * (one in first four bits and second in last four bits).
     * This method may not be supported by every implementation.
     *
     * @param taxon taxon
     * @param startSite starting site
     * @param newBases new diploid allele values
     */
    public void setBaseRange(int taxon, int startSite, byte[] newBases);

    public void addSite(int site);

    public void removeSite(int site);

    public void addTaxon(Identifier id);

    public void removeTaxon(int taxon);

    /**
     * Clean alignment including sorting sites by position.
     */
    public void clean();

    /**
     * True if changes since last clean().
     *
     * @return true if dirty.
     */
    public boolean isDirty();
}
