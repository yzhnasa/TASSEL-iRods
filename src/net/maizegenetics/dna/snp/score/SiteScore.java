/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.GenotypeTable;

/**
 *
 * @author Terry Casstevens
 */
public interface SiteScore {

    /**
     * Returns the site score of the given sequence and site.
     *
     * @param taxon taxon index
     * @param site site
     *
     * @return site score.
     */
    public float siteScore(int taxon, int site);

    /**
     * Returns the site scores.
     *
     * @return site scores.
     */
    public float[][] siteScores();

    /**
     * Returns true if this alignment has site scores.
     *
     * @return true if this alignment has site scores.
     */
    public boolean hasSiteScores();

    /**
     * Return what type of site scores this alignment has.
     *
     * @return site score type.
     */
    public GenotypeTable.SITE_SCORE_TYPE siteScoreType();
}
