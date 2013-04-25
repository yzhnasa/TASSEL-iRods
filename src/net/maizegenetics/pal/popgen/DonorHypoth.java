/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

/**
 *
 * @author edbuckler
 */
class DonorHypoth {
    int targetTaxon = -1;
    int donor1Taxon = -1;
    int donor2Taxon = -1;
    int startBlock = -1;
    int focusBlock = -1;
    int endBlock = -1;
    double pError = 1;
    double pHeterozygous = -1;
    double pHomoD1 = -1;
    double pHomoD2 = -11;
    int totalSites = 0;
    int mendelianErrors = 0;

    public DonorHypoth() {
    }

    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, int focusBlock, int endBlock, int totalSites, int mendelianErrors) {
        this(targetTaxon, donor1Taxon, donor2Taxon, startBlock, focusBlock, endBlock);
        this.totalSites = totalSites;
        this.mendelianErrors = mendelianErrors;
    }

    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, int focusBlock, int endBlock) {
        this.targetTaxon = targetTaxon;
        if (donor1Taxon < donor2Taxon) {
            this.donor1Taxon = donor1Taxon;
            this.donor2Taxon = donor2Taxon;
        } else {
            this.donor1Taxon = donor2Taxon;
            this.donor2Taxon = donor1Taxon;
        }
        this.startBlock = startBlock;
        this.focusBlock = focusBlock;
        this.endBlock = endBlock;
    }

    public double getErrorRate() {
        return (double) mendelianErrors / (double) totalSites;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final DonorHypoth other = (DonorHypoth) obj;
        if (this.targetTaxon != other.targetTaxon) {
            return false;
        }
        if (this.donor1Taxon != other.donor1Taxon) {
            return false;
        }
        if (this.donor2Taxon != other.donor2Taxon) {
            return false;
        }
        if (this.focusBlock != other.focusBlock) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 83 * hash + this.targetTaxon;
        hash = 83 * hash + this.donor1Taxon;
        hash = 83 * hash + this.donor2Taxon;
        hash = 83 * hash + this.focusBlock;
        return hash;
    }

    public String toString() {
        return String.format("FTx:%d D1Tx:%d D2Tx:%d SBk:%d FBk:%d EBk:%d TS:%d MS:%d ", targetTaxon, donor1Taxon, donor2Taxon, startBlock, focusBlock, endBlock, totalSites, mendelianErrors);
    }
    
}