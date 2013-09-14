package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.Taxon;

import java.util.ArrayList;
import java.util.List;

public class FilterPhenotype extends AbstractPhenotype {

    private Phenotype thePhenotype;
    private int[] taxaIndex = null;
    private int[] traitIndex = null;
    private int numberOfRows = 0;
    private int numberOfColumns = 0;

    private FilterPhenotype(Phenotype phenotype, int[] taxa, TaxaList taxaGroup, int[] traits, List<Trait> traitList) {
        super(taxaGroup, traitList);
        thePhenotype = phenotype;
        taxaIndex = taxa;
        traitIndex = traits;
        numberOfRows = taxa.length;
        numberOfColumns = traits.length;
    }

    /**
     * @param phenotype the input Phenotype
     * @param taxa	an index of taxa to include in the filtered subset. If null, all taxa will be included
     * @param traits an index of traits to include in the filtered subset. If null, all traits will be included.
     * @return a new FilterPhenotype 
     * the indices for taxa and traits must be valid indexes in the original phenotype. 
     * That is, only taxa and traits in the original phenotype can be used. This function 
     * is not appropriate for union joins.
     */ 
    public static FilterPhenotype getInstance(Phenotype phenotype, int[] taxa, int[] traits) {
        TaxaList taxaGroup;
        List<Trait> traitList;

        if (taxa == null) {
            taxaGroup = phenotype.getTaxa();
            int n = taxaGroup.getIdCount();
            taxa = new int[n];
            for (int i = 0; i < n; i++) {
                taxa[i] = i;
            }
        } else {
            int n = taxa.length;
            Taxon[] ids = new Taxon[n];
            for (int i = 0; i < n; i++) {
                ids[i] = phenotype.getTaxon(taxa[i]);
            }
            taxaGroup = new SimpleIdGroup(ids);
        }

        if (traits == null) {
            traitList = copyTraitsFromPhenotype(phenotype);
            int n = traitList.size();
            traits = new int[n];
            for (int i = 0; i < n; i++) {
                traits[i] = i;
            }
        } else {
            traitList = new ArrayList<Trait>();
            for (int t : traits) {
                traitList.add(Trait.getInstance(phenotype.getTrait(t)));
            }
        }

        return new FilterPhenotype(phenotype, taxa, taxaGroup, traits, traitList);
    }

    /**
     * @param phenotype the input Phenotype
     * @param taxa the taxa to be included in the output FilterPhenotype
     * @param traits the traits to be included in the output FilterPhenotyp
     * @return a FilterPhenotype built from the original Phenotype for the specified traits and taxa
     * Taxa and traits not included in the original Phenotype will return missing values for data. 
     * Data values cannot be set for these traits and taxa as the underlying data matrix will not contain storage for them.
     */
    public static FilterPhenotype getInstance(Phenotype phenotype, TaxaList taxa, List<Trait> traits) {
        TaxaList taxaGroup;
        List<Trait> traitList;
        int[] taxaIndex;
        int[] traitIndex;

        if (taxa == null) {
            taxaGroup = phenotype.getTaxa();
            int n = taxaGroup.getIdCount();
            taxaIndex = new int[n];
            for (int i = 0; i < n; i++) {
                taxaIndex[i] = i;
            }
        } else {
            taxaGroup = taxa;
            int n = taxaGroup.getIdCount();
            taxaIndex = new int[n];
            for (int i = 0; i < n; i++) {
                taxaIndex[i] = phenotype.whichTaxon(taxaGroup.getIdentifier(i));
            }
        }

        if (traits == null) {
            traitList = copyTraitsFromPhenotype(phenotype);
            int n = traitList.size();
            traitIndex = new int[n];
            for (int i = 0; i < n; i++) {
                traitIndex[i] = i;
            }
        } else {
            traitList = traits;
            int n = traitList.size();
            traitIndex = new int[n];
            for (int i = 0; i < n; i++) {
                traitIndex[i] = phenotype.whichTrait(traitList.get(i));
            }
        }

        return new FilterPhenotype(phenotype, taxaIndex, taxaGroup, traitIndex, traitList);
    }

    /**
     * @param phenotype the input Phenotype
     * @param taxa the taxa to be excluded in the output FilterPhenotype
     * @return a FilterPhenotype with excluded IDs
     */
    public static FilterPhenotype getInstanceRemoveIDs(Phenotype phenotype, TaxaList taxa) {
        List result = new ArrayList();
        TaxaList current = phenotype.getTaxa();
        for (int i = 0, n = current.getIdCount(); i < n; i++) {
            if (taxa.whichIdNumber(current.getIdentifier(i)) == -1) {
                result.add(current.getIdentifier(i));
            }
        }
        Taxon[] ids = new Taxon[result.size()];
        result.toArray(ids);
        return FilterPhenotype.getInstance(phenotype, new SimpleIdGroup(ids), null);
    }

    //implement Phenotype interface
    public double getData(int taxon, int trait) {
    	if (taxon == -1 || trait == -1) return Double.NaN;
    	int taxonNdx = taxaIndex[taxon];
    	int traitNdx = traitIndex[trait];
    	if (taxonNdx < 0 || taxonNdx >= thePhenotype.getNumberOfTaxa()) return Double.NaN;
    	if (traitNdx < 0 || traitNdx >= thePhenotype.getNumberOfTraits()) return Double.NaN;
        return thePhenotype.getData(taxaIndex[taxon], traitIndex[trait]);
    }

    public double getData(Taxon taxon, Trait trait) {
        return thePhenotype.getData(whichTaxon(taxon), whichTrait(trait));
    }

    public double[][] getData() {
        double[][] result = new double[numberOfRows][numberOfColumns];
        for (int r = 0; r < numberOfRows; r++) {
            for (int c = 0; c < numberOfColumns; c++) {
                result[r][c] = getData(r, c);
            }
        }
        return result;
    }

    public void setData(int taxon, int trait, double value) {
        thePhenotype.setData(taxaIndex[taxon], traitIndex[trait], value);
    }

    public void setData(Taxon taxon, Trait trait, double value) {
        setData(whichTaxon(taxon), whichTrait(trait), value);
    }

    public SimplePhenotype simpleCopy() {
        return new SimplePhenotype(getTaxa(), getTraits(), getData());
    }

}
