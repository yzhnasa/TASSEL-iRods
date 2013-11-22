package net.maizegenetics.pal.taxa;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.dna.snp.genotype.GenotypeBuilder;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * A builder for creating immutable {@link TaxaList} instances.
 * Example:<pre>   {@code
 *   TaxaListBuilder tlb=new TaxaListBuilder();
 *   for (int i = 0; i < 10; i++) {
 *      Taxon at= new Taxon.Builder("Z"+i+":Line:mays:Zea")
 *           .inbreedF(0.99f)
 *           .parents("B73","B97")
 *           .pedigree("(B73xB97)S6I1")
 *           .build();
 *       tlb.add(at);
 *       }
 *   TaxaList tl=tlb.build();}</pre>
 *   <p></p>
 *   If building from HDF5:<pre>
 *   {@code
 *   TaxaList tl=new TaxaListBuilder().buildFromHDF5(testMutFile);
 *   }</pre>
 *
 * @author Ed Buckler
 */
public class TaxaListBuilder {
    //TODO need to move union and intersection utils to the builder

    private List<Taxon> myTaxaList;

    public TaxaListBuilder() {
        myTaxaList = new ArrayList<Taxon>();
    }

    public TaxaListBuilder add(Taxon taxon) {
        myTaxaList.add(taxon);
        return this;
    }

    public TaxaListBuilder addAll(Collection<Taxon> taxa) {
        myTaxaList.addAll(taxa);
        return this;
    }

    public TaxaListBuilder addAll(Alignment a) {
        myTaxaList.addAll(a.getTaxaList());
        return this;
    }

    public TaxaListBuilder addAll(String[] taxa) {
        for (int i = 0, n = taxa.length; i < n; i++) {
            myTaxaList.add(new Taxon.Builder(taxa[i]).build());
        }
        return this;
    }

    public TaxaListBuilder addAll(Taxon[] taxa) {
        for (int i = 0, n = taxa.length; i < n; i++) {
            myTaxaList.add(new Taxon.Builder(taxa[i]).build());
        }
        return this;
    }

    public TaxaList build() {
        return new TaxaArrayList(this);
    }

    public TaxaList buildFromHDF5(IHDF5Reader reader) {
        //IHDF5Reader reader = HDF5Factory.openForReading(hdf5FileName);
        myTaxaList.clear();
        List<HDF5LinkInformation> fields = reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        for (HDF5LinkInformation is : fields) {
            if (is.isDataSet() == false) {
                continue;
            }
            myTaxaList.add(new Taxon.Builder(is.getName()).build());
        }
        return build();
    }

    //Default package private method to hand the list to the instance
    List<Taxon> getImmutableList() {
        return Collections.unmodifiableList(myTaxaList);
    }

    public TaxaListBuilder sortTaxaAlphabetically(GenotypeBuilder genotypes) {
        int numTaxa = myTaxaList.size();
        if (numTaxa != genotypes.getTaxaCount()) {
            throw new IllegalArgumentException("TaxaListBuilder: sortTaxaAlphabetically: taxa list size: " + numTaxa + " doesn't match genotypes num taxa: " + genotypes.getTaxaCount());
        }
        genotypes.reorderTaxa(sortAlphabetically());
        return this;
    }

    public TaxaListBuilder sortTaxaAlphabetically() {
        sortAlphabetically();
        return this;
    }

    private int[] sortAlphabetically() {

        int numTaxa = myTaxaList.size();

        final int indicesOfSortByTaxa[] = new int[numTaxa];
        for (int i = 0; i < indicesOfSortByTaxa.length; i++) {
            indicesOfSortByTaxa[i] = i;
        }

        Swapper swapTaxa = new Swapper() {
            @Override
            public void swap(int a, int b) {
                int temp = indicesOfSortByTaxa[a];
                indicesOfSortByTaxa[a] = indicesOfSortByTaxa[b];
                indicesOfSortByTaxa[b] = temp;
            }
        };

        IntComparator compTaxa = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return myTaxaList.get(indicesOfSortByTaxa[a]).compareTo(myTaxaList.get(indicesOfSortByTaxa[b]));
            }
        };

        GenericSorting.quickSort(0, indicesOfSortByTaxa.length, compTaxa, swapTaxa);

        List<Taxon> temp = new ArrayList<>(numTaxa);
        for (int t = 0; t < numTaxa; t++) {
            temp.add(myTaxaList.get(indicesOfSortByTaxa[t]));
        }

        myTaxaList = temp;

        return indicesOfSortByTaxa;

    }
}
