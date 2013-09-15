package net.maizegenetics.pal.taxa;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;

import java.util.*;

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

    private final List<Taxon> myTaxaList;

    public TaxaListBuilder() {
        myTaxaList = new ArrayList<Taxon>();
    }

    public TaxaListBuilder add(Taxon taxon) {
        myTaxaList.add(taxon);
        return this;
    }

//    public TaxaListBuilder addAll(Collection<Taxon> taxa) {
//        myTaxaList.addAll(taxa);
//        return this;
//    }

    public TaxaListBuilder addAll(Collection<Taxon> taxa) {
        for (Taxon t: taxa) {
            myTaxaList.add(new Taxon.Builder(t).build());
        }
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

    /*Sort the taxa by their natural order (alphabetically by name)*/
    public TaxaListBuilder sort() {
        Collections.sort(myTaxaList);
        return this;
    }

    public TaxaList build() {
        return new TaxaArrayList(this);
    }

    public TaxaList buildFromHDF5(String hdf5FileName) {
        IHDF5Reader reader = HDF5Factory.openForReading(hdf5FileName);
        myTaxaList.clear();
        List<HDF5LinkInformation> fields = reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        for (HDF5LinkInformation is : fields) {
            if (is.isDataSet() == false) {
                continue;
            }
            myTaxaList.add(new Taxon.Builder(is.getName()).build());
        }
        sort();
        return build();
    }

    //Default package private method to hand the list to the instance
    List<Taxon> getImmutableList() {
        return Collections.unmodifiableList(myTaxaList);
    }
}
